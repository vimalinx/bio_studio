from __future__ import annotations

import importlib
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
PLAYGROUND_PIPELINE = ROOT / "projects" / "ai_design_playground" / "scripts" / "pipeline.py"


def _runtime():
    return importlib.import_module("lib.template_runtime")


def _config(tmp_path: Path) -> SimpleNamespace:
    project_root = tmp_path / "ai_design_playground"
    raw_dir = project_root / "data" / "raw"
    processed_dir = project_root / "data" / "processed"
    results_dir = project_root / "data" / "results"
    references_dir = project_root / "data" / "references"
    logs_dir = project_root / "logs"

    for path in (raw_dir, processed_dir, results_dir, references_dir, logs_dir):
        path.mkdir(parents=True, exist_ok=True)

    return SimpleNamespace(
        PROJECT_NAME="ai_design_playground",
        PROJECT_TYPE="generic",
        PROJECT_ROOT=project_root,
        RAW_DIR=raw_dir,
        PROCESSED_DIR=processed_dir,
        RESULTS_DIR=results_dir,
        REFERENCES_DIR=references_dir,
        LOGS_DIR=logs_dir,
        PRODIGAL_MODE="meta",
    )


def test_modules_package_exports_sequence_analysis_module() -> None:
    modules = importlib.import_module("lib.modules")
    assert hasattr(modules, "sequence_analysis")


def test_template_runtime_supports_ai_design_playground_analysis(
    tmp_path: Path,
    monkeypatch,
) -> None:
    runtime = _runtime()
    config = _config(tmp_path)
    dna_fa = Path(config.RAW_DIR) / "toy_dna.fa"
    dna_fa.write_text(">toy_001\nATGCATGC\n", encoding="utf-8")

    calls: list[str] = []

    def fake_run_seqkit_stats(*args, **kwargs):
        calls.append("seqkit")
        Path(kwargs["output_txt"]).write_text("stats", encoding="utf-8")
        return SimpleNamespace(stdout="stats\n", returncode=0)

    def fake_predict_genes(*args, **kwargs):
        calls.append("prodigal")
        Path(kwargs["output_faa"]).write_text(">orf1\nMST*\n", encoding="utf-8")
        Path(kwargs["output_fna"]).write_text(">orf1\nATGTCTACT\n", encoding="utf-8")
        Path(kwargs["output_gff"]).write_text("seq\tProdigal\tCDS\t1\t9\t.\t+\t0\tID=orf1\n", encoding="utf-8")
        return SimpleNamespace(returncode=0)

    def fake_run_rnafold(*args, **kwargs):
        calls.append("rnafold")
        return {"structure": "....", "mfe": -3.2}

    monkeypatch.setattr(runtime.sequence_analysis, "run_seqkit_stats", fake_run_seqkit_stats)
    monkeypatch.setattr(runtime.gene_prediction, "predict_genes", fake_predict_genes)
    monkeypatch.setattr(runtime.sequence_analysis, "run_rnafold", fake_run_rnafold)

    result = runtime.run_ai_design_playground_analysis(
        config,
        dna_fa,
        [("toy_001", "AUGCAUGC"), ("toy_002", "AUGG")],
    )

    assert calls == ["seqkit", "prodigal", "rnafold", "rnafold"]
    assert Path(result["seqkit_stats"]).name == "seqkit_stats.txt"
    assert Path(result["prodigal"]["proteins"]).name == "prodigal_proteins.faa"
    assert len(result["rnafold"]) == 2
    assert result["rnafold"][0]["id"] == "toy_001"


def test_ai_design_playground_pipeline_prefers_shared_runtime() -> None:
    pipeline_text = PLAYGROUND_PIPELINE.read_text(encoding="utf-8")

    assert "load_template_runtime" in pipeline_text
    assert "run_shared_ai_design_analysis" in pipeline_text
