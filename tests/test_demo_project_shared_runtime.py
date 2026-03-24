from __future__ import annotations

import importlib
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
YEAST_PIPELINE = ROOT / "projects" / "yeast_rnaseq_demo" / "scripts" / "pipeline.py"


def _runtime():
    return importlib.import_module("lib.template_runtime")


def _rnaseq_config(tmp_path: Path) -> SimpleNamespace:
    project_root = tmp_path / "yeast_like_demo"
    raw_dir = project_root / "data" / "raw"
    processed_dir = project_root / "data" / "processed"
    results_dir = project_root / "data" / "results"
    references_dir = project_root / "data" / "references"
    logs_dir = project_root / "logs"

    for path in (raw_dir, processed_dir, results_dir, references_dir, logs_dir):
        path.mkdir(parents=True, exist_ok=True)

    star_index_dir = references_dir / "star_index"

    return SimpleNamespace(
        PROJECT_NAME="yeast_like_demo",
        PROJECT_TYPE="rnaseq",
        PROJECT_ROOT=project_root,
        RAW_DIR=raw_dir,
        PROCESSED_DIR=processed_dir,
        RESULTS_DIR=results_dir,
        REFERENCES_DIR=references_dir,
        LOGS_DIR=logs_dir,
        THREADS=8,
        SAMPLES=["sample_a"],
        REFERENCE_GENOME=references_dir / "genome.fa",
        ANNOTATION_GTF=references_dir / "genome.gff",
        GENOME_FASTA=references_dir / "genome.fa",
        GENOME_GTF=references_dir / "genome.gff",
        STAR_INDEX_DIR=star_index_dir,
        READ1_PATTERN="*_R1.fastq.gz",
        READ2_PATTERN="*_R2.fastq.gz",
        STAR_INDEX_EXTRA_ARGS=["--genomeSAindexNbases", "10"],
        STAR_ALIGN_EXTRA_ARGS=["--outSAMtype", "BAM", "SortedByCoordinate"],
        FEATURECOUNTS_FEATURE_TYPE="gene",
        FEATURECOUNTS_ATTRIBUTE_TYPE="ID",
    )


def test_template_runtime_supports_rnaseq_build_index_step(
    tmp_path: Path,
    monkeypatch,
) -> None:
    runtime = _runtime()
    config = _rnaseq_config(tmp_path)
    Path(config.GENOME_FASTA).write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(config.GENOME_GTF).write_text("chr1\ttest\tgene\t1\t4\t.\t+\t.\tID=gene1\n", encoding="utf-8")

    calls: list[tuple[tuple, dict]] = []

    def fake_build_star_index(*args, **kwargs):
        calls.append((args, kwargs))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(runtime.rnaseq, "build_star_index", fake_build_star_index)

    assert runtime.run_shared_step(config, "build_index") is True
    assert calls
    assert Path(calls[0][0][1]).name == "star_index"
    assert calls[0][1]["extra_args"] == ["--genomeSAindexNbases", "10"]


def test_template_runtime_supports_rnaseq_alignment_step(
    tmp_path: Path,
    monkeypatch,
) -> None:
    runtime = _runtime()
    config = _rnaseq_config(tmp_path)
    config.STAR_INDEX_DIR.mkdir(parents=True, exist_ok=True)
    read1 = Path(config.RAW_DIR) / "sample_a_R1.fastq.gz"
    read2 = Path(config.RAW_DIR) / "sample_a_R2.fastq.gz"
    read1.write_text("r1", encoding="utf-8")
    read2.write_text("r2", encoding="utf-8")

    called: list[str] = []

    def fake_align_star(*args, **kwargs):
        called.append("align_star")
        bam_path = Path(config.PROCESSED_DIR) / "sample_a_Aligned.sortedByCoord.out.bam"
        bam_path.write_text("bam", encoding="utf-8")
        return SimpleNamespace(returncode=0)

    def fake_index_bam(*args, **kwargs):
        called.append("index_bam")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(runtime.rnaseq, "align_star", fake_align_star)
    monkeypatch.setattr(runtime.samtools_wrapper, "index_bam", fake_index_bam)

    assert runtime.run_shared_step(config, "alignment") is True
    assert called == ["align_star", "index_bam"]


def test_template_runtime_supports_rnaseq_quantification_step(
    tmp_path: Path,
    monkeypatch,
) -> None:
    runtime = _runtime()
    config = _rnaseq_config(tmp_path)
    Path(config.GENOME_GTF).write_text("chr1\ttest\tgene\t1\t4\t.\t+\t.\tID=gene1\n", encoding="utf-8")
    bam_path = Path(config.PROCESSED_DIR) / "sample_a_Aligned.sortedByCoord.out.bam"
    bam_path.write_text("bam", encoding="utf-8")

    calls: list[tuple[tuple, dict]] = []

    def fake_run_featurecounts(*args, **kwargs):
        calls.append((args, kwargs))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(runtime.rnaseq, "run_featurecounts", fake_run_featurecounts)

    assert runtime.run_shared_step(config, "quantification") is True
    assert calls
    assert calls[0][1]["feature_type"] == "gene"
    assert calls[0][1]["attribute_type"] == "ID"


def test_yeast_demo_pipeline_prefers_shared_runtime() -> None:
    pipeline_text = YEAST_PIPELINE.read_text(encoding="utf-8")

    assert "load_template_runtime" in pipeline_text
    assert 'run_shared_step("build_index")' in pipeline_text
    assert 'run_shared_step("alignment")' in pipeline_text
    assert 'run_shared_step("quantification")' in pipeline_text
