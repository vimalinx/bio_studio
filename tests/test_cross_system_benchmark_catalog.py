from __future__ import annotations

import importlib.util
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
FETCH_SCRIPT = ROOT / "projects" / "cross_system_benchmark" / "scripts" / "fetch_benchmarks.py"


def _load_module(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_write_catalog_materializes_expected_files(tmp_path: Path) -> None:
    module = _load_module(FETCH_SCRIPT, "cross_system_fetch")
    written = module.write_catalog(tmp_path)

    datasets_path = written["datasets"]
    samples_path = written["samples"]
    comparisons_path = written["comparisons"]
    source_links_path = written["source_links"]

    assert datasets_path.exists()
    assert samples_path.exists()
    assert comparisons_path.exists()
    assert source_links_path.exists()

    payload = yaml.safe_load(datasets_path.read_text(encoding="utf-8"))
    assert sorted(payload) == ["ebola", "sars_cov_2", "yeast"]
    assert payload["ebola"]["geo_accession"] == "GSE114905"
    assert payload["sars_cov_2"]["geo_accession"] == "GSE147507"
    assert payload["yeast"]["geo_accession"] == "GSE67149"
