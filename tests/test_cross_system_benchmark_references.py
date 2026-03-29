from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REFERENCE_SCRIPT = ROOT / "projects" / "cross_system_benchmark" / "scripts" / "prepare_references.py"


def _load_module(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_write_reference_manifest_records_small_default_reference_policy(tmp_path: Path) -> None:
    module = _load_module(REFERENCE_SCRIPT, "cross_system_prepare_refs")
    manifest_path = module.write_reference_manifest(tmp_path)

    assert manifest_path.exists()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert sorted(manifest["references"]) == ["ebola", "sars_cov_2", "yeast"]
    assert manifest["references"]["ebola"]["download_by_default"] is True
    assert manifest["references"]["sars_cov_2"]["download_by_default"] is True
    assert manifest["references"]["yeast"]["download_by_default"] is False
    assert manifest["host_references"]["human"]["download_by_default"] is False
