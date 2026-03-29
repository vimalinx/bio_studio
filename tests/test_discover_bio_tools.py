from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "skills" / "discover_bio_tools.py"


def _load_module(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_discovery_returns_tools_from_known_sources() -> None:
    module = _load_module(SCRIPT_PATH, "discover_bio_tools")
    discovered = module.discover_tools(ROOT)
    names = {item["name"] for item in discovered}

    assert "samtools" in names
    assert "fastp" in names
    assert "mafft" in names


def test_discovery_marks_existing_domain_coverage() -> None:
    module = _load_module(SCRIPT_PATH, "discover_bio_tools")
    discovered = module.discover_tools(ROOT)
    by_name = {item["name"]: item for item in discovered}

    assert by_name["samtools"]["sources"]
    assert "tool_reference" in by_name["samtools"]["sources"] or "path" in by_name["samtools"]["sources"]
