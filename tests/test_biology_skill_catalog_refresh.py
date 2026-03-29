from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "skills" / "generate_bio_skills.py"


def _load_module(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_catalog_refresh_lists_existing_and_generated_skills(tmp_path: Path) -> None:
    module = _load_module(SCRIPT_PATH, "generate_bio_skills")
    catalog_path = tmp_path / "biology-skill-catalog.md"
    module.write_catalog_doc(
        catalog_path,
        existing_domain_skills=["bioinformatics-toolkit"],
        generated_tools=[{"name": "samtools", "sources": ["path"]}],
    )

    text = catalog_path.read_text(encoding="utf-8")
    assert "bioinformatics-toolkit" in text
    assert "samtools" in text
    assert "generated-tool" in text
