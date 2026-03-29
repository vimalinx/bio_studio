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


def test_render_command_reference_markdown_contains_help_and_man_sections() -> None:
    module = _load_module(SCRIPT_PATH, "generate_bio_skills")
    text = module.render_command_reference_markdown(
        {
            "name": "samtools",
            "command": "samtools",
            "sources": ["conda_bioconda", "path"],
            "path": "/usr/bin/samtools",
            "summary": "Utilities for SAM/BAM/CRAM files.",
            "version_text": "samtools 1.21",
            "help_text": "samtools view\nsamtools sort",
            "man_text": "SAMTOOLS(1)\nUtilities for alignments",
        }
    )

    assert "## Captured Version" in text
    assert "## Captured Help" in text
    assert "## Captured Man Page" in text
    assert "samtools 1.21" in text
    assert "SAMTOOLS(1)" in text
