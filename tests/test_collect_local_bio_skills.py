from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "skills" / "collect_local_bio_skills.py"


def _load_module(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_bio_skill_filter_matches_real_bio_descriptions() -> None:
    module = _load_module(SCRIPT_PATH, "collect_local_bio_skills")

    keep = module.classify_skill(
        "scanpy",
        "Single-cell RNA-seq analysis with QC, PCA, UMAP, Leiden clustering, and marker genes.",
    )
    drop = module.classify_skill(
        "code-reviewer",
        "Comprehensive code review skill for TypeScript, JavaScript, Python, Swift, Kotlin, Go.",
    )

    assert keep["is_bio"] is True
    assert keep["category"] == "single-cell-omics"
    assert drop["is_bio"] is False


def test_inventory_doc_writer_emits_scope_sections(tmp_path: Path) -> None:
    module = _load_module(SCRIPT_PATH, "collect_local_bio_skills")
    output_path = tmp_path / "local-bio-skill-inventory.md"

    module.write_inventory_doc(
        output_path,
        [
            {
                "name": "samtools",
                "scope": "project-local",
                "category": "core-bioinformatics",
                "path": "/tmp/project/.claude/skills/samtools/SKILL.md",
                "description": "Utilities for SAM/BAM/CRAM alignment files.",
            },
            {
                "name": "scanpy",
                "scope": "global-codex",
                "category": "single-cell-omics",
                "path": "/tmp/home/.codex/skills/scanpy/SKILL.md",
                "description": "Single-cell RNA-seq analysis.",
            },
        ],
    )

    text = output_path.read_text(encoding="utf-8")
    assert "project-local" in text
    assert "global-codex" in text
    assert "samtools" in text
    assert "scanpy" in text
