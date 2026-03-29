from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_biology_skill_docs_exist() -> None:
    spec_path = ROOT / "docs" / "skills" / "biology-skill-spec.md"
    catalog_path = ROOT / "docs" / "skills" / "biology-skill-catalog.md"

    assert spec_path.exists()
    assert catalog_path.exists()


def test_biology_skill_spec_mentions_required_frontmatter_and_layout() -> None:
    text = (ROOT / "docs" / "skills" / "biology-skill-spec.md").read_text(encoding="utf-8")

    assert "SKILL.md" in text
    assert "disable-model-invocation: true" in text
    assert ".claude/skills/<tool-name>/" in text


def test_biology_skill_catalog_mentions_existing_domain_skills() -> None:
    text = (ROOT / "docs" / "skills" / "biology-skill-catalog.md").read_text(encoding="utf-8")

    for skill_name in [
        "bioinformatics-toolkit",
        "biomni",
        "evo2",
        "protein-structure",
        "rnaseq-pipeline",
    ]:
        assert skill_name in text
