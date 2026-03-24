from __future__ import annotations

from pathlib import Path

from lib.create_project import (
    build_config_content,
    build_pipeline_content,
    build_readme_content,
)


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_ROOT = ROOT / "projects" / "example_rnaseq"


def test_example_rnaseq_matches_generated_rnaseq_template() -> None:
    assert (EXAMPLE_ROOT / "README.md").read_text(encoding="utf-8") == build_readme_content(
        "example_rnaseq",
        "rnaseq",
        "示例RNA-seq分析项目",
    )
    assert (
        EXAMPLE_ROOT / "scripts" / "config.py"
    ).read_text(encoding="utf-8") == build_config_content(
        "example_rnaseq",
        "rnaseq",
        "示例RNA-seq分析项目",
    )
    assert (
        EXAMPLE_ROOT / "scripts" / "pipeline.py"
    ).read_text(encoding="utf-8") == build_pipeline_content("example_rnaseq", "rnaseq")
