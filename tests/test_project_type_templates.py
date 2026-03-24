from __future__ import annotations

import os
from pathlib import Path

from lib.create_project import create_project


def _generate_project(tmp_path: Path, name: str, project_type: str) -> Path:
    original_cwd = Path.cwd()
    try:
        workspace_root = tmp_path.resolve()
        workspace_root.mkdir(parents=True, exist_ok=True)
        os.chdir(workspace_root)
        create_project(name, project_type, f"{project_type} template test")
        return workspace_root / "projects" / name
    finally:
        os.chdir(original_cwd)


def test_rnaseq_template_includes_rnaseq_specific_skeleton(tmp_path: Path) -> None:
    project_dir = _generate_project(tmp_path, "rnaseq_demo", "rnaseq")
    config_text = (project_dir / "scripts" / "config.py").read_text(encoding="utf-8")
    pipeline_text = (project_dir / "scripts" / "pipeline.py").read_text(encoding="utf-8")
    readme_text = (project_dir / "README.md").read_text(encoding="utf-8")

    assert "ANNOTATION_GTF" in config_text
    assert 'ALIGNER = "STAR"' in config_text
    assert "featureCounts" in pipeline_text
    assert "STAR" in pipeline_text
    assert "RNA-seq" in readme_text


def test_variant_template_includes_variant_specific_skeleton(tmp_path: Path) -> None:
    project_dir = _generate_project(tmp_path, "variant_demo", "variant")
    config_text = (project_dir / "scripts" / "config.py").read_text(encoding="utf-8")
    pipeline_text = (project_dir / "scripts" / "pipeline.py").read_text(encoding="utf-8")
    readme_text = (project_dir / "README.md").read_text(encoding="utf-8")

    assert "VARIANT_CALLER" in config_text
    assert 'ALIGNER = "bwa-mem"' in config_text
    assert "bcftools" in pipeline_text
    assert "bwa" in pipeline_text.lower()
    assert "变异检测" in readme_text


def test_phylogeny_template_includes_phylogeny_specific_skeleton(tmp_path: Path) -> None:
    project_dir = _generate_project(tmp_path, "phylo_demo", "phylogeny")
    config_text = (project_dir / "scripts" / "config.py").read_text(encoding="utf-8")
    pipeline_text = (project_dir / "scripts" / "pipeline.py").read_text(encoding="utf-8")
    readme_text = (project_dir / "README.md").read_text(encoding="utf-8")

    assert "ALIGNMENT_INPUT" in config_text
    assert 'TREE_BUILDER = "iqtree2"' in config_text
    assert "mafft" in pipeline_text.lower()
    assert "iqtree2" in pipeline_text.lower()
    assert "系统发育" in readme_text
