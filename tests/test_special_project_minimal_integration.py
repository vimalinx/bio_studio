from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _project_paths(project_name: str) -> tuple[Path, Path]:
    project_root = ROOT / "projects" / project_name
    scripts_dir = project_root / "scripts"
    return project_root, scripts_dir


def test_special_projects_are_classified_in_readme() -> None:
    specialist_readme = (
        ROOT / "projects" / "test_rnaseq_analysis" / "README.md"
    ).read_text(encoding="utf-8")
    benchmark_readme = (
        ROOT / "projects" / "cross_system_benchmark" / "README.md"
    ).read_text(encoding="utf-8")
    learning_readme = (
        ROOT / "projects" / "yeast_genome_learning" / "README.md"
    ).read_text(encoding="utf-8")

    assert "专项项目" in specialist_readme
    assert "专项项目" in benchmark_readme
    assert "学习项目" in learning_readme


def test_special_projects_have_validation_entrypoints() -> None:
    for project_name in ("test_rnaseq_analysis", "cross_system_benchmark", "yeast_genome_learning"):
        _, scripts_dir = _project_paths(project_name)
        assert (scripts_dir / "validate_project.py").exists(), project_name
        assert (scripts_dir / "pipeline.py").exists(), project_name


def test_special_projects_pipeline_validate_runs_and_writes_report() -> None:
    for project_name in ("test_rnaseq_analysis", "cross_system_benchmark", "yeast_genome_learning"):
        project_root, scripts_dir = _project_paths(project_name)

        result = subprocess.run(
            [sys.executable, "pipeline.py", "--validate"],
            cwd=scripts_dir,
            capture_output=True,
            text=True,
            check=True,
        )

        assert "validation_report.json" in result.stdout, project_name
        assert (project_root / "logs" / "validation_report.json").exists(), project_name
        assert (project_root / "logs" / "validation_report.md").exists(), project_name


def test_special_projects_pipeline_steps_runs() -> None:
    for project_name in ("test_rnaseq_analysis", "cross_system_benchmark", "yeast_genome_learning"):
        _, scripts_dir = _project_paths(project_name)

        result = subprocess.run(
            [sys.executable, "pipeline.py", "--steps"],
            cwd=scripts_dir,
            capture_output=True,
            text=True,
            check=True,
        )

        assert "step" in result.stdout.lower() or "步骤" in result.stdout, project_name
