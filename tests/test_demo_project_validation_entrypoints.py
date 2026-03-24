from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _project_paths(project_name: str) -> tuple[Path, Path]:
    project_root = ROOT / "projects" / project_name
    scripts_dir = project_root / "scripts"
    return project_root, scripts_dir


def test_demo_projects_have_validation_entrypoints() -> None:
    for project_name in ("ai_design_playground", "yeast_rnaseq_demo"):
        _, scripts_dir = _project_paths(project_name)
        validate_script = scripts_dir / "validate_project.py"
        pipeline_script = scripts_dir / "pipeline.py"

        assert validate_script.exists(), project_name
        assert "--validate" in pipeline_script.read_text(encoding="utf-8"), project_name
        assert "--steps" in pipeline_script.read_text(encoding="utf-8"), project_name


def test_demo_projects_pipeline_validate_runs_and_writes_report() -> None:
    for project_name in ("ai_design_playground", "yeast_rnaseq_demo"):
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


def test_demo_projects_pipeline_steps_runs() -> None:
    for project_name in ("ai_design_playground", "yeast_rnaseq_demo"):
        _, scripts_dir = _project_paths(project_name)

        result = subprocess.run(
            [sys.executable, "pipeline.py", "--steps"],
            cwd=scripts_dir,
            capture_output=True,
            text=True,
            check=True,
        )

        assert "step" in result.stdout.lower() or "步骤" in result.stdout, project_name
