from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from lib.create_project import create_project


def _generate_project(tmp_path: Path, name: str = "demo_project") -> Path:
    original_cwd = Path.cwd()
    try:
        project_root = tmp_path.resolve()
        project_root.mkdir(parents=True, exist_ok=True)
        os.chdir(project_root)
        create_project(name, "rnaseq", "test project")
        return project_root / "projects" / name
    finally:
        os.chdir(original_cwd)


def test_create_project_generates_validation_entrypoints(tmp_path: Path) -> None:
    project_dir = _generate_project(tmp_path)

    validate_script = project_dir / "scripts" / "validate_project.py"
    pipeline_script = project_dir / "scripts" / "pipeline.py"
    readme_file = project_dir / "README.md"

    assert validate_script.exists()
    assert "--validate" in pipeline_script.read_text(encoding="utf-8")
    assert "validate_project.py" in readme_file.read_text(encoding="utf-8")


def test_generated_project_validation_runs_and_writes_report(tmp_path: Path) -> None:
    project_dir = _generate_project(tmp_path, "validated_project")
    scripts_dir = project_dir / "scripts"

    result = subprocess.run(
        [sys.executable, "validate_project.py"],
        cwd=scripts_dir,
        capture_output=True,
        text=True,
        check=True,
    )

    assert "PASS" in result.stdout
    assert (project_dir / "logs" / "validation_report.json").exists()
    assert (project_dir / "logs" / "validation_report.md").exists()


def test_generated_pipeline_validate_flag_runs_validation(tmp_path: Path) -> None:
    project_dir = _generate_project(tmp_path, "pipeline_project")
    scripts_dir = project_dir / "scripts"

    result = subprocess.run(
        [sys.executable, "pipeline.py", "--validate"],
        cwd=scripts_dir,
        capture_output=True,
        text=True,
        check=True,
    )

    assert "validation_report.json" in result.stdout
