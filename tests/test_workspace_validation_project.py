from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT / "projects" / "test_env_validation" / "scripts"


def test_workspace_validation_project_has_project_validate_entrypoint() -> None:
    validate_script = SCRIPTS_DIR / "validate_project.py"
    pipeline_script = SCRIPTS_DIR / "pipeline.py"

    assert validate_script.exists()
    assert "--validate" in pipeline_script.read_text(encoding="utf-8")


def test_workspace_validation_project_validate_runs_and_writes_report() -> None:
    result = subprocess.run(
        [sys.executable, "validate_project.py"],
        cwd=SCRIPTS_DIR,
        capture_output=True,
        text=True,
        check=True,
    )

    assert "PASS" in result.stdout
    assert (
        ROOT / "projects" / "test_env_validation" / "logs" / "validation_report.json"
    ).exists()
    assert (
        ROOT / "projects" / "test_env_validation" / "logs" / "validation_report.md"
    ).exists()


def test_workspace_validation_project_pipeline_validate_runs() -> None:
    result = subprocess.run(
        [sys.executable, "pipeline.py", "--validate"],
        cwd=SCRIPTS_DIR,
        capture_output=True,
        text=True,
        check=True,
    )

    assert "validation_report.json" in result.stdout
