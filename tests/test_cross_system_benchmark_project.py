from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = ROOT / "projects" / "cross_system_benchmark"
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
CLI = ROOT / "scripts" / "project.py"


def test_cross_system_benchmark_project_exists() -> None:
    assert PROJECT_ROOT.exists()
    assert (PROJECT_ROOT / "README.md").exists()
    assert (SCRIPTS_DIR / "pipeline.py").exists()
    assert (SCRIPTS_DIR / "validate_project.py").exists()


def test_cross_system_benchmark_steps_list_fetch_prepare_baseline_report() -> None:
    result = subprocess.run(
        [sys.executable, "pipeline.py", "--steps"],
        cwd=SCRIPTS_DIR,
        capture_output=True,
        text=True,
        check=True,
    )

    assert "fetch" in result.stdout
    assert "prepare" in result.stdout
    assert "baseline" in result.stdout
    assert "report" in result.stdout


def test_workspace_cli_can_list_cross_system_benchmark_steps() -> None:
    result = subprocess.run(
        [sys.executable, str(CLI), "steps", "cross_system_benchmark"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )

    assert "fetch" in result.stdout
    assert "report" in result.stdout
