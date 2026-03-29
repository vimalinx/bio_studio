from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "scripts" / "project.py"


def test_project_cli_create_generates_project(tmp_path: Path) -> None:
    result = subprocess.run(
        [
            sys.executable,
            str(CLI),
            "create",
            "cli_demo",
            "--type",
            "rnaseq",
            "--description",
            "cli generated project",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=True,
    )

    project_dir = tmp_path / "projects" / "cli_demo"
    assert project_dir.exists()
    assert (project_dir / "scripts" / "validate_project.py").exists()
    assert "cli_demo" in result.stdout


def test_project_cli_validate_runs_project_validation(tmp_path: Path) -> None:
    subprocess.run(
        [
            sys.executable,
            str(CLI),
            "create",
            "cli_validate_demo",
            "--type",
            "rnaseq",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=True,
    )

    result = subprocess.run(
        [sys.executable, str(CLI), "validate", "cli_validate_demo"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=True,
    )

    assert "validation_report.json" in result.stdout
    assert (
        tmp_path / "projects" / "cli_validate_demo" / "logs" / "validation_report.json"
    ).exists()


def test_project_cli_steps_lists_generated_pipeline_steps(tmp_path: Path) -> None:
    subprocess.run(
        [sys.executable, str(CLI), "create", "cli_steps_demo", "--type", "rnaseq"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=True,
    )

    result = subprocess.run(
        [sys.executable, str(CLI), "steps", "cli_steps_demo"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=True,
    )

    assert "data_preparation" in result.stdout
    assert "main_analysis" in result.stdout


def test_project_cli_workspace_validate_targets_workspace_smoke_test() -> None:
    result = subprocess.run(
        [sys.executable, str(CLI), "workspace-validate", "--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )

    assert "test_env_validation" in result.stdout


def test_project_cli_validate_runs_special_project_validation() -> None:
    for project_name in ("test_rnaseq_analysis", "cross_system_benchmark"):
        result = subprocess.run(
            [sys.executable, str(CLI), "validate", project_name],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=True,
        )

        assert "Validation summary: PASS" in result.stdout
        assert "validation_report.json" in result.stdout


def test_project_cli_steps_lists_special_project_steps() -> None:
    expected_tokens = {
        "test_rnaseq_analysis": ("data_preparation", "quality_control"),
        "cross_system_benchmark": ("fetch", "report"),
    }
    for project_name, tokens in expected_tokens.items():
        result = subprocess.run(
            [sys.executable, str(CLI), "steps", project_name],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=True,
        )

        for token in tokens:
            assert token in result.stdout


def test_project_cli_run_learning_project_prints_workspace_level_guidance() -> None:
    result = subprocess.run(
        [sys.executable, str(CLI), "run", "yeast_genome_learning"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )

    assert f"python {CLI} steps yeast_genome_learning" in result.stdout
    assert f"python {CLI} validate yeast_genome_learning" in result.stdout
