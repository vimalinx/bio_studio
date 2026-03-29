from __future__ import annotations

from pathlib import Path

from lib.service_api.execution import get_run, submit_run
from lib.service_api.models import RunCreateRequest


def test_submit_run_creates_project_shell_and_manifest(tmp_path: Path) -> None:
    request = RunCreateRequest(prompt="analyze ebola benchmarks", project_name="api_ebola_run")

    record, preview = submit_run(request, workspace_root=tmp_path)

    project_dir = tmp_path / "projects" / "api_ebola_run"
    run_dir = project_dir / "api"

    assert record.status == "ready"
    assert preview.brief.task_type == "analysis"
    assert project_dir.exists()
    assert (run_dir / "run_manifest.json").exists()
    assert (run_dir / "plan_preview.json").exists()

    loaded = get_run(record.run_id, workspace_root=tmp_path)
    assert loaded.run_id == record.run_id
    assert loaded.project_name == "api_ebola_run"
    assert any(item.label == "run_manifest.json" for item in loaded.artifacts)


def test_submit_run_executes_selected_cli_capability(tmp_path: Path) -> None:
    request = RunCreateRequest(
        prompt="validate the generated project shell",
        project_name="api_validate_run",
        requested_capabilities=["workspace.project.validate"],
    )

    record, _ = submit_run(request, workspace_root=tmp_path)

    assert record.status == "succeeded"
    assert (
        tmp_path / "projects" / "api_validate_run" / "logs" / "validation_report.json"
    ).exists()
