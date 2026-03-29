from __future__ import annotations

import re
import uuid
from pathlib import Path

from .adapters.cli import WorkspaceCLIAdapter
from .models import PlanPreview, RunCreateRequest, RunRecord, utc_now_iso
from .planner import plan_request
from .registry import CapabilityRegistry, build_default_registry
from .storage import FileRunStore


REPO_ROOT = Path(__file__).resolve().parents[2]


def default_workspace_root() -> Path:
    return Path(REPO_ROOT).resolve()


def submit_run(
    request: RunCreateRequest,
    workspace_root: Path | None = None,
    registry: CapabilityRegistry | None = None,
) -> tuple[RunRecord, PlanPreview]:
    workspace_root = Path(workspace_root or default_workspace_root()).resolve()
    registry = registry or build_default_registry()

    preview = plan_request(request, registry)
    run_id = _generate_run_id()
    project_name = request.project_name or _build_project_name(preview.brief.prompt, run_id)
    store = FileRunStore(workspace_root)

    record = RunRecord(
        run_id=run_id,
        status="planning",
        project_name=project_name,
        brief=preview.brief,
        selected_capability_ids=[item.capability_id for item in preview.selected_capabilities],
        notes=list(preview.notes),
    )
    record = store.write_run(record)
    store.write_plan_preview(run_id, preview)
    store.append_event(run_id, {"stage": "planning", "status": "ok"})

    if request.execution_mode == "plan":
        record = record.model_copy(update={"status": "ready", "updated_at": utc_now_iso()})
        record = store.write_run(record)
        return _finalize_record(store, record), preview

    if request.create_project:
        adapter = WorkspaceCLIAdapter(
            workspace_root,
            cli_path=default_workspace_root() / "scripts" / "project.py",
        )
        result = adapter.invoke(
            [
                "create",
                project_name,
                "--type",
                _project_type_for_brief(preview.brief),
                "--description",
                preview.brief.intent_summary or preview.brief.prompt or "API-submitted run",
            ]
        )
        store.append_event(
            run_id,
            {
                "stage": "project_create",
                "status": "ok" if result.returncode == 0 else "error",
                "command": result.command,
                "returncode": result.returncode,
                "stdout": result.stdout,
                "stderr": result.stderr,
            },
        )
        if result.returncode != 0:
            failed = record.model_copy(
                update={
                    "status": "failed",
                    "notes": [*record.notes, "Workspace project creation failed."],
                    "updated_at": utc_now_iso(),
                }
            )
            failed = store.write_run(failed)
            return _finalize_record(store, failed), preview

    execution_result = _execute_supported_cli_capabilities(
        store=store,
        run_id=run_id,
        project_name=project_name,
        preview=preview,
        workspace_root=workspace_root,
    )
    if execution_result == "failed":
        failed = record.model_copy(
            update={
                "status": "failed",
                "notes": [*record.notes, "A selected CLI capability failed during execution."],
                "updated_at": utc_now_iso(),
            }
        )
        failed = store.write_run(failed)
        return _finalize_record(store, failed), preview

    terminal_status = "succeeded" if execution_result == "executed" else "ready"
    ready = record.model_copy(update={"status": terminal_status, "updated_at": utc_now_iso()})
    ready = store.write_run(ready)
    store.write_plan_preview(run_id, preview)
    store.append_event(run_id, {"stage": terminal_status, "status": "ok"})
    return _finalize_record(store, ready), preview


def get_run(run_id: str, workspace_root: Path | None = None) -> RunRecord:
    store = FileRunStore(Path(workspace_root or default_workspace_root()).resolve())
    return _finalize_record(store, store.read_run(run_id))


def _finalize_record(store: FileRunStore, record: RunRecord) -> RunRecord:
    return record.model_copy(update={"artifacts": store.list_artifacts(record.run_id)})


def _generate_run_id() -> str:
    return f"run_{uuid.uuid4().hex[:12]}"


def _build_project_name(prompt: str, run_id: str) -> str:
    words = re.sub(r"[^a-zA-Z0-9]+", "-", prompt.lower()).strip("-")
    prefix = "-".join([part for part in words.split("-") if part][:4]) or "api-run"
    return f"{prefix}-{run_id[-6:]}"


def _project_type_for_brief(brief) -> str:
    prompt = brief.prompt.lower()
    if "rna-seq" in prompt or "rnaseq" in prompt:
        return "rnaseq"
    return "generic"


def _execute_supported_cli_capabilities(
    *,
    store: FileRunStore,
    run_id: str,
    project_name: str,
    preview: PlanPreview,
    workspace_root: Path,
) -> str:
    command_map = {
        "workspace.project.validate": ["validate", project_name],
        "workspace.project.run": ["run", project_name],
        "workspace.validation.smoke": ["workspace-validate"],
    }

    executable = [
        capability
        for capability in preview.selected_capabilities
        if capability.transport == "cli"
        and capability.capability_id in command_map
        and (
            capability.capability_id in preview.brief.requested_capabilities
            or preview.brief.task_type == "workspace_validation"
        )
    ]

    if not executable:
        for capability in preview.selected_capabilities:
            if capability.transport != "cli":
                store.append_event(
                    run_id,
                    {
                        "stage": "capability_skipped",
                        "status": "skipped",
                        "capability_id": capability.capability_id,
                        "reason": "Transport is not wired for direct execution in phase 1.",
                    },
                )
        return "none"

    adapter = WorkspaceCLIAdapter(
        workspace_root,
        cli_path=default_workspace_root() / "scripts" / "project.py",
    )
    store.append_event(run_id, {"stage": "execution", "status": "running"})

    for capability in executable:
        result = adapter.invoke(command_map[capability.capability_id])
        store.append_event(
            run_id,
            {
                "stage": "capability_execute",
                "status": "ok" if result.returncode == 0 else "error",
                "capability_id": capability.capability_id,
                "command": result.command,
                "returncode": result.returncode,
                "stdout": result.stdout,
                "stderr": result.stderr,
            },
        )
        if result.returncode != 0:
            return "failed"

    for capability in preview.selected_capabilities:
        if capability.transport != "cli":
            store.append_event(
                run_id,
                {
                    "stage": "capability_skipped",
                    "status": "skipped",
                    "capability_id": capability.capability_id,
                    "reason": "Transport is not wired for direct execution in phase 1.",
                },
            )

    return "executed"
