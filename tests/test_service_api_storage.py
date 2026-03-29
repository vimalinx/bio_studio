from __future__ import annotations

from pathlib import Path

from lib.service_api.models import PlanPreview, RunRecord, TaskBrief
from lib.service_api.storage import FileRunStore


def test_file_run_store_writes_manifest_and_events(tmp_path: Path) -> None:
    store = FileRunStore(tmp_path)
    record = RunRecord(
        run_id="run_001",
        status="queued",
        project_name="api_run_001",
        brief=TaskBrief(prompt="analyze yeast expression strategy"),
    )

    record = store.write_run(record)
    store.write_plan_preview("run_001", PlanPreview(brief=record.brief))
    store.append_event("run_001", {"stage": "planning", "status": "ok"})

    run_dir = tmp_path / ".codex-runs" / "service_api" / "runs" / "run_001"
    assert (run_dir / "run_manifest.json").exists()
    assert (run_dir / "plan_preview.json").exists()
    assert (run_dir / "events.jsonl").exists()

    artifacts = store.list_artifacts("run_001")
    labels = {item.label for item in artifacts}
    events = store.read_events("run_001")

    assert labels == {"events.jsonl", "plan_preview.json", "run_manifest.json"}
    assert events[0]["stage"] == "planning"
