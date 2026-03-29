from __future__ import annotations

from lib.service_api.models import RunRecord, TaskBrief


def test_task_brief_defaults_are_stable() -> None:
    brief = TaskBrief(prompt="analyze this RNA-seq dataset")

    assert brief.task_type == "analysis"
    assert brief.network_policy == "allow_online"
    assert brief.output_granularity == "report"


def test_run_record_round_trips_to_dict() -> None:
    record = RunRecord(
        run_id="run_123",
        status="queued",
        project_name="api_run_123",
        brief=TaskBrief(prompt="analyze ebola sequence diversity"),
    )

    payload = record.model_dump()
    restored = RunRecord.model_validate(payload)

    assert restored.run_id == "run_123"
    assert restored.status == "queued"
