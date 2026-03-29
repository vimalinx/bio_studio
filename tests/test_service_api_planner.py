from __future__ import annotations

from lib.service_api.models import RunCreateRequest, TaskBrief
from lib.service_api.planner import normalize_request, plan_request


def test_plan_request_normalizes_prompt_into_brief_and_plan() -> None:
    request = RunCreateRequest(prompt="analyze SARS-CoV-2 literature and sequence evidence")
    preview = plan_request(request)

    capability_ids = [item.capability_id for item in preview.selected_capabilities]

    assert preview.brief.task_type == "analysis"
    assert "database.search.pubmed" in capability_ids
    assert "sequence.analysis.basic" in capability_ids
    assert any(step.kind == "capability" for step in preview.steps)


def test_normalize_request_respects_structured_brief_override() -> None:
    brief = TaskBrief(prompt="custom", task_type="design", output_granularity="candidate_table")
    request = RunCreateRequest(prompt="ignored", structured_brief=brief, network_policy="offline")

    normalized = normalize_request(request)

    assert normalized.task_type == "design"
    assert normalized.output_granularity == "candidate_table"
    assert normalized.network_policy == "offline"
