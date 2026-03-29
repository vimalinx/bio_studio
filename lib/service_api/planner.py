from __future__ import annotations

from .models import PlanPreview, PlanStep, RunCreateRequest, TaskBrief
from .registry import CapabilityRegistry, build_default_registry


def normalize_request(request: RunCreateRequest) -> TaskBrief:
    if request.structured_brief is not None:
        brief = request.structured_brief.model_copy(deep=True)
    else:
        prompt = request.prompt.strip()
        task_type = _infer_task_type(prompt)
        brief = TaskBrief(
            prompt=prompt,
            task_type=task_type,
            intent_summary=_summarize_prompt(prompt),
            molecule_type=_infer_molecule_type(prompt),
            output_granularity=_infer_output_granularity(prompt, task_type),
            requested_capabilities=list(request.requested_capabilities),
        )

    if request.network_policy:
        brief = brief.model_copy(update={"network_policy": request.network_policy})

    if request.requested_capabilities:
        merged = list(dict.fromkeys([*brief.requested_capabilities, *request.requested_capabilities]))
        brief = brief.model_copy(update={"requested_capabilities": merged})

    return brief


def plan_request(
    request: RunCreateRequest,
    registry: CapabilityRegistry | None = None,
) -> PlanPreview:
    registry = registry or build_default_registry()
    brief = normalize_request(request)

    scored = []
    for capability in registry.list_capabilities():
        score = _score_capability(capability, brief)
        if score > 0:
            scored.append((score, capability))

    scored.sort(key=lambda item: (-item[0], item[1].capability_id))
    selected_capabilities = [capability for _, capability in scored[:5]]

    if not selected_capabilities:
        fallback = registry.get("workspace.project.create")
        if fallback is not None:
            selected_capabilities = [fallback]

    steps: list[PlanStep] = []
    if request.create_project:
        steps.append(
            PlanStep(
                step_id="project-shell",
                kind="project",
                label="Create or select a workspace project shell",
                capability_id="workspace.project.create",
                backend_target="scripts/project.py create",
                rationale="Every remote run should land in a durable workspace boundary.",
            )
        )

    for index, capability in enumerate(selected_capabilities, start=1):
        steps.append(
            PlanStep(
                step_id=f"capability-{index}",
                kind="capability",
                label=capability.display_name,
                capability_id=capability.capability_id,
                backend_target=capability.backend_target,
                rationale=f"Matched to task type '{brief.task_type}'.",
            )
        )

    notes = [
        "Planner is rules-first in phase 1.",
        "Selected capabilities are ranked for preview before heavy execution is added.",
    ]
    if request.execution_mode == "plan":
        notes.append("Request is preview-only and will not execute project actions.")

    return PlanPreview(
        brief=brief,
        selected_capabilities=selected_capabilities,
        steps=steps,
        notes=notes,
    )


def _summarize_prompt(prompt: str) -> str | None:
    cleaned = " ".join(prompt.strip().split())
    if not cleaned:
        return None
    return cleaned[:140]


def _infer_task_type(prompt: str) -> str:
    lowered = prompt.lower()
    if any(token in lowered for token in ("workspace validate", "workspace-validate", "smoke test", "environment validate")):
        return "workspace_validation"
    if "create project" in lowered or "new project" in lowered:
        return "project_setup"
    if any(token in lowered for token in ("design", "optimize", "construct", "mrna")):
        return "design"
    if (
        any(token in lowered for token in ("sequence", "fasta", "blast", "genome", "alignment"))
        and not any(token in lowered for token in ("analyze", "analysis", "literature", "evidence"))
    ):
        return "sequence"
    return "analysis"


def _infer_molecule_type(prompt: str) -> str:
    lowered = prompt.lower()
    if "protein" in lowered or "peptide" in lowered:
        return "protein"
    if "rna" in lowered or "mrna" in lowered:
        return "rna"
    if "dna" in lowered:
        return "dna"
    return "mixed"


def _infer_output_granularity(prompt: str, task_type: str) -> str:
    lowered = prompt.lower()
    if task_type == "design":
        return "candidate_table"
    if any(token in lowered for token in ("sequence", "fasta", "alignment")):
        return "sequence_report"
    return "report"


def _score_capability(capability, brief: TaskBrief) -> int:
    score = 0

    if brief.task_type in capability.task_types:
        score += 8
    if "analysis" in capability.task_types and brief.task_type in {"analysis", "sequence"}:
        score += 2

    lowered_prompt = brief.prompt.lower()
    for keyword in capability.keywords:
        if keyword.lower() in lowered_prompt:
            score += 3

    if capability.capability_id in brief.requested_capabilities:
        score += 10

    if brief.task_type == "analysis" and capability.capability_id == "database.search.pubmed":
        score += 3
    if brief.task_type in {"analysis", "sequence"} and capability.capability_id == "sequence.analysis.basic":
        score += 3
    if brief.task_type == "workspace_validation" and capability.capability_id == "workspace.validation.smoke":
        score += 6

    return score
