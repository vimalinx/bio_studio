from __future__ import annotations

from datetime import datetime, timezone
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field


RunStatus = Literal[
    "queued",
    "planning",
    "ready",
    "running",
    "succeeded",
    "failed",
    "cancelled",
]


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


class TaskBrief(BaseModel):
    model_config = ConfigDict(extra="forbid")

    prompt: str = ""
    task_type: str = "analysis"
    intent_summary: str | None = None
    input_mode: str = "prompt"
    molecule_type: str = "mixed"
    candidate_hosts: list[str] = Field(default_factory=list)
    output_granularity: str = "report"
    network_policy: str = "allow_online"
    success_criteria: list[str] = Field(
        default_factory=lambda: ["reproducible_artifacts"]
    )
    risk_tolerance: str = "medium"
    requested_capabilities: list[str] = Field(default_factory=list)
    constraints: list[str] = Field(default_factory=list)


class CapabilityRecord(BaseModel):
    model_config = ConfigDict(extra="forbid")

    capability_id: str
    display_name: str
    task_types: list[str]
    keywords: list[str] = Field(default_factory=list)
    transport: str
    backend_target: str
    network_requirements: str = "optional"
    resource_class: str = "standard"
    supports_preview: bool = True
    supports_async: bool = False
    stability_tier: str = "stable"


class PlanStep(BaseModel):
    model_config = ConfigDict(extra="forbid")

    step_id: str
    kind: Literal["project", "capability", "artifact"]
    label: str
    capability_id: str | None = None
    backend_target: str | None = None
    rationale: str | None = None


class PlanPreview(BaseModel):
    model_config = ConfigDict(extra="forbid")

    brief: TaskBrief
    selected_capabilities: list[CapabilityRecord] = Field(default_factory=list)
    steps: list[PlanStep] = Field(default_factory=list)
    notes: list[str] = Field(default_factory=list)


class RunCreateRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    prompt: str = ""
    structured_brief: TaskBrief | None = None
    project_name: str | None = None
    execution_mode: Literal["plan", "sync"] = "sync"
    create_project: bool = True
    network_policy: str | None = None
    requested_capabilities: list[str] = Field(default_factory=list)


class RunArtifact(BaseModel):
    model_config = ConfigDict(extra="forbid")

    kind: str
    path: str
    label: str


class RunRecord(BaseModel):
    model_config = ConfigDict(extra="forbid")

    run_id: str
    status: RunStatus
    project_name: str
    brief: TaskBrief
    selected_capability_ids: list[str] = Field(default_factory=list)
    notes: list[str] = Field(default_factory=list)
    artifact_dir: str | None = None
    artifacts: list[RunArtifact] = Field(default_factory=list)
    created_at: str = Field(default_factory=utc_now_iso)
    updated_at: str = Field(default_factory=utc_now_iso)
