# Online API Control Plane Phase 1 Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a thin online API control plane that can accept a request, normalize it into a task brief, preview the execution plan, and submit a run that reuses the existing Bio Studio workspace CLI and MCP-adjacent execution surfaces.

**Architecture:** Build a focused `lib/service_api/` package for request models, capability registration, planning, execution, and run storage. Keep the first service layer intentionally thin: FastAPI handles transport, file-backed state handles persistence, and the executor reuses `scripts/project.py` rather than bypassing the current workspace spine. The API should expose planning and run submission first, while long-running orchestration, tenant isolation, and distributed workers stay out of phase 1.

**Tech Stack:** Python 3.10+, FastAPI, Uvicorn, Pydantic, stdlib `json`/`pathlib`/`subprocess`, existing `scripts/project.py`, existing `lib/` helpers, `pytest`

---

## File Structure

### New files

- `lib/service_api/__init__.py`
  - Package export surface for the online control plane.
- `lib/service_api/models.py`
  - Pydantic models for task briefs, capability metadata, plan steps, run requests, and run records.
- `lib/service_api/registry.py`
  - Built-in capability registry and lookup helpers.
- `lib/service_api/planner.py`
  - Rules-first request normalization and capability matching.
- `lib/service_api/execution.py`
  - Run creation, project naming, workspace CLI invocation, and lifecycle transitions.
- `lib/service_api/storage.py`
  - File-backed run manifests, plan snapshots, event logs, and artifact references.
- `lib/service_api/adapters/__init__.py`
  - Adapter package export surface.
- `lib/service_api/adapters/base.py`
  - Shared adapter protocol and result model.
- `lib/service_api/adapters/cli.py`
  - Workspace CLI adapter for `scripts/project.py`.
- `scripts/api_server.py`
  - FastAPI entrypoint and route wiring.
- `tests/test_service_api_models.py`
  - Contract tests for request and run schemas.
- `tests/test_service_api_registry.py`
  - Tests for the built-in capability registry.
- `tests/test_service_api_planner.py`
  - Tests for task normalization and capability matching.
- `tests/test_service_api_storage.py`
  - Tests for file-backed run persistence.
- `tests/test_service_api_execution.py`
  - Tests for run creation and workspace CLI execution orchestration.

### Modified files

- `requirements.txt`
  - Add service dependencies if the project chooses to keep them in the main requirements file.
- `docs/README.md`
  - Add the new spec and plan to the documentation index.
- `README.md`
  - Add a short section documenting the online API control plane entrypoint once phase 1 lands.
- `mcp-servers/bio-lab-mcp/lab_server.py`
  - Optionally add a thin MCP wrapper for plan preview after the HTTP control plane is stable.
- `tests/test_lab_mcp_server.py`
  - Add coverage only if the MCP wrapper is added in this phase.

## Delivery Notes

- Phase 1 should prefer determinism over LLM cleverness.
- The planner must support preview mode before execution mode.
- The executor should call the workspace CLI instead of duplicating project lifecycle logic.
- File-backed storage is acceptable in phase 1 because the repo already treats `projects/` as a durable artifact boundary.
- Auth, quotas, and multi-tenant isolation are explicitly deferred.

## Chunk 1: Models And Capability Registry

### Task 1: Add request, plan, and run models

**Files:**
- Create: `lib/service_api/__init__.py`
- Create: `lib/service_api/models.py`
- Test: `tests/test_service_api_models.py`

- [ ] **Step 1: Write the failing schema tests**

```python
from lib.service_api.models import (
    TaskBrief,
    CapabilityRecord,
    PlanPreview,
    RunCreateRequest,
    RunRecord,
)


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
```

- [ ] **Step 2: Run the focused tests to confirm failure**

Run: `python -m pytest tests/test_service_api_models.py -q`
Expected: FAIL because `lib.service_api` does not exist yet.

- [ ] **Step 3: Implement the minimal models**

Define:

- `TaskBrief`
- `CapabilityRecord`
- `PlanStep`
- `PlanPreview`
- `RunCreateRequest`
- `RunRecord`

Required run statuses:

- `queued`
- `planning`
- `ready`
- `running`
- `succeeded`
- `failed`
- `cancelled`

- [ ] **Step 4: Run the schema tests**

Run: `python -m pytest tests/test_service_api_models.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add lib/service_api/__init__.py lib/service_api/models.py tests/test_service_api_models.py
git commit -m "feat: add online api control plane models"
```

### Task 2: Add the built-in capability registry

**Files:**
- Create: `lib/service_api/registry.py`
- Test: `tests/test_service_api_registry.py`

- [ ] **Step 1: Write the failing registry tests**

```python
from lib.service_api.registry import build_default_registry


def test_default_registry_exposes_workspace_and_bio_capabilities() -> None:
    registry = build_default_registry()
    capability_ids = {item.capability_id for item in registry.list_capabilities()}

    assert "workspace.project.create" in capability_ids
    assert "workspace.project.validate" in capability_ids
    assert "workspace.project.run" in capability_ids
    assert "database.search.pubmed" in capability_ids
    assert "sequence.analysis.basic" in capability_ids
```

- [ ] **Step 2: Run the focused tests**

Run: `python -m pytest tests/test_service_api_registry.py -q`
Expected: FAIL because no registry exists yet.

- [ ] **Step 3: Implement the registry**

Create a registry class with:

- `list_capabilities()`
- `get(capability_id)`
- `find_by_task_type(task_type)`
- `find_by_keyword(keyword)`

Seed the default registry with entries that map to the current workspace surfaces:

- workspace CLI lifecycle commands
- sequence analysis
- structure analysis
- database querying
- design helpers

- [ ] **Step 4: Run the registry tests**

Run: `python -m pytest tests/test_service_api_registry.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add lib/service_api/registry.py tests/test_service_api_registry.py
git commit -m "feat: add built-in capability registry"
```

## Chunk 2: Planner And Storage

### Task 3: Add request normalization and plan preview

**Files:**
- Create: `lib/service_api/planner.py`
- Test: `tests/test_service_api_planner.py`

- [ ] **Step 1: Write the failing planner tests**

```python
from lib.service_api.models import RunCreateRequest
from lib.service_api.planner import plan_request


def test_plan_request_normalizes_prompt_into_brief_and_plan() -> None:
    request = RunCreateRequest(prompt="analyze SARS-CoV-2 literature and sequence evidence")
    preview = plan_request(request)

    assert preview.brief.task_type == "analysis"
    assert preview.selected_capabilities
    assert any(step.kind == "capability" for step in preview.steps)
```

- [ ] **Step 2: Run the planner tests**

Run: `python -m pytest tests/test_service_api_planner.py -q`
Expected: FAIL because the planner does not exist.

- [ ] **Step 3: Implement rules-first planning**

The planner should:

- infer `task_type`
- infer `output_granularity`
- build a normalized `TaskBrief`
- score capabilities from the registry
- emit `PlanPreview`

Keep heuristics explicit and inspectable. Do not call an LLM in phase 1.

- [ ] **Step 4: Run the planner tests**

Run: `python -m pytest tests/test_service_api_planner.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add lib/service_api/planner.py tests/test_service_api_planner.py
git commit -m "feat: add request planner and plan preview"
```

### Task 4: Add file-backed run storage

**Files:**
- Create: `lib/service_api/storage.py`
- Test: `tests/test_service_api_storage.py`

- [ ] **Step 1: Write the failing storage tests**

```python
from pathlib import Path

from lib.service_api.models import RunRecord, TaskBrief
from lib.service_api.storage import FileRunStore


def test_file_run_store_writes_manifest_and_events(tmp_path: Path) -> None:
    store = FileRunStore(tmp_path)
    record = RunRecord(
        run_id="run_001",
        status="queued",
        project_name="api_run_001",
        brief=TaskBrief(prompt="analyze yeast expression strategy"),
    )

    store.write_run(record)
    store.append_event("run_001", {"stage": "planning", "status": "ok"})

    assert (tmp_path / "run_001" / "run_manifest.json").exists()
    assert (tmp_path / "run_001" / "events.jsonl").exists()
```

- [ ] **Step 2: Run the storage tests**

Run: `python -m pytest tests/test_service_api_storage.py -q`
Expected: FAIL because storage has not been implemented.

- [ ] **Step 3: Implement file-backed storage**

The store should support:

- `write_run()`
- `read_run()`
- `append_event()`
- `write_plan_preview()`
- `list_artifacts()`

Phase 1 paths should be deterministic and human-readable.

- [ ] **Step 4: Run the storage tests**

Run: `python -m pytest tests/test_service_api_storage.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add lib/service_api/storage.py tests/test_service_api_storage.py
git commit -m "feat: add file-backed api run storage"
```

## Chunk 3: Execution Layer And HTTP Entry

### Task 5: Add workspace CLI execution adapter

**Files:**
- Create: `lib/service_api/adapters/__init__.py`
- Create: `lib/service_api/adapters/base.py`
- Create: `lib/service_api/adapters/cli.py`
- Test: `tests/test_service_api_execution.py`

- [ ] **Step 1: Write the failing execution tests**

```python
import subprocess

from lib.service_api.adapters.cli import WorkspaceCLIAdapter


def test_workspace_cli_adapter_invokes_project_cli(monkeypatch) -> None:
    calls = []

    def fake_run(cmd, cwd, env, capture_output, text):
        calls.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, stdout="ok\n", stderr="")

    adapter = WorkspaceCLIAdapter("/repo/scripts/project.py", "/repo")
    monkeypatch.setattr("lib.service_api.adapters.cli.subprocess.run", fake_run)

    result = adapter.invoke(["workspace-validate"])
    assert result.returncode == 0
    assert calls
```

- [ ] **Step 2: Run the execution tests**

Run: `python -m pytest tests/test_service_api_execution.py -q`
Expected: FAIL because adapters do not exist yet.

- [ ] **Step 3: Implement the adapter base and CLI adapter**

The adapter result should standardize:

- `command`
- `returncode`
- `stdout`
- `stderr`
- `started_at`
- `finished_at`

- [ ] **Step 4: Add orchestration helpers**

Inside `lib/service_api/execution.py`, add:

- `create_run_record()`
- `submit_run()`
- `execute_run_sync()`

The executor should:

- create a deterministic project name
- call planner preview first
- persist the plan
- invoke the workspace CLI as needed
- update run status

- [ ] **Step 5: Run the execution tests**

Run: `python -m pytest tests/test_service_api_execution.py -q`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add lib/service_api/adapters/__init__.py lib/service_api/adapters/base.py lib/service_api/adapters/cli.py lib/service_api/execution.py tests/test_service_api_execution.py
git commit -m "feat: add workspace cli adapter and sync run execution"
```

### Task 6: Add the HTTP server entrypoint

**Files:**
- Create: `scripts/api_server.py`
- Test: `tests/test_service_api_http.py`
- Modify: `requirements.txt`

- [ ] **Step 1: Write the failing HTTP tests**

```python
from fastapi.testclient import TestClient

from scripts.api_server import app


def test_healthz_returns_ok() -> None:
    client = TestClient(app)
    response = client.get("/healthz")
    assert response.status_code == 200
    assert response.json()["status"] == "ok"


def test_plan_preview_endpoint_returns_preview() -> None:
    client = TestClient(app)
    response = client.post("/v1/plans/preview", json={"prompt": "analyze ebola benchmarks"})
    assert response.status_code == 200
    assert "brief" in response.json()
```

- [ ] **Step 2: Run the HTTP tests**

Run: `python -m pytest tests/test_service_api_http.py -q`
Expected: FAIL because the server entrypoint does not exist yet.

- [ ] **Step 3: Add the minimal FastAPI app**

Expose:

- `GET /healthz`
- `GET /v1/capabilities`
- `POST /v1/plans/preview`
- `POST /v1/runs`
- `GET /v1/runs/{run_id}`

The first HTTP implementation can keep everything in one file if route count is still small.

- [ ] **Step 4: Add or isolate the web dependencies**

Either:

- add `fastapi` and `uvicorn` to `requirements.txt`

or:

- create a dedicated requirements file for the control plane

Choose one path and document it in the README updates.

- [ ] **Step 5: Run the HTTP tests**

Run: `python -m pytest tests/test_service_api_http.py -q`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add scripts/api_server.py tests/test_service_api_http.py requirements.txt
git commit -m "feat: expose online api control plane endpoints"
```

## Chunk 4: Docs And Integration Cleanup

### Task 7: Document the online API surface

**Files:**
- Modify: `README.md`
- Modify: `docs/README.md`

- [ ] **Step 1: Update the top-level README**

Add:

- what the online API control plane is
- how it differs from the current CLI and MCP layers
- how to start the API server locally

- [ ] **Step 2: Update the docs index**

Add links to:

- the online API design spec
- the online API phase-1 plan

- [ ] **Step 3: Verify docs references**

Run: `python -m pytest tests/test_mcp_readme_paths.py -q`
Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add README.md docs/README.md
git commit -m "docs: add online api control plane documentation"
```

### Task 8: Optional MCP bridge for planning preview

**Files:**
- Modify: `mcp-servers/bio-lab-mcp/lab_server.py`
- Modify: `tests/test_lab_mcp_server.py`

- [ ] **Step 1: Write the failing MCP test**

Add a test asserting that `bio-lab-mcp` exposes a new tool such as `preview_workspace_plan`.

- [ ] **Step 2: Run the MCP tests**

Run: `python -m pytest tests/test_lab_mcp_server.py -q`
Expected: FAIL because the tool does not exist yet.

- [ ] **Step 3: Implement the thin MCP wrapper**

The wrapper should call the same planner used by the HTTP API.

Do not duplicate planning logic inside the MCP server.

- [ ] **Step 4: Run the MCP tests**

Run: `python -m pytest tests/test_lab_mcp_server.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add mcp-servers/bio-lab-mcp/lab_server.py tests/test_lab_mcp_server.py
git commit -m "feat: expose plan preview through bio-lab-mcp"
```

## Exit Criteria

- [ ] A structured task brief model exists and is tested
- [ ] A default capability registry exists and is tested
- [ ] Plan preview works without execution
- [ ] File-backed run persistence exists and is tested
- [ ] The API can create and inspect runs
- [ ] The executor reuses `scripts/project.py` instead of bypassing the workspace spine
- [ ] Documentation explains how the online API fits with CLI and MCP

## Recommended Execution Order

1. Chunk 1
2. Chunk 2
3. Chunk 3
4. Chunk 4

This keeps the public API from being implemented before the internal contracts are stable.
