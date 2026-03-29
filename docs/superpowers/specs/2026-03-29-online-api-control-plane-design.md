# Online API Control Plane Design

## Overview

This document defines how Bio Studio should evolve from a local-first AI bio workspace into an online API control plane.

The key design choice is:

`do not replace the current workspace spine`

Instead:

`add a service layer above the existing workspace CLI, MCP servers, and project artifact model`

That means the future online system should still reuse:

- `projects/` as the durable artifact boundary
- `scripts/project.py` as an execution entrypoint
- `mcp-servers/` as the AI-facing tool surface
- `lib/` as the shared runtime and adapter layer

The online API is therefore a new control plane, not a rewrite of the current platform.

## Product Positioning

After this upgrade, Bio Studio should be described as:

`an AI-driven biological task orchestration platform with both local and online control surfaces`

The local workspace remains the execution substrate.

The online API becomes the remote orchestration surface for:

- browser clients
- external applications
- batch submissions
- agent-to-agent tool invocation
- future multi-user hosted workflows

## Goals

- Accept natural-language or structured task requests over HTTP
- Normalize requests into a machine-readable task brief
- Match requests to the right tools, profiles, and execution backends
- Create durable runs with IDs, status transitions, logs, and artifacts
- Preserve the existing `projects/` artifact model for reproducibility
- Keep the first online version thin and operationally simple

## Non-Goals For Phase 1

- Full multi-tenant SaaS isolation
- Distributed heavy compute scheduling across multiple machines
- End-user billing
- Arbitrary free-form LLM orchestration with no schema constraints
- Replacing all local scripts with network-native microservices

## What Already Exists And Should Be Reused

### 1. Workspace execution spine

Current reusable entrypoints already exist:

- `scripts/project.py`
- project-local `scripts/pipeline.py`
- project-local `scripts/validate_project.py`

This should remain the authoritative execution boundary for projectized runs.

### 2. AI tool surface

The repo already has multiple MCP servers:

- `bio-lab-mcp`
- `bio-design-mcp`
- `bio-sequence-mcp`
- `bio-structure-mcp`
- `bio-database-mcp`

These do not need to be discarded. They should be wrapped as capability providers inside the online control plane.

### 3. Shared runtime

Shared Python modules in `lib/` already provide:

- project creation
- environment resolution
- project validation helpers
- domain modules

The online API should call these shared modules directly where possible and fall back to CLI/MCP transport only where isolation or compatibility matters.

### 4. Artifact discipline

The current workspace already treats `projects/` as the durable landing zone for logs, reports, and analysis outputs. This is a strong foundation and should remain true after service-ification.

## Core Design Principle

The online API should not directly “pick a random script and run it.”

It should always follow this path:

`request -> task brief -> capability match -> execution plan -> run record -> artifacts`

If this contract is violated, the platform becomes hard to debug, hard to reproduce, and impossible to govern once tool count grows.

## Target Architecture

The target system should have eight layers.

### 1. API Gateway Layer

Responsibility:

- receive HTTP requests
- validate authentication
- assign request IDs
- expose run creation, status, logs, artifacts, and capability discovery APIs

Recommended initial stack:

- FastAPI
- Uvicorn
- Pydantic

### 2. Request Interpreter

Responsibility:

- normalize raw user requests into a structured task brief

Minimum task brief fields:

- `task_type`
- `intent_summary`
- `input_mode`
- `molecule_type`
- `candidate_hosts`
- `output_granularity`
- `network_policy`
- `success_criteria`
- `risk_tolerance`
- `requested_capabilities`
- `constraints`

Phase 1 should prefer deterministic parsing and explicit overrides.

LLM-assisted parsing can be added later, but only after the schema and fallback rules are stable.

### 3. Capability Registry

Responsibility:

- declare what the system can do
- standardize routing metadata for tools, profiles, and backends

Every capability entry should include:

- `capability_id`
- `display_name`
- `task_types`
- `keywords`
- `input_schema_ref`
- `output_schema_ref`
- `transport`
- `backend_target`
- `network_requirements`
- `resource_class`
- `supports_preview`
- `supports_async`
- `stability_tier`

Examples:

- `workspace.project.create`
- `workspace.project.validate`
- `sequence.analysis.basic`
- `structure.analysis.fold`
- `database.search.pubmed`
- `design.mrna.optimize`

### 4. Planner And Router

Responsibility:

- turn a task brief into an execution plan
- select the right profile and capabilities
- produce a reviewable plan before execution

The planner should support two modes:

- `preview`: return a dry-run execution plan without doing work
- `execute`: enqueue or run the plan

The first implementation should be rules-first:

- profile rules
- keyword matching
- task-type matching
- network-policy checks
- capability compatibility checks

Later versions can add LLM ranking as a secondary scorer, not the source of truth.

### 5. Execution Adapters

Responsibility:

- hide transport differences between local Python, CLI, and MCP

There should be three adapter classes:

- `DirectPythonAdapter`
- `WorkspaceCLIAdapter`
- `MCPToolAdapter`

This is how the control plane reuses the current codebase without flattening everything into one style.

### 6. Run Store

Responsibility:

- persist run state
- persist plan snapshots
- persist event logs
- persist artifact references

Phase 1 storage can be file-first:

- `projects/<run-project>/api/run_manifest.json`
- `projects/<run-project>/api/plan.json`
- `projects/<run-project>/api/events.jsonl`

Phase 2 can add SQLite or PostgreSQL without changing the public API.

### 7. Worker Layer

Responsibility:

- execute async runs
- stream status changes
- write structured events

Phase 1 can use an in-process background worker.

Phase 2 can move to:

- Celery
- RQ
- Dramatiq
- Arq

The worker boundary matters because some biological jobs are long-running and should not block HTTP request threads.

### 8. Artifact Index Layer

Responsibility:

- expose what a run produced
- map files back to URLs and logical artifact types

Artifact types should be standardized:

- `task_brief`
- `execution_plan`
- `log`
- `report`
- `candidate_table`
- `final_recommendation`
- `sequence_file`
- `structure_file`
- `auxiliary_data`

## Request Lifecycle

The online request path should be:

1. Client sends `POST /v1/runs` with either natural-language or structured request payload.
2. API validates payload and creates a run shell.
3. Request interpreter emits a structured task brief.
4. Planner loads capability registry and scores candidate execution paths.
5. Planner returns a normalized execution plan.
6. Executor creates or selects a workspace project directory.
7. Worker invokes CLI, direct Python modules, or MCP tools through adapters.
8. Events and artifact references are written continuously.
9. Client polls or subscribes to run state and artifact outputs.

## Initial API Surface

### Control Endpoints

- `GET /healthz`
- `GET /v1/capabilities`
- `POST /v1/plans/preview`
- `POST /v1/runs`
- `GET /v1/runs/{run_id}`
- `GET /v1/runs/{run_id}/events`
- `GET /v1/runs/{run_id}/artifacts`
- `POST /v1/runs/{run_id}/cancel`

### Suggested Request Contracts

#### `POST /v1/plans/preview`

Used when the caller wants to see decomposition and tool matching without running anything.

Input:

- natural-language prompt or structured task brief
- optional profile overrides
- optional capability allowlist

Output:

- normalized task brief
- selected profile
- ranked capabilities
- planned stages
- execution risk notes

#### `POST /v1/runs`

Used when the caller wants a durable run.

Input:

- request payload
- execution mode: `sync` or `async`
- optional project naming hint
- optional network policy override

Output:

- `run_id`
- `status`
- `project_name`
- `task_brief`
- `plan_summary`

## Routing Strategy

The routing system should be profile-first and capability-scored.

### Profile selection order

1. explicit profile override
2. exact task-type rule
3. heuristic score from task brief fields
4. fallback generic profile

### Capability scoring order

Each capability gets a score from:

- task-type compatibility
- input compatibility
- network-policy compatibility
- output-granularity compatibility
- stability tier
- runtime cost
- explicit user/tool constraints

The planner should always store score rationale in the plan output.

## Safety And Governance

The online control plane introduces risks that the local workspace does not fully face.

Required guardrails:

- never expose raw filesystem paths outside allowed artifact views
- never expose `repositories/active/` as a public writable surface
- restrict which capabilities are remotely invocable
- enforce network policy per run
- write provenance for every adapter invocation
- preserve a dry-run preview mode

## Recommended File Layout

The first service layer should live in the main repo and reuse existing patterns.

```text
lib/service_api/
├── __init__.py
├── models.py
├── registry.py
├── planner.py
├── execution.py
├── storage.py
└── adapters/
    ├── __init__.py
    ├── base.py
    ├── cli.py
    ├── mcp.py
    └── python.py

scripts/
└── api_server.py

tests/
├── test_service_api_models.py
├── test_service_api_registry.py
├── test_service_api_planner.py
├── test_service_api_execution.py
└── test_service_api_http.py
```

This keeps the service layer close to the existing workspace spine and avoids introducing a second disconnected codebase.

## Rollout Strategy

### Phase 0: Control Plane Contracts

Deliver:

- task brief schema
- capability registry
- planner preview
- run state model

No heavy execution required yet.

### Phase 1: Thin Online API

Deliver:

- FastAPI server
- `/v1/plans/preview`
- `/v1/runs`
- synchronous execution through workspace CLI adapters
- file-backed run persistence

This phase is enough to prove remote invocation and tool matching.

### Phase 2: Async Runs And Artifact APIs

Deliver:

- background workers
- event streaming or polling
- artifact listing and download
- cancellation

### Phase 3: Multi-User Hosted Mode

Deliver:

- auth
- rate limiting
- tenant-aware storage
- quota controls
- operational isolation

## Engineering Judgment

This upgrade is medium complexity because the reusable core already exists.

The difficult part is not invoking tools.

The difficult part is standardizing:

- request interpretation
- capability declaration
- plan generation
- run persistence
- operational governance

That is why the first implementation should optimize for explicit contracts and thin adapters instead of ambitious end-to-end autonomy.

## Decision Summary

The correct path is:

- keep Bio Studio as the execution substrate
- add a service control plane above it
- route all online requests through task briefs and capability plans
- preserve `projects/` as the durable artifact boundary

This produces an online system without losing the reproducibility and composability that already make the workspace valuable.
