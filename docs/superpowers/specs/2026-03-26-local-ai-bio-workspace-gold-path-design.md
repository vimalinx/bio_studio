# Local AI Bio Workspace Gold Path Design

## Overview

This document captures the current assessment and next-step design for Bio Studio as a local-first AI-driven bioinformatics workspace.

The platform goal is not to become a single-purpose yeast tool, RNA-seq tool, or DNA design script collection. The goal is to become a general local AI orchestration workspace for biological tasks, where AI interprets a task, creates a project, selects the right execution profile, runs local and online capabilities, and writes reproducible outputs into `projects/`.

The immediate product goal is narrower:

- Prove the workspace is not an empty shell
- Establish one credible end-to-end gold path
- Keep the platform general while validating it through one reference profile

## Current Assessment

### What Already Exists

Bio Studio already has a meaningful platform base:

- A workspace-level entrypoint in `scripts/project.py`
- A `projects/` model for task instances
- Shared runtime modules in `lib/`
- Five MCP servers covering design, lab orchestration, sequence, structure, and database access
- Bridges to external engines such as LiteFold
- Stable engineering guardrails with passing tests

The current stable test chain passed locally with `86 passed`, which indicates the workspace infrastructure is already beyond an ad hoc script collection.

### Current Completion Estimate

Two different layers should be scored separately:

1. General AI bio workspace foundation: about `70/100`
2. Credible AI-driven design gold path: about `20-30/100`

This gap is expected. The platform base exists, but the vertical path from natural language goal to structured project outputs has not been made real yet.

### Verified Gaps In The Current Design Layer

The existing design surface is not yet a true gold path:

- `bio-design-mcp` is currently a wrapper around point tools, not a decision-making orchestration path
- `tools/scripts/dna_design.py` does not currently support a yeast path correctly; a direct test in the `bio` environment raised a `TypeError`
- `tools/scripts/mrna_optimize.py` can run, but its defaults are still human-centered rather than profile-driven
- There is no single path that goes from natural language task to project materialization, candidate generation, ranking, and report output

## Product Positioning

Bio Studio should be treated as:

`a local-first, AI-driven biological task orchestration workspace`

It should not be defined by one organism or one workflow. Instead, it should provide a common backbone for multiple task families.

Examples of future task families:

- DNA or expression construct design
- mRNA design
- peptide or protein design
- target research and evidence synthesis
- analysis workflows such as RNA-seq

The platform stays general. Specific biological assumptions live in profiles.

## Design Principles

- Keep the backbone general and the biology-specific logic in profiles
- Materialize every run into a project, not just a chat answer
- Prefer multiple candidates plus ranking over single-answer generation
- Keep first-stage claims at the level of proxy evaluation, not wet-lab truth
- Use local execution as the default, with online resources allowed when needed
- Treat one fully working reference path as more valuable than many half-built capability surfaces

## Backbone Architecture

The platform backbone should be organized into five layers.

### 1. Task Interpreter

Convert natural language into a structured task brief.

Minimum fields:

- `task_type`
- `molecule_type`
- `candidate_hosts`
- `output_granularity`
- `network_policy`
- `success_criteria`
- `risk_tolerance`
- `reference_constraints`

This is the layer that turns a vague user request into a reproducible machine-readable task.

### 2. Profile Registry

Profiles hold domain-specific policy and assets.

Examples:

- `yeast_expression_profile`
- `ecoli_expression_profile`
- `mammalian_mrna_profile`
- `peptide_design_profile`
- `rnaseq_analysis_profile`

Each profile should declare:

- What task shapes it supports
- What inputs it expects
- What outputs it can generate
- What tools and MCP servers it needs
- How it scores candidates
- What failure modes it recognizes

### 3. Orchestrator

The orchestrator is the platform spine.

Responsibilities:

- Read the structured task brief
- Choose or rank profiles
- Create the project directory
- Plan stages
- Invoke MCP servers, local scripts, and external engines
- Retry or fail with structured reasons
- Track provenance across the run

This is distinct from the current workspace CLI. `scripts/project.py` is a good entrypoint foundation, but it is not yet an AI orchestrator.

### 4. Candidate And Evaluation Ledger

The platform should manage multiple candidate solutions rather than only one final answer.

For each candidate, store:

- Origin and generation rationale
- Profile and host assumptions
- CDS, mRNA, or construct draft outputs
- Scoring breakdown
- Evidence links
- Risk notes
- Uncertainty notes

This layer is critical. Without it, the platform remains a generator rather than a decision-support workspace.

### 5. Project Materializer

Every run should leave standard artifacts in `projects/`.

Recommended outputs:

- `design_brief.json`
- `run_manifest.json`
- `profile_decision.json`
- `candidate_table.json`
- `candidate_table.csv`
- `constructs/`
- `evaluation_report.md`
- `final_recommendation.md`

This is the step that makes the workspace reproducible and reviewable.

## Gold Path Definition

The first gold path should prove that the platform can reliably turn an open-ended design request into structured outputs.

Recommended MVP gold path:

`natural language goal -> structured task brief -> profile selection -> project creation -> candidate generation -> proxy evaluation -> ranking -> project artifacts`

### User-Facing Behavior

The user gives a natural-language design goal.

The system then:

1. Interprets the task
2. Creates a new project
3. Chooses an appropriate reference profile
4. Generates multiple candidate designs
5. Produces CDS, mRNA, or construct drafts as appropriate
6. Runs proxy evaluations
7. Ranks candidates
8. Writes all outputs into the project directory

### Internal Design Reality

Even if the platform outputs DNA- or construct-level artifacts, internal reasoning should not blindly search in raw DNA space first.

For design tasks, the more credible internal flow is:

`functional intent -> candidate molecule design -> construct back-translation and expression packaging`

The platform may output DNA-level artifacts, but should not confuse output granularity with the internal optimization layer.

## Reference Profile Strategy

The platform should remain general, but one reference profile must be implemented fully enough to validate the architecture.

This reference profile should be treated as:

- A proof path
- A contract test for the backbone
- A template for future profiles

It should not be treated as the permanent platform identity.

## Development Priorities

### Priority 0: Task Brief Schema

Define the structured task brief first.

Without this:

- AI interpretation remains inconsistent
- Runs are hard to compare
- Failures are hard to diagnose
- Later automation becomes brittle

### Priority 1: Profile Contract

Define the interface every profile must satisfy before building many profiles.

Without this:

- Logic leaks into prompts and ad hoc scripts
- Host-specific or task-specific assumptions become tangled
- The platform cannot scale beyond the first example

### Priority 2: One Fully Working Reference Profile

Do not try to support every host or task type at once.

Build one reference profile that proves:

- task interpretation works
- profile selection works
- candidate generation works
- evaluation works
- project materialization works

### Priority 3: Candidate Ledger And Ranking

The first real platform milestone is not single-answer generation. It is ranked multi-candidate output with traceable reasons.

### Priority 4: Project Materialization

Every gold-path run must leave standard artifacts in `projects/`.

This is required to make the workspace feel like a real platform instead of an agent demo.

## Recommended First Three Development Steps

1. Define the task brief schema and standard project outputs
2. Define the profile contract and create the registry shape
3. Implement one end-to-end orchestrated gold path using one reference profile

This order is deliberate.

If the orchestrator is built before schema and contract, the system will drift into tool-specific improvisation. If many profiles are built before the backbone contract, the system will fragment.

## Resource Reality

The current machine is strong enough for a serious local-first prototype:

- 16 CPU cores
- about 31 GB RAM
- one CUDA GPU with about 12 GB VRAM
- substantial remaining disk

This supports:

- local orchestration
- local candidate generation
- moderate local structure or proxy evaluation workloads
- online evidence augmentation

This does not justify pretending that the system can already do:

- large-scale exhaustive candidate search
- wet-lab truth validation
- highly reliable end-to-end therapeutic design decisions

The correct near-term claim level is:

`AI-assisted design and proxy evaluation platform`

not:

`fully validated autonomous biological R&D platform`

## Non-Goals For The First Stage

The first stage should explicitly avoid:

- claiming experimental success
- supporting all hosts or all task families equally
- producing vendor-ready construct orders for every case
- presenting proxy function scores as biological truth

## Success Criteria

The gold path should be considered successful when:

- The same style of natural language task can be turned into a structured brief reproducibly
- A project is created automatically with standard artifacts
- The system selects a reference profile and records why
- The system produces multiple candidates, not just one answer
- Each candidate has traceable rationale, outputs, and scores
- The platform can explain failures in structured terms

## Final Recommendation

Bio Studio should not expand by accumulating more one-off tools first.

The next stage should focus on creating the platform spine:

- task brief schema
- profile contract
- orchestrated gold path
- candidate ledger
- standardized project outputs

The foundation is already strong enough to support this move. The missing piece is not more raw capability. The missing piece is a coherent path that turns AI intent into durable project artifacts.
