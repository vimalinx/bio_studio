# Template Module Integration Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Connect Bio Studio project templates to reusable `lib/modules` wrappers without breaking self-contained template validation.

**Architecture:** Keep generated projects self-validating and runnable outside the main workspace, but let template pipelines prefer a shared workspace runtime when `lib/` is available. Fill the current RNA-seq wrapper gap in `lib/modules`, harden wrapper execution semantics, then route `rnaseq` / `variant` / `phylogeny` templates through a thin runtime adapter instead of duplicating tool calls in every project.

**Tech Stack:** Python 3, pathlib, importlib, pytest, subprocess-based CLI wrappers

---

### Task 1: Lock the module-integration contract with failing tests

**Files:**
- Create: `tests/test_project_template_module_integration.py`
- Modify: `tests/test_project_type_templates.py`

- [ ] Add focused tests for shared runtime fallback behavior and template wiring.
- [ ] Add tests that require RNA-seq wrappers to exist in `lib/modules`.
- [ ] Run the focused pytest target and confirm it fails for the missing integration.

### Task 2: Harden and extend reusable modules

**Files:**
- Modify: `lib/modules/utils.py`
- Modify: `lib/modules/qc.py`
- Modify: `lib/modules/__init__.py`
- Create: `lib/modules/rnaseq.py`

- [ ] Fix shell-command handling in the shared runner utility so wrapper modules can safely use pipes and redirection.
- [ ] Add a `fastp` wrapper and a dedicated RNA-seq module for STAR / featureCounts entrypoints.
- [ ] Export the new RNA-seq module through `lib/modules/__init__.py`.

### Task 3: Add a shared template runtime adapter

**Files:**
- Create: `lib/template_runtime.py`
- Test: `tests/test_project_template_module_integration.py`

- [ ] Implement workspace-root discovery and optional runtime loading helpers.
- [ ] Add small, type-specific execution helpers for `rnaseq`, `variant`, and `phylogeny`.
- [ ] Keep graceful fallback when generated projects are copied outside the main workspace.

### Task 4: Wire generated templates and example project to the runtime

**Files:**
- Modify: `lib/create_project.py`
- Modify: `projects/example_rnaseq/scripts/pipeline.py`
- Modify: `projects/example_rnaseq/README.md`

- [ ] Update generated pipeline content so templates prefer the shared runtime when available.
- [ ] Preserve `--validate`, `--steps`, and hint-only behavior when no workspace runtime is available.
- [ ] Sync `projects/example_rnaseq` back to the generated `rnaseq` template.

### Task 5: Verify and document

**Files:**
- Modify: `docs/AI_ANALYSIS_PROTOCOL.md`
- Modify: `docs/CHANGELOG.md`

- [ ] Run focused module-integration tests, then the existing regression suite.
- [ ] Run `python scripts/project.py workspace-validate` as final workspace-level evidence.
- [ ] Record the template-to-module integration change in docs and changelog.
