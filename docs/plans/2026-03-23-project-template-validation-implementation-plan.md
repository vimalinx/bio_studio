# Project Template Validation Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Upgrade the generated Bio Studio project template so each new project can run a built-in self-validation flow immediately after creation.

**Architecture:** Keep workspace-level validation in `projects/test_env_validation`, and add project-level validation to generated projects via a self-contained `scripts/validate_project.py` plus a `pipeline.py --validate` entrypoint. Update the example project and protocol docs so templates, samples, and documentation stay aligned.

**Tech Stack:** Python 3, pytest, pathlib, importlib, subprocess-free project validation scripts

---

### Task 1: Add failing tests for the new template contract

**Files:**
- Create: `tests/test_project_template_validation.py`

- [ ] **Step 1: Write failing tests for generated template files**
- [ ] **Step 2: Run the tests and verify they fail because the new validation artifacts do not exist yet**

### Task 2: Upgrade the generated project template

**Files:**
- Modify: `lib/create_project.py`

- [ ] **Step 1: Expand generated `config.py` with path objects and richer defaults**
- [ ] **Step 2: Generate `scripts/validate_project.py`**
- [ ] **Step 3: Update generated `scripts/pipeline.py` to support `--validate`**
- [ ] **Step 4: Update generated `README.md` with validation commands**

### Task 3: Refresh the sample project

**Files:**
- Modify: `projects/example_rnaseq/README.md`
- Modify: `projects/example_rnaseq/scripts/config.py`
- Modify: `projects/example_rnaseq/scripts/pipeline.py`
- Create: `projects/example_rnaseq/scripts/validate_project.py`

- [ ] **Step 1: Bring the sample project into parity with the generated template**
- [ ] **Step 2: Keep the example lightweight and self-validating**

### Task 4: Update protocol docs

**Files:**
- Modify: `docs/AI_ANALYSIS_PROTOCOL.md`
- Modify: `docs/CHANGELOG.md`

- [ ] **Step 1: Document the new project validation workflow**
- [ ] **Step 2: Record the template upgrade in the changelog**

### Task 5: Verify end-to-end

**Files:**
- Test: `tests/test_project_template_validation.py`

- [ ] **Step 1: Run the focused pytest suite**
- [ ] **Step 2: Confirm generated project validation writes report files**
- [ ] **Step 3: Summarize exact verification evidence**
