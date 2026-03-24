# Special Project Minimal Integration Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Give specialist and learning projects a minimal unified validation/steps surface without forcing them into the standard shared runtime template.

**Architecture:** Keep `projects/test_rnaseq_analysis` as a specialist research project and `projects/yeast_genome_learning` as a learning project. Add lightweight `validate_project.py` entrypoints and minimal `pipeline.py --validate/--steps` compatibility so workspace-level project commands can operate on them, while preserving their existing bespoke execution paths.

**Tech Stack:** Python 3, pathlib, argparse, pytest, subprocess, markdown docs

---

### Task 1: Lock the minimal integration contract with failing tests

**Files:**
- Create: `tests/test_special_project_minimal_integration.py`

- [ ] Step 1: Write tests requiring both projects to declare their special classification in README.
- [ ] Step 2: Write tests requiring project-level `validate_project.py` entrypoints.
- [ ] Step 3: Write tests requiring `pipeline.py --validate` and `pipeline.py --steps` compatibility for both projects.
- [ ] Step 4: Run the focused test file and verify it fails for the missing integration.

### Task 2: Add minimal compatibility for `test_rnaseq_analysis`

**Files:**
- Create: `projects/test_rnaseq_analysis/scripts/validate_project.py`
- Modify: `projects/test_rnaseq_analysis/scripts/pipeline.py`
- Modify: `projects/test_rnaseq_analysis/README.md`

- [ ] Step 1: Add a lightweight validation script tailored to the project’s existing structure.
- [ ] Step 2: Extend the existing Python pipeline with `--validate` and `--steps` only.
- [ ] Step 3: Mark the README as a specialist project that intentionally stays outside the generic template runtime.

### Task 3: Add minimal compatibility for `yeast_genome_learning`

**Files:**
- Create: `projects/yeast_genome_learning/scripts/validate_project.py`
- Create: `projects/yeast_genome_learning/scripts/pipeline.py`
- Modify: `projects/yeast_genome_learning/README.md`
- Modify: `projects/yeast_genome_learning/QUICKSTART.md`

- [ ] Step 1: Add a project-specific validation script focused on learning assets and shell workflows.
- [ ] Step 2: Add a small Python compatibility pipeline that exposes `--validate` and `--steps` and points to the shell-based lessons.
- [ ] Step 3: Mark the README and quickstart as a learning project outside the generic template runtime.

### Task 4: Document and verify

**Files:**
- Modify: `docs/WORKSPACE_ARCHITECTURE.md`
- Modify: `docs/CHANGELOG.md`

- [ ] Step 1: Document the distinction between template projects and special projects.
- [ ] Step 2: Run focused tests, then broader regression.
- [ ] Step 3: Run `python scripts/project.py workspace-validate` as final evidence.
