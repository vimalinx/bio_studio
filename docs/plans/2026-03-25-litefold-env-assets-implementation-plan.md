# LiteFold Env Assets Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add workspace-owned LiteFold environment assets and expose them through preflight without modifying the vendored LiteFold source tree.

**Architecture:** Keep the installation assets outside the LiteFold submodule. Extend the bridge to scan both `selfhosted` and `fold_models` imports, resolve workspace-level asset paths, and report the install command. Add a dry-run shell setup script that installs from the repo-owned LiteFold requirements file using an explicit Python executable.

**Tech Stack:** Python standard library, Bash, pytest

---

### Task 1: Add failing tests for env asset discovery and dry-run setup

**Files:**
- Modify: `tests/test_litefold_bridge.py`
- Test: `tests/test_litefold_bridge.py`

- [ ] **Step 1: Write the failing tests**

Add tests that assert:
- `preflight` includes workspace requirements/setup asset paths
- `preflight` scans `fold_models` imports and reports packages like `esm` and `einops`
- `scripts/setup/setup_litefold_env.sh --dry-run` prints the pip install plan for a target workspace

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_litefold_bridge.py -q`
Expected: FAIL because the bridge does not yet expose env assets and the setup script does not exist.

- [ ] **Step 3: Write minimal implementation**

Implement:
- workspace asset resolution in `lib/litefold_bridge.py`
- recursive import scan for `fold_models`
- repo-owned `requirements-litefold-selfhosted.txt`
- `scripts/setup/setup_litefold_env.sh`

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_litefold_bridge.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add tests/test_litefold_bridge.py lib/litefold_bridge.py requirements-litefold-selfhosted.txt scripts/setup/setup_litefold_env.sh
git commit -m "feat: add litefold env setup assets"
```

### Task 2: Expose install guidance in docs and preflight verification

**Files:**
- Modify: `README.md`
- Modify: `scripts/litefold.py`
- Test: `tests/test_litefold_bridge.py`

- [ ] **Step 1: Update preflight output and README**

Document and print:
- repo-owned LiteFold requirements file
- setup script path
- dry-run install command

- [ ] **Step 2: Run focused verification**

Run: `python scripts/litefold.py preflight --json`
Expected: output includes installation asset paths and missing module hints.

- [ ] **Step 3: Run full stable verification**

Run: `bash scripts/ci/run_stable_tests.sh`
Expected: PASS with zero failures.

- [ ] **Step 4: Commit**

```bash
git add README.md scripts/litefold.py
git commit -m "docs: expose litefold env setup guidance"
```
