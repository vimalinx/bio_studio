# LiteFold Launcher Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add LiteFold preflight and self-hosted launcher wrappers without modifying the vendored LiteFold source tree.

**Architecture:** Extend the existing workspace bridge in `lib/litefold_bridge.py` so it can derive the correct launch context from the vendored LiteFold layout, inspect self-hosted source imports, and report dependency candidates. Expose those capabilities through `scripts/litefold.py` as `preflight` and `start-selfhosted`, with `--dry-run` so the launcher stays testable and safe in CI.

**Tech Stack:** Python standard library (`ast`, `subprocess`, `json`, `shutil`, `pathlib`), pytest

---

### Task 1: Add failing tests for preflight and launcher behavior

**Files:**
- Modify: `tests/test_litefold_bridge.py`
- Test: `tests/test_litefold_bridge.py`

- [ ] **Step 1: Write the failing tests**

Add tests that assert:
- the recommended launch command carries the derived `PYTHONPATH`
- preflight reports external module roots and requirement candidates
- `start-selfhosted --dry-run --json` prints a structured launch plan

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_litefold_bridge.py -q`
Expected: FAIL because `preflight` and `start-selfhosted` do not exist yet and the command output is still incomplete.

- [ ] **Step 3: Write minimal implementation**

Implement in `lib/litefold_bridge.py`:
- launch-context derivation
- candidate requirements discovery
- external import discovery for self-hosted files
- Python module availability probing

Implement in `scripts/litefold.py`:
- `preflight`
- `start-selfhosted`

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_litefold_bridge.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add tests/test_litefold_bridge.py lib/litefold_bridge.py scripts/litefold.py
git commit -m "feat: add litefold launcher preflight"
```

### Task 2: Document and stabilize the new launcher path

**Files:**
- Modify: `README.md`
- Modify: `scripts/ci/run_stable_tests.sh`
- Test: `tests/test_litefold_bridge.py`

- [ ] **Step 1: Update the README**

Document the new commands:
- `python scripts/litefold.py preflight`
- `python scripts/litefold.py start-selfhosted --dry-run`

- [ ] **Step 2: Verify stable tests include LiteFold bridge coverage**

Check that `tests/test_litefold_bridge.py` remains in `scripts/ci/run_stable_tests.sh`.

- [ ] **Step 3: Run focused verification**

Run: `python scripts/litefold.py preflight --json`
Expected: exit 0 and JSON output describing launch context and dependency hints.

- [ ] **Step 4: Run full stable verification**

Run: `bash scripts/ci/run_stable_tests.sh`
Expected: PASS with zero failures.

- [ ] **Step 5: Commit**

```bash
git add README.md scripts/ci/run_stable_tests.sh
git commit -m "docs: describe litefold preflight launcher"
```
