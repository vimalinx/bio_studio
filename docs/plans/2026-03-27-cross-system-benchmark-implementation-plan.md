# Cross-System Benchmark Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a new `projects/cross_system_benchmark` special project that curates three public benchmarks and runs lightweight processed-data baselines for Ebola, SARS-CoV-2, and yeast.

**Architecture:** Create a project-local workflow instead of forcing the generic template to own all behavior. Keep benchmark metadata curated in project files, use small project scripts to normalize GEO processed files into a shared format, and generate per-dataset plus integrated markdown reports. Leave raw SRA reanalysis as a future extension point.

**Tech Stack:** Python 3, pathlib, csv/json/gzip/urllib, PyYAML, pandas, matplotlib, Bio Studio project CLI and validation helpers

---

### Task 1: Create The Special Project Skeleton

**Files:**
- Create: `projects/cross_system_benchmark/README.md`
- Create: `projects/cross_system_benchmark/scripts/config.py`
- Create: `projects/cross_system_benchmark/scripts/validate_project.py`
- Create: `projects/cross_system_benchmark/scripts/pipeline.py`
- Create: `projects/cross_system_benchmark/metadata/.gitkeep`
- Test: `tests/test_cross_system_benchmark_project.py`

**Step 1: Write the failing test**

```python
def test_cross_system_benchmark_project_exists():
    project_root = ROOT / "projects" / "cross_system_benchmark"
    assert project_root.exists()
    assert (project_root / "scripts" / "pipeline.py").exists()
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_cross_system_benchmark_project.py -v`
Expected: FAIL because project files do not exist yet

**Step 3: Write minimal implementation**

- Create the project with the Bio Studio CLI using `generic`
- Replace the generated README and pipeline with special-project wording
- Add a `metadata/` directory because this project is benchmark-driven

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_cross_system_benchmark_project.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add projects/cross_system_benchmark tests/test_cross_system_benchmark_project.py
git commit -m "feat: scaffold cross system benchmark project"
```

### Task 2: Materialize The Curated Benchmark Catalog

**Files:**
- Create: `projects/cross_system_benchmark/scripts/fetch_benchmarks.py`
- Modify: `projects/cross_system_benchmark/scripts/config.py`
- Modify: `projects/cross_system_benchmark/scripts/pipeline.py`
- Create: `tests/test_cross_system_benchmark_catalog.py`

**Step 1: Write the failing test**

```python
def test_fetch_benchmarks_writes_catalog_files(tmp_path):
    result = run_fetch(tmp_path)
    assert (tmp_path / "metadata" / "datasets.yaml").exists()
    assert (tmp_path / "metadata" / "samples.tsv").exists()
    assert (tmp_path / "metadata" / "comparisons.tsv").exists()
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_cross_system_benchmark_catalog.py -v`
Expected: FAIL because the fetch script does not exist yet

**Step 3: Write minimal implementation**

- Encode the three curated datasets in Python data structures
- Write deterministic outputs:
  - `datasets.yaml`
  - `samples.tsv`
  - `comparisons.tsv`
  - `source_links.md`
- Add a `fetch` stage to `pipeline.py`

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_cross_system_benchmark_catalog.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add projects/cross_system_benchmark/scripts tests/test_cross_system_benchmark_catalog.py
git commit -m "feat: add cross system benchmark catalog generation"
```

### Task 3: Add Reference Preparation With Small Default Downloads

**Files:**
- Create: `projects/cross_system_benchmark/scripts/prepare_references.py`
- Modify: `projects/cross_system_benchmark/scripts/config.py`
- Modify: `projects/cross_system_benchmark/scripts/pipeline.py`
- Create: `tests/test_cross_system_benchmark_references.py`

**Step 1: Write the failing test**

```python
def test_prepare_references_writes_manifest_without_large_host_downloads(tmp_path):
    result = run_prepare(tmp_path)
    manifest = tmp_path / "data" / "references" / "reference_manifest.json"
    assert manifest.exists()
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_cross_system_benchmark_references.py -v`
Expected: FAIL because the prepare script does not exist yet

**Step 3: Write minimal implementation**

- Download or stage only the small default references:
  - Ebola viral reference
  - SARS-CoV-2 viral reference
  - Yeast reference metadata or small reference assets
- Record host reference recommendations in `reference_manifest.json`
- Do not auto-download large human references by default

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_cross_system_benchmark_references.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add projects/cross_system_benchmark/scripts tests/test_cross_system_benchmark_references.py
git commit -m "feat: add cross system reference preparation"
```

### Task 4: Download And Normalize Processed Benchmark Data

**Files:**
- Create: `projects/cross_system_benchmark/scripts/run_baseline.py`
- Modify: `projects/cross_system_benchmark/scripts/config.py`
- Modify: `projects/cross_system_benchmark/scripts/pipeline.py`
- Create: `tests/test_cross_system_benchmark_baseline.py`

**Step 1: Write the failing test**

```python
def test_run_baseline_normalizes_each_dataset_to_common_outputs(tmp_path):
    outputs = run_baseline(tmp_path, offline=True)
    assert (tmp_path / "data" / "processed" / "ebola" / "expression.tsv").exists()
    assert (tmp_path / "data" / "processed" / "sars_cov_2" / "expression.tsv").exists()
    assert (tmp_path / "data" / "processed" / "yeast" / "expression.tsv").exists()
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_cross_system_benchmark_baseline.py -v`
Expected: FAIL because no baseline script exists

**Step 3: Write minimal implementation**

- Support two modes:
  - online: download GEO processed files from curated URLs
  - offline: read tiny local fixtures for tests
- Normalize each dataset into a common `expression.tsv`
- Write per-dataset summary stats and small result tables

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_cross_system_benchmark_baseline.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add projects/cross_system_benchmark/scripts tests/test_cross_system_benchmark_baseline.py
git commit -m "feat: add processed-data baseline workflow"
```

### Task 5: Generate Per-Dataset And Integrated Reports

**Files:**
- Create: `projects/cross_system_benchmark/scripts/summarize_cross_systems.py`
- Modify: `projects/cross_system_benchmark/scripts/run_baseline.py`
- Modify: `projects/cross_system_benchmark/scripts/pipeline.py`
- Create: `tests/test_cross_system_benchmark_reports.py`

**Step 1: Write the failing test**

```python
def test_report_stage_writes_dataset_and_integrated_markdown(tmp_path):
    run_report(tmp_path)
    assert (tmp_path / "data" / "results" / "per_dataset" / "ebola_summary.md").exists()
    assert (tmp_path / "data" / "results" / "integrated" / "cross_system_summary.md").exists()
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_cross_system_benchmark_reports.py -v`
Expected: FAIL because report generation does not exist yet

**Step 3: Write minimal implementation**

- Emit one markdown summary per dataset
- Emit one integrated markdown summary that compares only process-level themes
- Keep language explicit about comparability limits

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_cross_system_benchmark_reports.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add projects/cross_system_benchmark/scripts tests/test_cross_system_benchmark_reports.py
git commit -m "feat: add cross system benchmark reporting"
```

### Task 6: Tighten Validation And Workspace Integration

**Files:**
- Modify: `projects/cross_system_benchmark/scripts/validate_project.py`
- Modify: `projects/cross_system_benchmark/README.md`
- Modify: `tests/test_workspace_project_cli.py`
- Modify: `tests/test_special_project_minimal_integration.py`

**Step 1: Write the failing test**

```python
def test_workspace_cli_can_validate_and_list_steps_for_cross_system_benchmark():
    result = subprocess.run([...], check=True, capture_output=True, text=True)
    assert "fetch" in result.stdout
    assert "report" in result.stdout
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_workspace_project_cli.py tests/test_special_project_minimal_integration.py -v`
Expected: FAIL because the new project is not wired into existing expectations

**Step 3: Write minimal implementation**

- Extend project validation checks for `metadata/`
- Make sure `pipeline.py --steps` exposes `fetch`, `prepare`, `baseline`, `report`
- Update README wording so the project is clearly classified as a special project

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_workspace_project_cli.py tests/test_special_project_minimal_integration.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add projects/cross_system_benchmark tests/test_workspace_project_cli.py tests/test_special_project_minimal_integration.py
git commit -m "feat: wire cross system benchmark into workspace cli"
```

### Task 7: Run End-To-End Verification

**Files:**
- Modify: `projects/cross_system_benchmark/logs/` as generated output only
- Test: `tests/test_cross_system_benchmark_project.py`
- Test: `tests/test_cross_system_benchmark_catalog.py`
- Test: `tests/test_cross_system_benchmark_references.py`
- Test: `tests/test_cross_system_benchmark_baseline.py`
- Test: `tests/test_cross_system_benchmark_reports.py`

**Step 1: Run focused test suite**

Run:

```bash
pytest \
  tests/test_cross_system_benchmark_project.py \
  tests/test_cross_system_benchmark_catalog.py \
  tests/test_cross_system_benchmark_references.py \
  tests/test_cross_system_benchmark_baseline.py \
  tests/test_cross_system_benchmark_reports.py \
  tests/test_workspace_project_cli.py \
  tests/test_special_project_minimal_integration.py -v
```

Expected: PASS

**Step 2: Run project validation**

Run: `python scripts/project.py validate cross_system_benchmark`
Expected: PASS and `logs/validation_report.json` written

**Step 3: Run steps inspection**

Run: `python scripts/project.py steps cross_system_benchmark`
Expected: output lists `fetch`, `prepare`, `baseline`, `report`

**Step 4: Run the baseline pipeline**

Run:

```bash
python scripts/project.py run cross_system_benchmark -- fetch
python scripts/project.py run cross_system_benchmark -- prepare
python scripts/project.py run cross_system_benchmark -- baseline
python scripts/project.py run cross_system_benchmark -- report
```

Expected: metadata, processed outputs, per-dataset summaries, and integrated summary created

**Step 5: Commit**

```bash
git add projects/cross_system_benchmark docs/plans tests
git commit -m "feat: add cross system benchmark project"
```
