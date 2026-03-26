# Platform Spine Phase 1 Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the first executable platform spine for Bio Studio so a natural-language design request can be turned into a structured design project with a task brief, reference profile selection, candidate ledger, proxy evaluation outputs, and standard project artifacts.

**Architecture:** Add a focused `lib/design_runs/` package for task-brief modeling, profile contracts, candidate tracking, project materialization, and orchestration. Reuse the existing workspace CLI and MCP servers as front doors, but route them through the new package. Implement one reference expression-design profile as a contract test for the backbone while keeping the platform itself host-agnostic and profile-driven.

**Tech Stack:** Python 3.10+, stdlib `dataclasses`/`json`/`pathlib`, existing `scripts/project.py` CLI, existing MCP server stack, Biopython-backed design scripts, `pytest`

---

## File Structure

### New files

- `lib/design_runs/__init__.py`
  - Package entrypoint for phase-1 design-run spine.
- `lib/design_runs/brief.py`
  - Structured task-brief dataclasses, JSON serialization, and heuristic prompt-to-brief parsing.
- `lib/design_runs/candidates.py`
  - Candidate ledger model and score aggregation helpers.
- `lib/design_runs/materialize.py`
  - Standard artifact writers for design projects.
- `lib/design_runs/orchestrator.py`
  - Main flow: brief -> profile -> candidate generation -> artifact materialization.
- `lib/design_runs/profiles/__init__.py`
  - Profile package export surface.
- `lib/design_runs/profiles/base.py`
  - Base protocol / abstract class for all design profiles.
- `lib/design_runs/profiles/registry.py`
  - Profile registry, lookup, and scoring hooks.
- `lib/design_runs/profiles/yeast_expression.py`
  - Reference profile used only as a contract test for the backbone, not as platform identity.
- `tests/test_design_task_brief.py`
  - Contract tests for task-brief parsing and serialization.
- `tests/test_design_profile_registry.py`
  - Contract tests for profile registration and selection.
- `tests/test_design_project_materialization.py`
  - Artifact path and content tests.
- `tests/test_design_reference_profile.py`
  - Tests for the reference profile’s candidate and evaluation outputs.
- `tests/test_design_orchestrator.py`
  - End-to-end unit tests for orchestration without touching MCP transport.
- `tests/test_workspace_design_cli.py`
  - CLI tests for `scripts/project.py design`.

### Modified files

- `lib/create_project.py`
  - Add a `design` project template so design runs are not forced into a generic shell.
- `scripts/project.py`
  - Add a `design` subcommand and wire it to the new orchestrator.
- `tools/scripts/dna_design.py`
  - Replace the current broken non-`E.coli` path with explicit species codon preferences and deterministic fallback behavior.
- `tools/scripts/mrna_optimize.py`
  - Replace human-only defaults with species-profile defaults and explicit runtime warnings when a profile is incomplete.
- `mcp-servers/bio-lab-mcp/lab_server.py`
  - Expose the new design gold path through an MCP wrapper tool once CLI orchestration exists.
- `tests/test_project_type_templates.py`
  - Add coverage for the new `design` project template.
- `tests/test_workspace_project_cli.py`
  - Add coverage for the new CLI subcommand.
- `tests/test_design_mcp_server.py`
  - Keep existing point-tool behavior stable while new platform-level behavior is added elsewhere.
- `tests/test_lab_mcp_server.py`
  - Add tests for the new MCP wrapper operation.
- `scripts/ci/run_stable_tests.sh`
  - Include the new test files in the stable suite.
- `README.md`
  - Document the new platform-spine path and CLI entrypoint.
- `docs/README.md`
  - Add the plan and new design-run entrypoints to the docs index.
- `mcp-servers/bio-lab-mcp/README.md`
  - Document the new design-project tool.
- `docs/CHANGELOG.md`
  - Record the shipped platform-spine changes.

## Delivery Notes

- Phase 1 does **not** attempt arbitrary de novo biological generation quality.
- The reference profile exists to prove the platform spine and artifact model work.
- Natural-language interpretation in phase 1 should be deterministic and schema-first. It may be heuristic now and upgraded later.
- The platform remains general even if the first concrete profile is yeast-based.

## Chunk 1: Contracts And Project Skeleton

### Task 1: Add a `design` project template

**Files:**
- Modify: `lib/create_project.py`
- Test: `tests/test_project_type_templates.py`
- Test: `tests/test_workspace_project_cli.py`

- [ ] **Step 1: Write the failing template test**

```python
def test_design_template_includes_design_specific_skeleton(tmp_path: Path) -> None:
    project_dir = _generate_project(tmp_path, "design_demo", "design")
    config_text = (project_dir / "scripts" / "config.py").read_text(encoding="utf-8")
    pipeline_text = (project_dir / "scripts" / "pipeline.py").read_text(encoding="utf-8")
    readme_text = (project_dir / "README.md").read_text(encoding="utf-8")

    assert "DESIGN_MODE" in config_text
    assert "NETWORK_POLICY" in config_text
    assert "design_brief" in pipeline_text.lower()
    assert "设计" in readme_text
```

- [ ] **Step 2: Run the focused tests to confirm failure**

Run: `python -m pytest tests/test_project_type_templates.py tests/test_workspace_project_cli.py -q`
Expected: FAIL because `design` is not a supported project type yet.

- [ ] **Step 3: Extend the project template definitions**

Add a `design` entry to `PROJECT_TYPE_PROFILES` and extend CLI type choices:

```python
"design": {
    "readme_overview": "...设计项目骨架说明...",
    "config_extra": """DESIGN_MODE = "gold_path"
NETWORK_POLICY = "allow_online"
PRIMARY_GOAL = ""
REFERENCE_PROFILE = None
""",
    ...
}
```

- [ ] **Step 4: Run tests to confirm the new template works**

Run: `python -m pytest tests/test_project_type_templates.py tests/test_workspace_project_cli.py -q`
Expected: PASS for the new design-template coverage.

- [ ] **Step 5: Commit**

```bash
git add lib/create_project.py tests/test_project_type_templates.py tests/test_workspace_project_cli.py
git commit -m "feat: add design project template"
```

### Task 2: Add structured task-brief modeling

**Files:**
- Create: `lib/design_runs/__init__.py`
- Create: `lib/design_runs/brief.py`
- Test: `tests/test_design_task_brief.py`

- [ ] **Step 1: Write the failing tests for task brief serialization and heuristic parsing**

```python
from lib.design_runs.brief import DesignTaskBrief, brief_from_prompt


def test_design_task_brief_round_trips_to_json(tmp_path: Path) -> None:
    brief = DesignTaskBrief(
        goal="design a secreted stable construct",
        task_type="design",
        molecule_type="protein",
        candidate_hosts=["S_cerevisiae", "K_phaffii"],
        output_granularity="expression_construct",
        network_policy="allow_online",
        success_criteria=["multiple_ranked_candidates"],
        risk_tolerance="medium",
    )
    path = tmp_path / "brief.json"
    brief.write_json(path)

    restored = DesignTaskBrief.read_json(path)
    assert restored.goal == brief.goal
    assert restored.candidate_hosts == ["S_cerevisiae", "K_phaffii"]


def test_brief_from_prompt_populates_required_defaults() -> None:
    brief = brief_from_prompt("Design a stable secreted construct for yeast expression")
    assert brief.task_type == "design"
    assert brief.output_granularity == "expression_construct"
    assert brief.network_policy == "allow_online"
```

- [ ] **Step 2: Run the new test file**

Run: `python -m pytest tests/test_design_task_brief.py -q`
Expected: FAIL because the package and schema do not exist.

- [ ] **Step 3: Implement the minimal schema**

Create a focused dataclass module:

```python
@dataclass
class DesignTaskBrief:
    goal: str
    task_type: str
    molecule_type: str
    candidate_hosts: list[str]
    output_granularity: str
    network_policy: str
    success_criteria: list[str]
    risk_tolerance: str
    reference_constraints: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, object]: ...
    def write_json(self, path: Path) -> None: ...
    @classmethod
    def read_json(cls, path: Path) -> "DesignTaskBrief": ...


def brief_from_prompt(prompt: str) -> DesignTaskBrief:
    ...
```

Keep parsing heuristic and deterministic in phase 1. Do not add LLM-coupled code here.

- [ ] **Step 4: Run the schema tests**

Run: `python -m pytest tests/test_design_task_brief.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add lib/design_runs/__init__.py lib/design_runs/brief.py tests/test_design_task_brief.py
git commit -m "feat: add design task brief schema"
```

### Task 3: Add project artifact materialization helpers

**Files:**
- Create: `lib/design_runs/candidates.py`
- Create: `lib/design_runs/materialize.py`
- Test: `tests/test_design_project_materialization.py`

- [ ] **Step 1: Write failing tests for artifact paths and outputs**

```python
from lib.design_runs.brief import DesignTaskBrief
from lib.design_runs.candidates import CandidateRecord
from lib.design_runs.materialize import materialize_design_run


def test_materialize_design_run_writes_standard_artifacts(tmp_path: Path) -> None:
    project_root = tmp_path / "projects" / "demo"
    project_root.mkdir(parents=True)
    brief = DesignTaskBrief(...)
    candidates = [
        CandidateRecord(candidate_id="cand-1", rank=1, summary="top candidate", scores={"overall": 0.81})
    ]

    outputs = materialize_design_run(project_root, brief, {"profile": "yeast_expression"}, candidates)

    assert (project_root / "design" / "design_brief.json").exists()
    assert (project_root / "design" / "profile_decision.json").exists()
    assert (project_root / "design" / "candidate_table.json").exists()
    assert (project_root / "logs" / "run_manifest.json").exists()
    assert outputs["candidate_json"].endswith("candidate_table.json")
```

- [ ] **Step 2: Run the new tests**

Run: `python -m pytest tests/test_design_project_materialization.py -q`
Expected: FAIL because the candidate/materialization modules do not exist.

- [ ] **Step 3: Implement minimal ledger and materializer**

Use a small, explicit candidate record:

```python
@dataclass
class CandidateRecord:
    candidate_id: str
    rank: int
    summary: str
    scores: dict[str, float]
    host: str | None = None
    construct_draft: dict[str, object] = field(default_factory=dict)
    evidence: list[str] = field(default_factory=list)
    risks: list[str] = field(default_factory=list)
    uncertainty: list[str] = field(default_factory=list)
```

Write artifacts into:

- `project_root / "design" / "design_brief.json"`
- `project_root / "design" / "profile_decision.json"`
- `project_root / "design" / "candidate_table.json"`
- `project_root / "design" / "candidate_table.csv"`
- `project_root / "logs" / "run_manifest.json"`
- `project_root / "data" / "results" / "evaluation_report.md"`

- [ ] **Step 4: Run the artifact tests**

Run: `python -m pytest tests/test_design_project_materialization.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add lib/design_runs/candidates.py lib/design_runs/materialize.py tests/test_design_project_materialization.py
git commit -m "feat: add design artifact materialization"
```

## Chunk 2: Profile Contract And Reference Profile

### Task 4: Add profile contract and registry

**Files:**
- Create: `lib/design_runs/profiles/__init__.py`
- Create: `lib/design_runs/profiles/base.py`
- Create: `lib/design_runs/profiles/registry.py`
- Test: `tests/test_design_profile_registry.py`

- [ ] **Step 1: Write the failing registry tests**

```python
from lib.design_runs.profiles.base import DesignProfile
from lib.design_runs.profiles.registry import DesignProfileRegistry


class FakeProfile(DesignProfile):
    profile_name = "fake"
    supported_hosts = ["fake_host"]


def test_registry_registers_and_returns_profiles() -> None:
    registry = DesignProfileRegistry()
    registry.register(FakeProfile())
    assert registry.get("fake").profile_name == "fake"


def test_registry_rejects_duplicate_profile_names() -> None:
    registry = DesignProfileRegistry()
    registry.register(FakeProfile())
    with pytest.raises(ValueError):
        registry.register(FakeProfile())
```

- [ ] **Step 2: Run the registry tests**

Run: `python -m pytest tests/test_design_profile_registry.py -q`
Expected: FAIL because the contract and registry do not exist.

- [ ] **Step 3: Implement the base contract**

The contract should force each profile to define:

```python
class DesignProfile(Protocol):
    profile_name: str
    supported_hosts: list[str]

    def supports(self, brief: DesignTaskBrief) -> bool: ...
    def plan(self, brief: DesignTaskBrief) -> dict[str, object]: ...
    def generate_candidates(self, brief: DesignTaskBrief) -> list[CandidateRecord]: ...
    def evaluate_candidates(self, brief: DesignTaskBrief, candidates: list[CandidateRecord]) -> list[CandidateRecord]: ...
```

- [ ] **Step 4: Implement a minimal registry**

Support:

- registration
- lookup by name
- selection by `supports()`
- optional ranking hook for future host-choice logic

- [ ] **Step 5: Run the registry tests**

Run: `python -m pytest tests/test_design_profile_registry.py -q`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add lib/design_runs/profiles/__init__.py lib/design_runs/profiles/base.py lib/design_runs/profiles/registry.py tests/test_design_profile_registry.py
git commit -m "feat: add design profile contract and registry"
```

### Task 5: Fix species handling in design scripts before profile integration

**Files:**
- Modify: `tools/scripts/dna_design.py`
- Modify: `tools/scripts/mrna_optimize.py`
- Modify: `tests/test_design_tool_script_compat.py`
- Create: `tests/test_design_species_profiles.py`

- [ ] **Step 1: Write focused failing tests for non-human/non-`E.coli` support**

```python
from tools.scripts.dna_design import DNADesigner
from tools.scripts.mrna_optimize import MRNAOptimizer


def test_dna_designer_supports_yeast_without_type_error() -> None:
    designer = DNADesigner(target_species="yeast")
    assert designer.protein_to_dna("MKT").startswith("ATG")


def test_mrna_optimizer_uses_species_specific_defaults_for_yeast() -> None:
    optimizer = MRNAOptimizer(species="yeast")
    result = optimizer.optimize_for_expression("MKT", expression_level="medium")
    assert result["5_utr"] != "GCCACC"
```

- [ ] **Step 2: Run the script tests in the `bio` environment**

Run: `/home/vimalinx/miniforge3/envs/bio/bin/python -m pytest tests/test_design_tool_script_compat.py tests/test_design_species_profiles.py -q`
Expected: FAIL because `DNADesigner("yeast")` still breaks and mRNA defaults are still human-centered.

- [ ] **Step 3: Replace implicit fallback behavior with explicit species tables**

In `dna_design.py`, define:

- `ECOLI_PREF`
- `YEAST_PREF`
- `HUMAN_PREF`

and set `self.codon_table` from an explicit `CODON_PREFERENCE_TABLES` mapping rather than falling back to `REVERSE_CODE`.

In `mrna_optimize.py`, move UTR/codon/GC defaults into a species-profile mapping:

```python
SPECIES_DEFAULTS = {
    "human": {...},
    "yeast": {...},
}
```

If a requested species is not fully modeled, return a structured warning in the analysis payload instead of silently using human defaults.

- [ ] **Step 4: Run the focused tests**

Run: `/home/vimalinx/miniforge3/envs/bio/bin/python -m pytest tests/test_design_tool_script_compat.py tests/test_design_species_profiles.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add tools/scripts/dna_design.py tools/scripts/mrna_optimize.py tests/test_design_tool_script_compat.py tests/test_design_species_profiles.py
git commit -m "fix: harden species-aware design scripts"
```

### Task 6: Implement one reference profile as a contract test

**Files:**
- Create: `lib/design_runs/profiles/yeast_expression.py`
- Test: `tests/test_design_reference_profile.py`

- [ ] **Step 1: Write failing tests for the reference profile**

```python
from lib.design_runs.brief import DesignTaskBrief
from lib.design_runs.profiles.yeast_expression import YeastExpressionProfile


def test_reference_profile_generates_multiple_candidates() -> None:
    brief = DesignTaskBrief(
        goal="design a stable secreted construct",
        task_type="design",
        molecule_type="protein",
        candidate_hosts=["S_cerevisiae", "K_phaffii"],
        output_granularity="expression_construct",
        network_policy="allow_online",
        success_criteria=["multiple_ranked_candidates"],
        risk_tolerance="medium",
    )
    profile = YeastExpressionProfile()
    candidates = profile.generate_candidates(brief)
    assert len(candidates) >= 3


def test_reference_profile_evaluates_and_ranks_candidates() -> None:
    ...
```

- [ ] **Step 2: Run the reference-profile tests**

Run: `/home/vimalinx/miniforge3/envs/bio/bin/python -m pytest tests/test_design_reference_profile.py -q`
Expected: FAIL because the reference profile does not exist.

- [ ] **Step 3: Implement a deterministic contract-test profile**

Requirements:

- Keep this profile clearly labeled as a reference implementation
- Use the hardened design scripts
- Generate multiple candidates deterministically from a small scaffold/variation strategy
- Produce host choice reasoning and construct draft fields
- Emit proxy-only evaluation notes

Candidate generation does **not** need to solve arbitrary de novo biology in phase 1. It only needs to prove that the backbone can carry a multi-candidate design run.

- [ ] **Step 4: Run the reference-profile tests**

Run: `/home/vimalinx/miniforge3/envs/bio/bin/python -m pytest tests/test_design_reference_profile.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add lib/design_runs/profiles/yeast_expression.py tests/test_design_reference_profile.py
git commit -m "feat: add reference yeast expression profile"
```

## Chunk 3: Orchestrator, CLI, MCP, And Docs

### Task 7: Add the design-run orchestrator

**Files:**
- Create: `lib/design_runs/orchestrator.py`
- Test: `tests/test_design_orchestrator.py`

- [ ] **Step 1: Write the failing orchestrator tests**

```python
from lib.design_runs.orchestrator import run_design_goal


def test_run_design_goal_creates_project_and_returns_paths(tmp_path: Path) -> None:
    result = run_design_goal(
        workspace_root=tmp_path,
        goal="design a stable expression construct",
    )
    assert result.project_root.exists()
    assert (result.project_root / "design" / "design_brief.json").exists()
    assert len(result.candidates) >= 1
```

- [ ] **Step 2: Run the orchestrator tests**

Run: `/home/vimalinx/miniforge3/envs/bio/bin/python -m pytest tests/test_design_orchestrator.py -q`
Expected: FAIL because the orchestrator does not exist.

- [ ] **Step 3: Implement the main orchestration flow**

`run_design_goal()` should:

1. Build a `DesignTaskBrief` from the prompt
2. Select a reference profile from the registry
3. Create a `design` project skeleton
4. Generate and evaluate candidates
5. Materialize artifacts
6. Return a structured result object with paths and summary metadata

Keep the orchestrator free of transport concerns. No CLI parsing or MCP JSON handling in this module.

- [ ] **Step 4: Run the orchestrator tests**

Run: `/home/vimalinx/miniforge3/envs/bio/bin/python -m pytest tests/test_design_orchestrator.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add lib/design_runs/orchestrator.py tests/test_design_orchestrator.py
git commit -m "feat: add design run orchestrator"
```

### Task 8: Add a workspace CLI entrypoint for the gold path

**Files:**
- Modify: `scripts/project.py`
- Modify: `tests/test_workspace_project_cli.py`
- Create: `tests/test_workspace_design_cli.py`

- [ ] **Step 1: Write CLI tests for the new subcommand**

```python
def test_project_cli_design_creates_design_project(tmp_path: Path) -> None:
    result = subprocess.run(
        [sys.executable, str(CLI), "design", "design a stable construct"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=True,
    )

    assert "design_brief.json" in result.stdout
    assert any((tmp_path / "projects").iterdir())
```

- [ ] **Step 2: Run the CLI tests**

Run: `/home/vimalinx/miniforge3/envs/bio/bin/python -m pytest tests/test_workspace_project_cli.py tests/test_workspace_design_cli.py -q`
Expected: FAIL because the `design` subcommand does not exist.

- [ ] **Step 3: Add the new CLI subcommand**

Suggested command:

```bash
python scripts/project.py design "design a stable secreted construct"
```

Support optional flags:

- `--project-name`
- `--profile`
- `--network-policy`

Wire the subcommand to `lib.design_runs.orchestrator.run_design_goal`.

- [ ] **Step 4: Run the CLI tests**

Run: `/home/vimalinx/miniforge3/envs/bio/bin/python -m pytest tests/test_workspace_project_cli.py tests/test_workspace_design_cli.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add scripts/project.py tests/test_workspace_project_cli.py tests/test_workspace_design_cli.py
git commit -m "feat: add design gold path CLI"
```

### Task 9: Expose the design gold path through `bio-lab-mcp`

**Files:**
- Modify: `mcp-servers/bio-lab-mcp/lab_server.py`
- Modify: `tests/test_lab_mcp_server.py`

- [ ] **Step 1: Write MCP tests for the new design-run operation**

```python
def test_run_design_project_invokes_workspace_cli(monkeypatch) -> None:
    module = _load_module("bio_lab_server_design_test")
    ...
    result = module.run_design_project("design a stable construct")
    assert result["returncode"] == 0
```

- [ ] **Step 2: Run the MCP tests**

Run: `/home/vimalinx/miniforge3/envs/bio/bin/python -m pytest tests/test_lab_mcp_server.py -q`
Expected: FAIL because the design-run tool does not exist.

- [ ] **Step 3: Add the wrapper function and tool definition**

Add:

- `run_design_project(goal: str, project_name: str | None = None, profile: str | None = None)`
- a new MCP tool entry such as `run_design_project`

The MCP layer should remain a thin wrapper around the workspace CLI, matching the current `bio-lab-mcp` style.

- [ ] **Step 4: Run the MCP tests**

Run: `/home/vimalinx/miniforge3/envs/bio/bin/python -m pytest tests/test_lab_mcp_server.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add mcp-servers/bio-lab-mcp/lab_server.py tests/test_lab_mcp_server.py
git commit -m "feat: expose design gold path via lab MCP"
```

### Task 10: Update docs and stable test coverage

**Files:**
- Modify: `scripts/ci/run_stable_tests.sh`
- Modify: `README.md`
- Modify: `docs/README.md`
- Modify: `mcp-servers/bio-lab-mcp/README.md`
- Modify: `docs/CHANGELOG.md`

- [ ] **Step 1: Add the new tests to the stable suite**

Run: edit `scripts/ci/run_stable_tests.sh` to include:

- `tests/test_design_task_brief.py`
- `tests/test_design_profile_registry.py`
- `tests/test_design_project_materialization.py`
- `tests/test_design_reference_profile.py`
- `tests/test_design_orchestrator.py`
- `tests/test_workspace_design_cli.py`

- [ ] **Step 2: Update user-facing docs**

Document:

- the new `design` project type
- the new CLI command
- the new MCP tool
- the phase-1 limitation that this is a reference gold path with proxy evaluation, not wet-lab truth

- [ ] **Step 3: Run the full stable suite**

Run: `bash scripts/ci/run_stable_tests.sh`
Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add scripts/ci/run_stable_tests.sh README.md docs/README.md mcp-servers/bio-lab-mcp/README.md docs/CHANGELOG.md
git commit -m "docs: document platform spine phase 1"
```

## Implementation Notes

- Keep phase-1 parsing deterministic and schema-first. Avoid adding direct model-calling code inside the core package.
- Keep the reference profile visibly labeled as a contract test path.
- Prefer explicit structured warnings over silent fallbacks in species-sensitive code.
- Reuse `resolve_workspace_python()` and `build_subprocess_env()` whenever subprocess behavior is added.
- Do not move existing project template behavior around unnecessarily. Additive changes are preferred.

## Manual Review Checklist

- [ ] The plan keeps the platform general and the reference profile narrow
- [ ] The plan adds one clear CLI gold path rather than many partial entrypoints
- [ ] The plan preserves existing MCP point tools while adding the platform-level wrapper separately
- [ ] The plan writes artifacts into `projects/` deterministically
- [ ] The plan does not claim biological truth where only proxy evaluation exists

## Self-Review Notes

Because the current session policy does not permit spawning review subagents unless the user explicitly asks for delegated/sub-agent work, perform local review for this phase:

- Re-read the spec before execution
- Validate all planned file paths exist or are intentionally new
- Confirm every new task has at least one matching test target
- Confirm the reference profile is described as a contract test, not platform identity

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-03-26-platform-spine-phase-1.md`. Ready to execute?
