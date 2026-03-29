# Biology Skill Factory Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a repeatable skill-generation pipeline that discovers local bioinformatics tools, captures help/man content, and generates project-local `.claude/skills/<tool>/SKILL.md` assets plus a catalog and spec.

**Architecture:** Keep existing domain skills intact and layer a generator system alongside them. Use Python scripts to discover tools and collect structured evidence from local commands, then render skill drafts through a controlled template and optional LLM call. Keep secrets out of git and refresh a human-readable catalog from the discovered inventory.

**Tech Stack:** Python 3, pathlib, json, subprocess, argparse, urllib/http client, markdown text generation, existing `.claude/skills` tree

---

### Task 1: Add The Skill Spec And Catalog Docs

**Files:**
- Create: `docs/skills/biology-skill-spec.md`
- Create: `docs/skills/biology-skill-catalog.md`
- Test: `tests/test_biology_skill_docs.py`

**Step 1: Write the failing test**

```python
def test_biology_skill_docs_exist():
    assert (ROOT / "docs" / "skills" / "biology-skill-spec.md").exists()
    assert (ROOT / "docs" / "skills" / "biology-skill-catalog.md").exists()
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_biology_skill_docs.py -v`
Expected: FAIL because the docs do not exist yet

**Step 3: Write minimal implementation**

- Create the docs directory
- Write the spec doc with required frontmatter rules, naming rules, invocation defaults, and file layout
- Write the catalog doc with current known biology skills and generation workflow

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_biology_skill_docs.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add docs/skills tests/test_biology_skill_docs.py
git commit -m "docs: add biology skill factory spec and catalog"
```

### Task 2: Add Tool Discovery Script

**Files:**
- Create: `scripts/skills/discover_bio_tools.py`
- Create: `tests/test_discover_bio_tools.py`
- Modify: `docs/skills/biology-skill-catalog.md`

**Step 1: Write the failing test**

```python
def test_discovery_returns_tools_from_known_sources():
    tools = discover_tools(...)
    assert "samtools" in tools
    assert "STAR" in tools or "star" in tools
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_discover_bio_tools.py -v`
Expected: FAIL because the script does not exist yet

**Step 3: Write minimal implementation**

- Parse `scripts/maintenance/install_bio_tools.sh`
- Parse documented tool names from `.claude/skills/bioinformatics-toolkit/TOOLS.md`
- Probe `PATH` for matching executables
- Emit structured JSON and markdown-friendly rows

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_discover_bio_tools.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add scripts/skills/discover_bio_tools.py tests/test_discover_bio_tools.py
git commit -m "feat: add bio tool discovery"
```

### Task 3: Add Help Collection And Skill Rendering Primitives

**Files:**
- Create: `scripts/skills/render_skill_from_help.py`
- Create: `tests/test_render_skill_from_help.py`
- Create: `scripts/skills/config.example.json`

**Step 1: Write the failing test**

```python
def test_render_skill_from_help_outputs_valid_skill_markdown():
    text = render_skill(...)
    assert text.startswith("---")
    assert "name:" in text
    assert "description:" in text
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_render_skill_from_help.py -v`
Expected: FAIL because the renderer does not exist yet

**Step 3: Write minimal implementation**

- Add a deterministic prompt template builder
- Add LLM request config loading from env or local config
- Add offline fallback renderer for tests
- Keep secrets out of tracked files

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_render_skill_from_help.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add scripts/skills/render_skill_from_help.py scripts/skills/config.example.json tests/test_render_skill_from_help.py
git commit -m "feat: add skill renderer for bio tools"
```

### Task 4: Add Batch Skill Generation Script

**Files:**
- Create: `scripts/skills/generate_bio_skills.py`
- Create: `tests/test_generate_bio_skills.py`
- Modify: `.claude/skills/`

**Step 1: Write the failing test**

```python
def test_generate_bio_skills_creates_skill_directory(tmp_path):
    generate_skills(..., output_root=tmp_path)
    assert (tmp_path / "samtools" / "SKILL.md").exists()
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_generate_bio_skills.py -v`
Expected: FAIL because the generator does not exist yet

**Step 3: Write minimal implementation**

- Wire discovery + help collection + rendering together
- For each tool create:
  - `.claude/skills/<tool>/SKILL.md`
  - `.claude/skills/<tool>/references/help.md`
- Default generated tool skills to `disable-model-invocation: true`

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_generate_bio_skills.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add scripts/skills/generate_bio_skills.py tests/test_generate_bio_skills.py .claude/skills
git commit -m "feat: add biology skill batch generator"
```

### Task 5: Refresh Catalog From Real Generated Skills

**Files:**
- Modify: `docs/skills/biology-skill-catalog.md`
- Modify: `scripts/skills/generate_bio_skills.py`
- Create: `tests/test_biology_skill_catalog_refresh.py`

**Step 1: Write the failing test**

```python
def test_catalog_refresh_lists_existing_and_generated_skills():
    refresh_catalog(...)
    text = catalog_path.read_text()
    assert "bioinformatics-toolkit" in text
    assert "samtools" in text
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_biology_skill_catalog_refresh.py -v`
Expected: FAIL because catalog refresh is not wired yet

**Step 3: Write minimal implementation**

- List current domain skills
- List generated tool skills
- Include provenance such as `existing-domain` vs `generated-tool`

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_biology_skill_catalog_refresh.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add docs/skills/biology-skill-catalog.md scripts/skills/generate_bio_skills.py tests/test_biology_skill_catalog_refresh.py
git commit -m "feat: refresh biology skill catalog from generated skills"
```

### Task 6: Generate The First Batch Of Tool Skills

**Files:**
- Modify: `.claude/skills/<tool>/SKILL.md` for generated tools
- Modify: `.claude/skills/<tool>/references/help.md` for generated tools
- Test: `tests/test_generate_bio_skills.py`

**Step 1: Run generator against the workspace**

Run:

```bash
python scripts/skills/generate_bio_skills.py --source workspace
```

Expected: tool skill directories are created under `.claude/skills/`

**Step 2: Spot check generated outputs**

Run:

```bash
sed -n '1,80p' .claude/skills/samtools/SKILL.md
sed -n '1,80p' .claude/skills/mafft/SKILL.md
sed -n '1,80p' .claude/skills/fastp/SKILL.md
```

Expected: valid frontmatter + concise task guidance

**Step 3: Adjust generator if outputs are structurally wrong**

- Fix naming
- Fix duplicate descriptions
- Fix malformed frontmatter

**Step 4: Re-run generator and confirm stable output**

Run:

```bash
python scripts/skills/generate_bio_skills.py --source workspace
```

Expected: repeatable output, no uncontrolled duplication

**Step 5: Commit**

```bash
git add .claude/skills docs/skills scripts/skills
git commit -m "feat: generate initial biology tool skills"
```

### Task 7: Verify End-To-End Behavior

**Files:**
- Test: `tests/test_biology_skill_docs.py`
- Test: `tests/test_discover_bio_tools.py`
- Test: `tests/test_render_skill_from_help.py`
- Test: `tests/test_generate_bio_skills.py`
- Test: `tests/test_biology_skill_catalog_refresh.py`

**Step 1: Run focused tests**

Run:

```bash
pytest \
  tests/test_biology_skill_docs.py \
  tests/test_discover_bio_tools.py \
  tests/test_render_skill_from_help.py \
  tests/test_generate_bio_skills.py \
  tests/test_biology_skill_catalog_refresh.py -v
```

Expected: PASS

**Step 2: Run discovery command**

Run:

```bash
python scripts/skills/discover_bio_tools.py --format markdown
```

Expected: stable catalog output

**Step 3: Run generator in dry-run and real modes**

Run:

```bash
python scripts/skills/generate_bio_skills.py --dry-run
python scripts/skills/generate_bio_skills.py
```

Expected: dry-run previews actions, real run writes skills

**Step 4: Confirm docs and generated skills are aligned**

Run:

```bash
rg -n "generated-tool|existing-domain" docs/skills/biology-skill-catalog.md
find .claude/skills -maxdepth 2 -name SKILL.md | sort
```

Expected: catalog matches skill tree

**Step 5: Commit**

```bash
git add docs/skills scripts/skills tests .claude/skills
git commit -m "feat: add biology skill factory"
```
