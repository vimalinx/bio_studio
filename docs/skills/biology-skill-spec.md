# Biology Skill Factory Spec

This document is the local reference for every new biology or bioinformatics tool skill added to this repository.

## Scope

The factory keeps two layers of skills in the same project tree:

- Existing domain skills such as `bioinformatics-toolkit`, `rnaseq-pipeline`, and `protein-structure`
- Generated single-tool skills for concrete commands such as `samtools`, `fastp`, or `mafft`

The goal is to keep the domain skills intact while adding a repeatable path for `.claude/skills/<tool-name>/` generation.

## Required Layout

Every generated tool skill must use this file layout:

```text
.claude/skills/<tool-name>/
├── SKILL.md
└── references/
    └── help.md
```

- `SKILL.md` is the entrypoint used by Claude Code
- `references/help.md` stores captured `--help`, `-h`, or `man` output plus local command metadata

## Required Frontmatter

Generated biology tool skills default to manual invocation so the workspace does not flood the model context with dozens of auto-triggerable tool descriptions.

Use this frontmatter shape in `SKILL.md`:

```yaml
---
name: samtools
description: Use when you need samtools for BAM or CRAM inspection, conversion, sorting, indexing, or region extraction in this workspace.
disable-model-invocation: true
user-invocable: true
---
```

Required rules:

- `name` must be lowercase letters, numbers, and hyphens only
- `description` must say when to use the tool skill, not summarize the full workflow
- `disable-model-invocation: true` is required for generated single-tool skills
- `user-invocable: true` is required unless there is a clear reason to hide the skill from `/`

## Content Rules

Each generated `SKILL.md` should stay concise and operational:

- State what the tool is for
- Mention the local executable path when available
- Point the model to `references/help.md` for raw help text
- Give a short workflow for checking inputs, running the command, and validating outputs
- Avoid copying long help pages into `SKILL.md`

## Discovery Sources

The first factory version discovers biology tools from three local sources:

1. `scripts/maintenance/install_bio_tools.sh`
2. `.claude/skills/bioinformatics-toolkit/TOOLS.md`
3. Commands found on the workspace `PATH` after `lib.workspace_env.ensure_workspace_path()`

## Rendering Modes

The renderer supports two modes:

- `offline`: deterministic template rendering for tests and safe local fallback
- `online`: direct OpenAI-compatible API rendering driven by local config or environment variables

Tracked files must never contain a real API key. Keep live credentials in environment variables or an untracked local config file.

When `--render-mode online` is used for production regeneration, the factory should call the configured API directly and never silently switch the run into template-only rendering. API retries and reasoning-effort downgrades may happen inside the direct API client when the remote service is flaky, but the rendered result must still come from the API.

## Local Config

Factory scripts may read:

- `scripts/skills/config.example.json` as the tracked template
- `scripts/skills/config.local.json` as the untracked local override

The local config is intentionally ignored by git.

Recommended live config fields:

- `base_url`
- `model`
- `reasoning_effort`
- `max_completion_tokens`
- `request_retries`
- `api_key_env`

## Batch Workflow

For every new biology command, follow the same production line:

1. Discover the command from local sources
2. Capture `--help`, `-h`, or `man`
3. Render `SKILL.md`
4. Write `references/help.md`
5. Refresh `docs/skills/biology-skill-catalog.md`

## Validation

Before a large regeneration, probe the configured service directly:

```bash
BIO_SKILL_FACTORY_API_KEY=... python scripts/skills/validate_skill_factory_api.py
```

Then regenerate with direct API mode:

```bash
BIO_SKILL_FACTORY_API_KEY=... python scripts/skills/generate_bio_skills.py --source linux-all --render-mode online --continue-on-error
```
