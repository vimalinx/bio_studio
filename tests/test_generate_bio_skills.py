from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "skills" / "generate_bio_skills.py"


def _load_module(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_generate_bio_skills_creates_skill_directory(tmp_path: Path) -> None:
    module = _load_module(SCRIPT_PATH, "generate_bio_skills")
    tool = {
        "name": "samtools",
        "summary": "Utilities for SAM/BAM/CRAM files.",
        "help_text": "samtools view\nsamtools sort\nsamtools index",
        "path": "/usr/bin/samtools",
        "sources": ["path"],
    }

    written = module.generate_skill_for_tool(tool, tmp_path, render_mode="offline")

    assert (tmp_path / "samtools" / "SKILL.md").exists()
    assert (tmp_path / "samtools" / "references" / "help.md").exists()
    assert written["skill_dir"] == tmp_path / "samtools"


def test_generate_skills_can_continue_on_error(tmp_path: Path, monkeypatch) -> None:
    module = _load_module(SCRIPT_PATH, "generate_bio_skills")

    monkeypatch.setattr(
        module.discover_bio_tools,
        "discover_tools",
        lambda repo_root, source="combined": [
            {"name": "samtools", "command": "samtools", "sources": ["path"]},
            {"name": "bedtools", "command": "bedtools", "sources": ["path"]},
        ],
    )
    monkeypatch.setattr(module, "_should_skip_existing", lambda skill_dir: False)
    monkeypatch.setattr(
        module,
        "capture_command_reference",
        lambda tool, repo_root: {"version_text": "", "help_text": "", "man_text": ""},
    )

    def fake_generate_skill_for_tool(tool, output_root, render_mode="auto"):
        if tool["name"] == "samtools":
            raise ValueError("intentional failure")
        skill_dir = output_root / tool["name"]
        skill_dir.mkdir(parents=True, exist_ok=True)
        return {
            "skill_dir": skill_dir,
            "skill_path": skill_dir / "SKILL.md",
            "help_path": skill_dir / "references" / "help.md",
        }

    monkeypatch.setattr(module, "generate_skill_for_tool", fake_generate_skill_for_tool)

    failures: list[dict[str, str]] = []
    generated = module.generate_skills(
        repo_root=tmp_path,
        output_root=tmp_path / ".claude" / "skills",
        render_mode="online",
        continue_on_error=True,
        failures=failures,
    )

    assert [item["name"] for item in generated] == ["bedtools"]
    assert failures == [
        {
            "name": "samtools",
            "sources": ["path"],
            "error": "intentional failure",
        }
    ]
