"""
项目级验证辅助：统一声明模板项目的 CLI 工具需求并做 PATH 检查。
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path
from typing import Iterable

try:
    from .workspace_env import ensure_workspace_path
except ImportError:
    ROOT = Path(__file__).resolve().parent.parent
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    from lib.workspace_env import ensure_workspace_path


PROJECT_TYPE_REQUIRED_TOOLS: dict[str, list[str]] = {
    "rnaseq": ["fastp", "STAR", "featureCounts"],
    "variant": ["bwa", "samtools", "bcftools"],
    "phylogeny": ["mafft", "iqtree2"],
}

PROJECT_NAME_REQUIRED_TOOLS: dict[str, list[str]] = {
    "ai_design_playground": ["seqkit", "prodigal", "RNAfold"],
}

TOOL_CANDIDATES: dict[str, list[str]] = {
    "iqtree2": ["iqtree2", "iqtree"],
}


def _dedupe(values: Iterable[str]) -> list[str]:
    deduped: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = str(value).strip()
        if not text or text in seen:
            continue
        seen.add(text)
        deduped.append(text)
    return deduped


def required_tools_for_config(config) -> list[str]:
    explicit = getattr(config, "REQUIRED_TOOLS", None)
    if explicit:
        return _dedupe(explicit)

    project_name = str(getattr(config, "PROJECT_NAME", "")).strip()
    if project_name in PROJECT_NAME_REQUIRED_TOOLS:
        return list(PROJECT_NAME_REQUIRED_TOOLS[project_name])

    project_type = str(getattr(config, "PROJECT_TYPE", "generic")).strip()
    return list(PROJECT_TYPE_REQUIRED_TOOLS.get(project_type, []))


def _resolve_tool(tool_name: str) -> tuple[str | None, list[str]]:
    candidates = TOOL_CANDIDATES.get(tool_name, [tool_name])
    for candidate in candidates:
        resolved = shutil.which(candidate)
        if resolved:
            return candidate, candidates
    return None, candidates


def validate_cli_tools(context, config) -> None:
    ensure_workspace_path()

    required_tools = required_tools_for_config(config)
    project_type = str(getattr(config, "PROJECT_TYPE", "generic")).strip()
    if not required_tools:
        context.add(
            "required_tools",
            "WARN",
            f"no default CLI tool requirements declared for project type: {project_type}",
        )
        return

    resolved_tools: list[str] = []
    missing_tools: list[str] = []
    for tool_name in required_tools:
        matched, candidates = _resolve_tool(tool_name)
        if matched is None:
            if len(candidates) == 1:
                missing_tools.append(tool_name)
            else:
                missing_tools.append(f"{tool_name} ({'/'.join(candidates)})")
            continue
        resolved_tools.append(f"{tool_name}->{matched}" if matched != tool_name else tool_name)

    if missing_tools:
        context.add(
            "required_tools",
            "WARN",
            "missing on PATH: "
            + ", ".join(missing_tools)
            + "; available: "
            + (", ".join(resolved_tools) if resolved_tools else "none"),
        )
        return

    context.add(
        "required_tools",
        "PASS",
        "available on PATH: " + ", ".join(resolved_tools),
    )
