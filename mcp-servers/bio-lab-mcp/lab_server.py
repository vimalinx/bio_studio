#!/usr/bin/env python3
"""
Bio Lab MCP Server
为工作区项目发现、项目自检和环境验证提供 MCP 接口。
"""

from __future__ import annotations

import asyncio
import importlib.util
import json
import logging
import subprocess
import sys
from pathlib import Path
from typing import Any

from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import TextContent, Tool


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from lib.project_validation import required_tools_for_config
from lib.workspace_env import build_subprocess_env, resolve_workspace_python


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

server = Server("bio-lab-mcp")

PROJECTS_DIR = ROOT / "projects"
WORKSPACE_CLI = ROOT / "scripts" / "project.py"


def _project_root(project_name: str) -> Path:
    return PROJECTS_DIR / project_name


def _project_scripts_dir(project_name: str) -> Path:
    return _project_root(project_name) / "scripts"


def _project_entry(project_name: str) -> dict[str, Any]:
    scripts_dir = _project_scripts_dir(project_name)
    return {
        "name": project_name,
        "project_root": str(_project_root(project_name)),
        "scripts_dir": str(scripts_dir),
        "has_config": (scripts_dir / "config.py").exists(),
        "has_validate": (scripts_dir / "validate_project.py").exists(),
        "has_pipeline": (scripts_dir / "pipeline.py").exists(),
    }


def list_workspace_projects() -> list[dict[str, Any]]:
    projects: list[dict[str, Any]] = []
    if not PROJECTS_DIR.exists():
        return projects

    for child in sorted(PROJECTS_DIR.iterdir(), key=lambda item: item.name):
        if not child.is_dir():
            continue
        scripts_dir = child / "scripts"
        if not scripts_dir.exists():
            continue
        projects.append(_project_entry(child.name))
    return projects


def _load_project_config(project_name: str):
    config_path = _project_scripts_dir(project_name) / "config.py"
    if not config_path.exists():
        return None

    spec = importlib.util.spec_from_file_location(f"bio_lab_config_{project_name}", config_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to create import spec for {config_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _read_validation_summary(project_name: str) -> tuple[str | None, str | None]:
    report_path = _project_root(project_name) / "logs" / "validation_report.json"
    if not report_path.exists():
        return None, None

    payload = json.loads(report_path.read_text(encoding="utf-8"))
    return payload.get("summary"), str(report_path)


def get_project_summary(project_name: str) -> dict[str, Any]:
    entry = _project_entry(project_name)
    if not Path(entry["project_root"]).exists():
        raise FileNotFoundError(f"找不到项目 {project_name}")

    config = _load_project_config(project_name)
    validation_summary, validation_report = _read_validation_summary(project_name)

    project_root = _project_root(project_name)
    summary = {
        "project_name": project_name,
        "project_root": str(project_root),
        "scripts_dir": entry["scripts_dir"],
        "has_config": entry["has_config"],
        "has_validate": entry["has_validate"],
        "has_pipeline": entry["has_pipeline"],
        "latest_validation_summary": validation_summary,
        "validation_report_path": validation_report,
    }

    if config is None:
        summary.update(
            {
                "project_type": "unknown",
                "project_description": None,
                "samples_count": None,
                "threads": None,
                "required_tools": [],
            }
        )
        return summary

    samples = getattr(config, "SAMPLES", None) or []
    summary.update(
        {
            "project_type": getattr(config, "PROJECT_TYPE", "generic"),
            "project_description": getattr(config, "PROJECT_DESCRIPTION", None),
            "data_dir": str(getattr(config, "DATA_DIR", project_root / "data")),
            "logs_dir": str(getattr(config, "LOGS_DIR", project_root / "logs")),
            "samples_count": len(samples),
            "threads": getattr(config, "THREADS", None),
            "required_tools": required_tools_for_config(config),
        }
    )
    return summary


def _run_workspace_cli(arguments: list[str]) -> dict[str, Any]:
    command = [str(resolve_workspace_python()), str(WORKSPACE_CLI), *arguments]
    completed = subprocess.run(
        command,
        cwd=ROOT,
        env=build_subprocess_env(),
        capture_output=True,
        text=True,
    )
    return {
        "command": command,
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def get_project_steps(project_name: str) -> dict[str, Any]:
    result = _run_workspace_cli(["steps", project_name])
    result["project_name"] = project_name
    result["steps"] = [
        line.strip()
        for line in result["stdout"].splitlines()
        if line.strip()
    ]
    return result


def run_project_validation(project_name: str) -> dict[str, Any]:
    result = _run_workspace_cli(["validate", project_name])
    result["project_name"] = project_name
    return result


def run_workspace_validation() -> dict[str, Any]:
    return _run_workspace_cli(["workspace-validate"])


def get_server_capabilities() -> dict[str, Any]:
    return {
        "server": "bio-lab-mcp",
        "workspace_root": str(ROOT),
        "workspace_cli": str(WORKSPACE_CLI),
        "project_count": len(list_workspace_projects()),
        "tools": [
            "list_workspace_projects",
            "get_project_summary",
            "get_project_steps",
            "run_project_validation",
            "run_workspace_validation",
        ],
    }


@server.list_tools()
async def list_tools() -> list[Tool]:
    return [
        Tool(
            name="list_workspace_projects",
            description="列出当前工作区内可被统一入口识别的项目",
            inputSchema={"type": "object", "properties": {}},
        ),
        Tool(
            name="get_project_summary",
            description="读取项目配置、自检状态和所需工具摘要",
            inputSchema={
                "type": "object",
                "properties": {
                    "project_name": {
                        "type": "string",
                        "description": "项目名称，例如 ai_design_playground",
                    }
                },
                "required": ["project_name"],
            },
        ),
        Tool(
            name="get_project_steps",
            description="调用工作区统一入口，返回项目 pipeline 的步骤说明",
            inputSchema={
                "type": "object",
                "properties": {
                    "project_name": {
                        "type": "string",
                        "description": "项目名称",
                    }
                },
                "required": ["project_name"],
            },
        ),
        Tool(
            name="run_project_validation",
            description="运行项目级 validate_project.py 自检",
            inputSchema={
                "type": "object",
                "properties": {
                    "project_name": {
                        "type": "string",
                        "description": "项目名称",
                    }
                },
                "required": ["project_name"],
            },
        ),
        Tool(
            name="run_workspace_validation",
            description="运行工作区级 workspace-validate smoke test，执行时间可能较长",
            inputSchema={"type": "object", "properties": {}},
        ),
        Tool(
            name="get_server_capabilities",
            description="返回 bio-lab-mcp 的当前能力和工作区入口信息",
            inputSchema={"type": "object", "properties": {}},
        ),
    ]


@server.call_tool()
async def call_tool(name: str, arguments: Any) -> list[TextContent]:
    logger.info("调用工具: %s", name)

    try:
        arguments = arguments or {}

        if name == "list_workspace_projects":
            result = list_workspace_projects()
        elif name == "get_project_summary":
            result = get_project_summary(arguments["project_name"])
        elif name == "get_project_steps":
            result = get_project_steps(arguments["project_name"])
        elif name == "run_project_validation":
            result = run_project_validation(arguments["project_name"])
        elif name == "run_workspace_validation":
            result = run_workspace_validation()
        elif name == "get_server_capabilities":
            result = get_server_capabilities()
        else:
            result = {"error": f"未知工具: {name}"}

        return [
            TextContent(
                type="text",
                text=json.dumps(result, indent=2, ensure_ascii=False),
            )
        ]
    except Exception as exc:
        logger.error("工具执行错误: %s", exc)
        return [
            TextContent(
                type="text",
                text=json.dumps({"error": str(exc)}, indent=2, ensure_ascii=False),
            )
        ]


async def main_async():
    async with stdio_server() as (read_stream, write_stream):
        await server.run(
            read_stream,
            write_stream,
            server.create_initialization_options(),
        )


def main():
    asyncio.run(main_async())


if __name__ == "__main__":
    main()
