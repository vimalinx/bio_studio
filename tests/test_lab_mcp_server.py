from __future__ import annotations

import asyncio
import importlib.util
import json
import subprocess
import sys
import types
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SERVER_PATH = ROOT / "mcp-servers" / "bio-lab-mcp" / "lab_server.py"


def _install_fake_mcp_modules() -> None:
    mcp_module = types.ModuleType("mcp")
    mcp_server_module = types.ModuleType("mcp.server")
    mcp_stdio_module = types.ModuleType("mcp.server.stdio")
    mcp_types_module = types.ModuleType("mcp.types")

    class FakeServer:
        def __init__(self, name: str) -> None:
            self.name = name

        def list_tools(self):
            def decorator(func):
                return func

            return decorator

        def call_tool(self):
            def decorator(func):
                return func

            return decorator

        def create_initialization_options(self):
            return {}

        async def run(self, *args, **kwargs):
            return None

    async def fake_stdio_server():
        raise RuntimeError("not used in tests")

    class Tool:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

    class TextContent:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

    mcp_server_module.Server = FakeServer
    mcp_stdio_module.stdio_server = fake_stdio_server
    mcp_types_module.Tool = Tool
    mcp_types_module.TextContent = TextContent

    sys.modules["mcp"] = mcp_module
    sys.modules["mcp.server"] = mcp_server_module
    sys.modules["mcp.server.stdio"] = mcp_stdio_module
    sys.modules["mcp.types"] = mcp_types_module


def _load_module(name: str):
    _install_fake_mcp_modules()
    spec = importlib.util.spec_from_file_location(name, SERVER_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_list_tools_exposes_lab_workspace_operations() -> None:
    module = _load_module("bio_lab_server_tools_test")

    tools = asyncio.run(module.list_tools())
    names = [tool.kwargs["name"] for tool in tools]

    assert "list_workspace_projects" in names
    assert "get_project_summary" in names
    assert "get_project_steps" in names
    assert "run_project_validation" in names
    assert "run_workspace_validation" in names
    assert "get_server_capabilities" in names


def test_list_workspace_projects_discovers_projects(tmp_path: Path, monkeypatch) -> None:
    module = _load_module("bio_lab_server_discovery_test")
    projects_dir = tmp_path / "projects"
    alpha_scripts = projects_dir / "alpha" / "scripts"
    alpha_scripts.mkdir(parents=True)
    (alpha_scripts / "validate_project.py").write_text("print('ok')\n", encoding="utf-8")
    (alpha_scripts / "pipeline.py").write_text("print('ok')\n", encoding="utf-8")
    (alpha_scripts / "config.py").write_text("PROJECT_NAME = 'alpha'\n", encoding="utf-8")

    beta_scripts = projects_dir / "beta" / "scripts"
    beta_scripts.mkdir(parents=True)
    (beta_scripts / "pipeline.py").write_text("print('ok')\n", encoding="utf-8")

    monkeypatch.setattr(module, "PROJECTS_DIR", projects_dir)

    projects = module.list_workspace_projects()

    assert [item["name"] for item in projects] == ["alpha", "beta"]
    assert projects[0]["has_validate"] is True
    assert projects[1]["has_validate"] is False


def test_get_project_summary_reads_config_and_validation_report(
    tmp_path: Path, monkeypatch
) -> None:
    module = _load_module("bio_lab_server_summary_test")
    projects_dir = tmp_path / "projects"
    project_root = projects_dir / "alpha"
    scripts_dir = project_root / "scripts"
    logs_dir = project_root / "logs"
    scripts_dir.mkdir(parents=True)
    logs_dir.mkdir(parents=True)

    (scripts_dir / "config.py").write_text(
        "\n".join(
            [
                "from pathlib import Path",
                "PROJECT_NAME = 'alpha'",
                "PROJECT_TYPE = 'rnaseq'",
                "PROJECT_DESCRIPTION = 'alpha project'",
                "PROJECT_ROOT = Path(__file__).resolve().parent.parent",
                "DATA_DIR = PROJECT_ROOT / 'data'",
                "LOGS_DIR = PROJECT_ROOT / 'logs'",
                "SAMPLES = ['s1', 's2']",
                "THREADS = 8",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (scripts_dir / "validate_project.py").write_text("print('ok')\n", encoding="utf-8")
    (scripts_dir / "pipeline.py").write_text("print('ok')\n", encoding="utf-8")
    (logs_dir / "validation_report.json").write_text(
        json.dumps({"summary": "PASS"}) + "\n", encoding="utf-8"
    )

    monkeypatch.setattr(module, "PROJECTS_DIR", projects_dir)

    summary = module.get_project_summary("alpha")

    assert summary["project_name"] == "alpha"
    assert summary["project_type"] == "rnaseq"
    assert summary["samples_count"] == 2
    assert summary["required_tools"] == ["fastp", "STAR", "featureCounts"]
    assert summary["latest_validation_summary"] == "PASS"


def test_run_project_validation_invokes_workspace_cli(monkeypatch) -> None:
    module = _load_module("bio_lab_server_validation_test")
    calls: list[dict[str, object]] = []

    def fake_run(cmd, cwd, env, capture_output, text):
        calls.append(
            {
                "cmd": cmd,
                "cwd": cwd,
                "env": env,
                "capture_output": capture_output,
                "text": text,
            }
        )
        return subprocess.CompletedProcess(cmd, 0, stdout="Validation summary: PASS\n", stderr="")

    monkeypatch.setattr(module.subprocess, "run", fake_run)
    monkeypatch.setattr(module, "resolve_workspace_python", lambda: Path("/tmp/bio/bin/python"))
    monkeypatch.setattr(module, "build_subprocess_env", lambda extra_env=None: {"BIO_STUDIO": "1"})

    result = module.run_project_validation("ai_design_playground")

    assert result["returncode"] == 0
    assert "Validation summary: PASS" in result["stdout"]
    assert calls[0]["cmd"] == [
        "/tmp/bio/bin/python",
        str(ROOT / "scripts" / "project.py"),
        "validate",
        "ai_design_playground",
    ]


def test_run_workspace_validation_invokes_workspace_cli(monkeypatch) -> None:
    module = _load_module("bio_lab_server_workspace_validation_test")
    calls: list[list[str]] = []

    def fake_run(cmd, cwd, env, capture_output, text):
        calls.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, stdout="ALL TESTS PASSED\n", stderr="")

    monkeypatch.setattr(module.subprocess, "run", fake_run)
    monkeypatch.setattr(module, "resolve_workspace_python", lambda: Path("/tmp/bio/bin/python"))
    monkeypatch.setattr(module, "build_subprocess_env", lambda extra_env=None: {"BIO_STUDIO": "1"})

    result = module.run_workspace_validation()

    assert result["returncode"] == 0
    assert "ALL TESTS PASSED" in result["stdout"]
    assert calls[0] == [
        "/tmp/bio/bin/python",
        str(ROOT / "scripts" / "project.py"),
        "workspace-validate",
    ]


def test_main_wrapper_uses_asyncio_run(monkeypatch) -> None:
    module = _load_module("bio_lab_server_main_test")
    calls: list[object] = []

    monkeypatch.setattr(module.asyncio, "run", lambda coro: calls.append(coro))

    module.main()

    assert len(calls) == 1
    calls[0].close()
