"""
LiteFold 工作区桥接层：发现 vendored 源码、给出推荐入口、可选探测健康状态。
"""

from __future__ import annotations

import ast
import json
import os
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib import error, request


EXPECTED_ENDPOINTS = ["/health", "/predict", "/status/{job_id}"]
DEFAULT_LITEFOLD_URL = "http://localhost:7114"


@dataclass(frozen=True)
class LiteFoldLayout:
    workspace_root: Path
    mount_dir: Path
    source_dir: Path
    source_readme: Path
    litefold_dir: Path
    selfhosted_dir: Path
    selfhosted_script: Path
    selfhosted_dockerfile: Path
    selfhosted_requirements: Path
    managed_requirements: Path
    selfhosted_deployments_dir: Path


def resolve_litefold_layout(workspace_root: Path | str | None = None) -> LiteFoldLayout:
    root = Path(workspace_root or Path.cwd()).resolve()
    mount_dir = root / "repositories" / "active" / "litefold"
    source_dir = mount_dir / "source"
    litefold_dir = source_dir / "litefold"
    selfhosted_dir = litefold_dir / "selfhosted"
    return LiteFoldLayout(
        workspace_root=root,
        mount_dir=mount_dir,
        source_dir=source_dir,
        source_readme=source_dir / "README.md",
        litefold_dir=litefold_dir,
        selfhosted_dir=selfhosted_dir,
        selfhosted_script=selfhosted_dir / "selfhosted.py",
        selfhosted_dockerfile=selfhosted_dir / "Dockerfile",
        selfhosted_requirements=selfhosted_dir / "requirements.txt",
        managed_requirements=litefold_dir / "managed" / "requirements.txt",
        selfhosted_deployments_dir=litefold_dir / "selfhosted_deployments",
    )


def build_pythonpath_entries(
    workspace_root: Path | str | LiteFoldLayout | None = None,
) -> list[Path]:
    layout = workspace_root if isinstance(workspace_root, LiteFoldLayout) else resolve_litefold_layout(workspace_root)
    return [layout.source_dir, layout.litefold_dir]


def build_launch_context(
    workspace_root: Path | str | LiteFoldLayout | None = None,
    python_executable: str = "python",
) -> dict[str, Any]:
    layout = workspace_root if isinstance(workspace_root, LiteFoldLayout) else resolve_litefold_layout(workspace_root)
    return {
        "cwd": str(layout.litefold_dir),
        "script_path": str(layout.selfhosted_script),
        "python_executable": python_executable,
        "pythonpath_entries": [str(path) for path in build_pythonpath_entries(layout)],
    }


def build_launch_env(
    workspace_root: Path | str | LiteFoldLayout | None = None,
    base_env: dict[str, str] | None = None,
) -> dict[str, str]:
    env = dict(base_env) if base_env is not None else os.environ.copy()
    entries = [str(path) for path in build_pythonpath_entries(workspace_root)]
    existing = env.get("PYTHONPATH", "").strip()
    if existing:
        env["PYTHONPATH"] = os.pathsep.join([*entries, existing])
    else:
        env["PYTHONPATH"] = os.pathsep.join(entries)
    return env


def build_start_plan(
    workspace_root: Path | str | LiteFoldLayout | None = None,
    python_executable: str = "python",
) -> dict[str, Any]:
    layout = workspace_root if isinstance(workspace_root, LiteFoldLayout) else resolve_litefold_layout(workspace_root)
    launch_context = build_launch_context(layout, python_executable=python_executable)
    return {
        "launch_context": launch_context,
        "command": [python_executable, launch_context["script_path"]],
    }


def _normalize_health_url(base_url: str) -> str:
    cleaned = base_url.rstrip("/")
    if cleaned.endswith("/health"):
        return cleaned
    return f"{cleaned}/health"


def _decode_payload(raw_bytes: bytes) -> dict[str, Any] | list[Any] | str | None:
    raw_text = raw_bytes.decode("utf-8", errors="replace").strip()
    if not raw_text:
        return None

    try:
        return json.loads(raw_text)
    except json.JSONDecodeError:
        return raw_text


def probe_health(base_url: str, timeout: float = 3.0) -> dict[str, Any]:
    health_url = _normalize_health_url(base_url)
    req = request.Request(health_url, headers={"Accept": "application/json"})

    try:
        with request.urlopen(req, timeout=timeout) as response:
            return {
                "ok": True,
                "url": health_url,
                "status_code": getattr(response, "status", 200),
                "payload": _decode_payload(response.read()),
            }
    except error.HTTPError as exc:
        return {
            "ok": False,
            "url": health_url,
            "status_code": exc.code,
            "payload": _decode_payload(exc.read()),
            "error": str(exc),
        }
    except Exception as exc:
        return {
            "ok": False,
            "url": health_url,
            "error": str(exc),
        }


def build_selfhosted_command(
    workspace_root: Path | str | LiteFoldLayout | None = None,
    python_executable: str = "python",
) -> str:
    layout = workspace_root if isinstance(workspace_root, LiteFoldLayout) else resolve_litefold_layout(workspace_root)

    if not layout.selfhosted_script.exists():
        raise FileNotFoundError(f"找不到 LiteFold self-hosted 入口: {layout.selfhosted_script}")

    launch_context = build_launch_context(layout, python_executable=python_executable)
    pythonpath_prefix = ":".join(launch_context["pythonpath_entries"])
    return (
        f"cd {shlex.quote(launch_context['cwd'])} && "
        f"PYTHONPATH={shlex.quote(pythonpath_prefix)}${{PYTHONPATH:+:$PYTHONPATH}} "
        f"{shlex.quote(python_executable)} {shlex.quote(launch_context['script_path'])}"
    )


def _read_requirement_file(path: Path) -> list[str]:
    if not path.exists():
        return []

    return [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def collect_requirement_candidates(
    workspace_root: Path | str | LiteFoldLayout | None = None,
) -> list[dict[str, Any]]:
    layout = workspace_root if isinstance(workspace_root, LiteFoldLayout) else resolve_litefold_layout(workspace_root)
    candidate_paths = [
        layout.selfhosted_requirements,
        layout.managed_requirements,
        *sorted(layout.selfhosted_deployments_dir.glob("requirements*.txt")),
    ]

    candidates: list[dict[str, Any]] = []
    seen: set[Path] = set()
    for path in candidate_paths:
        if path in seen:
            continue
        seen.add(path)
        packages = _read_requirement_file(path)
        if not packages:
            continue
        candidates.append(
            {
                "path": str(path),
                "packages": packages,
            }
        )
    return candidates


def _iter_selfhosted_python_files(layout: LiteFoldLayout) -> list[Path]:
    return sorted(layout.selfhosted_dir.glob("*.py"))


def _local_module_roots(layout: LiteFoldLayout) -> set[str]:
    roots = {"litefold"}
    if layout.litefold_dir.is_dir():
        for child in layout.litefold_dir.iterdir():
            if child.is_dir():
                roots.add(child.name)
            elif child.suffix == ".py":
                roots.add(child.stem)
    return roots


def collect_external_python_modules(
    workspace_root: Path | str | LiteFoldLayout | None = None,
) -> list[str]:
    layout = workspace_root if isinstance(workspace_root, LiteFoldLayout) else resolve_litefold_layout(workspace_root)
    external_roots: set[str] = set()
    stdlib_modules = set(getattr(sys, "stdlib_module_names", set()))
    local_roots = _local_module_roots(layout)

    for path in _iter_selfhosted_python_files(layout):
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except SyntaxError:
            continue

        for node in ast.walk(tree):
            root_name: str | None = None
            if isinstance(node, ast.Import):
                for alias in node.names:
                    root_name = alias.name.split(".", 1)[0]
                    if root_name not in local_roots and root_name not in stdlib_modules:
                        external_roots.add(root_name)
            elif isinstance(node, ast.ImportFrom):
                if node.level != 0 or node.module is None:
                    continue
                root_name = node.module.split(".", 1)[0]
                if root_name not in local_roots and root_name not in stdlib_modules:
                    external_roots.add(root_name)

    return sorted(external_roots)


def probe_python_modules(
    python_executable: str,
    module_names: list[str],
) -> dict[str, Any]:
    resolved = shutil.which(python_executable)
    if resolved is None and not Path(python_executable).expanduser().exists():
        return {
            "python_executable": python_executable,
            "checked_modules": module_names,
            "missing_modules": module_names,
            "ok": False,
            "error": f"找不到 Python 可执行文件: {python_executable}",
        }

    script = (
        "import importlib.util, json, sys\n"
        "mods = sys.argv[1:]\n"
        "missing = [name for name in mods if importlib.util.find_spec(name) is None]\n"
        "print(json.dumps({'checked_modules': mods, 'missing_modules': missing, 'ok': not missing}))\n"
    )
    completed = subprocess.run(
        [python_executable, "-c", script, *module_names],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        return {
            "python_executable": python_executable,
            "checked_modules": module_names,
            "missing_modules": module_names,
            "ok": False,
            "error": completed.stderr.strip() or "模块探测失败",
        }

    payload = json.loads(completed.stdout)
    payload["python_executable"] = python_executable
    return payload


def collect_preflight(
    workspace_root: Path | str | None = None,
    python_executable: str = "python",
) -> dict[str, Any]:
    layout = resolve_litefold_layout(workspace_root)
    external_modules = collect_external_python_modules(layout)
    module_probe = probe_python_modules(python_executable, external_modules)
    requirements_candidates = collect_requirement_candidates(layout)
    start_plan = build_start_plan(layout, python_executable=python_executable)

    warnings = []
    if not requirements_candidates:
        warnings.append("没有发现可直接复用的 LiteFold requirements 候选文件。")
    if module_probe.get("missing_modules"):
        warnings.append("当前 Python 环境缺少 LiteFold self-hosted 所需模块。")

    return {
        "status": collect_status(workspace_root=layout.workspace_root),
        "launch_context": start_plan["launch_context"],
        "command": start_plan["command"],
        "external_python_modules": external_modules,
        "module_probe": module_probe,
        "requirements_candidates": requirements_candidates,
        "warnings": warnings,
        "ready_to_launch": bool(
            layout.selfhosted_script.exists() and not module_probe.get("missing_modules")
        ),
    }


def collect_status(
    workspace_root: Path | str | None = None,
    base_url: str | None = None,
) -> dict[str, Any]:
    layout = resolve_litefold_layout(workspace_root)
    docker_available = shutil.which("docker") is not None
    docker_build_ready = (
        layout.selfhosted_dockerfile.exists()
        and layout.selfhosted_requirements.exists()
    )
    selfhosted_ready = layout.selfhosted_script.exists()

    warnings: list[str] = []
    if layout.selfhosted_dockerfile.exists() and not layout.selfhosted_requirements.exists():
        warnings.append(
            "LiteFold self-hosted Dockerfile 已存在，但同目录下缺少 requirements.txt，当前不建议直接按 Docker 方式构建。"
        )

    recommended_launch_mode: str | None = None
    recommended_launch_command: str | None = None
    if selfhosted_ready:
        recommended_launch_mode = "python"
        recommended_launch_command = build_selfhosted_command(layout)

    status: dict[str, Any] = {
        "workspace_root": str(layout.workspace_root),
        "mount_dir": str(layout.mount_dir),
        "source_dir": str(layout.source_dir),
        "selfhosted_dir": str(layout.selfhosted_dir),
        "litefold_dir": str(layout.litefold_dir),
        "source_dir_exists": layout.source_dir.exists(),
        "source_readme_exists": layout.source_readme.exists(),
        "selfhosted_script_exists": layout.selfhosted_script.exists(),
        "selfhosted_dockerfile_exists": layout.selfhosted_dockerfile.exists(),
        "selfhosted_requirements_exists": layout.selfhosted_requirements.exists(),
        "docker_available": docker_available,
        "docker_build_ready": docker_build_ready,
        "recommended_launch_mode": recommended_launch_mode,
        "recommended_launch_command": recommended_launch_command,
        "default_url": DEFAULT_LITEFOLD_URL,
        "expected_endpoints": EXPECTED_ENDPOINTS,
        "warnings": warnings,
        "health": None,
    }

    if base_url:
        status["health"] = probe_health(base_url)

    return status
