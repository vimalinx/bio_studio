"""
LiteFold 工作区桥接层：发现 vendored 源码、给出推荐入口、可选探测健康状态。
"""

from __future__ import annotations

import json
import shlex
import shutil
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
    selfhosted_dir: Path
    selfhosted_script: Path
    selfhosted_dockerfile: Path
    selfhosted_requirements: Path


def resolve_litefold_layout(workspace_root: Path | str | None = None) -> LiteFoldLayout:
    root = Path(workspace_root or Path.cwd()).resolve()
    mount_dir = root / "repositories" / "active" / "litefold"
    source_dir = mount_dir / "source"
    selfhosted_dir = source_dir / "litefold" / "selfhosted"
    return LiteFoldLayout(
        workspace_root=root,
        mount_dir=mount_dir,
        source_dir=source_dir,
        source_readme=source_dir / "README.md",
        selfhosted_dir=selfhosted_dir,
        selfhosted_script=selfhosted_dir / "selfhosted.py",
        selfhosted_dockerfile=selfhosted_dir / "Dockerfile",
        selfhosted_requirements=selfhosted_dir / "requirements.txt",
    )


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

    return (
        f"cd {shlex.quote(str(layout.selfhosted_dir))} "
        f"&& {python_executable} selfhosted.py"
    )


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
