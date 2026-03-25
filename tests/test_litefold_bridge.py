from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "scripts" / "litefold.py"


def _create_fake_litefold_workspace(workspace_root: Path) -> Path:
    source_dir = workspace_root / "repositories" / "active" / "litefold" / "source"
    selfhosted_dir = source_dir / "litefold" / "selfhosted"
    selfhosted_dir.mkdir(parents=True)

    (source_dir / "README.md").write_text("# LiteFold\n", encoding="utf-8")
    (selfhosted_dir / "Dockerfile").write_text("FROM scratch\n", encoding="utf-8")
    (selfhosted_dir / "selfhosted.py").write_text("print('litefold')\n", encoding="utf-8")
    return selfhosted_dir


def test_collect_status_reports_vendored_litefold_layout(tmp_path: Path, monkeypatch) -> None:
    _create_fake_litefold_workspace(tmp_path)

    from lib import litefold_bridge

    monkeypatch.setattr(litefold_bridge.shutil, "which", lambda name: "/usr/bin/docker")

    status = litefold_bridge.collect_status(workspace_root=tmp_path)

    assert status["source_dir_exists"] is True
    assert status["source_readme_exists"] is True
    assert status["selfhosted_script_exists"] is True
    assert status["selfhosted_dockerfile_exists"] is True
    assert status["selfhosted_requirements_exists"] is False
    assert status["docker_available"] is True
    assert status["recommended_launch_mode"] == "python"
    assert status["recommended_launch_command"].endswith("python selfhosted.py")
    assert status["expected_endpoints"] == ["/health", "/predict", "/status/{job_id}"]


def test_collect_status_can_include_health_probe(tmp_path: Path, monkeypatch) -> None:
    _create_fake_litefold_workspace(tmp_path)

    from lib import litefold_bridge

    monkeypatch.setattr(
        litefold_bridge,
        "probe_health",
        lambda url, timeout=3.0: {
            "ok": True,
            "url": f"{url.rstrip('/')}/health",
            "status_code": 200,
            "payload": {"status": "healthy"},
        },
    )

    status = litefold_bridge.collect_status(
        workspace_root=tmp_path,
        base_url="http://localhost:7114",
    )

    assert status["health"]["ok"] is True
    assert status["health"]["payload"]["status"] == "healthy"


def test_probe_health_appends_health_path(monkeypatch) -> None:
    from lib import litefold_bridge

    captured: dict[str, str] = {}

    class FakeResponse:
        status = 200

        def read(self) -> bytes:
            return b'{"status":"healthy"}'

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb) -> None:
            return None

    def fake_urlopen(request, timeout=3.0):
        captured["url"] = request.full_url
        return FakeResponse()

    monkeypatch.setattr(litefold_bridge.request, "urlopen", fake_urlopen)

    result = litefold_bridge.probe_health("http://localhost:7114")

    assert captured["url"] == "http://localhost:7114/health"
    assert result["ok"] is True
    assert result["payload"]["status"] == "healthy"


def test_litefold_cli_status_json_reports_bridge_status(tmp_path: Path) -> None:
    _create_fake_litefold_workspace(tmp_path)

    result = subprocess.run(
        [
            sys.executable,
            str(CLI),
            "status",
            "--workspace-root",
            str(tmp_path),
            "--json",
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )

    payload = json.loads(result.stdout)
    assert payload["source_dir_exists"] is True
    assert payload["recommended_launch_mode"] == "python"
    assert "selfhosted.py" in payload["recommended_launch_command"]
