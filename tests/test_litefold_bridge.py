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
    managed_dir = source_dir / "litefold" / "managed"
    managed_dir.mkdir(parents=True)
    fold_models_dir = source_dir / "litefold" / "fold_models" / "esmfold"
    fold_models_dir.mkdir(parents=True)
    setup_dir = workspace_root / "scripts" / "setup"
    setup_dir.mkdir(parents=True)

    (source_dir / "README.md").write_text("# LiteFold\n", encoding="utf-8")
    (selfhosted_dir / "Dockerfile").write_text("FROM scratch\n", encoding="utf-8")
    (selfhosted_dir / "__init__.py").write_text("", encoding="utf-8")
    (selfhosted_dir / "selfhosted.py").write_text(
        "\n".join(
            [
                "from fastapi import FastAPI",
                "from sqlalchemy.orm import Session",
                "import torch",
                "from litefold.selfhosted.models import Job",
                "from fold_models import ESMFold",
                "from selfhosted.constants import CUDA_DEVICE",
                "from Bio.PDB import PDBParser",
                "import biotite.structure.io as bsio",
                "",
                "app = FastAPI()",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (selfhosted_dir / "models.py").write_text(
        "\n".join(
            [
                "from sqlalchemy import create_engine",
                "from selfhosted.constants import SQLALCHEMY_DATABASE_URL",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (selfhosted_dir / "schemas.py").write_text(
        "\n".join(
            [
                "from pydantic import BaseModel",
                "import numpy as np",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (selfhosted_dir / "constants.py").write_text(
        "CUDA_DEVICE = 'cuda:0'\nSQLALCHEMY_DATABASE_URL = 'sqlite:///db/jobs.db'\n",
        encoding="utf-8",
    )
    (source_dir / "litefold" / "fold_models" / "__init__.py").write_text("", encoding="utf-8")
    (fold_models_dir / "main.py").write_text(
        "\n".join(
            [
                "import torch",
                "import esm",
                "from scipy.stats import truncnorm",
                "from einops import rearrange",
                "import tree",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (managed_dir / "requirements.txt").write_text(
        "fastapi==0.104.1\nsqlalchemy==2.0.23\npydantic>=2.4.2\n",
        encoding="utf-8",
    )
    (workspace_root / "requirements-litefold-selfhosted.txt").write_text(
        "fastapi==0.104.1\nuvicorn==0.24.0\nbiopython==1.79\nbiotite>=0.38.0\n",
        encoding="utf-8",
    )
    (setup_dir / "setup_litefold_env.sh").write_text(
        "#!/bin/bash\nset -e\n",
        encoding="utf-8",
    )
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
    assert "PYTHONPATH=" in status["recommended_launch_command"]
    assert status["recommended_launch_command"].endswith("selfhosted.py")
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


def test_collect_preflight_reports_dependency_candidates(tmp_path: Path, monkeypatch) -> None:
    _create_fake_litefold_workspace(tmp_path)

    from lib import litefold_bridge

    monkeypatch.setattr(
        litefold_bridge,
        "probe_python_modules",
        lambda python_executable, module_names: {
            "python_executable": python_executable,
            "checked_modules": module_names,
            "missing_modules": ["torch"],
            "ok": False,
        },
    )

    preflight = litefold_bridge.collect_preflight(
        workspace_root=tmp_path,
        python_executable="python",
    )

    assert preflight["launch_context"]["cwd"].endswith("/source/litefold")
    assert preflight["launch_context"]["script_path"].endswith("/source/litefold/selfhosted/selfhosted.py")
    assert any(path.endswith("/source") for path in preflight["launch_context"]["pythonpath_entries"])
    assert any(path.endswith("/source/litefold") for path in preflight["launch_context"]["pythonpath_entries"])
    assert preflight["external_python_modules"] == [
        "Bio",
        "biotite",
        "einops",
        "esm",
        "fastapi",
        "numpy",
        "pydantic",
        "scipy",
        "sqlalchemy",
        "torch",
        "tree",
    ]
    assert preflight["module_probe"]["missing_modules"] == ["torch"]
    assert preflight["requirements_candidates"][0]["path"].endswith("/managed/requirements.txt")
    assert preflight["requirements_candidates"][0]["packages"] == [
        "fastapi==0.104.1",
        "sqlalchemy==2.0.23",
        "pydantic>=2.4.2",
    ]
    assert preflight["workspace_installation"]["requirements_file"].endswith("/requirements-litefold-selfhosted.txt")
    assert preflight["workspace_installation"]["setup_script"].endswith("/scripts/setup/setup_litefold_env.sh")
    assert preflight["workspace_installation"]["install_command"].startswith("bash ")


def test_litefold_cli_start_selfhosted_dry_run_json_prints_launch_plan(tmp_path: Path) -> None:
    _create_fake_litefold_workspace(tmp_path)

    result = subprocess.run(
        [
            sys.executable,
            str(CLI),
            "start-selfhosted",
            "--workspace-root",
            str(tmp_path),
            "--dry-run",
            "--json",
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )

    payload = json.loads(result.stdout)
    assert payload["dry_run"] is True
    assert payload["launch_context"]["cwd"].endswith("/source/litefold")
    assert payload["command"][-1].endswith("selfhosted.py")
    assert any(path.endswith("/source") for path in payload["launch_context"]["pythonpath_entries"])


def test_setup_litefold_env_script_supports_dry_run_for_target_workspace(tmp_path: Path) -> None:
    _create_fake_litefold_workspace(tmp_path)
    script_path = ROOT / "scripts" / "setup" / "setup_litefold_env.sh"

    result = subprocess.run(
        [
            "bash",
            str(script_path),
            "--workspace-root",
            str(tmp_path),
            "--python-executable",
            sys.executable,
            "--dry-run",
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )

    assert "requirements-litefold-selfhosted.txt" in result.stdout
    assert f"{sys.executable} -m pip install -r" in result.stdout
