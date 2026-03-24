from __future__ import annotations

import importlib.util
import os
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
PROJECT_CLI_PATH = ROOT / "scripts" / "project.py"
WORKSPACE_VALIDATE_PATH = (
    ROOT / "projects" / "test_env_validation" / "scripts" / "run_validation.py"
)


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_run_python_entrypoint_prefers_detected_bio_python(
    monkeypatch, tmp_path: Path
) -> None:
    project_cli = _load_module(PROJECT_CLI_PATH, "bio_studio_project_cli_test")

    fake_home = tmp_path / "home"
    fake_bio_python = fake_home / "miniforge3" / "envs" / "bio" / "bin" / "python"
    fake_bio_python.parent.mkdir(parents=True, exist_ok=True)
    fake_bio_python.write_text("#!/bin/sh\n", encoding="utf-8")
    fake_bio_python.chmod(0o755)

    script_path = tmp_path / "demo.py"
    script_path.write_text("print('demo')\n", encoding="utf-8")

    calls: list[tuple[list[str], Path | str | None]] = []

    class _Completed:
        returncode = 0

    def fake_run(args, cwd=None, **kwargs):
        calls.append((args, cwd))
        return _Completed()

    monkeypatch.setenv("HOME", str(fake_home))
    monkeypatch.setattr(project_cli.sys, "executable", "/usr/bin/python")
    monkeypatch.setattr(project_cli.subprocess, "run", fake_run)

    project_cli.run_python_entrypoint(script_path, ["--steps"])

    assert calls
    assert calls[0][0][0] == str(fake_bio_python)
    assert calls[0][0][1:] == [script_path.name, "--steps"]


def test_workspace_validation_treats_scanpy_as_optional(
    monkeypatch, capsys
) -> None:
    run_validation = _load_module(
        WORKSPACE_VALIDATE_PATH, "bio_studio_workspace_validation_test"
    )

    fake_dataframe = SimpleNamespace(shape=(1, 2))
    fake_pandas = SimpleNamespace(
        __version__="2.3.3",
        read_csv=lambda *args, **kwargs: fake_dataframe,
    )
    fake_numpy = SimpleNamespace(__version__="2.2.6")
    fake_scipy = SimpleNamespace(__version__="1.15.2")
    fake_torch = SimpleNamespace(
        __version__="2.10.0",
        cuda=SimpleNamespace(is_available=lambda: False),
    )
    fake_bio = SimpleNamespace(__version__="1.86")

    def fake_import_module(name: str):
        if name == "pandas":
            return fake_pandas
        if name == "numpy":
            return fake_numpy
        if name == "scipy":
            return fake_scipy
        if name == "torch":
            return fake_torch
        if name == "Bio":
            return fake_bio
        if name == "scanpy":
            raise ModuleNotFoundError("No module named 'scanpy'")
        raise AssertionError(f"unexpected import: {name}")

    monkeypatch.setattr(run_validation.importlib, "import_module", fake_import_module)

    run_validation.test_python_libs()

    output = capsys.readouterr().out.lower()
    assert "scanpy" in output
    assert "optional" in output or "可选" in output


def test_run_python_entrypoint_prepends_detected_bio_bin_to_path(
    monkeypatch, tmp_path: Path
) -> None:
    project_cli = _load_module(PROJECT_CLI_PATH, "bio_studio_project_cli_env_test")

    fake_home = tmp_path / "home"
    fake_bio_python = fake_home / "miniforge3" / "envs" / "bio" / "bin" / "python"
    fake_bio_python.parent.mkdir(parents=True, exist_ok=True)
    fake_bio_python.write_text("#!/bin/sh\n", encoding="utf-8")
    fake_bio_python.chmod(0o755)

    script_path = tmp_path / "demo.py"
    script_path.write_text("print('demo')\n", encoding="utf-8")

    captured_env: dict[str, str] = {}

    class _Completed:
        returncode = 0

    def fake_run(args, cwd=None, env=None, **kwargs):
        if env is not None:
            captured_env.update(env)
        return _Completed()

    monkeypatch.setenv("HOME", str(fake_home))
    monkeypatch.setattr(project_cli.sys, "executable", "/usr/bin/python")
    monkeypatch.setattr(project_cli.subprocess, "run", fake_run)

    project_cli.run_python_entrypoint(script_path)

    path_entries = captured_env["PATH"].split(os.pathsep)
    assert path_entries[0] == str(fake_bio_python.parent)
