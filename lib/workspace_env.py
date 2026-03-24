"""
工作区环境解析：统一定位 bio 环境的 python/bin，并为非交互执行补 PATH。
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Mapping


_COMMON_ENV_ROOTS = (
    "miniforge3",
    "mambaforge",
    "micromamba",
    "miniconda3",
    "anaconda3",
)


def _dedupe_paths(paths: list[Path]) -> list[Path]:
    deduped: list[Path] = []
    seen: set[Path] = set()
    for path in paths:
        if path in seen:
            continue
        seen.add(path)
        deduped.append(path)
    return deduped


def is_bio_python(python_path: Path) -> bool:
    parts = python_path.expanduser().resolve().parts
    return "envs" in parts and "bio" in parts


def _candidate_pythons() -> list[Path]:
    current_python = Path(sys.executable).expanduser().resolve()
    candidates: list[Path] = []

    override = os.environ.get("BIO_STUDIO_PYTHON", "").strip()
    if override:
        candidates.append(Path(override).expanduser())

    if is_bio_python(current_python):
        candidates.append(current_python)

    home = Path.home()
    for env_root in _COMMON_ENV_ROOTS:
        candidates.append(home / env_root / "envs" / "bio" / "bin" / "python")

    candidates.append(Path("/opt/conda/envs/bio/bin/python"))
    candidates.append(current_python)
    return _dedupe_paths(candidates)


def resolve_workspace_python() -> Path:
    for candidate in _candidate_pythons():
        if candidate.exists():
            return candidate
    return Path(sys.executable).expanduser().resolve()


def _candidate_bins() -> list[Path]:
    candidates: list[Path] = []

    override_bin = os.environ.get("BIO_STUDIO_CONDA_BIN", "").strip()
    if override_bin:
        candidates.append(Path(override_bin).expanduser())

    override_python = os.environ.get("BIO_STUDIO_PYTHON", "").strip()
    if override_python:
        candidates.append(Path(override_python).expanduser().parent)

    current_python = Path(sys.executable).expanduser().resolve()
    if is_bio_python(current_python):
        candidates.append(current_python.parent)

    candidates.append(resolve_workspace_python().parent)

    home = Path.home()
    for env_root in _COMMON_ENV_ROOTS:
        candidates.append(home / env_root / "envs" / "bio" / "bin")

    candidates.append(Path("/opt/conda/envs/bio/bin"))
    return _dedupe_paths(candidates)


def resolve_workspace_bin() -> Path | None:
    for candidate in _candidate_bins():
        if candidate.is_dir():
            return candidate
    return None


def _prepend_path_entry(env: dict[str, str], entry: str) -> None:
    entries = [item for item in env.get("PATH", "").split(os.pathsep) if item]
    if entries and entries[0] == entry:
        return
    entries = [entry, *[item for item in entries if item != entry]]
    env["PATH"] = os.pathsep.join(entries)


def build_subprocess_env(
    base_env: Mapping[str, str] | None = None,
    extra_env: Mapping[str, str] | None = None,
) -> dict[str, str]:
    env = dict(base_env) if base_env is not None else os.environ.copy()

    bin_dir = resolve_workspace_bin()
    if bin_dir is not None:
        _prepend_path_entry(env, str(bin_dir))
        env.setdefault("BIO_STUDIO_CONDA_BIN", str(bin_dir))

    env.setdefault("BIO_STUDIO_PYTHON", str(resolve_workspace_python()))

    if extra_env:
        for key, value in extra_env.items():
            env[str(key)] = str(value)

    return env


def ensure_workspace_path() -> Path | None:
    bin_dir = resolve_workspace_bin()
    if bin_dir is not None:
        _prepend_path_entry(os.environ, str(bin_dir))
        os.environ.setdefault("BIO_STUDIO_CONDA_BIN", str(bin_dir))

    os.environ.setdefault("BIO_STUDIO_PYTHON", str(resolve_workspace_python()))
    return bin_dir
