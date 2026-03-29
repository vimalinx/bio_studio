from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "skills" / "discover_bio_tools.py"


def _load_module(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_discover_conda_bio_tools_reads_bioconda_bin_entries(tmp_path: Path) -> None:
    module = _load_module(SCRIPT_PATH, "discover_bio_tools")

    conda_meta = tmp_path / "conda-meta"
    env_bin = tmp_path / "bin"
    conda_meta.mkdir()
    env_bin.mkdir()

    package_meta = {
        "name": "bedtools",
        "channel": "bioconda",
        "files": ["bin/bedtools", "bin/bamToBed", "bin/xfetch.ini", "share/doc/readme.txt"],
    }
    (conda_meta / "bedtools-1.0.0.json").write_text(json.dumps(package_meta), encoding="utf-8")

    for command_name in ["bedtools", "bamToBed", "xfetch.ini"]:
        command_path = env_bin / command_name
        command_path.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        command_path.chmod(0o755)

    discovered = module.discover_conda_bio_tools(conda_meta, env_bin)
    by_name = {item["name"]: item for item in discovered}

    assert "bedtools" in by_name
    assert "bam-to-bed" in by_name
    assert "xfetch-ini" not in by_name
    assert "conda_bioconda" in by_name["bedtools"]["sources"]
    assert "bedtools" in by_name["bedtools"]["package_names"]
