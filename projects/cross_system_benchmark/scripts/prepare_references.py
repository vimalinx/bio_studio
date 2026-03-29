#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import common


def build_reference_manifest() -> dict[str, object]:
    return {
        "references": {
            dataset: payload["reference"]
            for dataset, payload in common.CURATED_DATASETS.items()
        },
        "host_references": common.HOST_REFERENCES,
    }


def write_reference_manifest(project_root: Path) -> Path:
    project_root = Path(project_root).resolve()
    manifest_path = project_root / "data" / "references" / "reference_manifest.json"
    return common.write_json(manifest_path, build_reference_manifest())


def download_default_references(project_root: Path, force: bool = False) -> list[Path]:
    project_root = Path(project_root).resolve()
    references_root = common.ensure_dir(project_root / "data" / "references")
    downloaded: list[Path] = []
    for dataset, payload in common.CURATED_DATASETS.items():
        reference = payload["reference"]
        if not reference["download_by_default"] or not reference["filename"]:
            continue
        destination = references_root / dataset / reference["filename"]
        downloaded.append(common.download_file(reference["url"], destination, force=force))
    return downloaded


def prepare_references(project_root: Path, force: bool = False) -> Path:
    manifest_path = write_reference_manifest(project_root)
    download_default_references(project_root, force=force)
    return manifest_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare reference manifest and small default references")
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--force", action="store_true", help="re-download default references")
    args = parser.parse_args()
    manifest_path = prepare_references(args.project_root, force=args.force)
    print(manifest_path)


if __name__ == "__main__":
    main()
