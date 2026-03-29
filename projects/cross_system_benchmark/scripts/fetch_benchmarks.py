#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import common


def build_catalog() -> dict[str, dict[str, object]]:
    return common.CURATED_DATASETS


def write_catalog(project_root: Path) -> dict[str, Path]:
    project_root = Path(project_root).resolve()
    metadata_dir = common.ensure_dir(project_root / "metadata")

    datasets_path = metadata_dir / "datasets.yaml"
    samples_path = metadata_dir / "samples.tsv"
    comparisons_path = metadata_dir / "comparisons.tsv"
    source_links_path = metadata_dir / "source_links.md"

    datasets_path.write_text(
        yaml.safe_dump(build_catalog(), sort_keys=True, allow_unicode=True),
        encoding="utf-8",
    )
    common.write_tsv(samples_path, common.SAMPLE_ROWS)
    common.write_tsv(comparisons_path, common.COMPARISON_ROWS)

    lines = ["# Source Links", ""]
    for dataset, payload in common.CURATED_DATASETS.items():
        lines.append(f"## {dataset}")
        lines.append("")
        lines.append(f"- GEO: {payload['source_url']}")
        lines.append(f"- Processed file: {payload['processed_file']['url']}")
        lines.append("")
    source_links_path.write_text("\n".join(lines), encoding="utf-8")

    return {
        "datasets": datasets_path,
        "samples": samples_path,
        "comparisons": comparisons_path,
        "source_links": source_links_path,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Write cross-system benchmark catalog files")
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    written = write_catalog(args.project_root)
    for name, path in written.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
