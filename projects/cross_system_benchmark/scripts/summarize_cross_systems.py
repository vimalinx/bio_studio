#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import common


def write_dataset_summary(dataset_name: str, summary_json_path: Path, output_path: Path) -> Path:
    summary = json.loads(Path(summary_json_path).read_text(encoding="utf-8"))
    lines = [
        f"# {dataset_name} Summary",
        "",
        f"- Dataset: `{summary['dataset']}`",
        f"- Samples: `{summary['sample_count']}`",
        f"- Features: `{summary['feature_count']}`",
        f"- Top feature by mean signal: `{summary['top_feature']}`",
        "",
    ]
    common.ensure_dir(output_path.parent)
    output_path.write_text("\n".join(lines), encoding="utf-8")
    return output_path


def write_integrated_summary(per_dataset_dir: Path, output_path: Path) -> Path:
    summaries = []
    for summary_path in sorted(Path(per_dataset_dir).glob("*_summary.json")):
        summaries.append(json.loads(summary_path.read_text(encoding="utf-8")))

    lines = [
        "# Cross-System Summary",
        "",
        "This report compares dataset-level themes only. It does not claim direct one-to-one gene comparability across systems.",
        "",
    ]
    for summary in summaries:
        lines.append(f"## {summary['dataset']}")
        lines.append("")
        lines.append(f"- Samples: `{summary['sample_count']}`")
        lines.append(f"- Features: `{summary['feature_count']}`")
        lines.append(f"- Top feature: `{summary['top_feature']}`")
        lines.append("")

    common.ensure_dir(output_path.parent)
    output_path.write_text("\n".join(lines), encoding="utf-8")
    return output_path


def run_reporting(project_root: Path) -> Path:
    project_root = Path(project_root).resolve()
    per_dataset_dir = common.ensure_dir(project_root / "data" / "results" / "per_dataset")
    processed_root = project_root / "data" / "processed"

    for dataset in common.CURATED_DATASETS:
        source_summary = processed_root / dataset / "dataset_summary.json"
        target_summary = per_dataset_dir / f"{dataset}_summary.json"
        target_summary.write_text(source_summary.read_text(encoding="utf-8"), encoding="utf-8")
        write_dataset_summary(dataset, target_summary, per_dataset_dir / f"{dataset}_summary.md")

    output_path = project_root / "data" / "results" / "integrated" / "cross_system_summary.md"
    return write_integrated_summary(per_dataset_dir, output_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Write per-dataset and integrated summaries")
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    print(run_reporting(args.project_root))


if __name__ == "__main__":
    main()
