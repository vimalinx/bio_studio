#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import common


def _write_normalized_outputs(
    dataset: str,
    expression_df: pd.DataFrame,
    samples_df: pd.DataFrame,
    output_dir: Path,
) -> dict[str, Path]:
    output_dir = common.ensure_dir(output_dir)
    expression_path = output_dir / "expression.tsv"
    samples_path = output_dir / "samples.tsv"
    summary_path = output_dir / "dataset_summary.json"

    expression_df.to_csv(expression_path, sep="\t", index=False)
    samples_df.to_csv(samples_path, sep="\t", index=False)

    top_feature = (
        expression_df.groupby("feature", as_index=False)["value"]
        .mean()
        .sort_values("value", ascending=False)
        .iloc[0]["feature"]
    )
    summary = {
        "dataset": dataset,
        "sample_count": int(samples_df["sample"].nunique()),
        "feature_count": int(expression_df["feature"].nunique()),
        "top_feature": str(top_feature),
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    return {
        "expression": expression_path,
        "samples": samples_path,
        "summary": summary_path,
    }


def _ebola_sample_name(sheet_name: str) -> str:
    if "D1" in sheet_name:
        return "EBOV_D1"
    if "D2" in sheet_name:
        return "EBOV_D2"
    if "D3" in sheet_name:
        return "EBOV_D3"
    return sheet_name.replace(" ", "_")


def normalize_ebola_workbook(workbook_path: Path, output_dir: Path) -> dict[str, Path]:
    workbook = pd.ExcelFile(workbook_path)
    records: list[dict[str, object]] = []
    sample_rows: list[dict[str, str]] = []
    for sheet_name in workbook.sheet_names:
        frame = pd.read_excel(workbook_path, sheet_name=sheet_name)
        sample = _ebola_sample_name(sheet_name)
        sample_rows.append({"sample": sample, "condition": sample.lower(), "dataset": "ebola"})
        for row in frame[["Name", "TPM"]].itertuples(index=False):
            records.append(
                {
                    "dataset": "ebola",
                    "sample": sample,
                    "feature": row.Name,
                    "value": float(row.TPM),
                }
            )
    expression_df = pd.DataFrame(records)
    samples_df = pd.DataFrame(sample_rows).drop_duplicates()
    return _write_normalized_outputs("ebola", expression_df, samples_df, output_dir)


def normalize_sars_counts(counts_path: Path, output_dir: Path, sample_prefix: str = "Series1_NHBE") -> dict[str, Path]:
    frame = pd.read_csv(counts_path, sep="\t")
    feature_column = frame.columns[0]
    selected_columns = [column for column in frame.columns if column.startswith(sample_prefix)]
    subset = frame[[feature_column, *selected_columns]].rename(columns={feature_column: "feature"})
    expression_df = subset.melt(id_vars="feature", var_name="sample", value_name="value")
    expression_df["dataset"] = "sars_cov_2"
    expression_df = expression_df[["dataset", "sample", "feature", "value"]]

    sample_rows = []
    for sample in selected_columns:
        condition = "infected" if "SARS-CoV-2" in sample else "mock"
        sample_rows.append({"sample": sample, "condition": condition, "dataset": "sars_cov_2"})
    samples_df = pd.DataFrame(sample_rows)
    return _write_normalized_outputs("sars_cov_2", expression_df, samples_df, output_dir)


def normalize_yeast_workbook(workbook_path: Path, output_dir: Path) -> dict[str, Path]:
    frame = pd.read_excel(workbook_path, sheet_name=0)
    selected_columns = [column for column in ["WT_Log", "WT_Q"] if column in frame.columns]
    if not selected_columns:
        selected_columns = [column for column in frame.columns if str(column).startswith("WT_")]
    subset = frame[["Gene", *selected_columns]].rename(columns={"Gene": "feature"})
    expression_df = subset.melt(id_vars="feature", var_name="sample", value_name="value")
    expression_df["dataset"] = "yeast"
    expression_df = expression_df[["dataset", "sample", "feature", "value"]]

    sample_rows = []
    for sample in selected_columns:
        if sample == "WT_Log":
            condition = "log_phase"
        elif sample == "WT_Q":
            condition = "quiescence"
        else:
            condition = sample.lower()
        sample_rows.append({"sample": sample, "condition": condition, "dataset": "yeast"})
    samples_df = pd.DataFrame(sample_rows)
    return _write_normalized_outputs("yeast", expression_df, samples_df, output_dir)


def download_processed_inputs(project_root: Path, force: bool = False) -> dict[str, Path]:
    project_root = Path(project_root).resolve()
    downloaded = {}
    for dataset, payload in common.CURATED_DATASETS.items():
        raw_dir = common.ensure_dir(project_root / "data" / "raw" / dataset)
        processed_file = payload["processed_file"]
        destination = raw_dir / processed_file["filename"]
        downloaded[dataset] = common.download_file(processed_file["url"], destination, force=force)
    return downloaded


def run_project_baseline(project_root: Path, force_download: bool = False) -> dict[str, dict[str, Path]]:
    project_root = Path(project_root).resolve()
    downloaded = download_processed_inputs(project_root, force=force_download)
    processed_root = project_root / "data" / "processed"
    return {
        "ebola": normalize_ebola_workbook(downloaded["ebola"], processed_root / "ebola"),
        "sars_cov_2": normalize_sars_counts(downloaded["sars_cov_2"], processed_root / "sars_cov_2"),
        "yeast": normalize_yeast_workbook(downloaded["yeast"], processed_root / "yeast"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run processed-data baselines for curated GEO benchmarks")
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--force-download", action="store_true")
    args = parser.parse_args()
    outputs = run_project_baseline(args.project_root, force_download=args.force_download)
    for dataset, paths in outputs.items():
        print(dataset)
        for name, path in paths.items():
            print(f"  {name}: {path}")


if __name__ == "__main__":
    main()
