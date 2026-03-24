#!/usr/bin/env python3
"""
test_rnaseq_analysis 项目级自检
"""

from __future__ import annotations

import importlib.util
import json
import sys
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
CONFIG_PATH = SCRIPT_DIR / "config.py"


class ValidationContext:
    def __init__(self) -> None:
        self.results: list[dict[str, str]] = []

    def add(self, name: str, status: str, message: str) -> None:
        self.results.append({"name": name, "status": status, "message": message})
        print(f"[{status}] {name}: {message}")

    @property
    def failed(self) -> bool:
        return any(item["status"] == "FAIL" for item in self.results)


def load_config(context: ValidationContext):
    if not CONFIG_PATH.exists():
        context.add("config_file", "FAIL", f"missing {CONFIG_PATH.name}")
        return None

    spec = importlib.util.spec_from_file_location("special_project_config", CONFIG_PATH)
    if spec is None or spec.loader is None:
        context.add("config_import", "FAIL", "unable to create config import spec")
        return None

    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except Exception as exc:
        context.add("config_import", "FAIL", f"failed to import config: {exc}")
        return None

    context.add("config_import", "PASS", f"imported {CONFIG_PATH.name} successfully")
    return module


def validate_structure(context: ValidationContext) -> None:
    expected = [
        PROJECT_ROOT / "data" / "raw",
        PROJECT_ROOT / "data" / "processed",
        PROJECT_ROOT / "data" / "results",
        PROJECT_ROOT / "data" / "references",
        PROJECT_ROOT / "scripts",
    ]
    missing = [str(path.relative_to(PROJECT_ROOT)) for path in expected if not path.exists()]
    if missing:
        context.add("directory_structure", "FAIL", "missing: " + ", ".join(missing))
    else:
        context.add("directory_structure", "PASS", "specialist project directories exist")


def validate_project_assets(context: ValidationContext) -> None:
    required_scripts = [
        SCRIPT_DIR / "pipeline.py",
        SCRIPT_DIR / "evo2_gp_saturation.py",
        SCRIPT_DIR / "evo2_variant_analysis.py",
    ]
    missing = [path.name for path in required_scripts if not path.exists()]
    if missing:
        context.add("specialist_scripts", "FAIL", "missing: " + ", ".join(missing))
    else:
        context.add("specialist_scripts", "PASS", "core specialist analysis scripts exist")

    archived_reports = list((PROJECT_ROOT / "archive").glob("*.md"))
    if archived_reports:
        context.add("archived_reports", "PASS", f"found {len(archived_reports)} archived reports")
    else:
        context.add("archived_reports", "WARN", "no archived reports found")


def validate_data_presence(context: ValidationContext, config) -> None:
    if config is None:
        return

    raw_dir = PROJECT_ROOT / getattr(config, "RAW_DIR", "data/raw")
    processed_dir = PROJECT_ROOT / getattr(config, "PROCESSED_DIR", "data/processed")
    references_dir = PROJECT_ROOT / getattr(config, "REFERENCES_DIR", "data/references")

    raw_fastqs = list(raw_dir.glob("*.fastq.gz"))
    if raw_fastqs:
        context.add("raw_fastq", "PASS", f"found {len(raw_fastqs)} FASTQ files")
    else:
        context.add("raw_fastq", "WARN", "no FASTQ files found")

    reference_files = list(references_dir.glob("*"))
    if reference_files:
        context.add("references", "PASS", f"found {len(reference_files)} reference files")
    else:
        context.add("references", "WARN", "no reference files found")

    bam_files = list(processed_dir.glob("*.bam")) + list((processed_dir / "aligned").glob("*.bam"))
    if bam_files:
        context.add("processed_bams", "PASS", f"found {len(bam_files)} BAM files")
    else:
        context.add("processed_bams", "WARN", "no BAM outputs found yet")


def write_reports(context: ValidationContext) -> tuple[Path, Path]:
    logs_dir = PROJECT_ROOT / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    payload = {
        "project_root": str(PROJECT_ROOT),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "summary": "PASS" if not context.failed else "FAIL",
        "results": context.results,
    }

    json_path = logs_dir / "validation_report.json"
    md_path = logs_dir / "validation_report.md"
    json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md_lines = [
        f"# Validation Report: {payload['summary']}",
        "",
        f"- Project root: `{PROJECT_ROOT}`",
        f"- Generated at: `{payload['generated_at']}`",
        "",
        "## Checks",
        "",
    ]
    for item in context.results:
        md_lines.append(f"- **{item['status']}** `{item['name']}`: {item['message']}")
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    return json_path, md_path


def run_validation() -> int:
    context = ValidationContext()
    validate_structure(context)
    config = load_config(context)
    validate_project_assets(context)
    validate_data_presence(context, config)
    json_path, md_path = write_reports(context)
    summary = "PASS" if not context.failed else "FAIL"
    print(f"Validation summary: {summary}")
    print(f"Validation reports written to: {json_path} and {md_path}")
    return 0 if not context.failed else 1


if __name__ == "__main__":
    sys.exit(run_validation())
