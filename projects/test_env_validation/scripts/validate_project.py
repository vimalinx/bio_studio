#!/usr/bin/env python3
"""
test_env_validation 项目级自检
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


def load_validation_helpers():
    for candidate_root in SCRIPT_DIR.parents:
        helper_path = candidate_root / "lib" / "project_validation.py"
        if not helper_path.exists():
            continue
        spec = importlib.util.spec_from_file_location("bio_project_validation", helper_path)
        if spec is None or spec.loader is None:
            continue
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    return None


VALIDATION_HELPERS = load_validation_helpers()


class ValidationContext:
    def __init__(self) -> None:
        self.results = []
        self.config = None

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

    spec = importlib.util.spec_from_file_location("project_config", CONFIG_PATH)
    if spec is None or spec.loader is None:
        context.add("config_import", "FAIL", "unable to create config import spec")
        return None

    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except Exception as exc:  # pragma: no cover - defensive path for broken configs
        context.add("config_import", "FAIL", f"failed to import config: {exc}")
        return None

    context.config = module
    context.add("config_import", "PASS", f"imported {CONFIG_PATH.name} successfully")
    return module


def validate_directories(context: ValidationContext) -> None:
    expected_dirs = [
        PROJECT_ROOT / "data",
        PROJECT_ROOT / "data" / "raw",
        PROJECT_ROOT / "data" / "processed",
        PROJECT_ROOT / "data" / "results",
        PROJECT_ROOT / "data" / "references",
        PROJECT_ROOT / "scripts",
        PROJECT_ROOT / "notebooks",
        PROJECT_ROOT / "logs",
    ]

    missing = [str(path.relative_to(PROJECT_ROOT)) for path in expected_dirs if not path.exists()]
    if missing:
        context.add("directory_structure", "FAIL", "missing: " + ", ".join(missing))
    else:
        context.add("directory_structure", "PASS", "all expected directories exist")


def validate_config_fields(context: ValidationContext, config) -> None:
    required_fields = [
        "PROJECT_NAME",
        "PROJECT_TYPE",
        "PROJECT_ROOT",
        "DATA_DIR",
        "RAW_DIR",
        "PROCESSED_DIR",
        "RESULTS_DIR",
        "REFERENCES_DIR",
        "LOGS_DIR",
        "SCRIPTS_DIR",
        "NOTEBOOKS_DIR",
        "SAMPLES",
        "THREADS",
    ]

    missing = [name for name in required_fields if not hasattr(config, name)]
    if missing:
        context.add("config_fields", "FAIL", "missing fields: " + ", ".join(missing))
    else:
        context.add("config_fields", "PASS", "required config fields are present")


def validate_paths(context: ValidationContext, config) -> None:
    path_fields = [
        "PROJECT_ROOT",
        "DATA_DIR",
        "RAW_DIR",
        "PROCESSED_DIR",
        "RESULTS_DIR",
        "REFERENCES_DIR",
        "LOGS_DIR",
        "SCRIPTS_DIR",
        "NOTEBOOKS_DIR",
    ]

    invalid = []
    actual_project_root = PROJECT_ROOT.resolve()
    configured_project_root = Path(getattr(config, "PROJECT_ROOT", PROJECT_ROOT)).resolve()

    if configured_project_root != actual_project_root:
        invalid.append(
            f"PROJECT_ROOT={configured_project_root} (expected {actual_project_root})"
        )

    for field in path_fields:
        value = getattr(config, field, None)
        if value is None:
            continue
        path = Path(value).resolve()
        if path != actual_project_root and actual_project_root not in path.parents:
            invalid.append(f"{field}={path}")

    if invalid:
        context.add("path_scope", "FAIL", "paths outside project root: " + "; ".join(invalid))
    else:
        context.add("path_scope", "PASS", "configured paths stay within project root")


def validate_optional_inputs(context: ValidationContext, config) -> None:
    raw_files = list(Path(config.RAW_DIR).glob("*"))
    if raw_files:
        context.add("raw_inputs", "PASS", f"found {len(raw_files)} files in raw data directory")
    else:
        context.add("raw_inputs", "WARN", "no raw input files yet; template is ready for data ingestion")

    reference = getattr(config, "REFERENCE_GENOME", None)
    if reference:
        reference_path = Path(reference)
        if not reference_path.is_absolute():
            reference_path = Path(config.PROJECT_ROOT) / reference_path
        if reference_path.exists():
            context.add("reference_genome", "PASS", f"reference file exists: {reference_path}")
        else:
            context.add("reference_genome", "WARN", f"reference file not found yet: {reference_path}")
    else:
        context.add("reference_genome", "WARN", "REFERENCE_GENOME is not configured yet")


def write_reports(context: ValidationContext) -> tuple[Path, Path]:
    logs_dir = PROJECT_ROOT / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    summary = "PASS" if not context.failed else "FAIL"
    payload = {
        "project_root": str(PROJECT_ROOT),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
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
    validate_directories(context)
    config = load_config(context)
    if config is not None:
        validate_config_fields(context, config)
        validate_paths(context, config)
        validate_optional_inputs(context, config)
        if VALIDATION_HELPERS is None:
            context.add("required_tools", "WARN", "workspace validation helper unavailable; skipping CLI tool check")
        else:
            VALIDATION_HELPERS.validate_cli_tools(context, config)

    json_path, md_path = write_reports(context)
    summary = "PASS" if not context.failed else "FAIL"
    print(f"Validation summary: {summary}")
    print(f"Validation reports written to: {json_path} and {md_path}")
    return 0 if not context.failed else 1


if __name__ == "__main__":
    sys.exit(run_validation())
