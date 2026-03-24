#!/usr/bin/env python3
"""
yeast_genome_learning 项目级自检
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent


class ValidationContext:
    def __init__(self) -> None:
        self.results: list[dict[str, str]] = []

    def add(self, name: str, status: str, message: str) -> None:
        self.results.append({"name": name, "status": status, "message": message})
        print(f"[{status}] {name}: {message}")

    @property
    def failed(self) -> bool:
        return any(item["status"] == "FAIL" for item in self.results)


def validate_structure(context: ValidationContext) -> None:
    expected = [
        PROJECT_ROOT / "scripts",
        PROJECT_ROOT / "data",
        PROJECT_ROOT / "results",
        PROJECT_ROOT / "README.md",
        PROJECT_ROOT / "QUICKSTART.md",
    ]
    missing = [str(path.relative_to(PROJECT_ROOT)) for path in expected if not path.exists()]
    if missing:
        context.add("project_structure", "FAIL", "missing: " + ", ".join(missing))
    else:
        context.add("project_structure", "PASS", "learning project structure exists")


def validate_learning_assets(context: ValidationContext) -> None:
    lesson_scripts = [
        "01_setup_database.sh",
        "02_verify_install.sh",
        "03_extract_gene.sh",
        "04_sequence_analysis.sh",
        "05_blast_search.sh",
        "analyze_genes.sh",
    ]
    missing = [name for name in lesson_scripts if not (SCRIPT_DIR / name).exists()]
    if missing:
        context.add("lesson_scripts", "FAIL", "missing: " + ", ".join(missing))
    else:
        context.add("lesson_scripts", "PASS", "all lesson shell scripts are present")

    website_entry = PROJECT_ROOT / "www" / "index.html"
    if website_entry.exists():
        context.add("learning_site", "PASS", "found learning site entrypoint")
    else:
        context.add("learning_site", "WARN", "learning site entrypoint not found")


def validate_dataset_presence(context: ValidationContext) -> None:
    sequence_dir = PROJECT_ROOT / "data" / "sequence"
    annotation_dir = PROJECT_ROOT / "data" / "annotation"
    blastdb_dir = PROJECT_ROOT / "data" / "blastdb"

    if any(sequence_dir.glob("*")):
        context.add("sequence_data", "PASS", "sequence data is present")
    else:
        context.add("sequence_data", "WARN", "sequence data not downloaded yet")

    if any(annotation_dir.glob("*")):
        context.add("annotation_data", "PASS", "annotation data is present")
    else:
        context.add("annotation_data", "WARN", "annotation data not downloaded yet")

    if any(blastdb_dir.glob("*")):
        context.add("blastdb", "PASS", "BLAST databases are present")
    else:
        context.add("blastdb", "WARN", "BLAST databases not built yet")


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
    validate_learning_assets(context)
    validate_dataset_presence(context)
    json_path, md_path = write_reports(context)
    summary = "PASS" if not context.failed else "FAIL"
    print(f"Validation summary: {summary}")
    print(f"Validation reports written to: {json_path} and {md_path}")
    return 0 if not context.failed else 1


if __name__ == "__main__":
    sys.exit(run_validation())
