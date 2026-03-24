#!/usr/bin/env python3
"""
yeast_genome_learning 兼容入口
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
VALIDATOR_PATH = SCRIPT_DIR / "validate_project.py"


def run_validation() -> int:
    spec = importlib.util.spec_from_file_location("learning_project_validation", VALIDATOR_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"无法加载验证脚本: {VALIDATOR_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.run_validation()


def print_steps() -> None:
    print("可用步骤:")
    print("  setup_database: 运行 bash scripts/01_setup_database.sh")
    print("  verify_install: 运行 bash scripts/02_verify_install.sh")
    print("  extract_gene: 运行 bash scripts/03_extract_gene.sh <GENE>")
    print("  sequence_analysis: 运行 bash scripts/04_sequence_analysis.sh")
    print("  blast_search: 运行 bash scripts/05_blast_search.sh")
    print("  gene_summary: 运行 bash scripts/analyze_genes.sh")


def print_workspace_cli_guidance() -> None:
    if os.environ.get("BIO_STUDIO_WORKSPACE_CLI") != "1":
        return

    cli_path = os.environ.get("BIO_STUDIO_WORKSPACE_CLI_PATH")
    project_name = os.environ.get("BIO_STUDIO_PROJECT_NAME", PROJECT_ROOT.name)
    if not cli_path:
        return

    print("也可以继续使用工作区统一入口：")
    print(f"  python {cli_path} steps {project_name}")
    print(f"  python {cli_path} validate {project_name}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Yeast genome learning compatibility pipeline")
    parser.add_argument("--steps", action="store_true", help="只列出学习步骤，不执行")
    parser.add_argument("--validate", action="store_true", help="运行项目级自检并退出")
    args = parser.parse_args()

    if args.steps:
        print_steps()
        return

    if args.validate:
        exit_code = run_validation()
        report_path = PROJECT_ROOT / "logs" / "validation_report.json"
        print(f"项目验证完成，请查看: {report_path}")
        sys.exit(exit_code)

    print("这是学习型项目，默认主入口仍是 scripts/ 下的 Bash 教学脚本。")
    print("先运行以下命令查看步骤：")
    print("  python pipeline.py --steps")
    print_workspace_cli_guidance()
    print("或阅读：")
    print(f"  {PROJECT_ROOT / 'QUICKSTART.md'}")


if __name__ == "__main__":
    main()
