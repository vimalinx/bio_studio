#!/usr/bin/env python3
"""
cross_system_benchmark 主分析流程

使用方法:
    python pipeline.py [fetch|prepare|baseline|report] [--steps] [--validate]
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
CONFIG_PATH = SCRIPT_DIR / "config.py"
VALIDATOR_PATH = SCRIPT_DIR / "validate_project.py"

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import fetch_benchmarks
import prepare_references
import run_baseline
import summarize_cross_systems


def load_config():
    spec = importlib.util.spec_from_file_location("project_config", CONFIG_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"无法加载配置: {CONFIG_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CONFIG = load_config()


def step_fetch():
    """fetch: 写 benchmark catalog 与标准化元数据"""
    written = fetch_benchmarks.write_catalog(CONFIG.PROJECT_ROOT)
    print(f"  benchmark catalog 已写入: {written['datasets']}")
    return True


def step_prepare():
    """prepare: 生成 reference manifest 并下载默认小体量参考"""
    manifest_path = prepare_references.prepare_references(CONFIG.PROJECT_ROOT)
    print(f"  reference manifest 已写入: {manifest_path}")
    return True


def step_baseline():
    """baseline: 下载 processed benchmark 文件并生成统一表达表"""
    outputs = run_baseline.run_project_baseline(CONFIG.PROJECT_ROOT)
    print(f"  已生成 {len(outputs)} 套 baseline 输出")
    return True


def step_report():
    """report: 生成单数据集摘要与跨系统附加报告"""
    output_path = summarize_cross_systems.run_reporting(CONFIG.PROJECT_ROOT)
    print(f"  综合报告已写入: {output_path}")
    return True


def run_validation() -> int:
    spec = importlib.util.spec_from_file_location("project_validation", VALIDATOR_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"无法加载验证脚本: {VALIDATOR_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.run_validation()


def main() -> None:
    parser = argparse.ArgumentParser(description=f"{CONFIG.PROJECT_NAME} 分析流程")
    parser.add_argument("stage", nargs="?", choices=["fetch", "prepare", "baseline", "report", "all"])
    parser.add_argument("--steps", action="store_true", help="只列出可用步骤，不执行")
    parser.add_argument("--validate", action="store_true", help="运行项目级自检并退出")
    args = parser.parse_args()

    steps = [
        ("fetch", step_fetch),
        ("prepare", step_prepare),
        ("baseline", step_baseline),
        ("report", step_report),
    ]

    if args.steps:
        print("可用步骤:")
        for name, func in steps:
            print(f"  {name}: {func.__doc__}")
        return

    if args.validate:
        exit_code = run_validation()
        report_path = Path(CONFIG.LOGS_DIR) / "validation_report.json"
        print(f"项目验证完成，请查看: {report_path}")
        sys.exit(exit_code)

    print(f"项目: {CONFIG.PROJECT_NAME}")
    print(f"类型: {CONFIG.PROJECT_TYPE}")
    print("-" * 50)

    if args.stage and args.stage != "all":
        run_list = [item for item in steps if item[0] == args.stage]
    else:
        run_list = steps

    for name, func in run_list:
        print(f"\n执行: {name}")
        if not func():
            print(f"错误: 步骤 {name} 失败")
            sys.exit(1)

    print("\n" + "=" * 50)
    print("分析完成！")


if __name__ == "__main__":
    main()
