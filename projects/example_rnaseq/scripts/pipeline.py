#!/usr/bin/env python3
"""
example_rnaseq 主分析流程

使用方法:
    python pipeline.py [--step START_STEP] [--steps] [--validate]
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
CONFIG_PATH = SCRIPT_DIR / "config.py"
VALIDATOR_PATH = SCRIPT_DIR / "validate_project.py"


def load_config():
    spec = importlib.util.spec_from_file_location("project_config", CONFIG_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"无法加载配置: {CONFIG_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CONFIG = load_config()


def load_template_runtime():
    for candidate_root in SCRIPT_DIR.parents:
        runtime_path = candidate_root / "lib" / "template_runtime.py"
        if not runtime_path.exists():
            continue
        spec = importlib.util.spec_from_file_location("bio_template_runtime", runtime_path)
        if spec is None or spec.loader is None:
            continue
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    return None


TEMPLATE_RUNTIME = load_template_runtime()


def run_shared_step(step_name: str):
    if TEMPLATE_RUNTIME is None:
        return None
    return TEMPLATE_RUNTIME.run_shared_step(CONFIG, step_name)


def step_01_data_preparation():
    """步骤1: RNA-seq 输入与参考准备"""
    print("步骤1: RNA-seq 输入与参考准备")
    raw_dir = Path(CONFIG.RAW_DIR)
    files = list(raw_dir.glob("*"))
    if not files:
        print(f"  警告: {raw_dir} 为空，请先放入原始数据")
        print("  推荐准备双端 FASTQ、参考基因组和 annotation.gtf。")
        return True
    print(f"  找到 {len(files)} 个输入文件")
    print("  推荐准备双端 FASTQ、参考基因组和 annotation.gtf。")
    return True


def step_02_quality_control():
    """步骤2: 读段质控与清洗"""
    print("步骤2: 读段质控与清洗")
    shared_result = run_shared_step("quality_control")
    if shared_result is not None:
        print("  已通过共享 runtime 调用 lib/modules 执行该步骤。")
        return shared_result
    print("  推荐工具: FastQC / fastp，用于质控、过滤和接头清理。")
    return True


def step_03_main_analysis():
    """步骤3: 比对与定量"""
    print("步骤3: 比对与定量")
    shared_result = run_shared_step("main_analysis")
    if shared_result is not None:
        print("  已通过共享 runtime 调用 lib/modules 执行该步骤。")
        return shared_result
    print("  推荐主流程: STAR/HISAT2 比对 + featureCounts 定量。")
    return True


def step_04_results():
    """步骤4: 结果汇总与表达分析准备"""
    print("步骤4: 结果汇总与表达分析准备")
    print("  建议输出 counts matrix、MultiQC 报告和差异分析输入。")
    print(f"  结果目录: {CONFIG.RESULTS_DIR}")
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
    parser.add_argument("--step", help="从哪个步骤开始")
    parser.add_argument("--steps", action="store_true", help="只列出可用步骤，不执行")
    parser.add_argument("--validate", action="store_true", help="运行项目级自检并退出")
    args = parser.parse_args()

    steps = [
        ("data_preparation", step_01_data_preparation),
        ("quality_control", step_02_quality_control),
        ("main_analysis", step_03_main_analysis),
        ("results", step_04_results),
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

    start_step = args.step or "data_preparation"
    start_idx = next((i for i, (name, _) in enumerate(steps) if name == start_step), None)
    if start_idx is None:
        available = ", ".join(name for name, _ in steps)
        print(f"错误: 未知步骤 {start_step}。可用步骤: {available}")
        sys.exit(1)

    print(f"项目: {CONFIG.PROJECT_NAME}")
    print(f"类型: {CONFIG.PROJECT_TYPE}")
    print(f"从步骤 {start_step} 开始...")
    print("-" * 50)

    for name, func in steps[start_idx:]:
        print(f"\n执行: {name}")
        if not func():
            print(f"错误: 步骤 {name} 失败")
            sys.exit(1)

    print("\n" + "=" * 50)
    print("分析完成！")


if __name__ == "__main__":
    main()
