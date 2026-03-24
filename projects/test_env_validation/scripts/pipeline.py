#!/usr/bin/env python3
"""
test_env_validation 主分析流程

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


def step_01_data_preparation():
    """步骤1: 数据准备"""
    print("步骤1: 数据准备")
    raw_dir = Path(CONFIG.RAW_DIR)
    files = list(raw_dir.glob("*"))
    if not files:
        print(f"  警告: {raw_dir} 为空，请先放入原始数据")
        return True
    print(f"  找到 {len(files)} 个输入文件")
    return True


def step_02_quality_control():
    """步骤2: 质量控制"""
    print("步骤2: 质量控制")
    print("  工作区完整工具链验证请运行 run_validation.py")
    return True


def step_03_main_analysis():
    """步骤3: 主要分析"""
    print("步骤3: 主要分析")
    print("  这个项目的真实主入口是 run_validation.py，用于环境 smoke test")
    return True


def step_04_results():
    """步骤4: 结果整理"""
    print("步骤4: 结果整理")
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
