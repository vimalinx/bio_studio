#!/usr/bin/env python3
"""
创建新项目模板
"""

import argparse
from pathlib import Path
import sys


def create_project(project_name, project_type="generic", description=""):
    """
    创建新项目目录结构

    Args:
        project_name: 项目名称
        project_type: 项目类型 ('generic', 'rnaseq', 'variant', 'phylogeny', 'genome')
        description: 项目描述
    """
    project_dir = Path("projects") / project_name

    if project_dir.exists():
        print(f"错误: 项目 {project_name} 已存在")
        sys.exit(1)

    print(f"创建项目: {project_name}")
    print(f"类型: {project_type}")

    dirs = [
        project_dir / "data" / "raw",
        project_dir / "data" / "processed",
        project_dir / "data" / "results",
        project_dir / "data" / "references",
        project_dir / "scripts",
        project_dir / "notebooks",
        project_dir / "logs",
    ]

    for d in dirs:
        d.mkdir(parents=True)
        print(f"  创建目录: {d}")

    readme_content = f"""# {project_name}

## 项目类型
{project_type}

## 描述
{description}

## 项目结构

```
{project_name}/
├── data/                  # 项目数据（完全独立）
│   ├── raw/              # 原始数据
│   ├── processed/         # 中间结果
│   ├── results/          # 最终结果
│   └── references/       # 项目特定参考序列
├── scripts/              # 项目脚本（调用lib模块）
│   ├── pipeline.py       # 主要分析流程
│   ├── config.py         # 项目配置
│   └── analysis.py      # 项目特定分析
├── notebooks/           # Jupyter notebooks
├── logs/               # 日志文件
└── README.md           # 本文件
```

## 使用方法

### 1. 准备数据
将原始数据放入 `data/raw/`

### 2. 编辑配置
编辑 `scripts/config.py` 设置项目参数

### 3. 运行分析
```bash
cd scripts
python pipeline.py
```

## 项目状态

- [ ] 数据准备
- [ ] 质量控制
- [ ] 主要分析
- [ ] 结果验证

## 笔记
记录你的分析过程和发现
"""

    (project_dir / "README.md").write_text(readme_content)
    print(f"  创建文件: {project_dir}/README.md")

    config_content = f'''"""
{project_name} 项目配置
"""

# 项目配置
PROJECT_NAME = "{project_name}"
PROJECT_TYPE = "{project_type}"

# 数据路径
DATA_DIR = "data"
RAW_DIR = f"{{DATA_DIR}}/raw"
PROCESSED_DIR = f"{{DATA_DIR}}/processed"
RESULTS_DIR = f"{{DATA_DIR}}/results"
REFERENCES_DIR = f"{{DATA_DIR}}/references"

# 样本配置
SAMPLES = []

# 参考基因组
REFERENCE_GENOME = None

# 工具配置
THREADS = 4

# 分析参数
'''

    (project_dir / "scripts" / "config.py").write_text(config_content)
    print(f"  创建文件: {project_dir}/scripts/config.py")

    pipeline_content = f'''#!/usr/bin/env python3
"""
{project_name} 主分析流程

使用方法:
    python pipeline.py [--step START_STEP] [--steps]
"""

import sys
import argparse
from pathlib import Path

bio_studio_path = str(Path(__file__).resolve().parent.parent.parent.parent)
sys.path.insert(0, bio_studio_path)

from config import *
from lib.modules import utils


def step_01_data_preparation():
    """步骤1: 数据准备"""
    print("步骤1: 数据准备")
    print("  检查原始数据...")
    raw_dir = Path(RAW_DIR)

    if not raw_dir.exists():
        print(f"  错误: 找不到原始数据目录 {{RAW_DIR}}")
        return False

    files = list(raw_dir.glob('*'))
    if not files:
        print(f"  警告: {{RAW_DIR}} 为空")
        return False

    print(f"  找到 {{len(files)}} 个文件")
    return True


def step_02_quality_control():
    """步骤2: 质量控制"""
    print("步骤2: 质量控制")
    print("  TODO: 实现QC步骤")
    return True


def step_03_main_analysis():
    """步骤3: 主要分析"""
    print("步骤3: 主要分析")
    print("  TODO: 实现分析步骤")
    return True


def step_04_results():
    """步骤4: 结果整理"""
    print("步骤4: 结果整理")
    print("  TODO: 整理结果")
    return True


def main():
    parser = argparse.ArgumentParser(description='{{PROJECT_NAME}} 分析流程')
    parser.add_argument('--step', help='从哪个步骤开始')
    parser.add_argument('--steps', action='store_true',
                       help='只列出可用步骤，不执行')

    args = parser.parse_args()

    steps = [
        ('data_preparation', step_01_data_preparation),
        ('quality_control', step_02_quality_control),
        ('main_analysis', step_03_main_analysis),
        ('results', step_04_results),
    ]

    if args.steps:
        print("可用步骤:")
        for name, func in steps:
            print(f"  {{name}}: {{func.__doc__}}")
        return

    start_step = args.step or 'data_preparation'
    start_idx = next((i for i, (name, _) in enumerate(steps)
                     if name == start_step), 0)

    print(f"从步骤 {{start_step}} 开始...")
    print("-" * 50)

    for name, func in steps[start_idx:]:
        print(f"\\n执行: {{name}}")
        if not func():
            print(f"错误: 步骤 {{name}} 失败")
            sys.exit(1)

    print("\\n" + "=" * 50)
    print("分析完成！")


if __name__ == '__main__':
    main()
'''

    (project_dir / "scripts" / "pipeline.py").write_text(pipeline_content)
    (project_dir / "scripts" / "pipeline.py").chmod(0o755)
    print(f"  创建文件: {project_dir}/scripts/pipeline.py")

    print(f"\\n✓ 项目创建完成: {project_dir}")
    print(f"\\n下一步:")
    print(f"  1. cd {project_dir}")
    print(f"  2. 编辑 scripts/config.py")
    print(f"  3. 将数据放入 data/raw/")
    print(f"  4. cd scripts && python pipeline.py")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="创建新的生物信息学项目")
    parser.add_argument("project_name", help="项目名称")
    parser.add_argument(
        "--type",
        "-t",
        default="generic",
        choices=["generic", "rnaseq", "variant", "phylogeny", "genome"],
        help="项目类型",
    )
    parser.add_argument("--description", "-d", default="", help="项目描述")

    args = parser.parse_args()

    create_project(args.project_name, args.type, args.description)
