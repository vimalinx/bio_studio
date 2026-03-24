#!/usr/bin/env python3
"""
创建新项目模板
"""

import argparse
from pathlib import Path
import sys


PROJECT_TYPE_PROFILES = {
    "generic": {
        "readme_overview": """这个模板提供通用分析骨架，适合需要自行定义工具链的项目。

- 默认保留项目级自检、步骤入口和结果目录
- 适合后续按具体任务接入 `lib/modules` 或外部工具
""",
        "config_extra": """PRIMARY_INPUT = RAW_DIR
DELIVERABLE = RESULTS_DIR / "summary.md"
""",
        "step1_doc": "通用输入与参考准备",
        "step1_hint": "通用模板：先明确输入文件、参考文件和期望产出。",
        "step2_doc": "预处理或质量控制",
        "step2_hint": "在这里接入清洗、过滤、标准化或质量控制步骤。",
        "step3_doc": "核心分析",
        "step3_hint": "在这里实现你的核心分析逻辑或调用现有模块。",
        "step4_doc": "结果整理与报告",
        "step4_hint": "建议输出可复现报告、关键表格和最终图件。",
    },
    "rnaseq": {
        "readme_overview": """这个模板预置了 **RNA-seq** 最小骨架，默认围绕 `FastQC/Fastp`、`STAR`、`featureCounts` 和 `MultiQC` 组织。

- 推荐把参考基因组放在 `data/references/genome.fa`
- 推荐把注释放在 `data/references/annotation.gtf`
- 双端测序默认使用 `*_R1.fastq.gz` / `*_R2.fastq.gz`
""",
        "config_extra": """READ1_PATTERN = "*_R1.fastq.gz"
READ2_PATTERN = "*_R2.fastq.gz"
REFERENCE_GENOME = REFERENCES_DIR / "genome.fa"
ANNOTATION_GTF = REFERENCES_DIR / "annotation.gtf"
ALIGNER = "STAR"
QUANTIFIER = "featureCounts"
QC_TOOL = "fastp"
REPORT_TOOL = "MultiQC"
""",
        "step1_doc": "RNA-seq 输入与参考准备",
        "step1_hint": "推荐准备双端 FASTQ、参考基因组和 annotation.gtf。",
        "step2_doc": "读段质控与清洗",
        "step2_hint": "推荐工具: FastQC / fastp，用于质控、过滤和接头清理。",
        "step3_doc": "比对与定量",
        "step3_hint": "推荐主流程: STAR/HISAT2 比对 + featureCounts 定量。",
        "step4_doc": "结果汇总与表达分析准备",
        "step4_hint": "建议输出 counts matrix、MultiQC 报告和差异分析输入。",
    },
    "variant": {
        "readme_overview": """这个模板预置了 **变异检测** 最小骨架，默认围绕 `fastp`、`bwa-mem`、`samtools` 和 `bcftools` 组织。

- 推荐把参考基因组放在 `data/references/genome.fa`
- 适合病毒、小基因组或目标区域的快速变异检测原型
- 可在后续接入注释、过滤和可视化步骤
""",
        "config_extra": """READ1_PATTERN = "*_R1.fastq.gz"
READ2_PATTERN = "*_R2.fastq.gz"
REFERENCE_GENOME = REFERENCES_DIR / "genome.fa"
ALIGNER = "bwa-mem"
VARIANT_CALLER = "bcftools"
FILTER_EXPRESSION = "QUAL>30 && DP>10"
TARGET_REGIONS = None
""",
        "step1_doc": "测序输入与参考准备",
        "step1_hint": "推荐准备 FASTQ/BAM 输入和参考基因组 fasta。",
        "step2_doc": "读段质控与比对准备",
        "step2_hint": "推荐工具: fastp 质控，bwa-mem 或 Bowtie2 建立比对输入。",
        "step3_doc": "比对与变异检测",
        "step3_hint": "推荐主流程: bwa-mem 比对、samtools 排序、bcftools 调用变异。",
        "step4_doc": "VCF 过滤与报告",
        "step4_hint": "建议输出过滤后的 VCF、基础统计和可视化摘要。",
    },
    "phylogeny": {
        "readme_overview": """这个模板预置了 **系统发育** 最小骨架，默认围绕 `MAFFT` 多序列比对和 `IQ-TREE 2` 建树组织。

- 推荐把待分析序列放在 `data/raw/sequences.fasta`
- 适合快速搭建序列收集、比对、建树和树文件整理流程
- 可在后续接入注释、可视化和 clade 解释
""",
        "config_extra": """ALIGNMENT_INPUT = RAW_DIR / "sequences.fasta"
MSA_TOOL = "mafft"
TREE_BUILDER = "iqtree2"
TREE_MODEL = "MFP"
BOOTSTRAP = 1000
OUTGROUP = None
""",
        "step1_doc": "序列集与元数据准备",
        "step1_hint": "推荐准备序列 FASTA、样本元数据和必要的外群设置。",
        "step2_doc": "多序列比对",
        "step2_hint": "推荐工具: MAFFT / MUSCLE，用于生成一致的 alignment。",
        "step3_doc": "系统树构建",
        "step3_hint": "推荐主流程: MAFFT 比对后使用 iqtree2 建树与 bootstrap。",
        "step4_doc": "树文件与报告整理",
        "step4_hint": "建议输出 treefile、模型选择结果和可视化说明。",
    },
    "genome": {
        "readme_overview": """这个模板面向基因组分析任务，可作为组装、注释或比较基因组学的起点。

- 建议按任务补充 assembly、annotation 或 comparative genomics 子步骤
- 如需更明确流程，可从通用模板或其他类型模板继续裁剪
""",
        "config_extra": """REFERENCE_GENOME = REFERENCES_DIR / "genome.fa"
ASSEMBLY_INPUT = RAW_DIR
ANNOTATION_OUTPUT = RESULTS_DIR / "annotation"
""",
        "step1_doc": "基因组输入与参考准备",
        "step1_hint": "先确认原始读段、装配输入或参考基因组位置。",
        "step2_doc": "预处理与组装准备",
        "step2_hint": "在这里接入清洗、纠错或装配前准备步骤。",
        "step3_doc": "核心基因组分析",
        "step3_hint": "在这里实现组装、注释或比较基因组主流程。",
        "step4_doc": "结果归档与报告",
        "step4_hint": "建议输出最终装配、注释文件和关键统计摘要。",
    },
}


README_TEMPLATE = """# {project_name}

## 项目类型
{project_type}

## 描述
{description}

## 项目结构

```
{project_name}/
├── data/                  # 项目数据（完全独立）
│   ├── raw/               # 原始输入数据
│   ├── processed/         # 中间产物
│   ├── results/           # 最终结果
│   └── references/        # 项目内参考数据
├── scripts/               # 项目脚本
│   ├── config.py          # 项目配置与路径定义
│   ├── validate_project.py # 项目级自检入口
│   └── pipeline.py        # 主分析流程入口
├── notebooks/             # 探索性分析
├── logs/                  # 运行日志与验证报告
└── README.md              # 本文件
```

## 快速开始

### 1. 自检模板
```bash
cd scripts
python validate_project.py
python pipeline.py --validate
python pipeline.py --steps
```

### 2. 准备数据
将原始数据放入 `data/raw/`。

### 3. 编辑配置
按项目需要更新 `scripts/config.py` 中的样本、参考和线程数。

### 4. 运行流程
```bash
cd scripts
python pipeline.py
```

## 说明

- `validate_project.py` 负责项目级结构与配置自检
- `logs/validation_report.json` 和 `logs/validation_report.md` 会记录验证结果
- 缺少真实数据或参考文件时默认给出警告，不阻断模板自检

## 类型骨架

{type_overview}

## 项目状态

- [ ] 模板自检通过
- [ ] 数据准备
- [ ] 质量控制
- [ ] 主要分析
- [ ] 结果验证
"""


CONFIG_TEMPLATE = '''"""
{project_name} 项目配置
"""

from pathlib import Path

PROJECT_NAME = {project_name_literal}
PROJECT_TYPE = {project_type_literal}
PROJECT_DESCRIPTION = {description_literal}

SCRIPTS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPTS_DIR.parent
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = DATA_DIR / "results"
REFERENCES_DIR = DATA_DIR / "references"
LOGS_DIR = PROJECT_ROOT / "logs"
NOTEBOOKS_DIR = PROJECT_ROOT / "notebooks"

SAMPLES = []
REFERENCE_GENOME = None
THREADS = 4
ANALYSIS_PARAMETERS = {{}}

{type_config_extra}
'''


VALIDATE_TEMPLATE = '''#!/usr/bin/env python3
"""
{project_name} 项目级自检
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
        self.results.append({{"name": name, "status": status, "message": message}})
        print(f"[{{status}}] {{name}}: {{message}}")

    @property
    def failed(self) -> bool:
        return any(item["status"] == "FAIL" for item in self.results)


def load_config(context: ValidationContext):
    if not CONFIG_PATH.exists():
        context.add("config_file", "FAIL", f"missing {{CONFIG_PATH.name}}")
        return None

    spec = importlib.util.spec_from_file_location("project_config", CONFIG_PATH)
    if spec is None or spec.loader is None:
        context.add("config_import", "FAIL", "unable to create config import spec")
        return None

    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except Exception as exc:  # pragma: no cover - defensive path for broken configs
        context.add("config_import", "FAIL", f"failed to import config: {{exc}}")
        return None

    context.config = module
    context.add("config_import", "PASS", f"imported {{CONFIG_PATH.name}} successfully")
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
            f"PROJECT_ROOT={{configured_project_root}} (expected {{actual_project_root}})"
        )

    for field in path_fields:
        value = getattr(config, field, None)
        if value is None:
            continue
        path = Path(value).resolve()
        if path != actual_project_root and actual_project_root not in path.parents:
            invalid.append(f"{{field}}={{path}}")

    if invalid:
        context.add("path_scope", "FAIL", "paths outside project root: " + "; ".join(invalid))
    else:
        context.add("path_scope", "PASS", "configured paths stay within project root")


def validate_optional_inputs(context: ValidationContext, config) -> None:
    raw_files = list(Path(config.RAW_DIR).glob("*"))
    if raw_files:
        context.add("raw_inputs", "PASS", f"found {{len(raw_files)}} files in raw data directory")
    else:
        context.add("raw_inputs", "WARN", "no raw input files yet; template is ready for data ingestion")

    reference = getattr(config, "REFERENCE_GENOME", None)
    if reference:
        reference_path = Path(reference)
        if not reference_path.is_absolute():
            reference_path = Path(config.PROJECT_ROOT) / reference_path
        if reference_path.exists():
            context.add("reference_genome", "PASS", f"reference file exists: {{reference_path}}")
        else:
            context.add("reference_genome", "WARN", f"reference file not found yet: {{reference_path}}")
    else:
        context.add("reference_genome", "WARN", "REFERENCE_GENOME is not configured yet")


def write_reports(context: ValidationContext) -> tuple[Path, Path]:
    logs_dir = PROJECT_ROOT / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    summary = "PASS" if not context.failed else "FAIL"
    payload = {{
        "project_root": str(PROJECT_ROOT),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "results": context.results,
    }}

    json_path = logs_dir / "validation_report.json"
    md_path = logs_dir / "validation_report.md"
    json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\\n", encoding="utf-8")

    md_lines = [
        f"# Validation Report: {{payload['summary']}}",
        "",
        f"- Project root: `{{PROJECT_ROOT}}`",
        f"- Generated at: `{{payload['generated_at']}}`",
        "",
        "## Checks",
        "",
    ]
    for item in context.results:
        md_lines.append(f"- **{{item['status']}}** `{{item['name']}}`: {{item['message']}}")
    md_path.write_text("\\n".join(md_lines) + "\\n", encoding="utf-8")
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
    print(f"Validation summary: {{summary}}")
    print(f"Validation reports written to: {{json_path}} and {{md_path}}")
    return 0 if not context.failed else 1


if __name__ == "__main__":
    sys.exit(run_validation())
'''


PIPELINE_TEMPLATE = '''#!/usr/bin/env python3
"""
{project_name} 主分析流程

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
        raise RuntimeError(f"无法加载配置: {{CONFIG_PATH}}")
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
    """步骤1: {step1_doc}"""
    print("步骤1: {step1_doc}")
    raw_dir = Path(CONFIG.RAW_DIR)
    files = list(raw_dir.glob("*"))
    if not files:
        print(f"  警告: {{raw_dir}} 为空，请先放入原始数据")
        print("  {step1_hint}")
        return True
    print(f"  找到 {{len(files)}} 个输入文件")
    print("  {step1_hint}")
    return True


def step_02_quality_control():
    """步骤2: {step2_doc}"""
    print("步骤2: {step2_doc}")
    shared_result = run_shared_step("quality_control")
    if shared_result is not None:
        print("  已通过共享 runtime 调用 lib/modules 执行该步骤。")
        return shared_result
    print("  {step2_hint}")
    return True


def step_03_main_analysis():
    """步骤3: {step3_doc}"""
    print("步骤3: {step3_doc}")
    shared_result = run_shared_step("main_analysis")
    if shared_result is not None:
        print("  已通过共享 runtime 调用 lib/modules 执行该步骤。")
        return shared_result
    print("  {step3_hint}")
    return True


def step_04_results():
    """步骤4: {step4_doc}"""
    print("步骤4: {step4_doc}")
    print("  {step4_hint}")
    print(f"  结果目录: {{CONFIG.RESULTS_DIR}}")
    return True


def run_validation() -> int:
    spec = importlib.util.spec_from_file_location("project_validation", VALIDATOR_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"无法加载验证脚本: {{VALIDATOR_PATH}}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.run_validation()


def main() -> None:
    parser = argparse.ArgumentParser(description=f"{{CONFIG.PROJECT_NAME}} 分析流程")
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
            print(f"  {{name}}: {{func.__doc__}}")
        return

    if args.validate:
        exit_code = run_validation()
        report_path = Path(CONFIG.LOGS_DIR) / "validation_report.json"
        print(f"项目验证完成，请查看: {{report_path}}")
        sys.exit(exit_code)

    start_step = args.step or "data_preparation"
    start_idx = next((i for i, (name, _) in enumerate(steps) if name == start_step), None)
    if start_idx is None:
        available = ", ".join(name for name, _ in steps)
        print(f"错误: 未知步骤 {{start_step}}。可用步骤: {{available}}")
        sys.exit(1)

    print(f"项目: {{CONFIG.PROJECT_NAME}}")
    print(f"类型: {{CONFIG.PROJECT_TYPE}}")
    print(f"从步骤 {{start_step}} 开始...")
    print("-" * 50)

    for name, func in steps[start_idx:]:
        print(f"\\n执行: {{name}}")
        if not func():
            print(f"错误: 步骤 {{name}} 失败")
            sys.exit(1)

    print("\\n" + "=" * 50)
    print("分析完成！")


if __name__ == "__main__":
    main()
'''


def get_project_type_profile(project_type: str) -> dict[str, str]:
    return PROJECT_TYPE_PROFILES.get(project_type, PROJECT_TYPE_PROFILES["generic"])


def build_readme_content(project_name: str, project_type: str, description: str) -> str:
    profile = get_project_type_profile(project_type)
    return README_TEMPLATE.format(
        project_name=project_name,
        project_type=project_type,
        description=description or "待补充",
        type_overview=profile["readme_overview"],
    )


def build_config_content(project_name: str, project_type: str, description: str) -> str:
    profile = get_project_type_profile(project_type)
    return (
        CONFIG_TEMPLATE.format(
        project_name=project_name,
        project_name_literal=repr(project_name),
        project_type_literal=repr(project_type),
        description_literal=repr(description),
        type_config_extra=profile["config_extra"],
        ).rstrip()
        + "\n"
    )


def build_pipeline_content(project_name: str, project_type: str) -> str:
    profile = get_project_type_profile(project_type)
    return PIPELINE_TEMPLATE.format(
        project_name=project_name,
        step1_doc=profile["step1_doc"],
        step1_hint=profile["step1_hint"],
        step2_doc=profile["step2_doc"],
        step2_hint=profile["step2_hint"],
        step3_doc=profile["step3_doc"],
        step3_hint=profile["step3_hint"],
        step4_doc=profile["step4_doc"],
        step4_hint=profile["step4_hint"],
    )


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

    for directory in dirs:
        directory.mkdir(parents=True)
        print(f"  创建目录: {directory}")

    readme_content = build_readme_content(project_name, project_type, description)
    (project_dir / "README.md").write_text(readme_content, encoding="utf-8")
    print(f"  创建文件: {project_dir}/README.md")

    config_content = build_config_content(project_name, project_type, description)
    (project_dir / "scripts" / "config.py").write_text(config_content, encoding="utf-8")
    print(f"  创建文件: {project_dir}/scripts/config.py")

    validate_content = VALIDATE_TEMPLATE.format(project_name=project_name)
    validate_path = project_dir / "scripts" / "validate_project.py"
    validate_path.write_text(validate_content, encoding="utf-8")
    validate_path.chmod(0o755)
    print(f"  创建文件: {project_dir}/scripts/validate_project.py")

    pipeline_content = build_pipeline_content(project_name, project_type)
    pipeline_path = project_dir / "scripts" / "pipeline.py"
    pipeline_path.write_text(pipeline_content, encoding="utf-8")
    pipeline_path.chmod(0o755)
    print(f"  创建文件: {project_dir}/scripts/pipeline.py")

    print(f"\\n✓ 项目创建完成: {project_dir}")
    print("\\n下一步:")
    print(f"  1. cd {project_dir}/scripts")
    print("  2. python validate_project.py")
    print("  3. python pipeline.py --validate")
    print("  4. 准备 data/raw/ 后再运行 python pipeline.py")


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
