#!/usr/bin/env python3
"""
test_rnaseq_analysis 主分析流程 (Ebola Virus Variant Calling)
"""

import sys
import argparse
import subprocess
import os
import importlib.util
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
REPO_ROOT = PROJECT_DIR.parents[1]
for path in (REPO_ROOT, SCRIPT_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

try:
    from lib.workspace_env import ensure_workspace_path
except ModuleNotFoundError:
    ensure_workspace_path = None

if ensure_workspace_path is not None:
    ensure_workspace_path()

# 导入配置，但后续会重新定义为 Path 对象以避免类型混淆
try:
    import config
except ImportError:
    print("❌ Cannot import config.py")
    sys.exit(1)

# 将字符串路径转换为绝对路径 Path 对象
RAW_DIR = PROJECT_DIR / config.RAW_DIR
PROCESSED_DIR = PROJECT_DIR / config.PROCESSED_DIR
RESULTS_DIR = PROJECT_DIR / config.RESULTS_DIR
# config 中是 REFERENCES_DIR，这里统称为 REF_DIR
REF_DIR = PROJECT_DIR / config.REFERENCES_DIR

QC_DIR = PROCESSED_DIR / "qc"
ALIGN_DIR = PROCESSED_DIR / "aligned"

# 确保子目录存在
QC_DIR.mkdir(parents=True, exist_ok=True)
ALIGN_DIR.mkdir(parents=True, exist_ok=True)
VALIDATOR_PATH = SCRIPT_DIR / "validate_project.py"

def run_cmd(cmd, step_name):
    print(f"[{step_name}] 🚀 Executing: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
        print(f"[{step_name}] ✅ Success")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[{step_name}] ❌ Failed with code {e.returncode}")
        return False

def step_01_data_preparation():
    """步骤1: 数据检查与参考基因组准备"""
    print("\n=== 步骤1: 数据准备 ===")
    
    # 1. 检查原始数据
    r1 = list(RAW_DIR.glob("*_1.fastq.gz"))[0]
    r2 = list(RAW_DIR.glob("*_2.fastq.gz"))[0]
    print(f"  Raw Data: {r1.name}, {r2.name}")
    
    # 2. 检查参考基因组
    ref_fa = REF_DIR / "ebola.fa"
    if not ref_fa.exists():
        print(f"  ❌ Reference genome not found at {ref_fa}")
        return False
    print(f"  Reference: {ref_fa.name}")
    
    # 3. 构建 Bowtie2 索引
    bt2_idx = REF_DIR / "ebola_idx"
    if not (REF_DIR / "ebola_idx.1.bt2").exists():
        print("  Building Bowtie2 index...")
        cmd = f"bowtie2-build {ref_fa} {bt2_idx}"
        if not run_cmd(cmd, "Index"): return False
    else:
        print("  Bowtie2 index exists.")
        
    return True

def step_02_quality_control():
    """步骤2: 质控 (Fastp)"""
    print("\n=== 步骤2: 质量控制 (Fastp) ===")
    
    r1 = list(RAW_DIR.glob("*_1.fastq.gz"))[0]
    r2 = list(RAW_DIR.glob("*_2.fastq.gz"))[0]
    
    out1 = QC_DIR / "clean_R1.fq.gz"
    out2 = QC_DIR / "clean_R2.fq.gz"
    html = QC_DIR / "fastp.html"
    json = QC_DIR / "fastp.json"
    
    cmd = f"fastp -i {r1} -I {r2} -o {out1} -O {out2} -h {html} -j {json} --thread 4"
    
    if not out1.exists():
        return run_cmd(cmd, "Fastp")
    else:
        print("  Fastp results exist, skipping.")
        return True

def step_03_main_analysis():
    """步骤3: 比对与变异检测"""
    print("\n=== 步骤3: 比对与变异检测 ===")
    
    # 1. Bowtie2 Alignment
    idx = REF_DIR / "ebola_idx"
    r1 = QC_DIR / "clean_R1.fq.gz"
    r2 = QC_DIR / "clean_R2.fq.gz"
    bam = ALIGN_DIR / "aligned.bam"
    sorted_bam = ALIGN_DIR / "aligned.sorted.bam"
    
    # Pipe: bowtie2 -> samtools view -> samtools sort
    cmd_align = (
        f"bowtie2 -x {idx} -1 {r1} -2 {r2} -p 4 | "
        f"samtools view -bS - | "
        f"samtools sort -o {sorted_bam}"
    )
    
    if not sorted_bam.exists():
        if not run_cmd(cmd_align, "Alignment"): return False
        run_cmd(f"samtools index {sorted_bam}", "Index BAM")
    else:
        print("  BAM file exists, skipping alignment.")

    # 2. Variant Calling (bcftools)
    # mpileup -> call -> normalize
    ref = REF_DIR / "ebola.fa"
    vcf = RESULTS_DIR / "variants.vcf"
    
    cmd_call = (
        f"bcftools mpileup -Ou -f {ref} {sorted_bam} | "
        f"bcftools call -mv -Ob -o {RESULTS_DIR}/raw.bcf && "
        f"bcftools view {RESULTS_DIR}/raw.bcf > {vcf}"
    )
    
    if not vcf.exists():
        return run_cmd(cmd_call, "Variant Calling")
    else:
        print("  VCF exists, skipping calling.")
        return True

def step_04_results():
    """步骤4: 报告汇总"""
    print("\n=== 步骤4: 结果汇总 ===")
    
    # 1. Bam Stats
    bam = ALIGN_DIR / "aligned.sorted.bam"
    run_cmd(f"samtools flagstat {bam} > {RESULTS_DIR}/flagstat.txt", "Flagstat")
    run_cmd(f"samtools idxstats {bam} > {RESULTS_DIR}/idxstats.txt", "Idxstats")
    
    # 2. MultiQC
    # 扫描整个项目目录
    cmd_mqc = f"multiqc {PROJECT_DIR} -o {RESULTS_DIR}/multiqc_report --force"
    run_cmd(cmd_mqc, "MultiQC")
    
    print(f"\n✅ 分析完成! 报告位置: {RESULTS_DIR}/multiqc_report/multiqc_report.html")
    return True


def run_validation() -> int:
    spec = importlib.util.spec_from_file_location("special_project_validation", VALIDATOR_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"无法加载验证脚本: {VALIDATOR_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.run_validation()


def print_steps() -> None:
    steps = [
        ("data_preparation", step_01_data_preparation.__doc__),
        ("quality_control", step_02_quality_control.__doc__),
        ("main_analysis", step_03_main_analysis.__doc__),
        ("results", step_04_results.__doc__),
    ]
    print("可用步骤:")
    for name, doc in steps:
        print(f"  {name}: {doc}")

def main():
    parser = argparse.ArgumentParser(description='Ebola Virus Analysis Pipeline')
    parser.add_argument('--step', help='Start from specific step')
    parser.add_argument('--steps', action='store_true', help='List available steps only')
    parser.add_argument('--validate', action='store_true', help='Run project validation and exit')
    args = parser.parse_args()

    steps = [
        ('data_preparation', step_01_data_preparation),
        ('quality_control', step_02_quality_control),
        ('main_analysis', step_03_main_analysis),
        ('results', step_04_results),
    ]

    if args.steps:
        print_steps()
        return

    if args.validate:
        exit_code = run_validation()
        print(f"项目验证完成，请查看: {PROJECT_DIR / 'logs' / 'validation_report.json'}")
        sys.exit(exit_code)

    start_idx = 0
    if args.step:
        start_idx = next((i for i, (n, _) in enumerate(steps) if n == args.step), 0)

    for name, func in steps[start_idx:]:
        if not func():
            print(f"❌ Pipeline failed at step: {name}")
            sys.exit(1)

if __name__ == '__main__':
    main()
