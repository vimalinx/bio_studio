"""
变异检测模块 - bcftools包装器
"""

from pathlib import Path
from .utils import run_command, ensure_dir, check_file_exists


def call_variants(
    reference, input_bam, output_vcf, caller="mpileup", threads=4, extra_args=None
):
    """
    调用变异（使用bcftools mpileup + call）

    Args:
        reference: 参考序列FASTA文件
        input_bam: 输入BAM文件
        output_vcf: 输出VCF文件
        caller: 调用器类型 ('mpileup' 或 'mutect2')
        threads: 线程数
        extra_args: 额外参数

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(reference)
    check_file_exists(input_bam)
    ensure_dir(Path(output_vcf).parent)

    if caller == "mpileup":
        cmd = [
            "bcftools",
            "mpileup",
            "-f",
            reference,
            "-Ou",
            input_bam,
            "|",
            "bcftools",
            "call",
            "-mv",
            "-Ov",
            "-o",
            output_vcf,
        ]
        return run_command(" ".join(cmd), shell=True)

    return run_command([caller, reference, input_bam, "-o", output_vcf])


def filter_vcf(input_vcf, output_vcf, filters=None, extra_args=None):
    """
    过滤VCF文件

    Args:
        input_vcf: 输入VCF文件
        output_vcf: 输出VCF文件
        filters: 过滤条件（如 'QUAL>30 && DP>10'）
        extra_args: 额外参数

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(input_vcf)
    ensure_dir(Path(output_vcf).parent)

    cmd = ["bcftools", "filter"]

    if filters:
        cmd.extend(["-i", filters])

    if extra_args:
        cmd.extend(extra_args)

    cmd.extend(["-Ov", "-o", output_vcf, input_vcf])
    return run_command(cmd)


def normalize_vcf(input_vcf, output_vcf, reference=None):
    """
    规范化VCF文件（分解、合并多等位位点）

    Args:
        input_vcf: 输入VCF文件
        output_vcf: 输出VCF文件
        reference: 参考序列（可选，用于左对齐）

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(input_vcf)
    ensure_dir(Path(output_vcf).parent)

    cmd = ["bcftools", "norm", "-Ov", "-o", output_vcf, input_vcf]

    if reference:
        cmd.extend(["-f", reference])

    return run_command(cmd)


def get_vcf_stats(input_vcf):
    """
    获取VCF文件统计信息

    Args:
        input_vcf: 输入VCF文件

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(input_vcf)

    cmd = ["bcftools", "stats", input_vcf]
    return run_command(cmd)


def view_vcf(input_vcf, output_vcf=None, samples=None, min_qual=None, types=None):
    """
    查看和转换VCF文件

    Args:
        input_vcf: 输入VCF文件
        output_vcf: 输出文件（可选，默认标准输出）
        samples: 样本列表（只包含这些样本）
        min_qual: 最小质量分数
        types: 变异类型列表（如 'snps', 'indels'）

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(input_vcf)

    cmd = ["bcftools", "view"]

    if samples:
        for s in samples:
            cmd.extend(["-s", s])

    if min_qual:
        cmd.extend(["-i", f"QUAL>{min_qual}"])

    if types:
        cmd.extend(["-v", ",".join(types)])

    cmd.extend(["-Ov", "-o", output_vcf] if output_vcf else [input_vcf])

    return run_command(cmd)
