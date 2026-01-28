"""
基因预测模块 - Prodigal包装器
"""

from pathlib import Path
from .utils import run_command, ensure_dir, check_file_exists


def predict_genes(
    input_fasta,
    output_gff=None,
    output_fna=None,
    output_faa=None,
    mode="single",
    meta=False,
):
    """
    使用Prodigal预测原核生物基因

    Args:
        input_fasta: 输入基因组FASTA文件
        output_gff: 输出GFF文件
        output_fna: 输出核苷酸序列
        output_faa: 输出氨基酸序列
        mode: 预测模式 ('single', 'meta', 'train')
        meta: 是否使用元基因组模式

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(input_fasta)

    if output_gff:
        ensure_dir(Path(output_gff).parent)
    if output_fna:
        ensure_dir(Path(output_fna).parent)
    if output_faa:
        ensure_dir(Path(output_faa).parent)

    cmd = ["prodigal"]

    if mode == "single":
        cmd.extend(["-p", "single"])
    elif mode == "meta":
        cmd.extend(["-p", "meta"])
    elif mode == "train":
        cmd.extend(["-p", "train"])

    if meta:
        cmd.append("-m")

    if output_gff:
        cmd.extend(["-f", "gff", "-o", output_gff])
    if output_fna:
        cmd.extend(["-d", output_fna])
    if output_faa:
        cmd.extend(["-a", output_faa])

    cmd.extend(["-i", input_fasta])

    return run_command(cmd, capture_output=False)


def train_prodigal(input_fasta, training_file):
    """
    训练Prodigal模型

    Args:
        input_fasta: 输入基因组FASTA文件
        training_file: 输出训练文件

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(input_fasta)
    ensure_dir(Path(training_file).parent)

    cmd = ["prodigal", "-p", "train", "-i", input_fasta, "-o", training_file]
    return run_command(cmd, capture_output=False)


def use_training_file(input_fasta, training_file, output_gff):
    """
    使用训练好的文件预测

    Args:
        input_fasta: 输入基因组FASTA文件
        training_file: 训练文件
        output_gff: 输出GFF文件

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(input_fasta)
    check_file_exists(training_file)
    ensure_dir(Path(output_gff).parent)

    cmd = [
        "prodigal",
        "-p",
        "single",
        "-i",
        input_fasta,
        "-t",
        training_file,
        "-f",
        "gff",
        "-o",
        output_gff,
    ]
    return run_command(cmd, capture_output=False)
