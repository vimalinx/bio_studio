"""
多序列比对模块 - MUSCLE/MAFFT/ClustalW包装器
"""

from pathlib import Path
from .utils import run_command, ensure_dir, check_file_exists


def align_muscle(input_fasta, output_fasta, maxiters=16, clwstrict=False):
    """
    使用MUSCLE进行多序列比对

    Args:
        input_fasta: 输入FASTA文件
        output_fasta: 输出比对结果FASTA文件
        maxiters: 最大迭代次数
        clwstrict: 是否使用严格的CLUSTAL格式

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(input_fasta)
    ensure_dir(Path(output_fasta).parent)

    cmd = ["muscle", "-in", input_fasta, "-out", output_fasta]

    if maxiters:
        cmd.extend(["-maxiters", str(maxiters)])

    if clwstrict:
        cmd.append("-clwstrict")

    return run_command(cmd, capture_output=False)


def align_mafft(input_fasta, output_fasta, method="auto", threads=4):
    """
    使用MAFFT进行多序列比对

    Args:
        input_fasta: 输入FASTA文件
        output_fasta: 输出比对结果文件
        method: 比对方法 ('auto', 'linsi', 'ginsi', 'einsi', 'fft-ns-2')
        threads: 线程数

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(input_fasta)
    ensure_dir(Path(output_fasta).parent)

    cmd = ["mafft", "--thread", str(threads)]

    if method != "auto":
        cmd.append(method)

    cmd.extend([input_fasta, ">", output_fasta])
    return run_command(" ".join(cmd), shell=True, capture_output=False)


def align_clustalw(input_fasta, output_aln=None, output_dnd=None):
    """
    使用ClustalW进行多序列比对

    Args:
        input_fasta: 输入FASTA文件
        output_aln: 输出比对文件（默认.aln）
        output_dnd: 输出指导树文件（默认.dnd）

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(input_fasta)

    cmd = ["clustalw2", input_fasta]

    if output_aln:
        cmd.extend(["-OUTFILE", output_aln])

    if output_dnd:
        cmd.extend(["-OUTPUTTREE", output_dnd])

    return run_command(cmd, capture_output=False)


def convert_alignment(
    input_file, output_file, input_format="fasta", output_format="clustal"
):
    """
    转换比对格式

    Args:
        input_file: 输入文件
        output_file: 输出文件
        input_format: 输入格式
        output_format: 输出格式

    Returns:
        subprocess.CompletedProcess对象
    """
    from Bio import AlignIO

    check_file_exists(input_file)
    ensure_dir(Path(output_file).parent)

    alignments = AlignIO.read(input_file, input_format)
    AlignIO.write(alignments, output_file, output_format)

    return None
