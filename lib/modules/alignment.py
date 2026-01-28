"""
序列比对模块 - BWA/Bowtie2包装器
"""

from pathlib import Path
from .utils import run_command, ensure_dir, check_file_exists


def build_bwa_index(reference, prefix=None):
    """
    构建BWA索引

    Args:
        reference: 参考序列FASTA文件
        prefix: 输出前缀（默认与reference相同）

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(reference)

    if prefix is None:
        prefix = reference

    ensure_dir(Path(prefix).parent)

    cmd = ["bwa", "index", "-p", prefix, reference]
    return run_command(cmd)


def align_bwa_mem(
    reference, reads1, reads2=None, output_sam=None, threads=4, extra_args=None
):
    """
    使用BWA-MEM进行序列比对

    Args:
        reference: 参考序列（或BWA索引前缀）
        reads1: 第一个read文件（单端或双端的第一条）
        reads2: 第二个read文件（双端测序，可选）
        output_sam: 输出SAM文件路径
        threads: 线程数
        extra_args: 额外参数（列表）

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(reads1)

    if reads2:
        check_file_exists(reads2)

    cmd = ["bwa", "mem", "-t", str(threads)]

    if extra_args:
        cmd.extend(extra_args)

    cmd.append(reference)
    cmd.append(reads1)

    if reads2:
        cmd.append(reads2)

    if output_sam:
        ensure_dir(Path(output_sam).parent)
        result = run_command(cmd, capture_output=False)
        with open(output_sam, "w") as f:
            f.write(result.stdout)
        return result

    return run_command(cmd)


def build_bowtie2_index(reference, index_name=None):
    """
    构建Bowtie2索引

    Args:
        reference: 参考序列FASTA文件
        index_name: 索引名称（默认与reference相同）

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(reference)

    if index_name is None:
        index_name = reference

    ensure_dir(Path(index_name).parent)

    cmd = ["bowtie2-build", reference, index_name]
    return run_command(cmd)


def align_bowtie2(
    index, reads1, reads2=None, output_sam=None, threads=4, extra_args=None
):
    """
    使用Bowtie2进行序列比对

    Args:
        index: Bowtie2索引名称
        reads1: 第一个read文件
        reads2: 第二个read文件（可选）
        output_sam: 输出SAM文件路径
        threads: 线程数
        extra_args: 额外参数（列表）

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(reads1)

    if reads2:
        check_file_exists(reads2)

    cmd = ["bowtie2", "-p", str(threads)]

    if extra_args:
        cmd.extend(extra_args)

    cmd.extend(
        ["-x", index, "-U", reads1]
        if not reads2
        else ["-x", index, "-1", reads1, "-2", reads2]
    )

    if output_sam:
        ensure_dir(Path(output_sam).parent)
        cmd.extend(["-S", output_sam])

    return run_command(cmd, capture_output=False)
