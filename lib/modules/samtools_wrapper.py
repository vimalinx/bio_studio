"""
SAM/BAM文件处理模块 - SAMtools包装器
"""

from pathlib import Path
from .utils import run_command, ensure_dir, check_file_exists


def sam_to_bam(input_sam, output_bam=None, threads=4):
    """
    SAM转BAM格式

    Args:
        input_sam: 输入SAM文件
        output_bam: 输出BAM文件（默认与input_sam同名，扩展名改为.bam）
        threads: 线程数

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(input_sam)

    if output_bam is None:
        output_bam = Path(input_sam).with_suffix(".bam")

    ensure_dir(Path(output_bam).parent)

    cmd = ["samtools", "view", "-@", str(threads), "-bS", input_sam, "-o", output_bam]
    return run_command(cmd, capture_output=False)


def sort_bam(input_bam, output_bam=None, threads=4, by_name=False):
    """
    排序BAM文件

    Args:
        input_bam: 输入BAM文件
        output_bam: 输出BAM文件（默认添加.sorted后缀）
        threads: 线程数
        by_name: 是否按read名称排序（默认按坐标排序）

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(input_bam)

    if output_bam is None:
        output_bam = str(Path(input_bam).with_suffix(""))
        output_bam += ".sorted.bam"

    ensure_dir(Path(output_bam).parent)

    cmd = ["samtools", "sort", "-@", str(threads)]

    if by_name:
        cmd.append("-n")

    cmd.extend(["-o", output_bam, input_bam])
    return run_command(cmd, capture_output=False)


def index_bam(input_bam):
    """
    为BAM文件创建索引

    Args:
        input_bam: 输入BAM文件

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(input_bam)

    cmd = ["samtools", "index", input_bam]
    return run_command(cmd)


def filter_bam(input_bam, output_bam, filters=None, threads=4):
    """
    过滤BAM文件

    Args:
        input_bam: 输入BAM文件
        output_bam: 输出BAM文件
        filters: 过滤条件字符串（如 '-q 30 -F 4'）
        threads: 线程数

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(input_bam)
    ensure_dir(Path(output_bam).parent)

    cmd = ["samtools", "view", "-@", str(threads), "-b"]

    if filters:
        cmd.extend(filters.split())

    cmd.extend([input_bam, "-o", output_bam])
    return run_command(cmd, capture_output=False)


def get_coverage(input_bam, reference=None, min_maq=0):
    """
    获取覆盖度统计

    Args:
        input_bam: 输入BAM文件
        reference: 参考序列（可选）
        min_maq: 最小mapping quality

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(input_bam)

    cmd = ["samtools", "depth", "-a"]

    if min_maq > 0:
        cmd.extend(["-Q", str(min_maq)])

    if reference:
        cmd.extend(["-r", reference])

    cmd.append(input_bam)
    return run_command(cmd)


def flagstat(input_bam):
    """
    获取BAM文件统计信息

    Args:
        input_bam: 输入BAM文件

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(input_bam)

    cmd = ["samtools", "flagstat", input_bam]
    return run_command(cmd)
