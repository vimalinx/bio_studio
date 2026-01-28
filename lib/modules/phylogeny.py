"""
系统发育分析模块 - IQ-TREE 2包装器
"""

from pathlib import Path
from .utils import run_command, ensure_dir, check_file_exists


def build_tree(
    input_alignment,
    output_tree=None,
    model="MFP",
    bootstrap=1000,
    threads=4,
    extra_args=None,
):
    """
    使用IQ-TREE 2构建系统发育树

    Args:
        input_alignment: 输入比对文件（FASTA/PHYLIP/NEXUS）
        output_tree: 输出树文件（默认.treefile）
        model: 替换模型 ('MFP'自动选择, 'GTR+G', 'HKY+G'等)
        bootstrap: Bootstrap重复次数（0表示不计算）
        threads: 线程数
        extra_args: 额外参数

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(input_alignment)

    cmd = ["iqtree2", "-s", input_alignment, "-m", model, "-T", str(threads)]

    if bootstrap > 0:
        cmd.extend(["-b", str(bootstrap)])

    if extra_args:
        cmd.extend(extra_args)

    return run_command(cmd, capture_output=False)


def model_test(input_alignment, models=None, threads=4):
    """
    测试替换模型

    Args:
        input_alignment: 输入比对文件
        models: 模型列表（如 ['GTR', 'HKY', 'JC']）
        threads: 线程数

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(input_alignment)

    cmd = ["iqtree2", "-s", input_alignment, "-mset"]

    if models:
        cmd.append(",".join(models))
    else:
        cmd.append("ALL")

    cmd.extend(["-mtest", "-T", str(threads)])

    return run_command(cmd, capture_output=False)


def ancestral_reconstruction(input_alignment, tree_file, threads=4):
    """
    祖先序列重建

    Args:
        input_alignment: 输入比对文件
        tree_file: 输入树文件
        threads: 线程数

    Returns:
        subprocess.CompletedProcess对象
    """
    check_file_exists(input_alignment)
    check_file_exists(tree_file)

    cmd = [
        "iqtree2",
        "-s",
        input_alignment,
        "-t",
        tree_file,
        "-asr",
        "-T",
        str(threads),
    ]

    return run_command(cmd, capture_output=False)
