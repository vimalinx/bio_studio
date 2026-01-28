"""
通用工具函数模块
"""

from pathlib import Path
import subprocess
import sys


def run_command(cmd, check=True, capture_output=True, **kwargs):
    """
    执行Shell命令并返回结果

    Args:
        cmd: 命令字符串或列表
        check: 如果为True，命令失败时抛出异常
        capture_output: 是否捕获输出

    Returns:
        subprocess.CompletedProcess对象
    """
    if isinstance(cmd, str):
        cmd = cmd.split()

    try:
        result = subprocess.run(
            cmd, check=check, capture_output=capture_output, text=True, **kwargs
        )
        return result
    except subprocess.CalledProcessError as e:
        print(f"命令执行失败: {' '.join(e.cmd)}")
        print(f"错误信息: {e.stderr}")
        raise


def check_file_exists(filepath, raise_error=True):
    """
    检查文件是否存在

    Args:
        filepath: 文件路径
        raise_error: 如果为True，文件不存在时抛出异常

    Returns:
        bool: 文件是否存在
    """
    path = Path(filepath)
    if not path.exists():
        if raise_error:
            raise FileNotFoundError(f"文件不存在: {filepath}")
        return False
    return True


def ensure_dir(dirpath):
    """
    确保目录存在，不存在则创建

    Args:
        dirpath: 目录路径

    Returns:
        Path: 目录路径对象
    """
    path = Path(dirpath)
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_file_extension(filepath):
    """
    获取文件扩展名（小写）

    Args:
        filepath: 文件路径

    Returns:
        str: 文件扩展名（包含点）
    """
    return Path(filepath).suffix.lower()


def is_compressed(filepath):
    """
    检查文件是否是压缩文件

    Args:
        filepath: 文件路径

    Returns:
        bool: 是否是压缩文件
    """
    compressed_extensions = [".gz", ".bz2", ".zip", ".rar", ".7z"]
    ext = get_file_extension(filepath)
    return ext in compressed_extensions


def parse_fasta(filepath):
    """
    简单的FASTA文件解析器

    Args:
        filepath: FASTA文件路径

    Returns:
        dict: {sequence_id: sequence}
    """
    sequences = {}
    current_id = None
    current_seq = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_id is not None:
            sequences[current_id] = "".join(current_seq)

    return sequences


def write_fasta(sequences, filepath):
    """
    写入FASTA文件

    Args:
        sequences: dict {sequence_id: sequence} 或 SeqRecord列表
        filepath: 输出文件路径
    """
    ensure_dir(Path(filepath).parent)

    if isinstance(sequences, dict):
        with open(filepath, "w") as f:
            for seq_id, seq in sequences.items():
                f.write(f">{seq_id}\n")
                for i in range(0, len(seq), 60):
                    f.write(seq[i : i + 60] + "\n")
    else:
        from Bio import SeqIO

        SeqIO.write(sequences, filepath, "fasta")


def calculate_gc_content(sequence):
    """
    计算GC含量

    Args:
        sequence: 序列字符串

    Returns:
        float: GC百分比
    """
    sequence = sequence.upper()
    gc_count = sequence.count("G") + sequence.count("C")
    if len(sequence) == 0:
        return 0.0
    return (gc_count / len(sequence)) * 100


def reverse_complement(dna_seq):
    """
    计算DNA反向互补序列

    Args:
        dna_seq: DNA序列

    Returns:
        str: 反向互补序列
    """
    complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "a": "t",
        "t": "a",
        "g": "c",
        "c": "g",
        "N": "N",
        "n": "n",
    }
    rev_comp = "".join(complement.get(base, "N") for base in reversed(dna_seq))
    return rev_comp


def validate_sequence(sequence, seq_type="dna"):
    """
    验证序列是否有效

    Args:
        sequence: 序列字符串
        seq_type: 序列类型 ('dna', 'rna', 'protein')

    Returns:
        tuple: (is_valid, error_message)
    """
    valid_chars = {
        "dna": set("ATGCN"),
        "rna": set("AUGCN"),
        "protein": set("ACDEFGHIKLMNPQRSTVWY*-"),
    }

    allowed = valid_chars.get(seq_type.lower(), set())
    sequence = sequence.upper()

    invalid = [base for base in sequence if base not in allowed]

    if invalid:
        return False, f"发现非法字符: {set(invalid)}"

    return True, ""
