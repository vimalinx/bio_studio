"""
序列分析模块 - seqkit / RNAfold 包装器
"""

from __future__ import annotations

from pathlib import Path

from .utils import check_file_exists, ensure_dir, run_command


def run_seqkit_stats(input_fasta, output_txt=None, extra_args=None):
    """
    使用 seqkit stats 生成序列统计
    """
    check_file_exists(input_fasta)

    cmd = ["seqkit", "stats", str(input_fasta)]
    if extra_args:
        cmd.extend(extra_args)

    result = run_command(cmd)

    if output_txt:
        output_path = Path(output_txt)
        ensure_dir(output_path.parent)
        output_path.write_text(result.stdout or "", encoding="utf-8")

    return result


def run_rnafold(rna_sequence, sequence_id="seq", no_ps=True):
    """
    使用 RNAfold 计算单条 RNA 序列的结构和 MFE
    """
    cmd = ["RNAfold"]
    if no_ps:
        cmd.append("--noPS")

    result = run_command(
        cmd,
        input=f">{sequence_id}\n{rna_sequence}\n",
    )

    lines = [line.strip() for line in (result.stdout or "").splitlines() if line.strip()]
    best_parts: list[str] | None = None
    for line in lines:
        parts = line.split()
        if parts and parts[-1].startswith("(") and parts[-1].endswith(")"):
            best_parts = parts

    if best_parts is None:
        return {"structure": None, "mfe": None}

    structure = best_parts[0]
    mfe = None
    try:
        mfe = float(best_parts[-1].strip("()"))
    except ValueError:
        mfe = None

    return {"structure": structure, "mfe": mfe}
