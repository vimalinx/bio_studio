"""
质量控制模块 - FastQC包装器
"""

from pathlib import Path
from .utils import run_command, ensure_dir, check_file_exists


def run_fastqc(input_files, output_dir, threads=4, extra_args=None):
    """
    运行FastQC进行质量控制

    Args:
        input_files: 输入文件（字符串或列表）
        output_dir: 输出目录
        threads: 线程数
        extra_args: 额外参数（列表）

    Returns:
        subprocess.CompletedProcess对象
    """
    if isinstance(input_files, str):
        input_files = [input_files]

    for f in input_files:
        check_file_exists(f)

    ensure_dir(output_dir)

    cmd = ["fastqc", "-t", str(threads), "-o", output_dir] + input_files

    if extra_args:
        cmd.extend(extra_args)

    return run_command(cmd)


def parse_fastqc_report(report_file):
    """
    解析FastQC报告（简化版）

    Args:
        report_file: FastQC HTML报告文件路径

    Returns:
        dict: 包含关键QC指标的字典
    """
    import re
    from bs4 import BeautifulSoup

    with open(report_file, "r") as f:
        soup = BeautifulSoup(f, "html.parser")

    qc_data = {}

    per_base_quality = soup.find("h2", string="Per base sequence quality")
    if per_base_quality:
        qc_data["quality_status"] = per_base_quality.find_next("div").get_text(
            strip=True
        )

    gc_content = soup.find("h2", string="Per sequence GC content")
    if gc_content:
        gc_text = gc_content.find_next("div").get_text(strip=True)
        gc_match = re.search(r"(\d+)%", gc_text)
        if gc_match:
            qc_data["gc_content"] = int(gc_match.group(1))

    total_sequences = soup.find("h2", string="Basic Statistics")
    if total_sequences:
        table = total_sequences.find_next("table")
        rows = table.find_all("tr")
        for row in rows:
            cells = row.find_all("td")
            if len(cells) >= 2:
                label = cells[0].get_text(strip=True)
                value = cells[1].get_text(strip=True)
                qc_data[label] = value

    return qc_data
