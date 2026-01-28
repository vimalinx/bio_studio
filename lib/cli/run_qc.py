#!/usr/bin/env python3
"""
快速QC命令行工具
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from modules.qc import run_fastqc
from modules.utils import check_file_exists


def main():
    import argparse

    parser = argparse.ArgumentParser(description="快速质量控制")
    parser.add_argument("input_files", nargs="+", help="输入文件")
    parser.add_argument("-o", "--output", required=True, help="输出目录")
    parser.add_argument("-t", "--threads", type=int, default=4, help="线程数")

    args = parser.parse_args()

    for f in args.input_files:
        check_file_exists(f)

    result = run_fastqc(args.input_files, args.output, args.threads)

    if result.returncode == 0:
        print(f"FastQC完成: {args.output}")
    else:
        print(f"FastQC失败: {result.stderr}")
        sys.exit(1)


if __name__ == "__main__":
    main()
