#!/bin/bash
# 04_sequence_analysis.sh - 序列深度分析
# 使用 seqkit 进行统计和分析

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="$PROJECT_ROOT/data"
RESULTS_DIR="$PROJECT_ROOT/results"
mkdir -p "$RESULTS_DIR"

cd "$DATA_DIR" || exit 1

echo "🧬 酵母菌序列分析 (使用 seqkit)"
echo "================================="

if ! command -v seqkit &> /dev/null; then
    echo "❌ seqkit 未安装，请先安装: conda install -c bioconda seqkit"
    exit 1
fi

echo "1. 总体统计 (seqkit stats)"
echo "--------------------------"
seqkit stats sequence/genomic.fna proteins/protein.faa -a
echo ""

echo "2. 染色体 GC 含量分析"
echo "--------------------------"
echo "正在计算每条染色体的 GC 含量..."
# 提取 ID, 长度, GC%
seqkit fx2tab --name --length --gc sequence/genomic.fna | sort -k2,2nr > "$RESULTS_DIR/chromosome_stats.txt"
head -n 5 "$RESULTS_DIR/chromosome_stats.txt"
echo "... 完整结果已保存至 $RESULTS_DIR/chromosome_stats.txt"
echo ""

echo "3. 蛋白质长度分布"
echo "--------------------------"
echo "Top 10 最长蛋白质:"
seqkit fx2tab --name --length proteins/protein.faa | sort -k2,2nr | head -n 10
echo ""

echo "4. 寻找富含半胱氨酸(Cys)的蛋白 (可能涉及二硫键)"
echo "--------------------------"
# 计算每个蛋白中 C 的数量和比例
seqkit fx2tab -H -n -l -C C proteins/protein.faa | \
    sort -k4,4nr | head -n 5
echo ""

echo "✅ 分析完成"
