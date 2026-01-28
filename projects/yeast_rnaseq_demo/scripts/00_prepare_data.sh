#!/bin/bash
# 00_prepare_data.sh - 准备/模拟 RNA-seq 数据
# 我们将模拟两个样本: WT (野生型) 和 MUT (突变体)
# 使用 seqkit sample 从基因组中随机抽取 reads 来模拟 (快速且无需联网下载大文件)

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="$PROJECT_ROOT/data"
RAW_DIR="$DATA_DIR/raw"
REF_DIR="$DATA_DIR/references"

mkdir -p "$RAW_DIR"

echo "🧬 准备酵母菌 RNA-seq 模拟数据"
echo "================================="

if ! command -v seqkit &> /dev/null; then
    echo "❌ seqkit 未安装"
    exit 1
fi

# 检查参考基因组
if [ ! -f "$REF_DIR/genome.fa" ]; then
    echo "❌ 参考基因组未找到: $REF_DIR/genome.fa"
    exit 1
fi

# 模拟 Reads 生成函数
# 简单地从基因组切片模拟 SE (Single End) reads
# 注意: 这只是为了跑通流程的"假"数据，不具备真实生物学表达量分布
simulate_reads() {
    sample_name=$1
    seed=$2
    num_reads=10000  # 1万条 reads 用于快速测试
    
    echo "   正在生成样本: $sample_name ($num_reads reads)..."
    
    # 使用 seqkit sliding 生成碎片，然后随机抽样
    # 模拟 100bp reads
    seqkit sliding -s 50 -W 100 "$REF_DIR/genome.fa" 2>/dev/null | \
        seqkit sample -n "$num_reads" -s "$seed" 2>/dev/null | \
        seqkit seq -n -g 2>/dev/null | \
        awk -v n="$sample_name" '{print "@"n"_"NR"\n"$0"\n+\n"gensub(/./, "I", "g", $0)}' | \
        gzip > "$RAW_DIR/${sample_name}_R1.fastq.gz"
        
    echo "   ✅ 已生成: $RAW_DIR/${sample_name}_R1.fastq.gz"
}

echo "1. 生成模拟测序数据"
simulate_reads "WT_1" 100
simulate_reads "WT_2" 101
simulate_reads "MUT_1" 200
simulate_reads "MUT_2" 201

echo ""
echo "2. 数据清单"
ls -lh "$RAW_DIR"

echo ""
echo "✅ 数据准备完成"
