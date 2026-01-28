#!/bin/bash
# 05_blast_search.sh - BLAST 搜索示例

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="$PROJECT_ROOT/data"
RESULTS_DIR="$PROJECT_ROOT/results"
mkdir -p "$RESULTS_DIR"

cd "$DATA_DIR" || exit 1

echo "🔍 BLAST 搜索实战"
echo "================="

if ! command -v blastn &> /dev/null; then
    echo "❌ BLAST+ 未安装"
    exit 1
fi

# 1. 准备查询序列 (提取 ACT1 的一部分作为 query)
echo "1. 准备查询序列 (ACT1)"
# 假设我们已知 ACT1 序列，这里模拟创建一个 query 文件
# ACT1 序列片段
cat > "$RESULTS_DIR/query.fa" << EOF
>Query_ACT1_fragment
ATGGATTCTGAGGTTGCTGCTTTGGTTATTGATAACGGTTCTGGTATGTGTAAAGCCGGTTTTGCCGGTGACGATGCCCCCCGTGCCGTGTTTCCATCA
EOF
echo "   已创建查询文件: $RESULTS_DIR/query.fa"
echo ""

# 2. 运行 BLASTN
echo "2. 运行 BLASTN (核酸 vs 核酸数据库)"
echo "-----------------------------------"
# -outfmt 6 是表格格式
blastn -query "$RESULTS_DIR/query.fa" \
       -db blastdb/yeast_genome \
       -out "$RESULTS_DIR/blastn_results.txt" \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -evalue 1e-5

echo "   BLASTN 完成。前 5 条结果:"
head -n 5 "$RESULTS_DIR/blastn_results.txt"
echo ""

# 3. 结果解读
echo "3. 结果解读"
echo "-----------------------------------"
if [ -s "$RESULTS_DIR/blastn_results.txt" ]; then
    best_hit=$(head -n 1 "$RESULTS_DIR/blastn_results.txt" | awk '{print $2}')
    identity=$(head -n 1 "$RESULTS_DIR/blastn_results.txt" | awk '{print $3}')
    echo "   最佳匹配: $best_hit (一致性: $identity%)"
    echo "   这应该就是 ACT1 基因所在的染色体位置。"
else
    echo "   ⚠️ 未找到匹配结果"
fi

echo ""
echo "✅ 搜索演示完成"
echo "   结果文件: $RESULTS_DIR/blastn_results.txt"
