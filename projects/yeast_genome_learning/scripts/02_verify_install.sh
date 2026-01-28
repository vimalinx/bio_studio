#!/bin/bash
# 验证酵母菌数据库安装

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="$PROJECT_ROOT/data"

echo "🔍 验证酵母菌数据库安装..."
echo ""

if [ ! -d "$DATA_DIR" ]; then
    echo "❌ 数据目录不存在: $DATA_DIR"
    echo "请先运行: bash scripts/01_setup_database.sh"
    exit 1
fi

cd "$DATA_DIR"

# 检查必需文件
echo "📋 检查必需文件..."
required_files=(
    "sequence/genomic.fna"
    "annotation/genomic.gff"
    "proteins/protein.faa"
    "annotation/SGD_features.tab"
)

all_ok=true
for file in "${required_files[@]}"; do
    if [ -f "$file" ]; then
        size=$(ls -lh "$file" | awk '{print $5}')
        echo "   ✅ $file ($size)"
    else
        echo "   ❌ $file (缺失)"
        all_ok=false
    fi
done

echo ""

if [ "$all_ok" = false ]; then
    echo "❌ 安装不完整，请重新运行安装脚本"
    exit 1
fi

# 统计信息
echo "📊 基因组统计..."
echo ""

if command -v samtools &> /dev/null; then
    if [ ! -f "sequence/genomic.fna.fai" ]; then
        samtools faidx sequence/genomic.fna
    fi
    echo "染色体数量: $(cut -f1 sequence/genomic.fna.fai | wc -l)"
    echo "基因组大小: $(awk '{sum+=$2} END {printf "%.2f Mb\n", sum/1000000}' sequence/genomic.fna.fai)"
else
    echo "   ⚠️  samtools 未安装"
fi

echo "基因数量: $(grep -w 'gene' annotation/genomic.gff | wc -l)"
echo "蛋白数量: $(grep -c '^>' proteins/protein.faa)"
echo ""

# 查找 ACT1 基因
echo "🔍 查找 ACT1 基因（验证数据完整性）..."
if grep -q 'ACT1' annotation/SGD_features.tab; then
    echo "   ✅ 找到 ACT1 基因"
    grep 'ACT1' annotation/SGD_features.tab | head -1 | \
        awk -F'\t' '{printf "   标准名称: %s\n   系统名称: %s\n   染色体: %s\n", $3, $2, $10}'
else
    echo "   ⚠️  未找到 ACT1"
fi

echo ""

# 检查 BLAST 数据库
echo "💣 检查 BLAST 数据库..."
if [ -f "blastdb/yeast_genome.nhr" ]; then
    echo "   ✅ 核苷酸 BLAST 数据库已创建"
else
    echo "   ⚠️  BLAST 数据库未创建"
fi

echo ""
echo "================================"
echo "✅ 验证完成！数据库可以正常使用"
echo "================================"
echo ""
echo "💡 下一步:"
echo "   bash scripts/03_extract_gene.sh ACT1"
echo ""
