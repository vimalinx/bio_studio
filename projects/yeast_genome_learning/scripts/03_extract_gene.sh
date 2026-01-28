#!/bin/bash
# 03_extract_gene.sh - 提取特定基因序列
# 修复版: 正确处理 SGD_features.tab 列和 ID 映射

set -e

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="$PROJECT_ROOT/data"
RESULTS_DIR="$PROJECT_ROOT/results"
mkdir -p "$RESULTS_DIR"

to_roman() {
    case $1 in
        1) echo "I" ;;
        2) echo "II" ;;
        3) echo "III" ;;
        4) echo "IV" ;;
        5) echo "V" ;;
        6) echo "VI" ;;
        7) echo "VII" ;;
        8) echo "VIII" ;;
        9) echo "IX" ;;
        10) echo "X" ;;
        11) echo "XI" ;;
        12) echo "XII" ;;
        13) echo "XIII" ;;
        14) echo "XIV" ;;
        15) echo "XV" ;;
        16) echo "XVI" ;;
        17) echo "Mito" ;;
        *) echo "$1" ;;
    esac
}

if [ -z "$1" ]; then
    echo "用法: bash $0 <基因名>"
    echo "示例: bash $0 ACT1"
    exit 1
fi

GENE_NAME="$1"
echo "🔍 提取基因: $GENE_NAME"

# 1. 查找基因信息 (只匹配标准名或系统名)
# SGD_features.tab 列: 4=系统名, 5=标准名
# 使用 awk 精确匹配列
GENE_INFO=$(awk -F'\t' -v q="$GENE_NAME" '$4==q || $5==q {print $0; exit}' "$DATA_DIR/annotation/SGD_features.tab")

if [ -z "$GENE_INFO" ]; then
    echo "❌ 未找到基因: $GENE_NAME"
    exit 1
fi

# 解析列
# 4:Syst, 5:Std, 9:Chrom, 10:Start, 11:End, 12:Strand
SYST_NAME=$(echo "$GENE_INFO" | cut -f4)
STD_NAME=$(echo "$GENE_INFO" | cut -f5)
CHROM_NUM=$(echo "$GENE_INFO" | cut -f9)
START=$(echo "$GENE_INFO" | cut -f10)
END=$(echo "$GENE_INFO" | cut -f11)
STRAND=$(echo "$GENE_INFO" | cut -f12)

# 确保 Start < End (samtools 要求)
if [ "$START" -gt "$END" ]; then
    TEMP=$START
    START=$END
    END=$TEMP
fi

# 转换染色体 ID
CHROM_ID=$(to_roman "$CHROM_NUM")

echo "📋 基因定位:"
echo "   名称: ${STD_NAME:-$SYST_NAME} ($SYST_NAME)"
echo "   位置: 染色体 $CHROM_ID : $START - $END ($STRAND 链)"
echo ""

# 2. 提取 DNA 序列
echo "📥 提取 DNA 序列..."
if command -v samtools &> /dev/null; then
    # samtools faidx 需要 ID:Start-End
    # 注意: 如果是负链，samtools 默认提取正链序列。如果要反向互补，需处理。
    # 这里我们只提取基因组对应区域的序列。
    
    REGION="$CHROM_ID:$START-$END"
    OUT_FILE="$RESULTS_DIR/${GENE_NAME}.fna"
    
    samtools faidx "$DATA_DIR/sequence/genomic.fna" "$REGION" > "$OUT_FILE"
    
    if [ -s "$OUT_FILE" ]; then
        echo "   ✅ DNA 序列已保存: $OUT_FILE"
        seqkit stats "$OUT_FILE" 2>/dev/null | tail -n 1
        
        # 预览
        grep -v "^>" "$OUT_FILE" | head -n 1 | cut -c1-50 | awk '{print "   预览: " $0 "..."}'
    else
        echo "   ❌ 提取失败 (可能是索引问题或 ID 不匹配)"
        echo "   尝试使用的 Region ID: $REGION"
    fi
else
    echo "   ⚠️ samtools 未安装"
fi
echo ""

# 3. 提取蛋白质序列
echo "📥 提取蛋白质序列..."
# 蛋白文件通常用系统名作为 ID
PROT_FILE="$RESULTS_DIR/${GENE_NAME}.faa"
if grep -q "$SYST_NAME" "$DATA_DIR/proteins/protein.faa"; then
    # 使用 seqkit grep 更加稳健，或者简单的 awk/grep
    # 这里用 seqkit 如果有，否则用 awk
    if command -v seqkit &> /dev/null; then
        seqkit grep -p "$SYST_NAME" "$DATA_DIR/proteins/protein.faa" > "$PROT_FILE"
    else
        # 简单的 awk 提取 fasta 记录
        awk -v RS=">" -v q="$SYST_NAME" '$0 ~ q {print ">"$0}' "$DATA_DIR/proteins/protein.faa" > "$PROT_FILE"
    fi
    
    echo "   ✅ 蛋白序列已保存: $PROT_FILE"
else
    echo "   ⚠️ 未找到对应蛋白序列 ($SYST_NAME)"
fi

echo ""
echo "✅ 完成"
