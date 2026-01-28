#!/bin/bash
# 酵母菌基因组数据库安装脚本
# 按照 Bio Studio 规范：项目隔离，原子化模块

set -e

# 获取项目根目录
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="$PROJECT_ROOT/data"
LOG_DIR="$PROJECT_ROOT/logs"

echo "🍺 酵母菌基因组数据库安装"
echo "================================"
echo "项目根目录: $PROJECT_ROOT"
echo "数据目录: $DATA_DIR"
echo ""

# 创建目录
mkdir -p "$DATA_DIR"/{sequence,annotation,proteins,transcripts,blastdb}
mkdir -p "$LOG_DIR"

# ========================================
# 1. 下载参考基因组
# ========================================
echo "1️⃣  下载参考基因组 (R64-1-1)..."
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/reference/GCF_000146045.2_R64-1-1/genomic/GCF_000146045.2_R64-1-1_genomic.fna.gz -O "$DATA_DIR/sequence/genomic.fna.gz" 2>&1 | tee -a "$LOG_DIR/download.log"

echo "   解压中..."
gunzip -kf "$DATA_DIR/sequence/genomic.fna.gz"

echo "   创建索引..."
cd "$DATA_DIR/sequence"
if ! command -v samtools &> /dev/null; then
    echo "   ⚠️  samtools 未安装，跳过索引创建"
else
    samtools faidx genomic.fna
fi
cd "$PROJECT_ROOT"

echo "   ✅ 基因组序列完成"
echo ""

# ========================================
# 2. 下载基因注释
# ========================================
echo "2️⃣  下载基因注释 (GFF)..."
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/reference/GCF_000146045.2_R64-1-1/genomic/GCF_000146045.2_R64-1-1_genomic.gff.gz -O "$DATA_DIR/annotation/genomic.gff.gz" 2>&1 | tee -a "$LOG_DIR/download.log"

echo "   解压中..."
gunzip -kf "$DATA_DIR/annotation/genomic.gff.gz"

echo "   ✅ 基因注释完成"
echo ""

# ========================================
# 3. 下载蛋白质序列
# ========================================
echo "3️⃣  下载蛋白质序列..."
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Saccharomyces_cerevisiae/reference/GCF_000146045.2_R64-1-1/genomic/GCF_000146045.2_R64-1-1_protein.faa.gz -O "$DATA_DIR/proteins/protein.faa.gz" 2>&1 | tee -a "$LOG_DIR/download.log"

echo "   解压中..."
gunzip -kf "$DATA_DIR/proteins/protein.faa.gz"

echo "   ✅ 蛋白质序列完成"
echo ""

# ========================================
# 4. 下载 SGD 功能注释
# ========================================
echo "4️⃣  下载 SGD 功能注释..."
wget -c https://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff -O "$DATA_DIR/annotation/SGD.gff" 2>&1 | tee -a "$LOG_DIR/download.log"
wget -c https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab -O "$DATA_DIR/annotation/SGD_features.tab" 2>&1 | tee -a "$LOG_DIR/download.log"

echo "   ✅ SGD 注释完成"
echo ""

# ========================================
# 5. 创建 BLAST 数据库
# ========================================
echo "5️⃣  创建 BLAST 数据库..."
cd "$DATA_DIR"

if ! command -v makeblastdb &> /dev/null; then
    echo "   ⚠️  makeblastdb 未安装，跳过 BLAST 数据库创建"
else
    echo "   创建核苷酸数据库..."
    makeblastdb -in sequence/genomic.fna -dbtype nucl -parse_seqids \
               -title "Saccharomyces cerevisiae R64-1-1" \
               -out blastdb/yeast_genome 2>&1 | tee -a "$LOG_DIR/blast.log"

    echo "   创建蛋白质数据库..."
    makeblastdb -in proteins/protein.faa -dbtype prot -parse_seqids \
               -title "Saccharomyces cerevisiae proteins" \
               -out blastdb/yeast_protein 2>&1 | tee -a "$LOG_DIR/blast.log"

    echo "   ✅ BLAST 数据库完成"
fi

cd "$PROJECT_ROOT"
echo ""

# ========================================
# 完成信息
# ========================================
echo "================================"
echo "✅ 酵母菌数据库安装完成！"
echo "================================"
echo ""
echo "📊 数据统计:"
du -sh "$DATA_DIR"/* 2>/dev/null
echo ""
echo "📁 数据位置: $DATA_DIR"
echo "📝 日志位置: $LOG_DIR"
echo ""
echo "💡 下一步:"
echo "   bash scripts/02_verify_install.sh"
echo ""
