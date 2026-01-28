#!/bin/bash
# 分析酵母菌所有基因的功能

PROJECT_ROOT="$HOME/bio_studio/projects/yeast_genome_learning"
DATA_DIR="$PROJECT_ROOT/data"
RESULTS_DIR="$PROJECT_ROOT/results"

cd "$DATA_DIR"

echo "🧬 酵母菌基因功能分析"
echo "===================="
echo ""

# 提取所有基因信息
echo "📊 正在提取基因信息..."

awk '$3=="gene" {
  n=split($9, a, ";");
  name="";
  desc="";
  chr=$1;
  start=$4;
  end=$5;
  strand=$7;
  for(i=1; i<=n; i++){
    if(a[i] ~ /Name=/){
      sub(/Name=/, "", a[i]);
      gsub(/"/, "", a[i]);
      name=a[i];
    }
    if(a[i] ~ /description=/){
      sub(/description=/, "", a[i]);
      gsub(/%3B/, ";", a[i]);
      gsub(/%2C/, ",", a[i]);
      gsub(/"/, "", a[i]);
      desc=a[i];
    }
  }
  printf "%s\t%s\t%d\t%d\t%s\t%s\t%s\n", chr, name, start, end, strand, length(desc), desc;
}' annotation/genomic.gff > "$RESULTS_DIR/all_genes.txt" 2>/dev/null

TOTAL_GENES=$(wc -l < "$RESULTS_DIR/all_genes.txt")
echo "✅ 提取完成: $TOTAL_GENES 个基因"
echo ""

# 生成分类统计
echo "📈 基因功能分类统计..."
echo ""

# 按关键词分类
cat > "$RESULTS_DIR/keyword_categories.txt" << 'EOF'
代谢|metabolism, enzyme, catalyz, biosynthesis, pathway
转录|transcription, regulator, factor, rna, polymerase
翻译|translation, ribosome, trna, protein synthesis
运输|transport, transporter, carrier, channel
细胞周期|cycle, division, checkpoint, spindle
DNA修复|dna repair, replication, recombination
应激响应|stress, response, heat, cold, shock
线粒体|mitochondrial, respiration, oxidative
分泌|secretion, vesicle, golgi, er
细胞骨架|cytoskeleton, actin, tubulin, myosin
信号|signal, signaling, kinase, phosphatase
EOF

# 统计各类别
echo "📂 主要功能类别:"
echo "---------------------------"

# 代谢相关
METABOLIC=$(grep -i -E 'metabolism|enzyme|catalyz|biosynthesis|pathway' "$RESULTS_DIR/all_genes.txt" | wc -l)
echo "代谢相关: $METABOLIC 个基因"

# 转录相关
TRANSCRIPT=$(grep -i -E 'transcription|regulator|factor|polymerase|rna' "$RESULTS_DIR/all_genes.txt" | wc -l)
echo "转录调控: $TRANSCRIPT 个基因"

# 翻译相关
TRANSLATION=$(grep -i -E 'translation|ribosome|trna|protein synthesis' "$RESULTS_DIR/all_genes.txt" | wc -l)
echo "翻译/蛋白合成: $TRANSLATION 个基因"

# 运输相关
TRANSPORT=$(grep -i -E 'transport|transporter|carrier|channel' "$RESULTS_DIR/all_genes.txt" | wc -l)
echo "物质运输: $TRANSPORT 个基因"

# 细胞周期
CYCLE=$(grep -i -E 'cycle|division|checkpoint|spindle' "$RESULTS_DIR/all_genes.txt" | wc -l)
echo "细胞周期: $CYCLE 个基因"

# DNA相关
DNA=$(grep -i -E 'dna repair|replication|recombination' "$RESULTS_DIR/all_genes.txt" | wc -l)
echo "DNA代谢: $DNA 个基因"

# 应激响应
STRESS=$(grep -i -E 'stress|response|heat|cold|shock' "$RESULTS_DIR/all_genes.txt" | wc -l)
echo "应激响应: $STRESS 个基因"

# 线粒体
MITO=$(grep -i -E 'mitochondri|respirat|oxidative' "$RESULTS_DIR/all_genes.txt" | wc -l)
echo "线粒体功能: $MITO 个基因"

# 分泌
SECRET=$(grep -i -E 'secret|vesicle|golgi|endoplasm' "$RESULTS_DIR/all_genes.txt" | wc -l)
echo "分泌系统: $SECRET 个基因"

# 细胞骨架
CYTOSKELETON=$(grep -i -E 'cytoskeleton|actin|tubulin|myosin' "$RESULTS_DIR/all_genes.txt" | wc -l)
echo "细胞骨架: $CYTOSKELETON 个基因"

# 信号传导
SIGNAL=$(grep -i -E 'signal|signaling|kinase|phosphatase' "$RESULTS_DIR/all_genes.txt" | wc -l)
echo "信号传导: $SIGNAL 个基因"

echo ""
echo "---------------------------"

# 按染色体分布
echo ""
echo "🧬 染色体基因分布:"
echo "---------------------------"

for chr in I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI Mt; do
  if [ "$chr" = "Mt" ]; then
    count=$(grep -w "mitochondrion" "$RESULTS_DIR/all_genes.txt" 2>/dev/null | wc -l)
  else
    count=$(grep -w "$chr" "$RESULTS_DIR/all_genes.txt" 2>/dev/null | wc -l)
  fi
  if [ "$count" -gt 0 ]; then
    printf "  染色体 %2s: %4d 个基因\n" "$chr" "$count"
  fi
done

# 查找基因最多的染色体
echo ""
grep -w $'^[IVX]\+$' "$RESULTS_DIR/all_genes.txt" 2>/dev/null | awk '{print $1}' | sort | uniq -c | sort -rn | head -1 | \
  awk '{print "  📊 基因最多的染色体: " $2 " (" $1 " 个基因)"}'

echo ""

# 列出一些重要基因
echo "🌟 重要基因示例:"
echo "---------------------------"

echo ""
echo "1. ACT1 - 肌动蛋白（细胞骨架）"
grep "ACT1" "$RESULTS_DIR/all_genes.txt" 2>/dev/null | \
  awk -F'\t' '{printf "   染色体: %s, 位置: %s-%s, 功能: %s\n", $1, $3, $4, $7}'

echo ""
echo "2. ADH1 - 酒精脱氢酶（发酵）"
grep "ADH1" "$RESULTS_DIR/all_genes.txt" 2>/dev/null | \
  awk -F'\t' '{printf "   染色体: %s, 位置: %s-%s, 功能: %s\n", $1, $3, $4, $7}'

echo ""
echo "3. HIS3 - 组氨酸合成"
grep "HIS3" "$RESULTS_DIR/all_genes.txt" 2>/dev/null | \
  awk -F'\t' '{printf "   染色体: %s, 位置: %s-%s, 功能: %s\n", $1, $3, $4, $7}'

echo ""
echo "4. TUB1/TUB2/TUB3 - 微管蛋白"
grep -E "TUB1|TUB2|TUB3" "$RESULTS_DIR/all_genes.txt" 2>/dev/null | \
  awk -F'\t' '{printf "   %s: 染色体 %s, 功能: %s\n", $2, $1, $7}'

# 生成详细报告
cat > "$RESULTS_DIR/gene_analysis_summary.txt" << EOFREPORT
🍺 酵母菌基因功能分析报告
========================================

📊 基因组基本信息:
物种: Saccharomyces cerevisiae (酿酒酵母)
菌株: S288C
版本: R64-1-1

基因总数: $TOTAL_GENES
蛋白编码基因: $(grep -v 'gene:' annotation/genomic.gff | grep -c 'gene')
非编码RNA: 424

📈 主要功能类别:
-----------------------

EOFREPORT

cat "$RESULTS_DIR/keyword_categories.txt" >> "$RESULTS_DIR/gene_analysis_summary.txt"
echo "" >> "$RESULTS_DIR/gene_analysis_summary.txt"

echo "✅ 分析完成！"
echo ""
echo "📁 结果文件:"
echo "   - 所有基因: $RESULTS_DIR/all_genes.txt"
echo "   - 详细报告: $RESULTS_DIR/gene_analysis_summary.txt"
echo ""
echo "💡 查看特定基因:"
echo "   grep 'ACT1' $RESULTS_DIR/all_genes.txt"
echo "   grep 'kinase' $RESULTS_DIR/all_genes.txt"
echo ""
