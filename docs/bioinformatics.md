# 生物信息学基础技能

---

## 🧬 核心概念

### 什么是生物信息学？

生物信息学 = 生物学 + 计算机科学 + 统计学

**目标**: 从生物数据中提取生物学洞察

### 数据类型

| 数据类型 | 描述 | 常用格式 |
|---------|------|---------|
| **DNA序列** | 基因组、基因 | FASTA, FASTQ, GENBANK |
| **RNA序列** | 转录组、非编码RNA | FASTA, BAM |
| **蛋白质序列** | 氨基酸序列 | FASTA, PDB |
| **结构数据** | 3D结构 | PDB, MMCIF |
| **表达数据** | 基因表达 | CSV, TSV, HDF5 |
| **变异数据** | SNP、Indel | VCF, BCF |

---

## 🔄 标准工作流

### 1. 序列获取

```bash
# 从NCBI下载
efetch -db nucleotide -id NM_000546 -format fasta > p53.fa

# 从Ensembl下载
wget http://ensembl.org/homo_sapiens/Fasta

# 使用Biopython
from Bio import Entrez
Entrez.email = "your@email.com"
handle = Entrez.efetch(db="nucleotide", id="NM_000546", rettype="fasta")
```

### 2. 序列比对

**为什么比对？**
- 找相似序列
- 推断功能
- 研究进化

**工具选择**:

| 任务 | 工具 | 命令 |
|------|------|------|
| 局部比对 | BLAST | `blastn` |
| 全局比对 | Needle | `needle` |
| 多序列比对 | MAFFT | `mafft` |
| 短序列比对 | Bowtie2 | `bowtie2` |
| 基因组比对 | BWA | `bwa mem` |

**示例**:
```bash
# BLAST搜索
blastn -query gene.fa -db nt -out results.txt -evalue 1e-5

# 多序列比对
mafft input.fa > aligned.fa
```

### 3. 基因预测

**目标**: 从DNA序列中找到基因

**工具**:
- **Prodigal**: 原核生物
- **Glimmer**: 原核生物
- **Augustus**: 真核生物
- **GENSCAN**: 真核脊椎动物

```bash
# 原核生物基因预测
prodigal -i genome.fa -a proteins.fa -d genes.fa -f gff

# 真核生物基因预测
augustus --species=human genome.fa > predictions.gff
```

### 4. 功能注释

**数据库**:
- **GO** (Gene Ontology): 分子功能、生物过程、细胞组分
- **KEGG**: 代谢通路
- **Pfam**: 蛋白质结构域
- **InterPro**: 蛋白质家族

```bash
# BLAST2GO
# InterProScan
interproscan.sh -i proteins.fa -f tsv -o annotations.tsv
```

---

## 📊 统计概念

### p值和FDR

- **p值**: 假阳性的概率
- **FDR** (False Discovery Rate): 错误发现率
- **Benjamini-Hochberg**: 常用多重检验校正

### 富集分析

**目的**: 找过表达的基因/通路

```r
# DESeq2 (R)
dds <- DESeqDataSetFromMatrix(countData, colData, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
```

---

## 🧪 常用文件格式

### FASTA格式
```
>sequence_name description
ATGCGTACGTACGTACGT...
ATGCGTACGTACGTACGT...
```

### FASTQ格式 (测序数据)
```
@read_id
ATGCGTACGT
+
IIIIIIIII
```

### SAM/BAM格式 (比对)
```
@HD VN:1.0 SO:unsorted
r001 99 ref1 7 60 8M2I4M1D3M = 37 39 TTAGATAAAGGATACTG *
```

### VCF格式 (变异)
```
##fileformat=VCFv4.2
#CHROM POS ID REF ALT QUAL FILTER INFO
20 14370 . G A 29 PASS NS=3;DP=14
```

### GFF/GTF格式 (注释)
```
##gff-version 3
chr1  ensembl  gene  1000  2000  .  +  .  ID=gene1
```

---

## 🛠️ 常用操作

### 序列处理

```python
from Bio import SeqIO
from Bio.Seq import Seq

# 读取序列
record = SeqIO.read("file.fa", "fasta")
sequence = record.seq

# 反向互补
rev_comp = sequence.reverse_complement()

# 翻译
protein = sequence.translate()

# GC含量
from Bio.SeqUtils import GC
gc_content = GC(sequence)

# 分子量
from Bio.SeqUtils import molecular_weight
mw = molecular_weight(sequence)
```

### 批量处理

```bash
# 批量BLAST
cat *.fa | blastx -db swissprot -out results.txt

# 使用xargs
ls *.fa | xargs -I {} blastn -query {} -db nt
```

---

## 🎯 质量控制

### 序列质量

```bash
# 检查测序质量
FastQC reads.fastq

# 过滤低质量
fastp -i reads.fastq -o clean_reads.fastq
```

### 比对质量

```bash
# 统计比对率
samtools flagstat aligned.bam

# 覆盖度
samtools depth aligned.bam > coverage.txt
```

---

## 📈 可视化

### 基因组浏览器

- **IGV** (Integrative Genomics Viewer)
- **JBrowse**
- **UCSC Genome Browser**

### 蛋白质可视化

- **PyMOL**: 3D结构
- **ChimeraX**: 高级可视化
- **NGL Viewer**: 网页版

### 数据可视化

```python
import matplotlib.pyplot as plt
import seaborn as sns

# GC含量分布
plt.hist(gc_contents, bins=50)
plt.xlabel('GC %')
plt.ylabel('Frequency')

# 热图
sns.heatmap(expression_matrix)
```

---

## 🔍 故障排除

### 常见问题

**Q: BLAST太慢？**
A: 使用本地数据库或`-num_threads`参数

**Q: 内存不足？**
A: 使用流式处理或降低线程数

**Q: 格式不兼容？**
A: 使用`seqtk`或BioPython转换

---

## 📚 推荐资源

### 书籍
- "Bioinformatics Algorithms" by Compeau & Pevzner
- "Bioinformatics Data Skills" by Vince Buffalo
- "Practical Computing for Biologists"

### 在线资源
- Rosalind: 生信编程练习
- Coursera Bioinformatics Specialization
- NCBI Learning Resources

---

**最后更新**: 2025-01-18
