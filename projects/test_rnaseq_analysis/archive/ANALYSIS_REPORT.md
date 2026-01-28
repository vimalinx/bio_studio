# SRR1972739 RNA-seq 分析报告

## 项目信息

- **项目名称**: test_rnaseq_analysis
- **分析日期**: 2026-01-22
- **任务**: 分析 SRR1972739 RNA-seq 数据，完成质控和比对到 E. coli K12 基因组，报告比对率

---

## 数据下载

### 1. 参考基因组
- **来源**: NCBI (GCF_000005845.2_ASM584v2)
- **物种**: Escherichia coli str. K-12 substr. MG1655
- **文件**: ecoli_K12.fa
- **大小**: 4.5 MB
- **行数**: 58,022
- **位置**: shared_data/references/ecoli_K12.fa
- **状态**: ✓ 已下载

### 2. 测序数据 (SRR1972739)
- **来源**: ENA (ftp.sra.ebi.ac.uk)
- **Run ID**: SRR1972739
- **样本标题**: G4252.1
- **科学名称**: Zaire ebolavirus (扎伊尔埃博拉病毒)
- **文库来源**: TRANSCRIPTOMIC
- **测序策略**: RNA-Seq
- **文件**:
  - SRR1972739_1.fastq.gz (55 MB)
  - SRR1972739_2.fastq.gz (51 MB)
- **位置**: projects/test_rnaseq_analysis/data/raw/
- **状态**: ✓ 已下载

---

## 质量控制 (FastQC)

### 执行命令
```bash
fastqc -t 4 -o qc ../raw/SRR1972739_1.fastq.gz ../raw/SRR1972739_2.fastq.gz
```

### 结果文件
- `SRR1972739_1_fastqc.html` (597 KB)
- `SRR1972739_2_fastqc.html` (576 KB)
- `SRR1972739_1_fastqc.zip` (390 KB)
- `SRR1972739_2_fastqc.zip` (383 KB)
- **位置**: projects/test_rnaseq_analysis/data/processed/qc/

### 质控状态
- ✓ FastQC 报告已生成
- HTML 报告可在浏览器中查看详细质量评估

---

## 序列比对 (Bowtie2)

### 1. 索引建立
```bash
cd shared_data/references
bowtie2-build ecoli_K12.fa ecoli_K12_index
```

**状态**: ✓ 索引已成功建立
**索引文件**:
- ecoli_K12_index.1.bt2
- ecoli_K12_index.2.bt2
- ecoli_K12_index.3.bt2
- ecoli_K12_index.4.bt2
- ecoli_K12_index.rev.1.bt2
- ecoli_K12_index.rev.2.bt2

### 2. 序列比对
```bash
bowtie2 -x /media/vimalinx/Data/bio_studio/shared_data/references/ecoli_K12_index \
        -1 ../raw/SRR1972739_1.fastq.gz \
        -2 ../raw/SRR1972739_2.fastq.gz \
        -S SRR1972739_aligned.sam -p 4
```

### 比对统计

#### Bowtie2 输出摘要
```
758337 reads; of these:
  758337 (100.00%) were paired; of these:
    758311 (100.00%) aligned concordantly 0 times
    10 (0.00%) aligned concordantly exactly 1 time
    16 (0.00%) aligned concordantly >1 times
    ----
    758311 pairs aligned concordantly 0 times; of these:
      2 (0.00%) aligned discordantly 1 time
    ----
    758309 pairs aligned 0 times concordantly or discordantly; of these:
      1516618 mates make up the pairs; of these:
        1516513 (99.99%) aligned 0 times
        14 (0.00%) aligned exactly 1 time
        91 (0.01%) aligned >1 times
0.01% overall alignment rate
```

#### SAMtools 统计 (flagstat)
```
1516674 + 0 in total (QC-passed reads + QC-failed reads)
1516674 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
161 + 0 mapped (0.01% : N/A)
161 + 0 primary mapped (0.01% : N/A)
1516674 + 0 paired in sequencing
758337 + 0 read1
758337 + 0 read2
52 + 0 properly paired (0.00% : N/A)
58 + 0 with itself and mate mapped
103 + 0 singletons (0.01% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

### 比对结果文件
- `SRR1972739_aligned.sam` (358 MB)
- `SRR1972739_aligned.bam`
- `SRR1972739_aligned.sorted.bam`
- `SRR1972739_aligned.sorted.bam.bai` (索引)
- `alignment_stats.txt` (统计结果)
- **位置**: projects/test_rnaseq_analysis/data/processed/

---

## 结果分析

### ⚠️ 重要发现：物种不匹配

**问题**: 比对率极低（0.01%），这是因为：

| 项目 | 使用的参考基因组 | 实际样本来源 |
|------|---------------|------------|
| **参考基因组** | Escherichia coli str. K-12 substr. MG1655 (大肠杆菌) | - |
| **样本 SRR1972739** | - | Zaire ebolavirus (扎伊尔埃博拉病毒) |

**原因**: SRR1972739 是埃博拉病毒的 RNA-seq 数据，而不是大肠杆菌数据。

### 比对率总结

| 指标 | 数值 | 百分比 |
|------|------|--------|
| **总 reads 数** | 1,516,674 | 100% |
| **已比对 reads** | 161 | **0.01%** |
| **未比对 reads** | 1,516,513 | 99.99% |
| **配对 reads 数** | 758,337 | - |
| **正确配对** | 52 | 0.00% |

### 分析结论

1. **数据与参考基因组不匹配**: 埃博拉病毒序列无法有效比对到大肠杆菌基因组
2. **比对率符合预期**: 对于完全不同的物种，0.01%的比对率主要是随机匹配造成的假阳性
3. **质控成功**: FastQC 报告已生成，可用于评估原始数据质量

---

## 建议与后续分析

### 如果目标是分析埃博拉病毒：

**应该使用**:
- 扎伊尔埃博拉病毒参考基因组 (例如 NC_002549.1)
- 或者从 ENA/NCBI 搜索埃博拉病毒的完整基因组

**建议命令**:
```bash
# 下载埃博拉病毒参考基因组
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000848505.1_ViralProj14004/GCF_000848505.1_ViralProj14004_genomic.fna.gz

# 重新建立索引和比对
bowtie2-build ebov_genome.fa ebov_index
bowtie2 -x ebov_index -1 ../raw/SRR1972739_1.fastq.gz -2 ../raw/SRR1972739_2.fastq.gz -S ebov_aligned.sam
```

### 如果目标是大肠杆菌分析：

**建议**: 使用大肠杆菌的 RNA-seq 数据，例如：
- SRR1972739 (不是，这是埃博拉病毒)
- SRR1972740 (请验证)
- 或搜索 NCBI SRA: "Escherichia coli K12 RNA-seq"

---

## 文件位置总结

### 项目结构
```
/media/vimalinx/Data/bio_studio/
├── shared_data/
│   └── references/
│       ├── ecoli_K12.fa              # 大肠杆菌参考基因组 (4.5 MB)
│       └── ecoli_K12_index.*        # Bowtie2 索引文件
│
└── projects/test_rnaseq_analysis/
    ├── data/
    │   ├── raw/
    │   │   ├── SRR1972739_1.fastq.gz  # Read 1 (55 MB)
    │   │   └── SRR1972739_2.fastq.gz  # Read 2 (51 MB)
    │   │
    │   ├── processed/
    │   │   ├── qc/
    │   │   │   ├── SRR1972739_1_fastqc.html
    │   │   │   ├── SRR1972739_2_fastqc.html
    │   │   │   ├── SRR1972739_1_fastqc.zip
    │   │   │   └── SRR1972739_2_fastqc.zip
    │   │   │
    │   │   ├── SRR1972739_aligned.sam            # SAM 文件 (358 MB)
    │   │   ├── SRR1972739_aligned.bam            # BAM 文件
    │   │   ├── SRR1972739_aligned.sorted.bam     # 排序后的 BAM
    │   │   ├── SRR1972739_aligned.sorted.bam.bai # BAM 索引
    │   │   ├── bowtie2_summary.txt               # Bowtie2 摘要
    │   │   └── alignment_stats.txt               # SAMtools 统计
    │   │
    │   └── results/
    │
    ├── scripts/
    ├── notebooks/
    └── logs/
```

---

## 执行命令记录

### 完整命令列表

1. **创建项目**
```bash
python3 lib/create_project.py test_rnaseq_analysis --type generic --description "RNA-seq analysis for SRR1972739"
```

2. **下载参考基因组**
```bash
cd shared_data/references
curl -L "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz" -o ecoli_K12.fna.gz
gunzip -c ecoli_K12.fna.gz > ecoli_K12.fa
```

3. **下载测序数据**
```bash
cd projects/test_rnaseq_analysis/data/raw
curl -L -o SRR1972739_1.fastq.gz "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/009/SRR1972739/SRR1972739_1.fastq.gz"
curl -L -o SRR1972739_2.fastq.gz "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/009/SRR1972739/SRR1972739_2.fastq.gz"
```

4. **质量控制**
```bash
cd projects/test_rnaseq_analysis/data/processed
mkdir -p qc
fastqc -t 4 -o qc ../raw/SRR1972739_1.fastq.gz ../raw/SRR1972739_2.fastq.gz
```

5. **建立索引**
```bash
cd shared_data/references
bowtie2-build ecoli_K12.fa ecoli_K12_index
```

6. **序列比对**
```bash
cd projects/test_rnaseq_analysis/data/processed
bowtie2 -x /media/vimalinx/Data/bio_studio/shared_data/references/ecoli_K12_index \
        -1 ../raw/SRR1972739_1.fastq.gz \
        -2 ../raw/SRR1972739_2.fastq.gz \
        -S SRR1972739_aligned.sam -p 4
```

7. **转换和统计**
```bash
cd projects/test_rnaseq_analysis/data/processed
samtools view -bS SRR1972739_aligned.sam -o SRR1972739_aligned.bam
samtools sort SRR1972739_aligned.bam -o SRR1972739_aligned.sorted.bam
samtools index SRR1972739_aligned.sorted.bam
samtools flagstat SRR1972739_aligned.sorted.bam > alignment_stats.txt
```

---

## 总结

### 任务完成情况
- [x] 下载参考基因组 (E. coli K12)
- [x] 下载测序数据 (SRR1972739)
- [x] 执行 FastQC 质控
- [x] 建立 Bowtie2 索引
- [x] 序列比对
- [x] 转换 SAM 到 BAM 并排序
- [x] 生成比对统计

### 比对率
- **0.01%** (161/1,516,674 reads)

### ⚠️ 关键发现
- SRR1972739 是**埃博拉病毒**数据，而非大肠杆菌
- 因此与大杆菌基因组的比对率极低是**预期的**
- 建议使用正确的参考基因组重新分析

---

**报告生成时间**: 2026-01-22 17:20
**工作区**: /media/vimalinx/Data/bio_studio
**项目路径**: projects/test_rnaseq_analysis/
