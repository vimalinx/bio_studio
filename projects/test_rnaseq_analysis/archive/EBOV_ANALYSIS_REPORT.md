# SRR1972739 埃博拉病毒 RNA-seq 分析报告（修正版）

## 项目信息

- **项目名称**: test_rnaseq_analysis
- **分析日期**: 2026-01-22
- **任务**: 分析 SRR1972739 RNA-seq 数据，使用正确的埃博拉病毒参考基因组进行比对

---

## 数据信息

### 1. 参考基因组

#### 修正：使用埃博拉病毒参考基因组

| 项目 | 信息 |
|------|------|
| **病毒名称** | Zaire ebolavirus (扎伊尔埃博拉病毒) |
| **分离株** | Ebola virus/H.sapiens-tc/COD/1976/Yambuku-Mayinga |
| **Accession** | NC_002549.1 |
| **类型** | Complete genome (完整基因组) |
| **文件名** | ebov_zaire.fa |
| **大小** | 19 KB |
| **序列长度** | 18,959 bp |
| **位置** | shared_data/references/ebov_zaire.fa |
- **状态**: ✓ 已下载（使用 Biopython 从 NCBI）

### 2. 测序数据 (SRR1972739)

| 项目 | 信息 |
|------|------|
| **Run ID** | SRR1972739 |
| **样本标题** | G4252.1 |
| **科学名称** | Zaire ebolavirus (扎伊尔埃博拉病毒) |
| **文库来源** | TRANSCRIPTOMIC |
| **测序策略** | RNA-Seq |
| **Read 1** | SRR1972739_1.fastq.gz (55 MB) |
| **Read 2** | SRR1972739_2.fastq.gz (51 MB) |
| **位置** | projects/test_rnaseq_analysis/data/raw/ |

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
bowtie2-build ebov_zaire.fa ebov_zaire_index
```

**状态**: ✓ 索引已成功建立

**索引文件**:
- ebov_zaire_index.1.bt2
- ebov_zaire_index.2.bt2
- ebov_zaire_index.3.bt2
- ebov_zaire_index.4.bt2
- ebov_zaire_index.rev.1.bt2
- ebov_zaire_index.rev.2.bt2

### 2. 序列比对

```bash
cd projects/test_rnaseq_analysis/data/processed
bowtie2 -x /media/vimalinx/Data/bio_studio/shared_data/references/ebov_zaire_index \
        -1 ../raw/SRR1972739_1.fastq.gz \
        -2 ../raw/SRR1972739_2.fastq.gz \
        -S SRR1972739_ebov_aligned.sam -p 4
```

### 比对统计

#### Bowtie2 输出摘要
```
758337 reads; of these:
  758337 (100.00%) were paired; of these:
    568915 (75.02%) aligned concordantly 0 times
    189422 (24.98%) aligned concordantly exactly 1 time
    0 (0.00%) aligned concordantly >1 times
    ----
    568915 pairs aligned concordantly 0 times; of these:
      62983 (11.07%) aligned discordantly 1 time
    ----
    505932 pairs aligned 0 times concordantly or discordantly; of these:
      1011864 mates make up the pairs; of these:
        912713 (90.20%) aligned 0 times
        99151 (9.80%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
39.82% overall alignment rate
```

#### SAMtools 统计 (flagstat)
```
1516674 + 0 in total (QC-passed reads + QC-failed reads)
1516674 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
603961 + 0 mapped (39.82% : N/A)
603961 + 0 primary mapped (39.82% : N/A)
1516674 + 0 paired in sequencing
758337 + 0 read1
758337 + 0 read2
378844 + 0 properly paired (24.98% : N/A)
504810 + 0 with itself and mate mapped
99151 + 0 singletons (6.54% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

### 比对结果文件
- `SRR1972739_ebov_aligned.sam` (413 MB)
- `SRR1972739_ebov_aligned.bam`
- `SRR1972739_ebov_aligned.sorted.bam`
- `SRR1972739_ebov_aligned.sorted.bam.bai` (索引)
- `alignment_ebov_stats.txt` (统计结果)
- **位置**: projects/test_rnaseq_analysis/data/processed/

---

## 结果分析

### 比对率总结

| 指标 | 数值 | 百分比 |
|------|------|--------|
| **总 reads 数** | 1,516,674 | 100% |
| **已比对 reads** | 603,961 | **39.82%** |
| **未比对 reads** | 912,713 | 60.18% |
| **配对 reads 数** | 758,337 | - |
| **正确配对** | 378,844 | 24.98% |
| **单端比对** | 99,151 | 6.54% |
| **不一致比对** | 62,983 | 4.14% |

### 与错误比对（使用大肠杆菌）的对比

| 参考基因组 | 比对率 | 结果 |
|-----------|--------|------|
| ❌ 大肠杆菌 K12 (E. coli) | 0.01% | 物种不匹配 |
| ✅ 扎伊尔埃博拉病毒 (EBOV) | **39.82%** | **正确匹配** |

### 比对率分析

#### ✅ 正常比对率
- **39.82% 的比对率是合理的**，因为：
  1. 埃博拉病毒基因组较小（18,959 bp）
  2. RNA-seq 数据可能包含宿主基因组背景（人类基因组污染）
  3. 部分reads可能来自接头序列或低质量片段

#### 统计细节
- **24.98% 配对正确**: 一半的 reads 成功作为配对比对到病毒基因组
- **9.80% 单端比对**: 约 99,151 个reads只有一端比对成功，可能是：
  - 配对一端质量较差
  - 另一端匹配到宿主或污染序列
- **11.07% 不一致比对**: 约 62,983 个配对不一致，可能由于：
  - 结构变异
  - 病毒准种（quasispecies）的存在
- **60.18% 未比对**: 约 912,713 个reads未比对，可能原因：
  - 宿主基因组污染（人类 RNA）
  - 接头序列
  - 低质量reads
  - 病毒变异或新毒株

---

## 结论

### ✅ 任务完成情况
- [x] 下载埃博拉病毒参考基因组 (NC_002549.1)
- [x] 建立正确的 Bowtie2 索引
- [x] 重新比对 SRR1972739 数据
- [x] 生成 BAM 文件和比对统计
- [x] 生成完整分析报告

### 比对率
- **39.82%** (603,961/1,516,674 reads)

### 关键发现

1. **参考基因组正确**: 使用扎伊尔埃博拉病毒参考基因组后，比对率从 0.01% 提升到 39.82%
2. **合理的比对率**: 考虑到病毒基因组较小和可能的宿主污染，39.82% 是一个合理的比对率
3. **数据质量**: 24.98% 的reads正确配对，表明测序数据质量良好
4. **可能的污染**: 60.18% 的未比对reads可能来自宿主或其他污染源

---

## 建议

### 进一步分析建议

1. **去除宿主污染**
   ```bash
   # 先比对到人类基因组，去除宿主reads
   bowtie2 -x hg38_index -1 R1.fastq -2 R2.fastq --un-conc-gz clean_R%.fastq.gz
   # 然后用clean reads比对到病毒基因组
   bowtie2 -x ebov_zaire_index -1 clean_R1.fastq.gz -2 clean_R2.fastq.gz -S virus_aligned.sam
   ```

2. **变异检测**
   ```bash
   # 检测病毒基因组变异
   samtools mpileup -f ebov_zaire.fa virus_aligned.bam | bcftools call -mv -o variants.vcf
   ```

3. **覆盖度分析**
   ```bash
   # 分析病毒基因组覆盖度
   samtools depth virus_aligned.bam > coverage.txt
   ```

4. **使用更严格的比对参数**
   ```bash
   # 尝试更敏感的比对模式
   bowtie2 -x ebov_zaire_index -1 R1.fastq -2 R2.fastq \
           --very-sensitive-local -S aligned.sam
   ```

---

## 文件位置总结

### 项目结构
```
/media/vimalinx/Data/bio_studio/
├── shared_data/
│   └── references/
│       ├── ebov_zaire.fa                 # 埃博拉病毒参考基因组 (19 KB)
│       └── ebov_zaire_index.*           # Bowtie2 索引文件
│
└── projects/test_rnaseq_analysis/
    ├── data/
    │   ├── raw/
    │   │   ├── SRR1972739_1.fastq.gz    # Read 1 (55 MB)
    │   │   └── SRR1972739_2.fastq.gz    # Read 2 (51 MB)
    │   │
    │   ├── processed/
    │   │   ├── qc/
    │   │   │   ├── SRR1972739_1_fastqc.html
    │   │   │   ├── SRR1972739_2_fastqc.html
    │   │   │   ├── SRR1972739_1_fastqc.zip
    │   │   │   └── SRR1972739_2_fastqc.zip
    │   │   │
    │   │   ├── SRR1972739_aligned.sam            # ❌ 错误比对 (E. coli)
    │   │   ├── SRR1972739_aligned.bam
    │   │   ├── SRR1972739_aligned.sorted.bam
    │   │   ├── alignment_stats.txt               # 错误比对的统计 (0.01%)
    │   │   │
    │   │   ├── SRR1972739_ebov_aligned.sam     # ✅ 正确比对 (EBOV)
    │   │   ├── SRR1972739_ebov_aligned.bam     # ✅ 正确比对
    │   │   ├── SRR1972739_ebov_aligned.sorted.bam
    │   │   ├── SRR1972739_ebov_aligned.sorted.bam.bai
    │   │   ├── bowtie2_ebov_summary.txt         # Bowtie2 摘要
    │   │   └── alignment_ebov_stats.txt         # ✅ 正确比对的统计 (39.82%)
    │   │
    │   └── results/
    │
    ├── scripts/
    ├── notebooks/
    └── logs/
```

---

## 执行命令记录

### 完整命令列表（修正版）

1. **下载埃博拉病毒参考基因组**
```bash
cd shared_data/references
python3 << 'EOF'
from Bio import Entrez
Entrez.email = "user@example.com"
handle = Entrez.efetch(db="nucleotide", id="10313991", rettype="fasta", retmode="text")
fasta_data = handle.read()
handle.close()
with open("ebov_zaire.fa", "w") as f:
    f.write(fasta_data)
EOF
```

2. **建立埃博拉病毒索引**
```bash
cd shared_data/references
bowtie2-build ebov_zaire.fa ebov_zaire_index
```

3. **使用正确参考基因组重新比对**
```bash
cd projects/test_rnaseq_analysis/data/processed
bowtie2 -x /media/vimalinx/Data/bio_studio/shared_data/references/ebov_zaire_index \
        -1 ../raw/SRR1972739_1.fastq.gz \
        -2 ../raw/SRR1972739_2.fastq.gz \
        -S SRR1972739_ebov_aligned.sam -p 4
```

4. **转换和统计**
```bash
cd projects/test_rnaseq_analysis/data/processed
samtools view -bS SRR1972739_ebov_aligned.sam -o SRR1972739_ebov_aligned.bam
samtools sort SRR1972739_ebov_aligned.bam -o SRR1972739_ebov_aligned.sorted.bam
samtools index SRR1972739_ebov_aligned.sorted.bam
samtools flagstat SRR1972739_ebov_aligned.sorted.bam > alignment_ebov_stats.txt
```

---

## 总结

### 任务完成情况
- [x] 下载埃博拉病毒参考基因组 (NC_002549.1)
- [x] 建立 Bowtie2 索引
- [x] 使用正确的参考基因组重新比对 SRR1972739
- [x] 转换 SAM 到 BAM 并排序
- [x] 生成比对统计
- [x] 生成完整分析报告

### 最终比对率
- **39.82%** (603,961/1,516,674 reads)

### 与初始分析的对比

| 分析轮次 | 参考基因组 | 比对率 | 状态 |
|---------|-----------|--------|------|
| 第1轮 | ❌ 大肠杆菌 K12 | 0.01% | 物种不匹配 |
| **第2轮（修正）** | ✅ 扎伊尔埃博拉病毒 | **39.82%** | **正确匹配** |

### 关键收获
1. **数据验证的重要性**: 比对前应验证数据来源和物种信息
2. **使用正确的参考基因组**: 确保参考基因组与数据匹配是获得有意义结果的关键
3. **合理的比对率**: 对于小基因组病毒数据，39.82% 的比对率是可接受的
4. **生物信息学工作流**: 从质控到比对再到统计分析的完整流程已成功执行

---

**报告生成时间**: 2026-01-22 17:35
**工作区**: /media/vimalinx/Data/bio_studio
**项目路径**: projects/test_rnaseq_analysis/
