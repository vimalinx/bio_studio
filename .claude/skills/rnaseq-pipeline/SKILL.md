---
name: rnaseq-pipeline
description: Complete RNA-seq data analysis pipeline. Process raw sequencing reads, align to genome, quantify expression, perform differential expression analysis, and visualize results. Use when user mentions RNA-seq, transcriptomics, gene expression, or differential expression.
allowed-tools:
  - Read
  - Write
  - Bash(fastqc:*)
  - Bash(star:*)
  - Bash(hisat2:*)
  - Bash(salmon:*)
  - Bash(featurecounts:*)
  - Bash(r:*)
  - Bash(python:*)
context: fork
agent: bio-expert
---

# RNA-seq Analysis Pipeline Skill

## Overview

Complete RNA-seq analysis from raw reads to differential expression and visualization. Supports:
- Quality control (FastQC, MultiQC)
- Read trimming (Trimmomatic, fastp)
- Alignment (STAR, HISAT2)
- Quantification (featureCounts, HTSeq, Salmon)
- Differential expression (DESeq2, edgeR, limma)
- Visualization (PCA, heatmaps, volcano plots)
- Functional enrichment (GO, KEGG)

## Complete Workflow

### 1. Quality Control

```bash
# FastQC on raw reads
fastqc -t 8 raw_reads/*.fastq.gz -o qc/raw_fastqc/

# Aggregate reports
multiqc qc/raw_fastqc/ -o qc/multiqc/
```

### 2. Read Trimming

```bash
# Trimmomatic
trimmomatic PE -threads 8 \
  input_R1.fastq.gz input_R2.fastq.gz \
  output_R1_paired.fastq.gz output_R1_unpaired.fastq.gz \
  output_R2_paired.fastq.gz output_R2_unpaired.fastq.gz \
  ILLUMINACLIP:adapters.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Or fastp (faster)
fastp -i input_R1.fastq.gz -I input_R2.fastq.gz \
  -o output_R1.fastq.gz -O output_R2.fastq.gz \
  --detect_adapter_for_pe \
  --trim_front1 5 --trim_front2 5 \
  --cut_mean_quality 20 \
  --length_required 36 \
  --thread 8 \
  --html qc/fastp.html \
  --json qc/fastp.json
```

### 3. Post-trimming QC

```bash
fastqc trimmed/*.fastq.gz -o qc/trimmed_fastqc/
multiqc qc/trimmed_fastqc/ -o qc/multiqc_trimmed/
```

### 4. Alignment

#### Option A: STAR (Spliced alignment)

```bash
# Build genome index (one-time)
STAR --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir genome_index/ \
  --genomeFastaFiles genome.fa \
  --sjdbGTFfile annotation.gtf \
  --sjdbOverhang 100

# Align reads
STAR --runThreadN 8 \
  --genomeDir genome_index/ \
  --readFilesIn read_R1.fastq.gz read_R2.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix sample1_ \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts \
  --twopassMode Basic
```

#### Option B: HISAT2 (Faster, less memory)

```bash
# Build index
hisat2-build genome.fa genome_index

# Align
hisat2 -p 8 -x genome_index \
  -1 read_R1.fastq.gz -2 read_R2.fastq.gz \
  -S sample1.sam

# Convert to BAM
samtools view -bS sample1.sam | samtools sort -o sample1.sorted.bam
samtools index sample1.sorted.bam
```

### 5. Quantification

#### Option A: featureCounts (Alignment-based)

```bash
featureCounts -T 8 -p -B -C \
  -a annotation.gtf \
  -o counts.txt \
  *.sorted.bam
```

#### Option B: Salmon (Alignment-free, faster)

```bash
# Build transcriptome index
salmon index -t transcripts.fa -i transcriptome_index --type quasi -k 31

# Quantify
salmon quant -i transcriptome_index \
  -l A \
  -1 read_R1.fastq.gz -2 read_R2.fastq.gz \
  -p 8 \
  --validateMappings \
  -o sample1_quant
```

### 6. Differential Expression (R/DESeq2)

```r
# scripts/run_deseq2.R
library(DESeq2)
library(apeglm)
library(pheatmap)

# Read count data
counts <- read.table("counts.txt", header=TRUE, row.names=1)
counts <- counts[,6:ncol(counts)]  # Remove first 5 columns from featureCounts

# Sample information
coldata <- data.frame(
  sample = colnames(counts),
  condition = c("control", "control", "treated", "treated")
)
rownames(coldata) <- coldata$sample

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

# Pre-filter low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq2
dds <- DESeq(dds)

# Results
res <- results(dds, contrast=c("condition", "treated", "control"))
res <- lfcShrink(dds, coef="condition_treated_vs_control", type="apeglm")

# Save results
write.csv(as.data.frame(res), "deseq2_results.csv")

# Significant genes
sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(sig), "significant_genes.csv")

# PCA
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition")

# Heatmap
topGenes <- head(order(rowMeans(counts(dds)), decreasing=TRUE), 50)
mat <- assay(vsd)[topGenes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col=coldata)

# MA plot
plotMA(res, main="DESeq2 MA Plot")

# Volcano plot
# See scripts/volcano_plot.R
```

### 7. Functional Enrichment

```r
# GO enrichment
library(clusterProfiler)
library(org.Hs.eg.db)

gene_list <- res$log2FoldChange
names(gene_list) <- rownames(res)
gene_list <- sort(gene_list, decreasing=TRUE)

go <- gseGO(geneList=gene_list,
            OrgDb=org.Hs.eg.db,
            ont="BP",
            nPerm=1000,
            minGSSize=10,
            maxGSSize=500,
            pvalueCutoff=0.05,
            verbose=FALSE)

# KEGG pathway
kegg <- gseKEGG(geneList=gene_list,
                organism='hsa',
                nPerm=1000,
                minGSSize=10,
                maxGSSize=500,
                pvalueCutoff=0.05)

# Save results
write.csv(as.data.frame(go), "go_enrichment.csv")
write.csv(as.data.frame(kegg), "kegg_enrichment.csv")
```

### 8. Downstream Variant Analysis (AI-Enhanced)

After calling variants with `bcftools` from RNA-seq data, use Evo 2 to predict functional effects:

```bash
# 1. Call variants
bcftools mpileup -f genome.fa aligned.sorted.bam | bcftools call -mv -O v -o variants.vcf

# 2. Analyze with Evo 2 (see bioinformatics-toolkit for setup)
# Predict variant effects (zero-shot)
docker exec evo2-analysis python3 evo2_variant_analysis.py \
  genome.fa variants.vcf results.json
```

## Utility Scripts

### Generate Sample Sheet

Create sample sheet manually or use a template:
```csv
sample_id,condition,fastq_r1,fastq_r2
sample1,control,data/raw/sample1_R1.fastq.gz,data/raw/sample1_R2.fastq.gz
sample2,treated,data/raw/sample2_R1.fastq.gz,data/raw/sample2_R2.fastq.gz
```

### Quality Summary

```bash
# Aggregate QC reports
multiqc qc/ -o qc/multiqc_summary/

# Manual QC summary required
# Aggregate metrics from:
# - FastQC reports
# - MultiQC summary
# - Alignment statistics
```

### Generate Report

```r
# Manual report generation required
# Tools: Rmarkdown, Jupyter, or custom scripts
# Include: QC metrics, DEG tables, plots, enrichment results
```

## Project Structure

```
project/
├── data/
│   ├── raw/              # Raw FASTQ files
│   ├── trimmed/          # Trimmed reads
│   ├── aligned/          # BAM files
│   └── quant/            # Quantification files
├── qc/
│   ├── raw_fastqc/       # FastQC reports
│   ├── trimmed_fastqc/
│   └── multiqc/          # MultiQC reports
├── results/
│   ├── counts.txt
│   ├── deseq2_results.csv
│   ├── significant_genes.csv
│   └── figures/
│       ├── pca.png
│       ├── heatmap.png
│       ├── volcano.png
│       └── ma_plot.png
├── scripts/              # Analysis scripts
├── genome/               # Reference genome
│   ├── genome.fa
│   ├── annotation.gtf
│   └── index/            # Alignment indices
└── report.html           # Final report
```

## Input Requirements

### Sequencing Data

- **Format**: FASTQ (gzipped)
- **Read length**: ≥ 50 bp single-end, ≥ 75 bp paired-end recommended
- **Depth**: 20-30 million reads for standard analysis
- **Quality**: Q30 ≥ 80%

### Reference Files

- **Genome**: FASTA format
- **Annotation**: GTF/GFF format (Ensembl, RefSeq, GENCODE)
- **Consistent**: Use matching genome and annotation versions

## Common Parameters

### STAR

```bash
--runThreadN 8              # Threads
--outFilterMultimapNmax 10  # Max multimap
--outSAMtype BAM SortedByCoordinate
--quantMode GeneCounts      # Output gene counts
--twopassMode Basic         # Two-pass mapping
```

### DESeq2

```r
# Design formula
design = ~ batch + condition  # Include batch effects

# Shrinkage
type="apeglm"  # OR "normal", "ashr"

# Significance
padj < 0.05           # FDR threshold
|log2FoldChange| > 1  # Fold change threshold
```

## Best Practices

### Experimental Design

1. **Biological replicates**: Minimum 3 per condition
2. **Randomization**: Randomize library prep
3. **Batch effects**: Process samples together when possible
4. **Controls**: Include appropriate controls
5. **Power analysis**: Determine required replicates

### Quality Control

1. **Check raw data**: FastQC before processing
2. **Trim adapters**: Always trim sequencing adapters
3. **Filter low quality**: Remove low-quality reads
4. **Check alignment**: Assess mapping rates
5. **Verify results**: Check known markers

### Statistical Analysis

1. **Normalize data**: Use appropriate methods (TMM, RLE, etc.)
2. **Model batch effects**: Include in design if present
3. **Correct for testing**: Use FDR, not p-values alone
4. **Validate findings**: Use independent method
5. **Report everything**: Include all parameters

## Troubleshooting

### Low Mapping Rate

**Issue**: < 70% reads mapped

**Causes**:
- Poor quality reads
- Wrong reference genome
- Species mismatch
- Contamination

**Solutions**:
- Check quality scores
- Verify reference
- Screen for contamination (FastQ Screen)
- Check adapter trimming

### Too Few Differentially Expressed Genes

**Issue**: Only a few genes significant

**Causes**:
- Low statistical power
- Similar conditions
- Technical issues

**Solutions**:
- Check sample size (power)
- Review experimental design
- Verify condition differences
- Check batch effects

### Batch Effects

**Issue**: Samples cluster by batch, not condition

**Solutions**:
- Include batch in DESeq2 design
- Use sva/Combat for correction
- Re-analyze with batch correction

## Integration with Other Skills

- **sequence-analysis**: Process raw sequences
- **genome-annotation**: Get gene annotations
- **pathway-analysis**: Enrichment analysis
- **visualization**: Create publication figures
- **single-cell**: Single-cell RNA-seq

## Advanced Features

### Alternative Splicing

```r
# Use DEXSeq or rMATS
library(DEXSeq)
```

### Fusion Detection

```bash
# STAR-Fusion
STAR-Fusion --genome_lib_dir /path/to/StarFusion/ \
  --left_fq read_R1.fastq.gz \
  --right_fq read_R2.fastq.gz \
  --output_dir fusion_out
```

### Allele-Specific Expression

```bash
# GATK ASEReadCounter
gatk ASEReadCounter \
  -R genome.fa \
  -I sample.bam \
  -V variants.vcf \
  -O asereadcounter.table
```

## For More Details

- Bioinformatics tools: See [bioinformatics-toolkit/SKILL.md](../bioinformatics-toolkit/SKILL.md)
- DESeq2: See https://bioconductor.org/packages/release/bioc/html/DESeq2.html
- STAR: See https://github.com/alexdobin/STAR
- HISAT2: See http://daehwankimlab.github.io/hisat2/
- featureCounts: See http://subread.sourceforge.net/
- MultiQC: See https://multiqc.info/
