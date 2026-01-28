# Bioinformatics Skills & Tools Summary

> **Last Updated**: 2026-01-25
> **Environment**: Bio Studio v2.1 (Conda `bio` env + Docker AI Models)
> **Versions**: 以 `docs/ENVIRONMENT.md` 为准

## 🛠️ Installed Tools & Skills

This system is equipped with a comprehensive suite of bioinformatics tools, all accessible via CLI.

### 1. Quality Control & Preprocessing
| Tool | Description | Skill Trigger |
|---|---|---|
| **Fastp** | All-in-one FASTQ preprocessor (QC + Trimming) | "qc", "filter reads" |
| **FastQC** | Classic QC metrics visualization | "qc report" |
| **MultiQC** | Aggregate reports from FastQC, STAR, etc. | "merge reports" |
| **SeqKit** | Cross-platform FASTA/Q swiss army knife | "fasta stats", "convert" |

### 2. Sequence Alignment
| Tool | Application | Skill Trigger |
|---|---|---|
| **BWA-MEM** | DNA alignment (Genome) | "align dna" |
| **Bowtie2** | DNA alignment (Short reads) | "align bowtie" |
| **STAR** | RNA-seq alignment (Spliced) | "align rna" |
| **HISAT2** | RNA-seq alignment (Fast) | "align rna fast" |

### 3. Variant Calling & Processing
| Tool | Description | Skill Trigger |
|---|---|---|
| **Samtools** | BAM/SAM manipulation | "sort bam", "index bam" |
| **Bcftools** | Variant calling & VCF utils | "call variants" |
| **Bedtools** | Genomic arithmetic | "intersect", "merge bed" |

### 4. Phylogenetics & MSA
| Tool | Description | Skill Trigger |
|---|---|---|
| **MAFFT** | High-speed Multiple Sequence Alignment | "align multiple" |
| **Muscle** | MSA (Classic) | "align multiple" |
| **IQ-TREE 2** | Maximum Likelihood Phylogeny | "build tree" |

### 5. AI & Deep Learning Models
| Tool | Source | Description | Skill Trigger |
|---|---|---|---|
| **Evo 2 (1B)** | Docker | Genomic foundation model (Variant effects) | "predict variant effect" |
| **ESM** | Python | Protein language model (Structure/Function) | "analyze protein" |
| **Biomni** | Agent | Biomedical Research Agent | "search paper", "biomni" |

---

## 🤖 Skill Usage Patterns

Bio Studio uses an "Intent-Based" skill system. You don't need to remember flags.

### Example: RNA-seq Pipeline
**User**: "Run QC on these files, then align to hg38 and call variants."

**AI Execution Plan**:
1.  **QC**: Invokes `fastp` (Skill: `qc`)
    *   `fastp -i R1.fq -I R2.fq -o clean_R1.fq -O clean_R2.fq`
2.  **Align**: Invokes `STAR` (Skill: `align_rna`)
    *   `STAR --genomeDir hg38 --readFilesIn clean_R1.fq ...`
3.  **Process**: Invokes `samtools` (Skill: `process_bam`)
    *   `samtools sort -o aligned.sorted.bam`
4.  **Variant**: Invokes `bcftools` (Skill: `variant_call`)
    *   `bcftools mpileup ... | bcftools call ...`

### Example: Evolutionary Analysis
**User**: "Predict the effect of these mutations using Evo 2."

**AI Execution Plan**:
1.  **Model**: Invokes `evo2` container (Skill: `evo2`)
2.  **Input**: Parses VCF and Reference Genome
3.  **Output**: JSON report with log-likelihood scores

---

## 📂 Environment Paths

- **Conda Env**: `/home/vimalinx/miniforge3/envs/bio/bin/` (Add to PATH in scripts!)
- **AI Models**: `docker run evo2:latest`
- **Reference Data**: `shared_data/references/`
