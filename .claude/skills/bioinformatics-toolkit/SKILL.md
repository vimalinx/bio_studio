---
name: bioinformatics-toolkit
description: Execute bioinformatics analysis including sequence alignment (BLAST, bowtie2, bwa), NGS data processing (samtools, bcftools), variant calling, gene prediction (prodigal), RNA structure analysis (ViennaRNA), and Python-based genomics workflows (Biopython). Automatically triggered when working with sequencing data, genome analysis, molecular biology data, or using bioinformatics tools.
allowed-tools: Bash, Read, Write, Grep, Glob, Edit
user-invocable: true
---

# Bioinformatics Toolkit

Comprehensive toolkit for molecular biology and genomics analysis. This Skill provides command-line tools and best practices for common bioinformatics workflows.

## Quick Reference

### Sequence Analysis
- **BLAST** (`blastn`): Nucleotide sequence alignment and similarity search
- **HMMER** (`hmmer`): Protein sequence profile searches
- **Biopython**: Python library for sequence manipulation

### Alignment & Mapping
- **bowtie2**: Fast and sensitive read alignment
- **bwa**: Burrows-Wheeler aligner for Illumina reads
- **samtools**: SAM/BAM file processing and analysis

### Genome Analysis
- **prodigal**: Prokaryotic gene prediction
- **ViennaRNA**: RNA secondary structure prediction
- **bcftools**: Variant calling and VCF file processing
- **bedtools**: Genome arithmetic toolkit

### Quality Control
- **fastp**: All-in-one FASTQ preprocessor (QC + Trimming)
- **FastQC**: Classic quality control tool
- **MultiQC**: Aggregate results from bioinformatics analyses

### Phylogenetics
- **iqtree**: Efficient phylogenomic inference (v3.x)

### AI-Powered Analysis
- **Biomni**: Specialized biological agents and tools (ORF finding, restriction sites, protocol design)
- **Evo 2**: Genomic foundation model for zero-shot variant effect prediction and sequence generation

## Environment Setup

Always use conda environments for bioinformatics tools:

```bash
# Activate bio environment (has torch and bio packages)
conda activate bio

# Or use base environment for system tools
# System has: blastn, bowtie2, bwa, samtools
```

## Common Workflows

### 1. BLAST Sequence Search

```bash
# Basic nucleotide BLAST
blastn -query input.fasta -db nt -out results.txt -outfmt 6

# With parameters
blastn -query sequence.fasta -db nt \
  -evalue 1e-10 -num_threads 4 \
  -out blast_results.txt -outfmt 6
```

### 2. Read Alignment with bowtie2

```bash
# Build index
bowtie2-build reference.fasta ref_index

# Align reads
bowtie2 -x ref_index -U reads.fq -S aligned.sam

# Convert to BAM and sort
samtools view -bS aligned.sam | samtools sort -o aligned_sorted.bam
samtools index aligned_sorted.bam
```

### 3. BAM File Processing

```bash
# View alignment statistics
samtools flagstat aligned.bam
samtools idxstats aligned.bam

# Extract mapped reads
samtools view -b -F 4 aligned.bam > mapped.bam

# Get coverage
samtools depth aligned.bam > coverage.txt
```

### 4. Gene Prediction with Prodigal

```bash
# Predict protein-coding genes
prodigal -i genome.fasta -a proteins.faa -d genes.fna \
  -f gff -o genes.gff -p single
```

### 5. RNA Secondary Structure (ViennaRNA)

```bash
# Predict minimum free energy structure
RNAfold < sequence.fasta

# Calculate base pairing probabilities
RNAplfold < sequence.fasta
```

## Python with Biopython

```python
from Bio import SeqIO
from Bio.Blast import NCBIWWW

# Read sequences
for record in SeqIO.parse("input.fasta", "fasta"):
    print(f"{record.id}: {len(record.seq)} bp")

# Online BLAST
result = NCBIWWW.qblast("blastn", "nt", sequence)
```

## Python with Biomni (No API Key Required)

Biomni provides specialized tools that can be imported directly:

```python
import sys
import subprocess
from pathlib import Path

# Locate Bio Studio root
root = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip()
sys.path.append(str(Path(root) / "repositories/active/Biomni"))

from biomni.tool.molecular_biology import annotate_open_reading_frames, find_restriction_sites

# Find ORFs
orfs = annotate_open_reading_frames(sequence, min_length=300)

# Find Restriction Sites
sites = find_restriction_sites(sequence, ["EcoRI", "BamHI"])
```

## Deep Learning with Evo 2 (Docker)

For advanced variant effect prediction:

```bash
# Start Evo 2 container
docker compose -f repositories/active/evo2/docker-compose.override.yml run --rm -d --name evo2-analysis evo2 sleep infinity

# Run analysis script
docker exec evo2-analysis python3 /workdir/projects/my_project/scripts/evo2_analysis.py
```

## Best Practices

1. **Always check tool versions**: `tool --version` or `conda list tool-name`
2. **Use appropriate output formats**: `-outfmt 6` for BLAST (tabular), BAM for alignments
3. **Parallel processing**: Use `-num_threads` or `-t` for multi-core tools
4. **File compression**: Use `.gz` for large files, pipe with `gzip`
5. **Verify results**: Check output files, use `samtools quickcheck` for BAM files

## Additional Resources

- For detailed tool usage, see [TOOLS.md](TOOLS.md)
- For complete analysis workflows, see [WORKFLOWS.md](WORKFLOWS.md)
- Utility scripts are in the `scripts/` directory

## Troubleshooting

### Tool not found
- Check conda environment: `conda list | grep tool-name`
- Switch to bio environment: `conda activate bio`

### File format issues
- Validate FASTA: Check for proper headers and sequence lines
- Check SAM/BAM: `samtools quickcheck file.bam`
- Convert formats: Use `samtools view` or Biopython

### Memory errors
- Reduce threads: `-num_threads 2`
- Process in chunks: Split input files
- Use streaming: Pipe commands instead of intermediate files

---

## Getting Help

For each tool, check built-in help:
```bash
tool --help
man tool  # for detailed manual
```

For complex analyses, Claude can also invoke specialized Skills:
- `sequence-analysis`: Advanced sequence operations
- `rnaseq-pipeline`: Complete RNA-seq workflows
- `protein-structure`: Protein 3D structure prediction
