---
name: sequence-analysis
description: Analyze DNA, RNA, and protein sequences. Use BLAST for similarity searches, translate sequences, predict ORFs, calculate GC content, and perform sequence alignments. Automatically triggered when user mentions sequences, BLAST, alignment, translation, or sequence analysis.
allowed-tools:
  - Read
  - Bash(blast*:*)
  - Bash(hmmer*:*)
  - Bash(mafft:*)
  - Bash(bowtie2:*)
  - Bash(bwa:*)
  - Write
---

# Sequence Analysis Skill

## Overview

This Skill provides comprehensive sequence analysis capabilities including:
- BLAST similarity searches
- Sequence translation (DNA/RNA → Protein)
- ORF prediction
- GC content calculation
- Multiple sequence alignments
- Sequence format conversions

## Quick Reference

### BLAST Searches

```bash
# Protein BLAST against Swiss-Prot
blastp -query input.fa -db swissprot -out results.txt -evalue 1e-5

# Nucleotide BLAST
blastn -query gene.fa -db nt -out results.txt

# Translate and search
tblastn -query protein.fa -db nt -out results.txt
```

### Sequence Translation

```python
from Bio.Seq import Seq
dna_seq = Seq("ATGCGTACGT")
protein_seq = dna_seq.translate()
print(protein_seq)
```

### ORF Prediction

```bash
# Prokaryotic genomes (CLI)
prodigal -i genome.fa -a proteins.faa -d genes.fna -f gff

# Eukaryotic genomes (CLI)
augustus --species=human genome.fa > predictions.gff
```

```python
# Python (using Biomni)
from biomni.tool.molecular_biology import annotate_open_reading_frames
orfs = annotate_open_reading_frames(sequence, min_length=300)
```

### Multiple Sequence Alignment

```bash
# MAFFT (recommended for accuracy)
mafft input.fa > aligned.fa

# ClustalW (faster)
clustalw input.fa

# MUSCLE (balanced)
muscle -in input.fa -out aligned.fa
```

### HMMER for Protein Domains

```bash
# Search Pfam
hmmscan --domtblout pfam_results.txt Pfam-A.hmm proteins.fa

# Build custom HMM
hmmbuild my_family.hmm alignment.fa
hmmsearch my_family.hmm database.fa
```

## Common Workflows

### 1. Identify Unknown Sequence

```bash
# Step 1: BLAST search
blastp -query unknown.fa -db nr -out blast_results.txt

# Step 2: If no hit, try structure prediction
# See: protein-structure Skill

# Step 3: Check for domains
hmmscan --domtblout domains.txt Pfam-A.hmm unknown.fa
```

### 2. Annotate a Genome

```bash
# Predict genes
prodigal -i genome.fa -a proteins.fa -d genes.fa -o coordinates.gff

# Annotate proteins
blastp -query proteins.fa -db swissprot -out annotations.txt
hmmscan --domtblout pfam.txt Pfam-A.hmm proteins.fa

# Combine results
# See: genome-annotation Skill
```

### 3. Analyze Gene Expression

```bash
# Align RNA-seq reads
hisat2 -x genome_index -1 reads_1.fq -2 reads_2.fq -S aligned.sam

# Convert to BAM
samtools view -bS aligned.sam > aligned.bam
samtools sort aligned.bam -o aligned.sorted.bam
samtools index aligned.sorted.bam

# Count reads
featureCounts -a annotation.gff -o counts.txt aligned.sorted.bam
```

## Utility Scripts

### Validate FASTA Format

```bash
# Manual validation or use Biopython
python -c "from Bio import SeqIO; SeqIO.read('input.fa', 'fasta')"

# Check:
# - Proper FASTA header (>name)
# - Valid nucleotide/amino acid characters
# - No empty sequences
```

### Extract Sequences by ID

```python
# Use Biopython
from Bio import SeqIO

with open('id_list.txt') as ids:
    id_set = set(line.strip() for line in ids)

records = [r for r in SeqIO.parse('input.fa', 'fasta') if r.id in id_set]
SeqIO.write(records, 'output.fa', 'fasta')
```

### Calculate Sequence Statistics

```python
# Use Biopython
from Bio import SeqIO
from Bio.SeqUtils import GC

for record in SeqIO.parse('input.fa', 'fasta'):
    print(f"{record.id}: {len(record.seq)} bp, GC={GC(record.seq):.1f}%")
```

## Troubleshooting

### BLAST Too Slow

- Use local database instead of online
- Add `-num_threads` option
- Reduce database size with `-db` specific database

### Format Not Recognized

Check file format:
```bash
head -1 input.fa
```

Should start with `>` for FASTA. Use `seqtk` to convert:
```bash
seqtk format input.fa > output.fa
```

### Alignment Fails

- Check sequence quality (no ambiguous bases)
- Remove duplicate sequences
- Try different aligner (MAFFT more tolerant)

## Integration with Other Skills

- **protein-structure**: Predict 3D structure of proteins
- **genome-annotation**: Full genome annotation pipeline
- **rnaseq-analysis**: RNA-seq differential expression
- **variant-calling**: SNP and indel detection

## Best Practices

1. **Always validate input** - Check FASTA format before analysis
2. **Use appropriate e-value** - 1e-5 for proteins, 1e-10 for nucleotides
3. **Save intermediate results** - Don't redo expensive computations
4. **Document parameters** - Keep track of which parameters were used
5. **Validate results** - Cross-check with multiple methods

## For More Details

- Bioinformatics tools: See [bioinformatics-toolkit/SKILL.md](../bioinformatics-toolkit/SKILL.md)
- BLAST: See https://blast.ncbi.nlm.nih.gov/Blast.cgi
- HMMER: See http://hmmer.org/
- Biopython: See https://biopython.org/
- MAFFT: See https://mafft.cbrc.jp/alignment/software/
- MUSCLE: See https://www.drive5.com/muscle/
- Prodigal: See https://github.com/hyattpd/Prodigal
