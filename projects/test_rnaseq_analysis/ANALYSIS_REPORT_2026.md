# Ebola Virus Analysis Report (2026 Update)

**Project**: `test_rnaseq_analysis`  
**Date**: 2026-01-24  
**Sample**: SRR1972739 (Zaire ebolavirus, 2014 Outbreak)  
**Pipeline**: Bio Studio v2.0 (Fastp -> Bowtie2 -> Bcftools -> Evo2 AI)

---

## 1. Executive Summary

This project re-analyzed the RNA-seq data from the 2014 Ebola outbreak using a modernized bioinformatics pipeline integrated with the **Evo 2 genomic foundation model**. 

**Key Findings:**
- **High Mutation Load**: Identified **563 high-quality variants** compared to the 1976 reference genome (NC_002549.1), reflecting decades of viral evolution.
- **AI-Predicted Functional Hotspots**: Evo 2 identified critical mutations in the **L gene (Polymerase)** and **Leader region**, suggesting potential functional shifts in viral replication.
- **Data Quality**: The sample showed a **53.3% alignment rate** to the Ebola genome, confirming high viral load in the clinical sample.

---

## 2. Methodology

### 2.1 Bioinformatics Pipeline
The analysis was automated using `scripts/pipeline.py` with the following workflow:

1.  **Quality Control**: `fastp` (v0.23.4) for adapter trimming and quality filtering.
2.  **Alignment**: `bowtie2` (v2.5.4) mapping against *Zaire ebolavirus* (NC_002549.1).
3.  **Processing**: `samtools` (v1.19) for sorting and indexing.
4.  **Variant Calling**: `bcftools` (v1.19) using mpileup/call model.

### 2.2 AI Impact Analysis
Variants were analyzed using **Evo 2 (1B parameters)**, a genomic foundation model, to predict the functional impact of mutations based on evolutionary likelihood.
- **Model**: `evo2_1b_base`
- **Context Window**: 100bp flanking sequence
- **Metric**: Effect Score (Log-likelihood difference between ALT and REF alleles)

---

## 3. Results

### 3.1 Sequencing Statistics
- **Raw Reads**: 758,337 pairs
- **Clean Reads**: 482,396 pairs (63.6% retention after strict filtering)
- **Alignment Rate**: 53.27% (High viral content)
- **Duplication Rate**: 1.78% (Low PCR bias)

### 3.2 Variant Profile
- **Total Variants**: 563
- **Type**: 100% SNPs (No Indels detected in consensus)
- **Dominant Patterns**: Strong bias towards `C->T` and `A->G` transitions, consistent with viral evolution signatures and potential host RNA editing (ADAR/APOBEC).

### 3.3 Evo 2 AI Analysis
The AI model evaluated the "fitness" of each mutation. Strongly negative scores indicate mutations that deviate from the model's learned evolutionary constraints (potentially deleterious or functionally altering).

**Top 10 High-Impact Mutations:**

| Position | Ref | Alt | Region | Effect Score | Quality | Interpretation |
|----------|-----|-----|--------|--------------|---------|----------------|
| **17,535** | T | G | **L Gene** | **-0.809** | 225 | High impact on Polymerase; potential replication modulation |
| **13,856** | A | G | **L Gene** | **-0.609** | 225 | Polymerase mutation |
| **14,154** | G | A | **L Gene** | **-0.570** | 228 | Polymerase mutation |
| **6,476** | G | A | **GP Gene** | **-0.568** | 228 | Surface Glycoprotein mutation |
| **381** | G | A | **Leader** | **-0.508** | 73 | Non-coding regulatory region; may affect transcription |
| **15,816** | T | C | **L Gene** | **-0.508** | 225 | Polymerase mutation |
| **15,387** | C | T | **L Gene** | **-0.496** | 225 | Polymerase mutation |
| **11,079** | T | C | **VP24** | **+0.480** | 225 | Matrix protein mutation (Positive score = High likelihood) |
| **11,982** | G | A | **VP24** | **-0.469** | 225 | Matrix protein mutation |
| **17,568** | A | G | **L Gene** | **-0.465** | 225 | Polymerase mutation |

*Note: Region mapping is approximate based on NC_002549.1 coordinates.*

---

## 4. Conclusion

This analysis successfully reconstructed the mutational landscape of the 2014 Ebola virus sample using a fully automated Bio Studio pipeline. The integration of **Evo 2** provided unique insights, highlighting that while most mutations are neutral (mean score ~0), specific hotspots in the **Polymerase (L gene)** and **Regulatory regions** show significant evolutionary divergence.

These high-impact variants are prime candidates for further experimental validation regarding viral replication efficiency and pathogenicity.
