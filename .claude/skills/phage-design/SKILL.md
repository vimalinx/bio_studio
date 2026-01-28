---
name: phage-design
description: Design synthetic bacteriophages for specific bacterial hosts. Analyze host receptors, design tail fibers, assemble genomes, optimize codons, and assess safety. Automatically triggered when user mentions phage, bacteriophage, phage therapy, or viral design.
allowed-tools:
  - Read
  - Write
  - Bash(blast*:*)
  - Bash(python:*)
  - Bash(prodigal:*)
context: fork
agent: bio-expert
---

# Synthetic Phage Design Skill

## Overview

Design custom bacteriophages targeting specific bacteria. This Skill integrates:
- Host receptor analysis
- Tail fiber protein design
- Genome assembly and optimization
- Safety and biosecurity assessment
- In silico validation

## Design Workflow

### Complete Pipeline (Evo2-based)

```bash
# Use the genome design filtering pipeline
cd repositories/active/evo2/phage_gen
python pipelines/genome_design_filtering_pipeline.py config.yaml

# This pipeline includes:
# - Host genome analysis
# - Receptor identification
# - Tail fiber design
# - Genome assembly
# - Optimization
# - Safety filtering
```

### Phase 1: Host Analysis

Analyze the host bacterial genome to identify:
- Surface proteins (potential receptors)
- Restriction-modification systems
- CRISPR-Cas arrays
- Antiphage defense systems

### Phase 2: Receptor Identification

Use structural prediction to identify potential receptor proteins:
```bash
# Use protein-structure Skill for receptor modeling
colabfold_batch receptor_sequences/ structures/
```

### Phase 3: Tail Fiber Design

Design tail fiber proteins to bind specific host receptors.

### Phase 4: Genome Assembly

Assemble complete phage genome using Evo2:
```bash
cd repositories/active/evo2/phage_gen
python pipelines/genome_design_filtering_pipeline.py config.yaml
```

Includes modules:
- Structural proteins (capsid, tail, fiber)
- DNA packaging (terminase, portal)
- Lysis (holin, lysin, spanin)
- Replication (DNA pol, helicase)
- Regulation (promoters, repressors)

### Phase 5: Optimization

Optimize codon usage for target host:
```bash
# Use DNA design tools
python tools/scripts/dna_design.py --protein protein.fa --species E.coli
```

Optimizes:
- Codon adaptation index (CAI)
- GC content
- mRNA stability
- Removal of restriction sites

### Phase 6: Safety Assessment

Perform comprehensive safety checks:
- Virulence genes (BLAST against VFDB)
- Antibiotic resistance (CARD)
- Toxin genes
- Horizontal transfer risk
- Environmental impact

## Utility Scripts

### Validate Genome

```bash
# Use Prodigal for gene prediction
prodigal -i phage.fa -a proteins.fa -d genes.fa -f gff

# Check completeness
# Manual validation required for:
# - Complete genome (no gaps)
# - All essential genes present
# - No frame-shift mutations
# - Proper gene order
```

### Predict Proteomics

```bash
# Predict all proteins
prodigal -i phage.fa -a proteins.fa -d genes.fa -f gff

# Annotate proteins
blastp -query proteins.fa -db phage_db -out annotations.txt
hmmscan --domtblout domains.txt Pfam-A.hmm proteins.fa
```

### Estimate Titer

```bash
# Titer prediction requires experimental data
# Use growth modeling tools from:
# evo2/phage_gen/pipelines/genome_design_filtering_pipeline.py
```

## Reference Data

### Essential Phage Genes

| Gene | Function | Required |
|------|----------|----------|
| **terminase** | DNA packaging | Yes |
| **portal** | DNA entry/exit | Yes |
| **capsid** | Head protein | Yes |
| **tail** | Tail proteins | Yes |
| **fiber** | Host recognition | Yes* |
| **holin** | Membrane disruption | Yes |
| **lysin** | Peptidoglycan degradation | Yes |
| **integrase** | Integration | Optional† |
| **repressor** | Lysogeny control | Optional† |

* Required for virulent phages
† For temperate phages

### Common Modules

```
[Regulatory]
  ├── Early promoter (constitutive)
  ├── Late promoter (σ-dependent)
  └── Terminator (Rho-independent)

[Replication]
  ├── DNA polymerase
  ├── Helicase
  ├── Primase
  └── Single-strand binding

[Structural]
  ├── Capsid (scaffold + major)
  ├── Portal (connector)
  ├── Tail sheath
  ├── Tail tube
  └── Tail fibers

[Packaging]
  ├── Terminate large subunit
  ├── Terminate small subunit
  └── Portal

[Lysis]
  ├── Holin (timing)
  ├── Endolysin (degradation)
  └── Spanin (outer membrane)
```

## Example Projects

### Project 1: E. coli Phage

```bash
# Target: Uropathogenic E. coli
# Receptor: FimH

python scripts/design_phage.py \
    --target_fimh fimH_sequence.fa \
    --host E_coli \
    --type virulent \
    --output ecophage.fa
```

### Project 2: MRSA Phage

```bash
# Target: Methicillin-resistant S. aureus
# Receptor: Wall teichoic acid

python scripts/design_phage.py \
    --target_wta WTA_biosynthesis.fa \
    --host S_aureus \
    --type virulent \
    --output mrsaphage.fa
```

### Project 3: Phage Cocktail

```bash
# Design multiple phages for same host
python scripts/design_cocktail.py \
    --host E_coli \
    --num_phages 3 \
    --diversity high \
    --output cocktail_genomes/
```

## Validation Workflow

### In Silico

```bash
# 1. Genome validation
python scripts/validate_phage_genome.py designed.fa

# 2. Protein structure prediction
colabfold_batch proteins/ structures/

# 3. Docking validation
python/scripts/validate_binding.py tail_fiber.pdb receptor.pdb

# 4. Safety check
python scripts/safety_assessment.py designed.fa
```

### In Vitro (Experimental)

```bash
# Order synthetic genome
# (See: ordering_genome.md)

# Assemble genome
# Gibson assembly or yeast assembly

# Package genome
# In vitro packaging extract

# Titer determination
# Plaque assay

# Host range testing
# Spot test against bacterial panel
```

## Safety and Biosecurity

### Mandatory Checks

1. **Virulence factors** - BLAST against VFDB
2. **Antibiotic resistance** - Check against CARD
3. **Toxin genes** - Search known toxins
4. **Lysogeny potential** - Check for integrase
5. **Transduction risk** - Assess packaging site
6. **Environmental impact** - Model ecosystem effects

### Safety Score

```python
# Calculate safety score (0-100)
safety_score = calculate_safety(phage_genome)

# Components:
# - Absence of virulence genes (40%)
# - Absence of resistance genes (30%)
# - Limited host range (15%)
# - No lysogeny capability (10%)
# - Environmental stability (5%)

# Threshold: ≥80 to proceed
```

### Required Documentation

- [ ] Safety assessment report
- [ ] Genome annotation
- [ ] Predicted proteomics
- [ ] Structural models
- [ ] Experimental plan
- [ ] Risk assessment
- [ ] Containment strategy

## Integration with Other Skills

- **sequence-analysis**: Analyze phage genome
- **protein-structure**: Predict structural proteins
- **protein-design**: Optimize tail fibers
- **genome-annotation**: Annotate phage genes
- **lab-automation**: Automate experiments

## Troubleshooting

### Design Fails

**No receptor identified**
- Expand surface protein search
- Use transcriptomics data
- Experimental validation needed

**Tail fiber unstable**
- Check folding predictions
- Try different scaffolds
- Add stabilizing mutations

**Genome too large**
- Remove non-essential genes
- Use compact promoters
- Optimize overlap regions

### Safety Concerns

**Virulence genes found**
- Remove or replace genes
- Redesign to eliminate function
- Consider alternative design

**Broad host range**
- Add targeting specificity
- Use narrow-range fibers
- Experimental validation required

## Best Practices

1. **Start with known phage** - Modify existing genome
2. **Modular design** - Use interchangeable parts
3. **Validate computationally** - Before synthesis
4. **Test iteratively** - Build and test modules
5. **Document everything** - Track all modifications
6. **Safety first** - Never skip safety checks

## For More Details

- Evo2 phage design: See [repositories/active/evo2/phage_gen/README.md](../../repositories/active/evo2/phage_gen/README.md)
- Complete pipeline: See [repositories/active/evo2/phage_gen/pipelines/genome_design_filtering_pipeline.py](../../repositories/active/evo2/phage_gen/pipelines/genome_design_filtering_pipeline.py)
- Related tools:
  - `repositories/active/evo2/phage_gen/pipelines/genetic_architecture.py` - Genome architecture analysis
  - `repositories/active/evo2/phage_gen/pipelines/genome_gibson_assembly.py` - Gibson assembly design
  - `tools/scripts/dna_design.py` - DNA design and optimization
