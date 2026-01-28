---
name: protein-structure
description: Predict and analyze protein 3D structures using AlphaFold, ESM, and molecular visualization. Use when user mentions protein structure, folding, AlphaFold, PyMOL, structure prediction, or molecular visualization.
allowed-tools:
  - Read
  - Write
  - Bash(colabfold:*)
  - Bash(pymol:*)
  - Bash(chimera:*)
  - Bash(python:*)
context: fork
agent: general-purpose
---

# Protein Structure Prediction Skill

## Overview

Predict protein 3D structures, analyze structural features, and visualize molecular interactions. Supports:
- AlphaFold2 high-accuracy prediction
- ESM fast protein structure prediction
- Structure visualization (PyMOL, ChimeraX)
- Structural analysis (domains, binding sites)
- Molecular docking preparation

## Quick Start

### Option 1: AlphaFold (High Accuracy)

```bash
# Using ColabFold (recommended for most users)
colabfold_batch input/ output/

# Input: FASTA file with protein sequences
# Output: PDB files with confidence scores (pLDDT)
```

### Option 2: ESM (Fast Prediction)

```python
import esm
import torch

# Load model
model, alphabet = esm.pretrained.load_model("esm2_t33_650M_UR50D")
model = model.eval()  # Disable dropout

# Predict structure
# (See reference.md for complete example)
```

### Option 3: Local AlphaFold (Requires GPU)

```bash
# Requires setup (see installation guide)
run_alphafold.sh --fasta_dir input/ --output_dir output/
```

## Input Requirements

### FASTA Format

```
>protein_name
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFDSGDLSTPDAVMGNPK
VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG
KEFTPPVQAAYQKVVAGVANALAHKYH
```

**Requirements**:
- Standard amino acids only (no X, B, Z)
- Length: 50-1000 residues optimal
- Remove signal peptides if present
- Check for unusual sequences (low complexity regions)

## Understanding Results

### AlphaFold Output

- **ranked_0.pdb**: Top prediction
- **ranking_debug.json**: Confidence scores
- **pLDDT**: Per-residue confidence (0-100)
  - >90: Very high confidence
  - 70-90: Good confidence
  - 50-70: Low confidence
  - <50: Very low confidence

### ESM Output

- Structure in PDB format
- Per-residue confidence scores
- Predicted aligned error (PAE) matrix

## Visualization

### PyMOL

```bash
# Load structure
pymol predicted_structure.pdb

# Color by pLDDT (if available)
spectrum b, plddt

# Show cartoon
show cartoon

# Identify domains
# See reference.md for domain identification
```

### ChimeraX

```bash
# Open structure
chimerax predicted_structure.pdb

# Color by confidence
color bfactor palette

# Measure distances
distance :1.A :1.B

# Create publication-quality images
lighting full
copy image.png 2000,2000 supersample 3
```

## Structural Analysis

### Domain Identification

```python
# Using HMMER
from subprocess import run
run(["hmmscan", "--domtblout", "domains.txt",
      "Pfam-A.hmm", "predicted_structure.pdb"])
```

### Binding Site Prediction

```bash
# Use fpocket
fpocket -f predicted_structure.pdb

# Output in predicted_structure_pocket/
```

### Surface Analysis

```python
# Calculate solvent accessible surface area
from Bio.PDB import *
parser = PDBParser()
structure = parser.get_structure('protein', 'predicted.pdb')

# See reference.md for complete example
```

## Common Workflows

### 1. Predict Structure for Unknown Protein

```bash
# Step 1: Check sequence quality
python scripts/check_sequence.py input.fa

# Step 2: Run AlphaFold
colabfold_batch input/ output/

# Step 3: Evaluate confidence
# Check pLDDT scores in ranking_debug.json

# Step 4: If low confidence, try:
# - Remove N-terminal region
# - Split into domains
# - Use ESM as alternative
```

### 2. Compare Multiple Structures

```bash
# Predict all structures
colabfold_batch multimer_input/ output/

# Align structures
pymol structure1.pdb structure2.pdb
# In PyMOL: align structure2, structure1

# Calculate RMSD
# See scripts/calculate_rmsd.py
```

### 3. Prepare for Docking

```bash
# Step 1: Predict structure
colabfold_batch input/ output/

# Step 2: Prepare PDB
# Remove water, add hydrogens, optimize
python scripts/prepare_for_docking.py ranked_0.pdb

# Step 3: Define binding site
# Use fpocket or manual selection

# Step 4: Run docking
# See: molecular-docking Skill
```

## Utility Scripts

### Protein Prediction Tool

```bash
# Use the main protein prediction script
python tools/scripts/protein_predict.py --sequence "MVHLTPEEKSA..."

# Features:
# - ESM-based structure prediction
# - Confidence scoring
# - Structure validation
# - Domain identification
```

### Check Sequence Quality

```bash
# Manual quality checks:
# - Verify standard amino acids only
# - Check for low-complexity regions
# - Ensure length: 50-1000 residues optimal
# - Remove signal peptides if present
```

### Add Confidence Color Scheme to PDB

```bash
# Manual PDB editing required
# Tools: PyMOL, ChimeraX, or Biopython
```

## Troubleshooting

### Prediction Fails

**Error: "Sequence too long"**
- ESM: Max 1024 residues, split into domains
- AlphaFold: Max 2000 residues, may need more memory

**Error: "Low confidence"**
- Check if protein is intrinsically disordered
- Try removing low-complexity regions
- Use multiple sequence alignment if possible

**Error: "Out of memory"**
- Use smaller model (esm2_t30_150M_UR50D)
- Reduce batch size
- Close other programs

### Poor Quality Structure

1. **Check confidence scores** - If all <50, structure may be wrong
2. **Verify sequence** - Mutations or errors cause problems
3. **Check for domains** - Predict domains separately
4. **Use template** - AlphaFold finds templates automatically

### Visualization Issues

**PyMOL crashes**
- Structure too large? Use ChimeraX
- Not enough memory? Close other programs

**Colors not showing**
- Load pLDDT as B-factor first
- Check file format (should be PDB)

## Integration with Other Skills

- **sequence-analysis**: Prepare sequences, analyze domains
- **molecular-docking**: Prepare structures for docking
- **molecular-dynamics**: Set up MD simulations
- **protein-design**: Design variants based on structure
- **mutagenesis**: Predict effects of mutations

## Best Practices

1. **Always check confidence** - Don't trust low-confidence regions
2. **Validate experimentally** - Confirm key features
3. **Use multiple methods** - Compare AlphaFold and ESM
4. **Document parameters** - Save which model/version was used
5. **Consider templates** - Known structures improve accuracy

## Performance Tips

### Speed

- **ESM**: ~1 min for 300 residues (CPU)
- **ColabFold**: ~20 min for 500 residues (GPU)
- **AlphaFold**: ~1-2 hours for 500 residues (GPU)

### Accuracy

- **AlphaFold**: Highest accuracy, especially with templates
- **ESM**: Good for novel folds, faster
- **ColabFold**: Balanced, easy to use

### Resource Requirements

| Method | RAM | GPU | Time |
|--------|-----|-----|------|
| ESM-650M | 4 GB | Optional | 1-5 min |
| ESM-3B | 16 GB | Optional | 5-20 min |
| ColabFold | 16 GB | Required | 10-30 min |
| AlphaFold | 32 GB | Required | 1-2 hours |

## For More Details

- Main prediction tool: See [tools/scripts/protein_predict.py](../../tools/scripts/protein_predict.py)
- Evo2 documentation: See [evo2/README.md](../../evo2/README.md)
- ColabFold: See https://github.com/sokrypton/ColabFold
- ESM: See https://github.com/facebookresearch/esm
- PyMOL: See https://pymol.org/2/
- ChimeraX: See https://www.cgl.ucsf.edu/chimerax/
