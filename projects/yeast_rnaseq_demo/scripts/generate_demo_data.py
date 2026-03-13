import pandas as pd
import numpy as np
import os

# Configuration
OUTPUT_FILE = "projects/yeast_rnaseq_demo/data/results/counts_matrix_demo.csv"
GENE_COUNT = 6000
SAMPLES = ["WT_1", "WT_2", "MUT_1", "MUT_2"]

np.random.seed(42)

print("🧪 Generating synthetic RNA-seq data for demonstration...")

# 1. Generate Gene IDs (Yeast-like)
genes = [f"Y{c}{str(i).zfill(3)}W" for c in "ABCDEF" for i in range(1, 1001)][:GENE_COUNT]

# 2. Base expression levels (Log-normal distribution)
base_expr = np.random.lognormal(mean=2, sigma=1, size=GENE_COUNT)

# 3. Create Sample Data
data = {}
for sample in SAMPLES:
    # Add random technical noise
    noise = np.random.normal(loc=1, scale=0.1, size=GENE_COUNT)
    expr = base_expr * noise
    data[sample] = expr

# 4. Add Biological Signal (Differential Expression in MUT)
# Select 200 genes to be Up-regulated in MUT
up_genes_idx = np.random.choice(GENE_COUNT, 200, replace=False)
for idx in up_genes_idx:
    fold_change = np.random.uniform(2, 5) # 2x to 5x upregulation
    data["MUT_1"][idx] *= fold_change
    data["MUT_2"][idx] *= fold_change

# Select 200 genes to be Down-regulated in MUT
remaining_idx = list(set(range(GENE_COUNT)) - set(up_genes_idx))
down_genes_idx = np.random.choice(remaining_idx, 200, replace=False)
for idx in down_genes_idx:
    fold_change = np.random.uniform(0.2, 0.5) # 0.2x to 0.5x downregulation
    data["MUT_1"][idx] *= fold_change
    data["MUT_2"][idx] *= fold_change

# 5. Convert to Integer Counts
df = pd.DataFrame(data, index=genes)
df = df.round().astype(int)
df.index.name = "GeneID"

# Save
df.to_csv(OUTPUT_FILE)
print(f"✅ Generated demo data: {OUTPUT_FILE}")
print(f"   Genes: {GENE_COUNT}")
print(f"   Samples: {SAMPLES}")
print(f"   Simulated DE genes: ~400")
