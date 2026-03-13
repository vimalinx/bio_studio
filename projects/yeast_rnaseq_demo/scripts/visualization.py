import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy import stats
import os

# Set style
sns.set_theme(style="whitegrid")
plt.rcParams['figure.figsize'] = (10, 8)

# Paths
INPUT_FILE = "projects/yeast_rnaseq_demo/data/results/counts_matrix_demo.csv"
OUTPUT_DIR = "projects/yeast_rnaseq_demo/data/results/plots"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("📊 Starting Data Visualization Analysis...")

# 1. Load Data
print(f"   Loading data from {INPUT_FILE}...")
df = pd.read_csv(INPUT_FILE, index_col=0)
print(f"   Shape: {df.shape}")

# Filter low counts (simple filtering)
df = df[df.sum(axis=1) > 10]
print(f"   After filtering: {df.shape}")

# Normalize (Simple CPM-like: Count Per Million)
# Add 1 to avoid log(0)
norm_df = df.div(df.sum()) * 1e6
log_df = np.log2(norm_df + 1)

# 2. PCA Analysis
print("   Generating PCA Plot...")
pca = PCA(n_components=2)
pca_result = pca.fit_transform(log_df.T) # Transpose for samples as rows

pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'], index=log_df.columns)
pca_df['Condition'] = ['WT', 'WT', 'MUT', 'MUT'] # Extract from names in real scenario

plt.figure(figsize=(8, 6))
sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Condition', style='Condition', s=200)
plt.title('PCA of Gene Expression')
plt.savefig(f"{OUTPUT_DIR}/01_pca_plot.png")
plt.close()
print(f"   ✅ Saved: {OUTPUT_DIR}/01_pca_plot.png")

# 3. Differential Expression (Simple T-test)
print("   Calculating Differential Expression...")
wt_cols = [c for c in df.columns if "WT" in c]
mut_cols = [c for c in df.columns if "MUT" in c]

results = []
for gene in log_df.index:
    wt_vals = log_df.loc[gene, wt_cols]
    mut_vals = log_df.loc[gene, mut_cols]
    
    # Calculate Log2 Fold Change (Mean MUT - Mean WT)
    log2fc = mut_vals.mean() - wt_vals.mean()
    
    # T-test
    t_stat, p_val = stats.ttest_ind(mut_vals, wt_vals)
    
    results.append({
        'GeneID': gene,
        'Log2FC': log2fc,
        'PValue': p_val
    })

res_df = pd.DataFrame(results).set_index('GeneID')
res_df['-Log10P'] = -np.log10(res_df['PValue'] + 1e-300) # Avoid inf

# 4. Volcano Plot
print("   Generating Volcano Plot...")
plt.figure(figsize=(10, 8))
# Significant genes
sig = (res_df['PValue'] < 0.05) & (abs(res_df['Log2FC']) > 1)
res_df['Significant'] = sig

sns.scatterplot(data=res_df, x='Log2FC', y='-Log10P', hue='Significant', palette={True:'red', False:'grey'}, alpha=0.6)
plt.title('Volcano Plot (WT vs MUT)')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 P-value')
plt.axvline(x=1, color='k', linestyle='--')
plt.axvline(x=-1, color='k', linestyle='--')
plt.axhline(y=-np.log10(0.05), color='k', linestyle='--')
plt.savefig(f"{OUTPUT_DIR}/02_volcano_plot.png")
plt.close()
print(f"   ✅ Saved: {OUTPUT_DIR}/02_volcano_plot.png")

# 5. Heatmap (Top 50 Variable Genes)
print("   Generating Heatmap...")
# Calculate variance
var_genes = log_df.var(axis=1).sort_values(ascending=False).head(50).index
heatmap_data = log_df.loc[var_genes]

# Z-score normalization for heatmap
# Note: zscore computes along axis=0 by default (columns), we want rows (genes)
heatmap_data_z = stats.zscore(heatmap_data, axis=1)
# Convert back to DataFrame to keep index/columns
heatmap_data_z = pd.DataFrame(heatmap_data_z, index=heatmap_data.index, columns=heatmap_data.columns)

plt.figure(figsize=(10, 12))
sns.heatmap(heatmap_data_z, cmap='vlag', center=0)
plt.title('Top 50 Variable Genes')
plt.savefig(f"{OUTPUT_DIR}/03_heatmap.png")
plt.close()
print(f"   ✅ Saved: {OUTPUT_DIR}/03_heatmap.png")

print("\n✅ Visualization Analysis Complete!")
