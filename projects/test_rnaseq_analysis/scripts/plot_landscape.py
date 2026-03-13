#!/usr/bin/env python3
"""
Plot Ebola Evolutionary Landscape based on Evo2 AI Scores
"""

import json
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path

# Paths
PROJECT_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_DIR / "data" / "results"
VCF_FILE = RESULTS_DIR / "variants.vcf"
JSON_FILE = RESULTS_DIR / "evo2_variant_effects.json"
OUTPUT_PLOT = RESULTS_DIR / "ebola_evolutionary_landscape.png"

# Gene Coordinates (Zaire ebolavirus, approx)
GENES = {
    "NP": (470, 3140, "#1f77b4"),
    "VP35": (3134, 4151, "#ff7f0e"),
    "VP40": (4479, 5459, "#2ca02c"),
    "GP": (6039, 8068, "#d62728"),  # Red for GP
    "VP30": (8509, 9375, "#9467bd"),
    "VP24": (10345, 11100, "#8c564b"),
    "L": (11581, 18219, "#e377c2"),  # Pink for L
    "Intergenic": (0, 0, "#7f7f7f"), # Grey for non-coding regions
}

def get_gene_at_pos(pos):
    for gene, (start, end, color) in GENES.items():
        if start <= pos <= end:
            return gene
    return "Intergenic"

def main():
    print("📊 Loading data...")
    
    # 1. Load Evo2 Scores
    if not JSON_FILE.exists():
        print(f"❌ JSON file not found: {JSON_FILE}")
        return

    with open(JSON_FILE, "r") as f:
        ai_data_raw = json.load(f)
    
    # Extract variants list from dictionary
    # Structure: {"model": "...", "variants": [...]}
    if isinstance(ai_data_raw, dict) and "variants" in ai_data_raw:
        ai_data = ai_data_raw["variants"]
    else:
        print(f"⚠️ Unexpected JSON structure. Keys found: {list(ai_data_raw.keys())}")
        return
    
    data = []
    for item in ai_data:
        # Normalize position (1-based)
        pos = item.get("pos") or item.get("position")
        score = item.get("score") or item.get("effect_score")
        
        # Evo2 score can be None if model failed to predict
        if score is None or pos is None:
            continue
            
        # Determine gene
        gene = get_gene_at_pos(pos)
        
        data.append({
            "Position": pos,
            "Score": score,
            "Gene": gene,
            "Label": f"{item.get('ref')}{pos}{item.get('alt')}"
        })
    
    df = pd.DataFrame(data)
    print(f"  Loaded {len(df)} variants.")
    
    # 2. Plotting
    print("🎨 Generating plot...")
    plt.figure(figsize=(14, 8))
    sns.set_style("whitegrid")
    
    # Draw Gene Regions (Background Rectangles)
    ylim = (df["Score"].min() - 0.2, df["Score"].max() + 0.2)
    for gene, (start, end, color) in GENES.items():
        plt.axvspan(start, end, alpha=0.1, color=color)
        plt.text((start+end)/2, ylim[1] - 0.1, gene, ha='center', fontsize=10, fontweight='bold', color=color)

    # Scatter Plot
    sns.scatterplot(
        data=df, 
        x="Position", 
        y="Score", 
        hue="Gene", 
        palette={g: GENES[g][2] for g in GENES}, 
        s=60, 
        edgecolor='k',
        alpha=0.8
    )
    
    # Highlight Top Hits (Negative Score = High Impact)
    # Sort by score ascending (most negative first)
    top_hits = df.sort_values("Score").head(8)
    
    from adjustText import adjust_text
    texts = []
    for _, row in top_hits.iterrows():
        texts.append(plt.text(row["Position"], row["Score"], row["Label"], fontsize=9, fontweight='bold'))
    
    try:
        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray'))
    except ImportError:
        print("⚠️ 'adjustText' not installed, labels might overlap. (pip install adjustText)")

    # Formatting
    plt.title("Ebola Virus Evolutionary Landscape (AI Predicted Impact)", fontsize=16)
    plt.xlabel("Genome Position (bp)", fontsize=12)
    plt.ylabel("Evo2 Effect Score (Log-Likelihood Ratio)", fontsize=12)
    plt.axhline(0, color='gray', linestyle='--')
    
    # Save
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300)
    print(f"✅ Plot saved to: {OUTPUT_PLOT}")

if __name__ == "__main__":
    main()
