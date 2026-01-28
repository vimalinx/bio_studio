# RNA-seq pipeline configuration
import os
from pathlib import Path

# Project paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = DATA_DIR / "results"
REF_DIR = DATA_DIR / "references"
LOG_DIR = PROJECT_ROOT / "logs"

# Ensure directories exist
for d in [PROCESSED_DIR, RESULTS_DIR, LOG_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Reference files
GENOME_FASTA = REF_DIR / "genome.fa"
GENOME_GTF = REF_DIR / "genome.gff" # Note: featureCounts can handle GFF
STAR_INDEX_DIR = REF_DIR / "star_index"

# Sample information (Auto-detect)
SAMPLES = sorted(list(set([f.name.split("_R1")[0] for f in RAW_DIR.glob("*_R1.fastq.gz")])))

# Tools configuration
THREADS = 8 # Adjust based on CPU
