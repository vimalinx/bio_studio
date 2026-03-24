"""
yeast_rnaseq_demo 项目配置
"""

from pathlib import Path

PROJECT_NAME = "yeast_rnaseq_demo"
PROJECT_TYPE = "rnaseq"
PROJECT_DESCRIPTION = "Yeast RNA-seq demonstration pipeline using STAR and featureCounts."

SCRIPTS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPTS_DIR.parent
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = DATA_DIR / "results"
REFERENCES_DIR = DATA_DIR / "references"
LOGS_DIR = PROJECT_ROOT / "logs"
NOTEBOOKS_DIR = PROJECT_ROOT / "notebooks"

REF_DIR = REFERENCES_DIR
LOG_DIR = LOGS_DIR

for d in [PROCESSED_DIR, RESULTS_DIR, LOGS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

SAMPLES = sorted({f.name.split("_R1")[0] for f in RAW_DIR.glob("*_R1.fastq.gz")})
REFERENCE_GENOME = REFERENCES_DIR / "genome.fa"
THREADS = 8
ANALYSIS_PARAMETERS = {}

READ1_PATTERN = "*_R1.fastq.gz"
READ2_PATTERN = "*_R2.fastq.gz"
GENOME_FASTA = REFERENCE_GENOME
GENOME_GTF = REFERENCES_DIR / "genome.gff"
ANNOTATION_GTF = GENOME_GTF
STAR_INDEX_DIR = REF_DIR / "star_index"
ALIGNMENT_OUTPUT_DIR = PROCESSED_DIR
COUNTS_PATH = RESULTS_DIR / "counts.txt"
ALIGNER = "STAR"
QUANTIFIER = "featureCounts"
STAR_INDEX_EXTRA_ARGS = ["--genomeSAindexNbases", "10", "--sjdbGTFtagExonParentTranscript", "Parent"]
STAR_ALIGN_EXTRA_ARGS = ["--outSAMtype", "BAM", "SortedByCoordinate", "--outSAMunmapped", "Within", "--outSAMattributes", "Standard"]
FEATURECOUNTS_FEATURE_TYPE = "gene"
FEATURECOUNTS_ATTRIBUTE_TYPE = "ID"
