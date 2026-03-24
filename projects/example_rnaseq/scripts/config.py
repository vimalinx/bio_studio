"""
example_rnaseq 项目配置
"""

from pathlib import Path

PROJECT_NAME = 'example_rnaseq'
PROJECT_TYPE = 'rnaseq'
PROJECT_DESCRIPTION = '示例RNA-seq分析项目'

SCRIPTS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPTS_DIR.parent
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = DATA_DIR / "results"
REFERENCES_DIR = DATA_DIR / "references"
LOGS_DIR = PROJECT_ROOT / "logs"
NOTEBOOKS_DIR = PROJECT_ROOT / "notebooks"

SAMPLES = []
REFERENCE_GENOME = None
THREADS = 4
ANALYSIS_PARAMETERS = {}

READ1_PATTERN = "*_R1.fastq.gz"
READ2_PATTERN = "*_R2.fastq.gz"
REFERENCE_GENOME = REFERENCES_DIR / "genome.fa"
ANNOTATION_GTF = REFERENCES_DIR / "annotation.gtf"
ALIGNER = "STAR"
QUANTIFIER = "featureCounts"
QC_TOOL = "fastp"
REPORT_TOOL = "MultiQC"
