"""
cross_system_benchmark 项目配置
"""

from pathlib import Path

PROJECT_NAME = "cross_system_benchmark"
PROJECT_TYPE = "generic"
PROJECT_DESCRIPTION = "Cross-system benchmark panel for Ebola, SARS-CoV-2, and yeast"

SCRIPTS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPTS_DIR.parent
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = DATA_DIR / "results"
REFERENCES_DIR = DATA_DIR / "references"
METADATA_DIR = PROJECT_ROOT / "metadata"
LOGS_DIR = PROJECT_ROOT / "logs"
NOTEBOOKS_DIR = PROJECT_ROOT / "notebooks"

RAW_DATASET_DIRS = {
    "ebola": RAW_DIR / "ebola",
    "sars_cov_2": RAW_DIR / "sars_cov_2",
    "yeast": RAW_DIR / "yeast",
}
PROCESSED_DATASET_DIRS = {
    "ebola": PROCESSED_DIR / "ebola",
    "sars_cov_2": PROCESSED_DIR / "sars_cov_2",
    "yeast": PROCESSED_DIR / "yeast",
}
PER_DATASET_RESULTS_DIR = RESULTS_DIR / "per_dataset"
INTEGRATED_RESULTS_DIR = RESULTS_DIR / "integrated"

SAMPLES = []
REFERENCE_GENOME = None
THREADS = 4
ANALYSIS_PARAMETERS = {}
PRIMARY_INPUT = RAW_DIR
DELIVERABLE = INTEGRATED_RESULTS_DIR / "cross_system_summary.md"
REQUIRED_TOOLS = ["python"]
REQUIRED_PYTHON_PACKAGES = ["pandas", "yaml", "openpyxl"]
