"""
test_env_validation 项目配置
"""

from pathlib import Path

PROJECT_NAME = "test_env_validation"
PROJECT_TYPE = "rnaseq"
PROJECT_DESCRIPTION = "Environment validation test project"

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
