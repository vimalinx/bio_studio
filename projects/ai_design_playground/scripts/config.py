"""
ai_design_playground 项目配置
"""

from pathlib import Path

PROJECT_NAME = "ai_design_playground"
PROJECT_TYPE = "generic"
PROJECT_DESCRIPTION = "In-silico playground: generate toy sequences, analyze, and validate with local tools/AI models."

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

SAFETY_NOTICE = "In-silico only. Do not use outputs to create/modify organisms."

N_SEQUENCES = 8
DNA_LENGTH_BP = 900
GC_TARGET = 0.5
RANDOM_SEED = 1337

AVOID_RESTRICTION_SITES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "XhoI": "CTCGAG",
    "NdeI": "CATATG",
    "NotI": "GCGGCCGC",
}

PRODIGAL_MODE = "meta"

USE_ESM_CONTACTS_DEFAULT = False
ESM_MODEL_NAME_DEFAULT = "esm2_t6_8M_UR50D"
