"""
ai_design_playground 项目配置
"""

# 项目配置
PROJECT_NAME = "ai_design_playground"
PROJECT_TYPE = "generic"

# 数据路径
DATA_DIR = "data"
RAW_DIR = f"{DATA_DIR}/raw"
PROCESSED_DIR = f"{DATA_DIR}/processed"
RESULTS_DIR = f"{DATA_DIR}/results"
REFERENCES_DIR = f"{DATA_DIR}/references"

# 样本配置
SAMPLES = []

# 参考基因组
REFERENCE_GENOME = None

# 工具配置
THREADS = 4

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
