from __future__ import annotations

import csv
import json
import shutil
import urllib.request
from pathlib import Path


CURATED_DATASETS = {
    "ebola": {
        "label": "Ebola virus",
        "geo_accession": "GSE114905",
        "title": "Transcriptome in Huh7 cells infected with 4 Ebolaviruses",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114905",
        "processed_file": {
            "filename": "GSE114905_RNA-476-477-478_EBOV_cells_D1-D3_GE_.xlsx",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114905/suppl/GSE114905_RNA-476-477-478_EBOV_cells_D1-D3_GE_.xlsx",
        },
        "reference": {
            "accession": "NC_002549.1",
            "download_by_default": True,
            "filename": "NC_002549.1.fa",
            "url": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_002549.1&rettype=fasta&retmode=text",
        },
    },
    "sars_cov_2": {
        "label": "SARS-CoV-2",
        "geo_accession": "GSE147507",
        "title": "SARS-CoV-2 launches a unique transcriptional signature from in vitro, ex vivo, and in vivo systems",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507",
        "processed_file": {
            "filename": "GSE147507_RawReadCounts_Human.tsv.gz",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147507/suppl/GSE147507_RawReadCounts_Human.tsv.gz",
        },
        "reference": {
            "accession": "NC_045512.2",
            "download_by_default": True,
            "filename": "NC_045512.2.fa",
            "url": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&retmode=text",
        },
    },
    "yeast": {
        "label": "Saccharomyces cerevisiae",
        "geo_accession": "GSE67149",
        "title": "Rpd3 drives transcriptional quiescence",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67149",
        "processed_file": {
            "filename": "GSE67149_Normalized_FPKM_Values.xlsx",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE67nnn/GSE67149/suppl/GSE67149_Normalized_FPKM_Values.xlsx",
        },
        "reference": {
            "accession": "SGD:S288C",
            "download_by_default": False,
            "filename": None,
            "url": "https://www.yeastgenome.org/strain/S288C",
        },
    },
}

HOST_REFERENCES = {
    "human": {
        "accession": "GRCh38",
        "download_by_default": False,
        "url": "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/",
    }
}

SAMPLE_ROWS = [
    {"dataset": "ebola", "sample": "EBOV_D1", "condition": "ebov_day_1", "comparison": "ebov_timecourse"},
    {"dataset": "ebola", "sample": "EBOV_D2", "condition": "ebov_day_2", "comparison": "ebov_timecourse"},
    {"dataset": "ebola", "sample": "EBOV_D3", "condition": "ebov_day_3", "comparison": "ebov_timecourse"},
    {"dataset": "sars_cov_2", "sample": "Series1_NHBE_Mock_1", "condition": "mock", "comparison": "nhbe_mock_vs_sars"},
    {"dataset": "sars_cov_2", "sample": "Series1_NHBE_Mock_2", "condition": "mock", "comparison": "nhbe_mock_vs_sars"},
    {"dataset": "sars_cov_2", "sample": "Series1_NHBE_Mock_3", "condition": "mock", "comparison": "nhbe_mock_vs_sars"},
    {"dataset": "sars_cov_2", "sample": "Series1_NHBE_SARS-CoV-2_1", "condition": "infected", "comparison": "nhbe_mock_vs_sars"},
    {"dataset": "sars_cov_2", "sample": "Series1_NHBE_SARS-CoV-2_2", "condition": "infected", "comparison": "nhbe_mock_vs_sars"},
    {"dataset": "sars_cov_2", "sample": "Series1_NHBE_SARS-CoV-2_3", "condition": "infected", "comparison": "nhbe_mock_vs_sars"},
    {"dataset": "yeast", "sample": "WT_Log", "condition": "log_phase", "comparison": "wt_log_vs_q"},
    {"dataset": "yeast", "sample": "WT_Q", "condition": "quiescence", "comparison": "wt_log_vs_q"},
]

COMPARISON_ROWS = [
    {
        "dataset": "ebola",
        "comparison": "ebov_timecourse",
        "group_a": "EBOV_D1",
        "group_b": "EBOV_D3",
        "note": "viral expression shift across day 1-3",
    },
    {
        "dataset": "sars_cov_2",
        "comparison": "nhbe_mock_vs_sars",
        "group_a": "mock",
        "group_b": "infected",
        "note": "NHBE mock versus SARS-CoV-2 infected samples",
    },
    {
        "dataset": "yeast",
        "comparison": "wt_log_vs_q",
        "group_a": "log_phase",
        "group_b": "quiescence",
        "note": "wild-type log phase versus quiescence",
    },
]


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_tsv(path: Path, rows: list[dict[str, object]]) -> Path:
    ensure_dir(path.parent)
    if not rows:
        path.write_text("", encoding="utf-8")
        return path
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_json(path: Path, payload: dict[str, object]) -> Path:
    ensure_dir(path.parent)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return path


def download_file(url: str, destination: Path, force: bool = False) -> Path:
    ensure_dir(destination.parent)
    if destination.exists() and not force:
        return destination
    with urllib.request.urlopen(url) as response, destination.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    return destination
