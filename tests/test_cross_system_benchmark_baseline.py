from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BASELINE_SCRIPT = ROOT / "projects" / "cross_system_benchmark" / "scripts" / "run_baseline.py"


def _load_module(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_ebola_fixture(path: Path) -> None:
    with pd.ExcelWriter(path) as writer:
        for day, values in [
            ("RNA-476 EBOV cells D1 (GE)", [100.0, 70.0]),
            ("RNA-477 EBOV cells D2 (GE)", [120.0, 60.0]),
            ("RNA-478 EBOV cells D3 (GE)", [140.0, 50.0]),
        ]:
            pd.DataFrame(
                {
                    "Name": ["NP", "GP"],
                    "TPM": values,
                    "Expression value": [1000, 500],
                }
            ).to_excel(writer, sheet_name=day, index=False)


def _write_sars_fixture(path: Path) -> None:
    pd.DataFrame(
        {
            "Unnamed: 0": ["IFIT1", "CXCL10"],
            "Series1_NHBE_Mock_1": [10, 20],
            "Series1_NHBE_Mock_2": [11, 19],
            "Series1_NHBE_SARS-CoV-2_1": [100, 180],
            "Series1_NHBE_SARS-CoV-2_2": [95, 175],
        }
    ).to_csv(path, sep="\t", index=False)


def _write_yeast_fixture(path: Path) -> None:
    pd.DataFrame(
        {
            "Gene": ["AAC1", "AAC3"],
            "WT_Log": [4.0, 3.0],
            "WT_Q": [5.5, -2.0],
            "Rpd3_Log": [4.5, 2.5],
            "Rpd3_Q": [6.0, 1.5],
        }
    ).to_excel(path, sheet_name="Normalized_FPKM_Values", index=False)


def test_normalizers_write_common_expression_outputs(tmp_path: Path) -> None:
    module = _load_module(BASELINE_SCRIPT, "cross_system_baseline")

    ebola_input = tmp_path / "ebola.xlsx"
    sars_input = tmp_path / "sars.tsv"
    yeast_input = tmp_path / "yeast.xlsx"
    _write_ebola_fixture(ebola_input)
    _write_sars_fixture(sars_input)
    _write_yeast_fixture(yeast_input)

    ebola_output = module.normalize_ebola_workbook(ebola_input, tmp_path / "ebola")
    sars_output = module.normalize_sars_counts(sars_input, tmp_path / "sars_cov_2")
    yeast_output = module.normalize_yeast_workbook(yeast_input, tmp_path / "yeast")

    for output in (ebola_output, sars_output, yeast_output):
        expression_path = output["expression"]
        samples_path = output["samples"]
        summary_path = output["summary"]
        assert expression_path.exists()
        assert samples_path.exists()
        assert summary_path.exists()

    ebola_df = pd.read_csv(ebola_output["expression"], sep="\t")
    sars_df = pd.read_csv(sars_output["expression"], sep="\t")
    yeast_df = pd.read_csv(yeast_output["expression"], sep="\t")

    assert {"feature", "sample", "value", "dataset"} <= set(ebola_df.columns)
    assert {"feature", "sample", "value", "dataset"} <= set(sars_df.columns)
    assert {"feature", "sample", "value", "dataset"} <= set(yeast_df.columns)
