from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPORT_SCRIPT = ROOT / "projects" / "cross_system_benchmark" / "scripts" / "summarize_cross_systems.py"


def _load_module(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_integrated_summary_uses_dataset_summaries(tmp_path: Path) -> None:
    module = _load_module(REPORT_SCRIPT, "cross_system_reports")
    per_dataset_dir = tmp_path / "data" / "results" / "per_dataset"
    per_dataset_dir.mkdir(parents=True)

    for dataset, top_feature in {
        "ebola": "NP",
        "sars_cov_2": "IFIT1",
        "yeast": "AAC1",
    }.items():
        summary_json = per_dataset_dir / f"{dataset}_summary.json"
        summary_json.write_text(
            json.dumps(
                {
                    "dataset": dataset,
                    "top_feature": top_feature,
                    "sample_count": 2,
                    "feature_count": 2,
                }
            ),
            encoding="utf-8",
        )

    output_path = tmp_path / "data" / "results" / "integrated" / "cross_system_summary.md"
    module.write_integrated_summary(per_dataset_dir, output_path)

    assert output_path.exists()
    text = output_path.read_text(encoding="utf-8")
    assert "ebola" in text
    assert "sars_cov_2" in text
    assert "yeast" in text
