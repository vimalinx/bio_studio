from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_VALIDATE_PATH = (
    ROOT / "projects" / "test_env_validation" / "scripts" / "run_validation.py"
)


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_build_mock_read_pair_uses_reverse_complement_for_read2() -> None:
    run_validation = _load_module(
        WORKSPACE_VALIDATE_PATH, "bio_studio_workspace_validation_mock_data_test"
    )
    fragment = "ATGCCGTTAACCGGTTAACCGGTTAACCGGTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGG"
    read1, read2 = run_validation.build_mock_read_pair(fragment, read_length=20)

    expected_r1 = fragment[:20]
    expected_r2 = run_validation.reverse_complement(fragment[-20:])

    assert read1 == expected_r1
    assert read2 == expected_r2
    assert read1 != read2


def test_test_bio_tools_uses_paired_end_featurecounts(monkeypatch) -> None:
    run_validation = _load_module(
        WORKSPACE_VALIDATE_PATH, "bio_studio_workspace_validation_featurecounts_test"
    )
    commands: list[str] = []

    monkeypatch.setattr(run_validation, "run_cmd", lambda cmd: commands.append(cmd))

    run_validation.test_bio_tools()

    featurecounts_cmd = next(cmd for cmd in commands if cmd.startswith("featureCounts "))
    assert " -p " in f" {featurecounts_cmd} "


def test_test_bio_tools_ignores_existing_multiqc_output(monkeypatch) -> None:
    run_validation = _load_module(
        WORKSPACE_VALIDATE_PATH, "bio_studio_workspace_validation_multiqc_test"
    )
    commands: list[str] = []

    monkeypatch.setattr(run_validation, "run_cmd", lambda cmd: commands.append(cmd))

    run_validation.test_bio_tools()

    multiqc_cmd = next(cmd for cmd in commands if cmd.startswith("multiqc "))
    assert "--ignore" in multiqc_cmd
    assert "multiqc_report" in multiqc_cmd
