from __future__ import annotations

from types import SimpleNamespace

from lib import project_validation


class _Context:
    def __init__(self) -> None:
        self.results: list[dict[str, str]] = []

    def add(self, name: str, status: str, message: str) -> None:
        self.results.append({"name": name, "status": status, "message": message})


def test_required_tools_for_rnaseq_project_type() -> None:
    config = SimpleNamespace(PROJECT_TYPE="rnaseq", PROJECT_NAME="demo")
    assert project_validation.required_tools_for_config(config) == [
        "fastp",
        "STAR",
        "featureCounts",
    ]


def test_required_tools_for_known_project_name() -> None:
    config = SimpleNamespace(PROJECT_TYPE="generic", PROJECT_NAME="ai_design_playground")
    assert project_validation.required_tools_for_config(config) == [
        "seqkit",
        "prodigal",
        "RNAfold",
    ]


def test_validate_cli_tools_reports_missing_tools_as_warn(monkeypatch) -> None:
    context = _Context()
    config = SimpleNamespace(PROJECT_TYPE="variant", PROJECT_NAME="demo")

    monkeypatch.setattr(project_validation, "ensure_workspace_path", lambda: None)
    monkeypatch.setattr(project_validation.shutil, "which", lambda tool: None)

    project_validation.validate_cli_tools(context, config)

    assert context.results == [
        {
            "name": "required_tools",
            "status": "WARN",
            "message": "missing on PATH: bwa, samtools, bcftools; available: none",
        }
    ]


def test_validate_cli_tools_accepts_iqtree_alias(monkeypatch) -> None:
    context = _Context()
    config = SimpleNamespace(PROJECT_TYPE="phylogeny", PROJECT_NAME="demo")

    monkeypatch.setattr(project_validation, "ensure_workspace_path", lambda: None)

    def fake_which(tool: str):
        if tool == "mafft":
            return "/opt/bin/mafft"
        if tool == "iqtree":
            return "/opt/bin/iqtree"
        return None

    monkeypatch.setattr(project_validation.shutil, "which", fake_which)

    project_validation.validate_cli_tools(context, config)

    assert context.results == [
        {
            "name": "required_tools",
            "status": "PASS",
            "message": "available on PATH: mafft, iqtree2->iqtree",
        }
    ]
