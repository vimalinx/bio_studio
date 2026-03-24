from __future__ import annotations

import importlib
import importlib.util
import os
from pathlib import Path
from types import SimpleNamespace

from lib.create_project import build_pipeline_content
from lib.modules.utils import run_command


def _load_template_runtime():
    spec = importlib.util.find_spec("lib.template_runtime")
    assert spec is not None, "lib.template_runtime should exist"
    return importlib.import_module("lib.template_runtime")


def _config(tmp_path: Path, project_type: str) -> SimpleNamespace:
    project_root = tmp_path / project_type
    raw_dir = project_root / "data" / "raw"
    processed_dir = project_root / "data" / "processed"
    results_dir = project_root / "data" / "results"
    references_dir = project_root / "data" / "references"
    logs_dir = project_root / "logs"

    for path in (raw_dir, processed_dir, results_dir, references_dir, logs_dir):
        path.mkdir(parents=True, exist_ok=True)

    return SimpleNamespace(
        PROJECT_NAME=f"{project_type}_demo",
        PROJECT_TYPE=project_type,
        PROJECT_ROOT=project_root,
        RAW_DIR=raw_dir,
        PROCESSED_DIR=processed_dir,
        RESULTS_DIR=results_dir,
        REFERENCES_DIR=references_dir,
        LOGS_DIR=logs_dir,
        THREADS=4,
        READ1_PATTERN="*_R1.fastq.gz",
        READ2_PATTERN="*_R2.fastq.gz",
        REFERENCE_GENOME=references_dir / "genome.fa",
        ANNOTATION_GTF=references_dir / "annotation.gtf",
        FILTER_EXPRESSION="QUAL>30 && DP>10",
        ALIGNMENT_INPUT=raw_dir / "sequences.fasta",
        TREE_MODEL="MFP",
        BOOTSTRAP=1000,
    )


def test_run_command_supports_shell_redirection(tmp_path: Path) -> None:
    output_path = tmp_path / "shell.txt"

    result = run_command(f"printf 'bio' > {output_path}", shell=True, capture_output=False)

    assert result.returncode == 0
    assert output_path.read_text(encoding="utf-8") == "bio"


def test_run_command_prepends_detected_bio_bin_to_subprocess_env(
    monkeypatch, tmp_path: Path
) -> None:
    fake_home = tmp_path / "home"
    fake_bio_python = fake_home / "miniforge3" / "envs" / "bio" / "bin" / "python"
    fake_bio_python.parent.mkdir(parents=True, exist_ok=True)
    fake_bio_python.write_text("#!/bin/sh\n", encoding="utf-8")
    fake_bio_python.chmod(0o755)

    captured_kwargs: dict[str, object] = {}

    class _Completed:
        returncode = 0
        stdout = ""
        stderr = ""

    def fake_run(cmd, **kwargs):
        captured_kwargs.update(kwargs)
        return _Completed()

    monkeypatch.setenv("HOME", str(fake_home))
    monkeypatch.setattr(importlib.import_module("lib.modules.utils").subprocess, "run", fake_run)

    result = run_command(["fastp", "--version"])

    assert result.returncode == 0
    env = captured_kwargs["env"]
    assert isinstance(env, dict)
    assert env["PATH"].split(os.pathsep)[0] == str(fake_bio_python.parent)


def test_generated_rnaseq_pipeline_prefers_shared_runtime() -> None:
    pipeline_text = build_pipeline_content("rnaseq_demo", "rnaseq")

    assert "load_template_runtime" in pipeline_text
    assert "run_shared_step" in pipeline_text


def test_modules_package_exports_rnaseq_module() -> None:
    modules = importlib.import_module("lib.modules")
    assert hasattr(modules, "rnaseq")


def test_template_runtime_runs_rnaseq_quality_control_with_modules(
    tmp_path: Path,
    monkeypatch,
) -> None:
    runtime = _load_template_runtime()
    config = _config(tmp_path, "rnaseq")
    read1 = Path(config.RAW_DIR) / "sample_R1.fastq.gz"
    read2 = Path(config.RAW_DIR) / "sample_R2.fastq.gz"
    read1.write_text("r1", encoding="utf-8")
    read2.write_text("r2", encoding="utf-8")

    calls: list[tuple[str, tuple, dict]] = []

    def fake_run_fastp(*args, **kwargs):
        calls.append(("fastp", args, kwargs))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(runtime.qc, "run_fastp", fake_run_fastp)

    assert runtime.run_shared_step(config, "quality_control") is True
    assert calls
    assert calls[0][0] == "fastp"
    assert Path(calls[0][2]["output_read1"]).name == "sample_R1.clean.fastq.gz"


def test_template_runtime_runs_variant_main_analysis_with_modules(
    tmp_path: Path,
    monkeypatch,
) -> None:
    runtime = _load_template_runtime()
    config = _config(tmp_path, "variant")
    read1 = Path(config.RAW_DIR) / "sample_R1.fastq.gz"
    read2 = Path(config.RAW_DIR) / "sample_R2.fastq.gz"
    reference = Path(config.REFERENCES_DIR) / "genome.fa"
    read1.write_text("r1", encoding="utf-8")
    read2.write_text("r2", encoding="utf-8")
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")

    called: list[str] = []

    def fake_build_bwa_index(*args, **kwargs):
        called.append("build_bwa_index")
        return SimpleNamespace(returncode=0)

    def fake_align_bwa_mem(*args, **kwargs):
        called.append("align_bwa_mem")
        Path(kwargs["output_sam"]).write_text("@HD\tVN:1.6\n", encoding="utf-8")
        return SimpleNamespace(returncode=0)

    def fake_sam_to_bam(*args, **kwargs):
        called.append("sam_to_bam")
        Path(kwargs["output_bam"]).write_text("bam", encoding="utf-8")
        return SimpleNamespace(returncode=0)

    def fake_sort_bam(*args, **kwargs):
        called.append("sort_bam")
        Path(kwargs["output_bam"]).write_text("sorted", encoding="utf-8")
        return SimpleNamespace(returncode=0)

    def fake_index_bam(*args, **kwargs):
        called.append("index_bam")
        return SimpleNamespace(returncode=0)

    def fake_call_variants(*args, **kwargs):
        called.append("call_variants")
        Path(args[2]).write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
        return SimpleNamespace(returncode=0)

    def fake_filter_vcf(*args, **kwargs):
        called.append("filter_vcf")
        Path(args[1]).write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(runtime.alignment, "build_bwa_index", fake_build_bwa_index)
    monkeypatch.setattr(runtime.alignment, "align_bwa_mem", fake_align_bwa_mem)
    monkeypatch.setattr(runtime.samtools_wrapper, "sam_to_bam", fake_sam_to_bam)
    monkeypatch.setattr(runtime.samtools_wrapper, "sort_bam", fake_sort_bam)
    monkeypatch.setattr(runtime.samtools_wrapper, "index_bam", fake_index_bam)
    monkeypatch.setattr(runtime.variant, "call_variants", fake_call_variants)
    monkeypatch.setattr(runtime.variant, "filter_vcf", fake_filter_vcf)

    assert runtime.run_shared_step(config, "main_analysis") is True
    assert called == [
        "build_bwa_index",
        "align_bwa_mem",
        "sam_to_bam",
        "sort_bam",
        "index_bam",
        "call_variants",
        "filter_vcf",
    ]


def test_template_runtime_runs_phylogeny_main_analysis_with_modules(
    tmp_path: Path,
    monkeypatch,
) -> None:
    runtime = _load_template_runtime()
    config = _config(tmp_path, "phylogeny")
    alignment_input = Path(config.ALIGNMENT_INPUT)
    alignment_input.write_text(">a\nACGT\n>b\nACGA\n", encoding="utf-8")

    called: list[str] = []

    def fake_align_mafft(*args, **kwargs):
        called.append("align_mafft")
        Path(args[1]).write_text(">a\nACGT\n>b\nACGA\n", encoding="utf-8")
        return SimpleNamespace(returncode=0)

    def fake_build_tree(*args, **kwargs):
        called.append("build_tree")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(runtime.alignment_msa, "align_mafft", fake_align_mafft)
    monkeypatch.setattr(runtime.phylogeny, "build_tree", fake_build_tree)

    assert runtime.run_shared_step(config, "main_analysis") is True
    assert called == ["align_mafft", "build_tree"]
