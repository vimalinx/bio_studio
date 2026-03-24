"""
共享模板运行时：在工作区内为模板项目提供模块化步骤实现。
"""

from __future__ import annotations

from pathlib import Path
import sys


WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
if str(WORKSPACE_ROOT) not in sys.path:
    sys.path.insert(0, str(WORKSPACE_ROOT))

from lib.modules import (
    alignment,
    alignment_msa,
    gene_prediction,
    phylogeny,
    qc,
    rnaseq,
    samtools_wrapper,
    sequence_analysis,
    variant,
)


def _as_path(value) -> Path:
    return value if isinstance(value, Path) else Path(value)


def _find_paired_reads(raw_dir: Path, read1_pattern: str, read2_pattern: str) -> list[tuple[Path, Path | None]]:
    pairs: list[tuple[Path, Path | None]] = []
    for read1 in sorted(raw_dir.glob(read1_pattern)):
        mate_name = read1.name.replace("_R1", "_R2")
        read2 = read1.with_name(mate_name)
        if read2_pattern and not read2.exists():
            read2 = None
        pairs.append((read1, read2))
    return pairs


def _sample_name_from_read1(read1: Path) -> str:
    name = read1.name
    suffixes = [
        "_R1.clean.fastq.gz",
        "_R1.fastq.gz",
        "_R1.clean.fq.gz",
        "_R1.fq.gz",
        "_R1.clean.fastq",
        "_R1.fastq",
        "_R1.clean.fq",
        "_R1.fq",
    ]
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return read1.stem


def _get_rnaseq_reference(config) -> Path:
    reference = getattr(config, "REFERENCE_GENOME", None) or getattr(config, "GENOME_FASTA", None)
    return _as_path(reference) if reference is not None else Path("")


def _get_rnaseq_annotation(config) -> Path:
    annotation = getattr(config, "ANNOTATION_GTF", None) or getattr(config, "GENOME_GTF", None)
    return _as_path(annotation) if annotation is not None else Path("")


def _get_rnaseq_index_dir(config) -> Path:
    value = getattr(config, "STAR_INDEX_DIR", None)
    if value is not None:
        return _as_path(value)
    return _as_path(config.REFERENCES_DIR) / "star_index"


def _get_rnaseq_alignment_dir(config) -> Path:
    value = getattr(config, "ALIGNMENT_OUTPUT_DIR", None)
    if value is not None:
        return _as_path(value)
    return _as_path(config.PROCESSED_DIR) / "aligned"


def _get_rnaseq_input_pairs(config, prefer_cleaned: bool = True) -> list[tuple[Path, Path | None]]:
    if prefer_cleaned:
        qc_dir = _as_path(config.PROCESSED_DIR) / "qc"
        if qc_dir.exists():
            cleaned_pairs = _find_paired_reads(qc_dir, "*_R1.clean.fastq.gz", "*_R2.clean.fastq.gz")
            if cleaned_pairs:
                return cleaned_pairs
    return _find_paired_reads(
        _as_path(config.RAW_DIR),
        getattr(config, "READ1_PATTERN", "*_R1.fastq.gz"),
        getattr(config, "READ2_PATTERN", "*_R2.fastq.gz"),
    )


def _run_rnaseq_build_index(config) -> bool | None:
    reference = _get_rnaseq_reference(config)
    annotation = _get_rnaseq_annotation(config)
    if not reference.exists() or not annotation.exists():
        return None

    genome_dir = _get_rnaseq_index_dir(config)
    if genome_dir.exists() and any(genome_dir.iterdir()):
        return True

    rnaseq.build_star_index(
        str(reference),
        str(genome_dir),
        annotation_gtf=str(annotation),
        threads=getattr(config, "THREADS", 4),
        extra_args=getattr(config, "STAR_INDEX_EXTRA_ARGS", None),
    )
    return True


def _run_rnaseq_alignment(config) -> bool | None:
    genome_dir = _get_rnaseq_index_dir(config)
    if not genome_dir.exists():
        return None

    input_pairs = _get_rnaseq_input_pairs(config, prefer_cleaned=True)
    if not input_pairs:
        return None

    output_dir = _get_rnaseq_alignment_dir(config)
    output_dir.mkdir(parents=True, exist_ok=True)

    read_files_command = getattr(config, "STAR_READ_FILES_COMMAND", None)
    if read_files_command is None and any(str(read1).endswith(".gz") for read1, _ in input_pairs):
        read_files_command = "zcat"

    align_extra_args = getattr(
        config,
        "STAR_ALIGN_EXTRA_ARGS",
        ["--outSAMtype", "BAM", "SortedByCoordinate"],
    )

    for read1, read2 in input_pairs:
        sample_name = _sample_name_from_read1(read1)
        output_prefix = output_dir / f"{sample_name}_"
        bam_path = output_dir / f"{sample_name}_Aligned.sortedByCoord.out.bam"
        rnaseq.align_star(
            str(genome_dir),
            str(read1),
            read2=str(read2) if read2 else None,
            output_prefix=str(output_prefix),
            threads=getattr(config, "THREADS", 4),
            read_files_command=read_files_command,
            extra_args=align_extra_args,
        )
        samtools_wrapper.index_bam(str(bam_path))
    return True


def _run_rnaseq_quantification(config) -> bool | None:
    annotation = _get_rnaseq_annotation(config)
    if not annotation.exists():
        return None

    alignment_dir = _get_rnaseq_alignment_dir(config)
    bam_files = sorted(alignment_dir.glob("*_Aligned.sortedByCoord.out.bam"))
    if not bam_files:
        processed_dir = _as_path(config.PROCESSED_DIR)
        if alignment_dir != processed_dir:
            bam_files = sorted(processed_dir.glob("*_Aligned.sortedByCoord.out.bam"))
    if not bam_files:
        return None

    input_pairs = _get_rnaseq_input_pairs(config, prefer_cleaned=True)
    counts_path = _as_path(getattr(config, "COUNTS_PATH", _as_path(config.RESULTS_DIR) / "counts.txt"))
    rnaseq.run_featurecounts(
        str(annotation),
        [str(path) for path in bam_files],
        str(counts_path),
        threads=getattr(config, "THREADS", 4),
        paired_end=bool(input_pairs) and all(read2 is not None for _, read2 in input_pairs),
        feature_type=getattr(config, "FEATURECOUNTS_FEATURE_TYPE", None),
        attribute_type=getattr(config, "FEATURECOUNTS_ATTRIBUTE_TYPE", None),
        extra_args=getattr(config, "FEATURECOUNTS_EXTRA_ARGS", None),
    )
    return True


def _run_rnaseq_quality_control(config) -> bool | None:
    raw_dir = _as_path(config.RAW_DIR)
    pairs = _find_paired_reads(
        raw_dir,
        getattr(config, "READ1_PATTERN", "*_R1.fastq.gz"),
        getattr(config, "READ2_PATTERN", "*_R2.fastq.gz"),
    )
    if not pairs:
        return None

    qc_dir = _as_path(config.PROCESSED_DIR) / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    for read1, read2 in pairs:
        output_read1 = qc_dir / read1.name.replace("_R1", "_R1.clean")
        output_read2 = None
        if read2 is not None:
            output_read2 = qc_dir / read2.name.replace("_R2", "_R2.clean")
        stem = read1.name.replace("_R1.fastq.gz", "").replace("_R1.fastq", "")
        html_report = qc_dir / f"{stem}.fastp.html"
        json_report = qc_dir / f"{stem}.fastp.json"
        qc.run_fastp(
            str(read1),
            input_read2=str(read2) if read2 else None,
            output_read1=str(output_read1),
            output_read2=str(output_read2) if output_read2 else None,
            html_report=str(html_report),
            json_report=str(json_report),
            threads=getattr(config, "THREADS", 4),
        )
    return True


def _run_rnaseq_main_analysis(config) -> bool | None:
    if _run_rnaseq_build_index(config) is None:
        return None
    if _run_rnaseq_alignment(config) is None:
        return None
    return _run_rnaseq_quantification(config)


def _run_variant_main_analysis(config) -> bool | None:
    raw_dir = _as_path(config.RAW_DIR)
    reference = _as_path(getattr(config, "REFERENCE_GENOME", ""))
    if not reference.exists():
        return None

    pairs = _find_paired_reads(
        raw_dir,
        getattr(config, "READ1_PATTERN", "*_R1.fastq.gz"),
        getattr(config, "READ2_PATTERN", "*_R2.fastq.gz"),
    )
    if not pairs:
        return None

    processed_dir = _as_path(config.PROCESSED_DIR)
    results_dir = _as_path(config.RESULTS_DIR)
    processed_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    alignment.build_bwa_index(str(reference))

    bam_files: list[Path] = []
    for read1, read2 in pairs:
        sample_name = read1.name.replace("_R1.fastq.gz", "")
        sam_path = processed_dir / f"{sample_name}.sam"
        bam_path = processed_dir / f"{sample_name}.bam"
        sorted_bam = processed_dir / f"{sample_name}.sorted.bam"
        alignment.align_bwa_mem(
            str(reference),
            str(read1),
            reads2=str(read2) if read2 else None,
            output_sam=str(sam_path),
            threads=getattr(config, "THREADS", 4),
        )
        samtools_wrapper.sam_to_bam(str(sam_path), output_bam=str(bam_path), threads=getattr(config, "THREADS", 4))
        samtools_wrapper.sort_bam(str(bam_path), output_bam=str(sorted_bam), threads=getattr(config, "THREADS", 4))
        samtools_wrapper.index_bam(str(sorted_bam))
        bam_files.append(sorted_bam)

    raw_vcf = results_dir / "variants.raw.vcf"
    filtered_vcf = results_dir / "variants.filtered.vcf"
    variant.call_variants(str(reference), str(bam_files[0]), str(raw_vcf), threads=getattr(config, "THREADS", 4))
    variant.filter_vcf(str(raw_vcf), str(filtered_vcf), filters=getattr(config, "FILTER_EXPRESSION", None))
    return True


def _run_phylogeny_main_analysis(config) -> bool | None:
    input_fasta = _as_path(getattr(config, "ALIGNMENT_INPUT", ""))
    if not input_fasta.exists():
        return None

    results_dir = _as_path(config.RESULTS_DIR)
    results_dir.mkdir(parents=True, exist_ok=True)
    alignment_path = results_dir / "alignment.fasta"
    alignment_msa.align_mafft(
        str(input_fasta),
        str(alignment_path),
        threads=getattr(config, "THREADS", 4),
    )
    phylogeny.build_tree(
        str(alignment_path),
        model=getattr(config, "TREE_MODEL", "MFP"),
        bootstrap=getattr(config, "BOOTSTRAP", 1000),
        threads=getattr(config, "THREADS", 4),
    )
    return True


def run_ai_design_playground_analysis(config, dna_fa, rna_sequences):
    """
    复用共享 modules 执行 ai_design_playground 的本地工具链部分。
    """
    dna_fa = _as_path(dna_fa)
    if not dna_fa.exists():
        return None

    processed_dir = _as_path(config.PROCESSED_DIR)
    results_dir = _as_path(config.RESULTS_DIR)
    processed_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    seqkit_txt = results_dir / "seqkit_stats.txt"
    sequence_analysis.run_seqkit_stats(str(dna_fa), output_txt=str(seqkit_txt))

    proteins = processed_dir / "prodigal_proteins.faa"
    genes = processed_dir / "prodigal_genes.fna"
    out = processed_dir / "prodigal.out"
    gene_prediction.predict_genes(
        str(dna_fa),
        output_gff=str(out),
        output_fna=str(genes),
        output_faa=str(proteins),
        mode=getattr(config, "PRODIGAL_MODE", "single"),
    )

    rnafold_results = []
    for seq_id, rna_seq in rna_sequences:
        fold = sequence_analysis.run_rnafold(rna_seq, sequence_id=seq_id)
        structure = fold.get("structure")
        if not isinstance(structure, str):
            structure = None
        mfe = fold.get("mfe")
        if not isinstance(mfe, float):
            mfe = None
        rnafold_results.append(
            {
                "id": seq_id,
                "length_nt": len(rna_seq),
                "structure": structure,
                "mfe": mfe,
            }
        )

    return {
        "seqkit_stats": str(seqkit_txt),
        "prodigal": {
            "proteins": str(proteins),
            "genes": str(genes),
            "out": str(out),
        },
        "rnafold": rnafold_results,
    }


STEP_HANDLERS = {
    ("rnaseq", "build_index"): _run_rnaseq_build_index,
    ("rnaseq", "alignment"): _run_rnaseq_alignment,
    ("rnaseq", "quantification"): _run_rnaseq_quantification,
    ("rnaseq", "quality_control"): _run_rnaseq_quality_control,
    ("rnaseq", "main_analysis"): _run_rnaseq_main_analysis,
    ("variant", "main_analysis"): _run_variant_main_analysis,
    ("phylogeny", "main_analysis"): _run_phylogeny_main_analysis,
}


def run_shared_step(config, step_name: str) -> bool | None:
    """
    执行共享模块化步骤。返回 None 表示当前步骤仍应回退到模板内置提示逻辑。
    """
    project_type = getattr(config, "PROJECT_TYPE", "generic")
    handler = STEP_HANDLERS.get((project_type, step_name))
    if handler is None:
        return None
    return handler(config)
