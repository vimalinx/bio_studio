"""
RNA-seq模块 - STAR / featureCounts 包装器
"""

from pathlib import Path

from .utils import check_file_exists, ensure_dir, run_command


def build_star_index(
    reference_fasta,
    genome_dir,
    annotation_gtf=None,
    threads=4,
    sjdb_overhang=99,
    extra_args=None,
):
    """
    构建 STAR 索引
    """
    check_file_exists(reference_fasta)
    ensure_dir(genome_dir)

    cmd = [
        "STAR",
        "--runThreadN",
        str(threads),
        "--runMode",
        "genomeGenerate",
        "--genomeDir",
        str(genome_dir),
        "--genomeFastaFiles",
        str(reference_fasta),
    ]

    if annotation_gtf:
        check_file_exists(annotation_gtf)
        cmd.extend(["--sjdbGTFfile", str(annotation_gtf)])
        cmd.extend(["--sjdbOverhang", str(sjdb_overhang)])

    if extra_args:
        cmd.extend(extra_args)

    return run_command(cmd, capture_output=False)


def align_star(
    genome_dir,
    read1,
    read2=None,
    output_prefix=None,
    threads=4,
    read_files_command=None,
    extra_args=None,
):
    """
    使用 STAR 进行读段比对
    """
    check_file_exists(read1)
    if read2:
        check_file_exists(read2)

    cmd = [
        "STAR",
        "--runThreadN",
        str(threads),
        "--genomeDir",
        str(genome_dir),
        "--readFilesIn",
        str(read1),
    ]

    if read2:
        cmd.append(str(read2))

    if read_files_command:
        cmd.extend(["--readFilesCommand", str(read_files_command)])

    if output_prefix:
        ensure_dir(Path(output_prefix).parent)
        cmd.extend(["--outFileNamePrefix", str(output_prefix)])

    if extra_args:
        cmd.extend(extra_args)

    return run_command(cmd, capture_output=False)


def run_featurecounts(
    annotation_file,
    input_bams,
    output_counts,
    threads=4,
    paired_end=False,
    feature_type=None,
    attribute_type=None,
    extra_args=None,
):
    """
    使用 featureCounts 进行定量
    """
    check_file_exists(annotation_file)

    if isinstance(input_bams, (str, Path)):
        bam_list = [input_bams]
    else:
        bam_list = list(input_bams)

    for bam in bam_list:
        check_file_exists(bam)

    ensure_dir(Path(output_counts).parent)

    cmd = [
        "featureCounts",
        "-T",
        str(threads),
        "-a",
        str(annotation_file),
        "-o",
        str(output_counts),
    ]

    if paired_end:
        cmd.append("-p")

    if feature_type:
        cmd.extend(["-t", str(feature_type)])

    if attribute_type:
        cmd.extend(["-g", str(attribute_type)])

    if extra_args:
        cmd.extend(extra_args)

    cmd.extend(str(bam) for bam in bam_list)
    return run_command(cmd, capture_output=False)
