#!/usr/bin/env python3

from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
import time
from pathlib import Path

import config


VALIDATOR_PATH = Path(__file__).resolve().with_name("validate_project.py")


def load_template_runtime():
    script_dir = Path(__file__).resolve().parent
    for candidate_root in script_dir.parents:
        runtime_path = candidate_root / "lib" / "template_runtime.py"
        if not runtime_path.exists():
            continue
        spec = importlib.util.spec_from_file_location("bio_template_runtime", runtime_path)
        if spec is None or spec.loader is None:
            continue
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    return None


TEMPLATE_RUNTIME = load_template_runtime()


def run_command(cmd, log_file=None):
    """Run shell command with logging"""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] Running: {cmd}")
    
    if log_file:
        with open(log_file, "w") as f:
            process = subprocess.Popen(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT, executable="/bin/bash")
    else:
        process = subprocess.Popen(cmd, shell=True, executable="/bin/bash")
        
    process.wait()
    
    if process.returncode != 0:
        print(f"❌ Command failed with return code {process.returncode}")
        exit(1)
    return True


def run_shared_step(step_name: str):
    if TEMPLATE_RUNTIME is None:
        return None
    return TEMPLATE_RUNTIME.run_shared_step(config, step_name)


def simplify_counts_matrix(out_counts: Path, simple_counts: Path) -> None:
    cleanup_script = f"""
import pandas as pd
try:
    df = pd.read_csv('{out_counts}', sep='\\t', comment='#')
    cols = ['Geneid'] + [c for c in df.columns if 'bam' in c]
    df_simple = df[cols]
    df_simple.columns = ['GeneID'] + {config.SAMPLES}
    df_simple.to_csv('{simple_counts}', index=False)
    print("Matrix simplified.")
except Exception as e:
    print(f"Error simplifying matrix: {{e}}")
    raise
"""
    subprocess.run([sys.executable, "-c", cleanup_script], check=True)


def step_1_build_index():
    """Build STAR index"""
    print("\n>>> STEP 1: Building STAR Index")

    shared_result = run_shared_step("build_index")
    if shared_result is not None:
        print("   Shared runtime handled STAR index construction.")
        return shared_result

    if config.STAR_INDEX_DIR.exists() and any(config.STAR_INDEX_DIR.iterdir()):
        print("   Index already exists, skipping.")
        return True

    config.STAR_INDEX_DIR.mkdir(parents=True, exist_ok=True)
    
    # Calculate genome length for SAindexNbases
    # Formula: min(14, log2(GenomeLength)/2 - 1)
    # Yeast genome ~12Mb -> log2(12*10^6) ~23.5 -> 23.5/2 - 1 = 10.75 -> 10
    cmd = f"""
    STAR --runThreadN {config.THREADS} \\
         --runMode genomeGenerate \\
         --genomeDir {config.STAR_INDEX_DIR} \\
         --genomeFastaFiles {config.GENOME_FASTA} \\
         --genomeSAindexNbases 10 \\
         --sjdbGTFfile {config.GENOME_GTF} \\
         --sjdbGTFtagExonParentTranscript Parent \\
         --sjdbOverhang 99
    """
    run_command(cmd, log_file=config.LOG_DIR / "star_index.log")
    print("   ✅ STAR Index built successfully.")
    return True

def step_2_alignment():
    """Align reads using STAR"""
    print("\n>>> STEP 2: STAR Alignment")

    shared_result = run_shared_step("alignment")
    if shared_result is not None:
        print("   Shared runtime handled alignment and BAM indexing.")
        return shared_result

    for sample in config.SAMPLES:
        print(f"   Processing sample: {sample}")
        
        # Define I/O
        r1 = config.RAW_DIR / f"{sample}_R1.fastq.gz"
        # Note: Simulated data is SE (Single End) for simplicity
        
        out_prefix = config.PROCESSED_DIR / f"{sample}_"
        bam_file = config.PROCESSED_DIR / f"{sample}_Aligned.sortedByCoord.out.bam"
        
        if bam_file.exists():
            print(f"   BAM already exists for {sample}, skipping.")
            continue
            
        r2 = config.RAW_DIR / f"{sample}_R2.fastq.gz"
        read_files = f"{r1} {r2}" if r2.exists() else f"{r1}"
        
        cmd = f"""
        STAR --runThreadN {config.THREADS} \\
             --genomeDir {config.STAR_INDEX_DIR} \\
             --readFilesIn {read_files} \\
             --readFilesCommand zcat \\
             --outFileNamePrefix {out_prefix} \\
             --outSAMtype BAM SortedByCoordinate \\
             --outSAMunmapped Within \\
             --outSAMattributes Standard
        """
        run_command(cmd, log_file=config.LOG_DIR / f"star_align_{sample}.log")
        
        # Index BAM
        run_command(f"samtools index {bam_file}")
        print(f"   ✅ Aligned: {bam_file}")
    return True

def step_3_quantification():
    """Count reads using featureCounts"""
    print("\n>>> STEP 3: Quantification (featureCounts)")

    bam_files = [str(config.PROCESSED_DIR / f"{sample}_Aligned.sortedByCoord.out.bam") for sample in config.SAMPLES]
    bam_str = " ".join(bam_files)
    out_counts = config.COUNTS_PATH
    simple_counts = config.RESULTS_DIR / "counts_matrix.csv"

    shared_result = run_shared_step("quantification")
    if shared_result is not None:
        simplify_counts_matrix(out_counts, simple_counts)
        print(f"   ✅ Quantification complete: {simple_counts}")
        return shared_result
    
    # featureCounts parameters:
    # -t gene: Count reads in gene regions (GFF uses 'gene' type)
    # -g ID: Group by 'ID' attribute in GFF
    # -a: Annotation file
    # -T: Threads
    cmd = f"""
    featureCounts -T {config.THREADS} \\
                  -t gene \\
                  -g ID \\
                  -a {config.GENOME_GTF} \\
                  -o {out_counts} \\
                  {bam_str}
    """
    
    run_command(cmd, log_file=config.LOG_DIR / "featureCounts.log")
    
    # Simplify output (GeneID + Counts)
    simplify_counts_matrix(out_counts, simple_counts)

    print(f"   ✅ Quantification complete: {simple_counts}")
    return True


def run_validation() -> int:
    spec = importlib.util.spec_from_file_location("project_validation", VALIDATOR_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"无法加载验证脚本: {VALIDATOR_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.run_validation()


def print_steps() -> None:
    steps = [
        ("build_index", step_1_build_index.__doc__),
        ("alignment", step_2_alignment.__doc__),
        ("quantification", step_3_quantification.__doc__),
    ]
    print("可用步骤:")
    for name, doc in steps:
        print(f"  {name}: {doc}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Yeast RNA-seq demo pipeline")
    parser.add_argument("--steps", action="store_true", help="只列出可用步骤，不执行")
    parser.add_argument("--validate", action="store_true", help="运行项目级自检并退出")
    args = parser.parse_args()

    if args.steps:
        print_steps()
        return

    if args.validate:
        exit_code = run_validation()
        report_path = Path(config.LOGS_DIR) / "validation_report.json"
        print(f"项目验证完成，请查看: {report_path}")
        sys.exit(exit_code)

    print("🚀 Starting RNA-seq Pipeline (Yeast Demo)")
    print(f"   Samples: {config.SAMPLES}")

    step_1_build_index()
    step_2_alignment()
    step_3_quantification()

    print("\n🎉 Pipeline completed successfully!")
    print(f"   Results are in: {config.RESULTS_DIR}")

if __name__ == "__main__":
    main()
