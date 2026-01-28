import os
import subprocess
import time
from pathlib import Path
import config

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

def step_1_build_index():
    """Build STAR index"""
    print("\n>>> STEP 1: Building STAR Index")
    
    if config.STAR_INDEX_DIR.exists() and any(config.STAR_INDEX_DIR.iterdir()):
        print("   Index already exists, skipping.")
        return

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

def step_2_alignment():
    """Align reads using STAR"""
    print("\n>>> STEP 2: STAR Alignment")
    
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

def step_3_quantification():
    """Count reads using featureCounts"""
    print("\n>>> STEP 3: Quantification (featureCounts)")
    
    bam_files = [str(config.PROCESSED_DIR / f"{sample}_Aligned.sortedByCoord.out.bam") for sample in config.SAMPLES]
    bam_str = " ".join(bam_files)
    out_counts = config.RESULTS_DIR / "counts.txt"
    
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
    simple_counts = config.RESULTS_DIR / "counts_matrix.csv"
    # Extract columns: GeneID (1) + Counts (7,8,9,10...)
    # Skip header (2 lines)
    # Note: featureCounts header line starts with Geneid, Chr, Start...
    
    # Python script to clean up the matrix cleanly
    cleanup_script = f"""
import pandas as pd
try:
    df = pd.read_csv('{out_counts}', sep='\\t', comment='#')
    # Keep Geneid and sample columns (which correspond to bam filenames)
    cols = ['Geneid'] + [c for c in df.columns if 'bam' in c]
    df_simple = df[cols]
    # Rename columns to sample names
    df_simple.columns = ['GeneID'] + {config.SAMPLES}
    df_simple.to_csv('{simple_counts}', index=False)
    print("Matrix simplified.")
except Exception as e:
    print(f"Error simplifying matrix: {{e}}")
    """
    
    # Executing the cleanup python snippet
    subprocess.run(["python", "-c", cleanup_script])
    
    print(f"   ✅ Quantification complete: {simple_counts}")

def main():
    print("🚀 Starting RNA-seq Pipeline (Yeast Demo)")
    print(f"   Samples: {config.SAMPLES}")
    
    step_1_build_index()
    step_2_alignment()
    step_3_quantification()
    
    print("\n🎉 Pipeline completed successfully!")
    print(f"   Results are in: {config.RESULTS_DIR}")

if __name__ == "__main__":
    main()
