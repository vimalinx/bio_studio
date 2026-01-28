import os
import sys
import subprocess
import random
from pathlib import Path

# 配置路径
PROJECT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
REF_DIR = DATA_DIR / "references"
RES_DIR = DATA_DIR / "results"
QC_DIR = DATA_DIR / "processed" / "qc"
ALIGN_DIR = DATA_DIR / "processed" / "aligned"

# 确保目录存在
for d in [RAW_DIR, REF_DIR, RES_DIR, QC_DIR, ALIGN_DIR]:
    d.mkdir(parents=True, exist_ok=True)

def run_cmd(cmd):
    print(f"🚀 Running: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
        print("  ✅ Success")
    except subprocess.CalledProcessError as e:
        print(f"  ❌ Failed with error code {e.returncode}")
        sys.exit(1)

def generate_mock_data():
    print("\n[1/6] Generating mock data...")
    
    # 1. 生成微型参考基因组 (1条染色体, 2KB)
    ref_fa = REF_DIR / "ref.fa"
    chrom_seq = "".join(random.choices("ACGT", k=2000))
    with open(ref_fa, "w") as f:
        f.write(f">chr1\n{chrom_seq}\n")
    print(f"  Created {ref_fa}")

    # 2. 生成微型GTF (1个基因)
    gtf_file = REF_DIR / "genes.gtf"
    with open(gtf_file, "w") as f:
        # Gene from 100-500
        f.write('chr1\tMOCK\texon\t100\t500\t.\t+\t.\tgene_id "geneA"; transcript_id "txA";\n')
    print(f"  Created {gtf_file}")

    # 3. 生成微型FASTQ (Reads 来源于该基因)
    # 取基因片段 100-200
    gene_seq = chrom_seq[99:499] 
    read_seq = gene_seq[0:50] # 50bp read
    qual = "I" * 50 # High quality
    
    r1 = RAW_DIR / "sample_R1.fastq"
    r2 = RAW_DIR / "sample_R2.fastq"
    
    with open(r1, "w") as f1, open(r2, "w") as f2:
        for i in range(100): # 100 reads
            f1.write(f"@read{i}/1\n{read_seq}\n+\n{qual}\n")
            # Reverse complement for R2 (simplified: just same for validation test)
            f2.write(f"@read{i}/2\n{read_seq}\n+\n{qual}\n")
            
    print(f"  Created mock FASTQ files")

def test_bio_tools():
    print("\n[2/6] Testing Bioinformatics Tools...")
    
    # 1. Fastp
    run_cmd(f"fastp -i {RAW_DIR}/sample_R1.fastq -I {RAW_DIR}/sample_R2.fastq -o {QC_DIR}/clean_R1.fq -O {QC_DIR}/clean_R2.fq -h {QC_DIR}/fastp.html -j {QC_DIR}/fastp.json")
    
    # 2. STAR Index & Align
    # STAR index needs a directory
    star_idx_dir = REF_DIR / "star_index"
    star_idx_dir.mkdir(exist_ok=True)
    
    # 建索引
    run_cmd(f"STAR --runMode genomeGenerate --genomeDir {star_idx_dir} --genomeFastaFiles {REF_DIR}/ref.fa --genomeSAindexNbases 4")
    
    # 比对
    run_cmd(f"STAR --genomeDir {star_idx_dir} --readFilesIn {QC_DIR}/clean_R1.fq {QC_DIR}/clean_R2.fq --outFileNamePrefix {ALIGN_DIR}/sample_ --outSAMtype BAM SortedByCoordinate")
    
    # 3. FeatureCounts
    run_cmd(f"featureCounts -a {REF_DIR}/genes.gtf -o {RES_DIR}/counts.txt {ALIGN_DIR}/sample_Aligned.sortedByCoord.out.bam")
    
    # 4. MultiQC
    run_cmd(f"multiqc {DATA_DIR} -o {RES_DIR}/multiqc_report")

def test_python_libs():
    print("\n[3/6] Testing Python Libraries...")
    
    import pandas as pd
    import numpy as np
    import scipy
    import torch
    import Bio
    
    print(f"  ✅ Pandas {pd.__version__}")
    print(f"  ✅ Numpy {np.__version__}")
    print(f"  ✅ Scipy {scipy.__version__}")
    print(f"  ✅ Torch {torch.__version__} (CUDA: {torch.cuda.is_available()})")
    print(f"  ✅ Biopython {Bio.__version__}")
    
    # Test Scanpy import (usually heavy)
    print("  ⏳ Importing Scanpy...")
    import scanpy as sc
    print(f"  ✅ Scanpy {sc.__version__}")
    
    # Simple Pandas test
    df = pd.read_csv(f"{RES_DIR}/counts.txt", sep="\t", comment="#")
    print(f"  📊 FeatureCounts Output Shape: {df.shape}")

if __name__ == "__main__":
    print("🧬 Starting Bio Studio Environment Validation")
    print(f"📂 Project: {PROJECT_DIR}")
    
    # Add conda bin to PATH (crucial for non-interactive shells)
    os.environ["PATH"] = f"/home/vimalinx/miniforge3/envs/bio/bin:{os.environ['PATH']}"
    
    generate_mock_data()
    test_bio_tools()
    test_python_libs()
    
    print("\n✨ ALL TESTS PASSED! Your environment is rock solid.")
