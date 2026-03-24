import os
import sys
import subprocess
import random
import importlib
import shutil
from pathlib import Path

# 配置路径
PROJECT_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = PROJECT_DIR.parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from lib.workspace_env import ensure_workspace_path

DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
REF_DIR = DATA_DIR / "references"
RES_DIR = DATA_DIR / "results"
QC_DIR = DATA_DIR / "processed" / "qc"
ALIGN_DIR = DATA_DIR / "processed" / "aligned"

# 确保目录存在
for d in [RAW_DIR, REF_DIR, RES_DIR, QC_DIR, ALIGN_DIR]:
    d.mkdir(parents=True, exist_ok=True)


def reverse_complement(seq):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def build_mock_read_pair(fragment_seq, read_length=50):
    if len(fragment_seq) < read_length * 2:
        raise ValueError("fragment_seq must be at least twice the read length")
    read1 = fragment_seq[:read_length]
    read2 = reverse_complement(fragment_seq[-read_length:])
    return read1, read2

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
    # 构造一个最小 paired-end fragment，R2 使用末端反向互补，避免双端比对全空
    gene_seq = chrom_seq[99:499]
    fragment_seq = gene_seq[0:150]
    read1_seq, read2_seq = build_mock_read_pair(fragment_seq, read_length=50)
    qual = "I" * 50 # High quality
    
    r1 = RAW_DIR / "sample_R1.fastq"
    r2 = RAW_DIR / "sample_R2.fastq"
    
    with open(r1, "w") as f1, open(r2, "w") as f2:
        for i in range(100): # 100 reads
            f1.write(f"@read{i}/1\n{read1_seq}\n+\n{qual}\n")
            f2.write(f"@read{i}/2\n{read2_seq}\n+\n{qual}\n")
            
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
    run_cmd(f"featureCounts -p -a {REF_DIR}/genes.gtf -o {RES_DIR}/counts.txt {ALIGN_DIR}/sample_Aligned.sortedByCoord.out.bam")
    
    # 4. MultiQC
    multiqc_dir = RES_DIR / "multiqc_report"
    if multiqc_dir.exists():
        shutil.rmtree(multiqc_dir)
    run_cmd(f"multiqc {DATA_DIR} --ignore '*/multiqc_report/*' -o {multiqc_dir}")

def test_python_libs():
    print("\n[3/6] Testing Python Libraries...")

    try:
        pd = importlib.import_module("pandas")
        np = importlib.import_module("numpy")
        scipy = importlib.import_module("scipy")
        torch = importlib.import_module("torch")
        Bio = importlib.import_module("Bio")
    except ModuleNotFoundError as exc:
        print(f"  ❌ Missing required Python library: {exc.name}")
        print(f"  🐍 Current interpreter: {sys.executable}")
        print("  提示: 请先激活 conda `bio` 环境，或使用该环境内的 python 运行 workspace-validate")
        sys.exit(1)
    
    print(f"  ✅ Pandas {pd.__version__}")
    print(f"  ✅ Numpy {np.__version__}")
    print(f"  ✅ Scipy {scipy.__version__}")
    print(f"  ✅ Torch {torch.__version__} (CUDA: {torch.cuda.is_available()})")
    print(f"  ✅ Biopython {Bio.__version__}")
    
    print("  ⏳ Importing Scanpy...")
    try:
        sc = importlib.import_module("scanpy")
    except ModuleNotFoundError:
        print("  ⚠️ Scanpy not installed; optional for current smoke test, skipping.")
    else:
        print(f"  ✅ Scanpy {sc.__version__}")
    
    # Simple Pandas test
    df = pd.read_csv(f"{RES_DIR}/counts.txt", sep="\t", comment="#")
    print(f"  📊 FeatureCounts Output Shape: {df.shape}")

if __name__ == "__main__":
    print("🧬 Starting Bio Studio Environment Validation")
    print(f"📂 Project: {PROJECT_DIR}")

    bio_bin = ensure_workspace_path()
    if bio_bin is not None:
        print(f"🔧 Bio tool PATH ready: {bio_bin}")
    else:
        print("⚠️ Bio env bin not detected automatically; relying on current PATH")

    generate_mock_data()
    test_bio_tools()
    test_python_libs()
    
    print("\n✨ ALL TESTS PASSED! Your environment is rock solid.")
