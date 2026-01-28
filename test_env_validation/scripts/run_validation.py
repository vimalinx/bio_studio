import subprocess
import os
import sys
from pathlib import Path

# Color codes
GREEN = "\033[92m"
RED = "\033[91m"
RESET = "\033[0m"

def run_check(name, cmd, expected_file=None):
    print(f"Testing {name}...", end=" ")
    try:
        # Run command
        result = subprocess.run(
            cmd, 
            shell=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            executable="/bin/bash"
        )
        
        # Check return code
        if result.returncode != 0:
            print(f"{RED}FAILED{RESET}")
            print(f"  Command: {cmd}")
            print(f"  Error: {result.stderr.decode()}")
            return False
            
        # Check output file if specified
        if expected_file and not os.path.exists(expected_file):
            print(f"{RED}FAILED (No Output){RESET}")
            print(f"  Missing: {expected_file}")
            return False
            
        print(f"{GREEN}PASSED{RESET}")
        return True
    except Exception as e:
        print(f"{RED}ERROR{RESET}")
        print(f"  Exception: {e}")
        return False

def main():
    root = Path("test_env_validation")
    data = root / "data"
    refs = data / "references"
    raw = data / "raw"
    res = root / "results"
    
    print("=== Bio Studio Toolchain Validation ===")
    print(f"Workdir: {os.getcwd()}")
    
    # 1. QC Tools
    run_check("FastQC", 
              f"fastqc {raw}/test_reads_R1.fq -o {res}", 
              f"{res}/test_reads_R1_fastqc.html")
              
    run_check("SeqKit", 
              f"seqkit stats {raw}/test_reads_R1.fq")

    # 2. Alignment Tools (DNA)
    # Bowtie2
    run_check("Bowtie2 Index", 
              f"bowtie2-build {refs}/tiny_genome.fa {refs}/tiny_genome")
              
    run_check("Bowtie2 Align", 
              f"bowtie2 -x {refs}/tiny_genome -U {raw}/test_reads_R1.fq -S {res}/bowtie2.sam",
              f"{res}/bowtie2.sam")

    # BWA
    run_check("BWA Index", 
              f"bwa index {refs}/tiny_genome.fa")
              
    run_check("BWA Align", 
              f"bwa mem {refs}/tiny_genome.fa {raw}/test_reads_R1.fq > {res}/bwa.sam",
              f"{res}/bwa.sam")

    # 3. Alignment Tools (RNA)
    # STAR
    os.makedirs(f"{refs}/star_index", exist_ok=True)
    run_check("STAR Index", 
              f"STAR --runMode genomeGenerate --genomeDir {refs}/star_index --genomeFastaFiles {refs}/tiny_genome.fa --genomeSAindexNbases 2")
              
    run_check("STAR Align", 
              f"STAR --genomeDir {refs}/star_index --readFilesIn {raw}/test_reads_R1.fq --outFileNamePrefix {res}/star_",
              f"{res}/star_Aligned.out.sam")

    # 4. SAM/BAM Processing
    run_check("Samtools View", 
              f"samtools view -bS {res}/bowtie2.sam > {res}/bowtie2.bam",
              f"{res}/bowtie2.bam")
              
    run_check("Samtools Sort", 
              f"samtools sort {res}/bowtie2.bam -o {res}/bowtie2.sorted.bam",
              f"{res}/bowtie2.sorted.bam")
              
    run_check("Samtools Index", 
              f"samtools index {res}/bowtie2.sorted.bam",
              f"{res}/bowtie2.sorted.bam.bai")

    # 5. Quantification
    # Specify -F GFF because input is minimal GFF, not GTF
    run_check("featureCounts", 
              f"featureCounts -t gene -g ID -F GFF -a {refs}/tiny_genome.gff -o {res}/counts.txt {res}/bowtie2.sorted.bam",
              f"{res}/counts.txt")

    # 6. Variant Calling
    run_check("Bcftools Call", 
              f"bcftools mpileup -f {refs}/tiny_genome.fa {res}/bowtie2.sorted.bam | bcftools call -mv -o {res}/variants.vcf",
              f"{res}/variants.vcf")

    # 7. Bedtools
    run_check("Bedtools GenomeCov", 
              f"bedtools genomecov -ibam {res}/bowtie2.sorted.bam",
              None)

    print("\n=== Validation Complete ===")

if __name__ == "__main__":
    main()
