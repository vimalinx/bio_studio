#!/bin/bash
# 01_download_real_data.sh - Download real yeast RNA-seq data (SRR17226388)
# Source: ENA (European Nucleotide Archive) is usually faster/easier than NCBI SRA for direct FASTQ download.

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_DIR="$PROJECT_ROOT/data/raw"

mkdir -p "$RAW_DIR"

# SRR17226388 - Yeast RNA-seq
# URL pattern for ENA: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/088/SRR17226388/SRR17226388_1.fastq.gz
# We will download a subset (first 100k reads) if possible, or full file if curl supports range, 
# or just full file (it's small enough ~500MB).
# Actually, let's look for a smaller run or just take the full one.
# SRR17226388 is ~4.5M spots.

# ENA FTP URLs
R1_URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/088/SRR17226388/SRR17226388_1.fastq.gz"
# R2_URL="..." (This is Single End or Paired? ENA usually has _1.fastq.gz and _2.fastq.gz for paired)
# Let's check metadata first or assume Single End for now based on typical yeast RNA-seq demos. 
# Actually most modern are Paired.

echo "🌐 Downloading real RNA-seq data: SRR17226388"
echo "   Target: $RAW_DIR"
echo "   Source: ENA"

if ! command -v curl &> /dev/null; then
    echo "❌ curl not found"
    exit 1
fi

# Download Function with Retry
download_file() {
    url=$1
    out_file=$2
    
    echo "   ⬇️ Downloading $out_file ..."
    # -L: Follow redirects
    # -C -: Resume transfer
    # --retry 3: Retry 3 times
    curl -L -C - --retry 3 -o "$out_file" "$url"
    
    if [ $? -eq 0 ]; then
        echo "   ✅ Download complete: $out_file"
    else
        echo "   ❌ Download failed"
        rm -f "$out_file" # Clean up partial
        return 1
    fi
}

# Download R1
download_file "$R1_URL" "$RAW_DIR/SRR17226388_1.fastq.gz"

# Try downloading R2 (if paired)
# If it fails (404), it might be single end.
R2_URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR172/088/SRR17226388/SRR17226388_2.fastq.gz"
echo "   ⬇️ Checking for R2 (Paired-end)..."
# Check header first
if curl --head --silent --fail "$R2_URL" > /dev/null; then
    download_file "$R2_URL" "$RAW_DIR/SRR17226388_2.fastq.gz"
else
    echo "   ℹ️ R2 not found, assuming Single-End or interleaved."
fi

# Rename to match our pipeline convention (SampleName_R1.fastq.gz)
# Our config expects *_R1.fastq.gz
# Let's rename SRR17226388_1.fastq.gz -> Real_Yeast_R1.fastq.gz
if [ -f "$RAW_DIR/SRR17226388_1.fastq.gz" ]; then
    mv "$RAW_DIR/SRR17226388_1.fastq.gz" "$RAW_DIR/Real_Yeast_R1.fastq.gz"
    echo "   ✅ Renamed to: Real_Yeast_R1.fastq.gz"
fi

if [ -f "$RAW_DIR/SRR17226388_2.fastq.gz" ]; then
    mv "$RAW_DIR/SRR17226388_2.fastq.gz" "$RAW_DIR/Real_Yeast_R2.fastq.gz"
    echo "   ✅ Renamed to: Real_Yeast_R2.fastq.gz"
fi

echo ""
echo "✅ Real data download process finished."
