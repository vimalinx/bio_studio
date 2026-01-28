#!/bin/bash
# Bio Studio Environment Setup Script
# Detects Miniforge/Miniconda and installs 'bio' environment

set -e

# 1. Detect Conda
CONDA_EXE=""
if [ -f "$HOME/miniforge3/bin/conda" ]; then
    CONDA_EXE="$HOME/miniforge3/bin/conda"
    CONDA_BASE="$HOME/miniforge3"
elif [ -f "$HOME/miniconda3/bin/conda" ]; then
    CONDA_EXE="$HOME/miniconda3/bin/conda"
    CONDA_BASE="$HOME/miniconda3"
elif command -v conda &> /dev/null; then
    CONDA_EXE=$(command -v conda)
    CONDA_BASE=$(dirname $(dirname "$CONDA_EXE"))
else
    echo "❌ Error: Conda not found. Please install Miniforge or Miniconda."
    exit 1
fi

echo "✅ Found Conda at: $CONDA_EXE"
source "$CONDA_BASE/etc/profile.d/conda.sh"

# 2. Create/Update 'bio' environment
echo "🚀 Creating/Updating 'bio' conda environment..."
echo "   (This may take a few minutes depending on your internet connection)"

# Create environment with Python 3.10 and core bioconda tools
conda create -n bio -c conda-forge -c bioconda -y \
    python=3.10 \
    blast bowtie2 bwa samtools bcftools vcftools bedtools \
    hmmer muscle clustalw iqtree fastqc seqtk \
    cutadapt biopython pandas numpy matplotlib seaborn jupyterlab ipywidgets \
    pytorch

# 3. Activate and install Pip dependencies
echo "📦 Installing PyPI dependencies..."
conda activate bio

# Install remaining pip packages from requirements.txt
# excluding ones we already installed via conda to speed up
grep -vE "biopython|pandas|numpy|matplotlib|seaborn|jupyter|ipywidgets|cutadapt|torch" requirements.txt > requirements_pip.tmp

pip install -r requirements_pip.tmp
rm requirements_pip.tmp

echo ""
echo "🎉 Setup Complete!"
echo "To use the environment, run: conda activate bio"
