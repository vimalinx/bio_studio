#!/bin/bash
set -e

ENV_NAME="bio_evo2"
MAMBA_EXE="/home/vimalinx/miniforge3/bin/mamba"

echo "🧬 Setting up Conda environment: $ENV_NAME"
echo "=========================================="

# 1. Create Environment
if $MAMBA_EXE env list | grep -q "$ENV_NAME"; then
    echo "⚠️  Environment '$ENV_NAME' already exists."
    read -p "Recreate it? (This will delete existing env) [y/N] " -n 1 -r
    echo ""
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        $MAMBA_EXE env remove -n "$ENV_NAME" -y
    else
        echo "Skipping creation."
    fi
fi

if ! $MAMBA_EXE env list | grep -q "$ENV_NAME"; then
    echo "🚀 Creating environment..."
    $MAMBA_EXE create -n "$ENV_NAME" python=3.10 -y
fi

# 2. Install PyTorch & Core Libs
# 使用 PyTorch 官方源或清华源
echo "⬇️  Installing PyTorch (CUDA 12.4) and dependencies..."
# 注意：我们显式安装 pytorch-cuda=12.4 以匹配编译器，这对编译 flash-attn 很重要
$MAMBA_EXE install -n "$ENV_NAME" \
    pytorch torchvision torchaudio pytorch-cuda=12.4 \
    -c pytorch -c nvidia -y

# 3. Install Pip Packages
echo "⬇️  Installing HuggingFace & Evo2 ecosystem..."
# 使用清华源加速 pip
PIP_CMD="/home/vimalinx/miniforge3/envs/$ENV_NAME/bin/pip"
$PIP_CMD install \
    transformers \
    accelerate \
    bitsandbytes \
    scipy \
    biopython \
    einops \
    evo2 \
    -i https://pypi.tuna.tsinghua.edu.cn/simple

echo ""
echo "✨ Environment '$ENV_NAME' is ready!"
echo "➡️  Activate with: conda activate $ENV_NAME"
