# 使用 NVIDIA 最新 CUDA 13.0 镜像
FROM nvidia/cuda:13.0.0-cudnn-devel-ubuntu22.04

# 引入 uv (极速包管理器)
COPY --from=ghcr.io/astral-sh/uv:latest /uv /bin/uv

# 设置环境变量
ENV DEBIAN_FRONTEND=noninteractive
# 告诉 uv 安装到系统环境，而不是创建 venv
ENV UV_SYSTEM_PYTHON=1

# 1. 安装系统依赖
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    git \
    wget \
    curl \
    python-is-python3 \
    && rm -rf /var/lib/apt/lists/*

# 2. 使用 uv 安装 PyTorch (尝试 cu130 尝鲜版)
# uv 的下载稳定性通常优于 pip
RUN uv pip install --pre torch torchvision torchaudio --index-url https://download.pytorch.org/whl/nightly/cu130

# 3. 安装构建依赖
RUN uv pip install packaging ninja setuptools wheel psutil

# 4. 安装 Evo2 及相关依赖
RUN uv pip install \
    transformers \
    scipy \
    biopython \
    einops \
    flash-attn --no-build-isolation \
    stripedhyena

# 5. 设置工作目录
WORKDIR /workspace

# 6. 验证环境
CMD ["bash", "-c", "nvcc --version && python3 -c \"import torch; print(f'🔥 PyTorch: {torch.__version__}'); print(f'🚀 CUDA (PyTorch internal): {torch.version.cuda}'); print(f'✅ GPU: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'None'}\""]

