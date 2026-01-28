#!/bin/bash
set -e
cd "$(dirname "$0")"

echo "🚀 Building modern bio-AI environment (PyTorch Nightly + CUDA 12.6)..."
echo "   This will install PyTorch 2.6.0.dev (preview) compatible with RTX 5070."

docker build -t bio-modern .

echo ""
echo "✨ Environment built successfully!"
echo "🏃 Running container to verify GPU access..."
docker run --rm --gpus all bio-modern
