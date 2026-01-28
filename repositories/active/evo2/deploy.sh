#!/bin/bash
# Evo2 Docker Deployment Script
# For RTX 5070 (12GB VRAM) - Recommend starting with 1B model

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Evo2 Docker Deployment ===${NC}"

# Configuration
IMAGE_NAME="evo2"
CONTAINER_NAME="evo2-container"
HF_CACHE_DIR="$(pwd)/huggingface_cache"
MODEL_NAME="${1:-evo2_1b_base}"  # Default to 1B model for 12GB VRAM

# Create cache directory
echo -e "${YELLOW}Creating HuggingFace cache directory...${NC}"
mkdir -p "$HF_CACHE_DIR"

# Build Docker image
echo -e "${YELLOW}Building Docker image...${NC}"
docker build -t "$IMAGE_NAME" .

# Check if container already exists and remove it
if docker ps -a --format '{{.Names}}' | grep -q "^${CONTAINER_NAME}$"; then
    echo -e "${YELLOW}Removing existing container...${NC}"
    docker rm -f "$CONTAINER_NAME"
fi

# Run container
echo -e "${YELLOW}Starting container with GPU support...${NC}"
docker run -it --rm \
    --name "$CONTAINER_NAME" \
    --gpus all \
    --shm-size=16g \
    -v "$HF_CACHE_DIR:/root/.cache/huggingface" \
    -v "$(pwd):/workdir" \
    -e NVIDIA_VISIBLE_DEVICES=all \
    "$IMAGE_NAME" \
    bash -c "echo 'Container ready! Running test with model: $MODEL_NAME' && python -m evo2.test.test_evo2_generation --model_name $MODEL_NAME"
