---
name: evo2
description: Genome modeling and design using Evo 2 DNA language model. Supports generation, scoring, and embeddings.
allowed-tools: Bash
---

Evo 2 is a DNA language model for long-context modeling and design.
This skill uses the `evo2` Docker container which is already available in the environment.

## Usage

**Generation**:
To generate DNA sequences:
```bash
docker run --rm --gpus all -v $(pwd):/workspace -w /workspace evo2 python -m evo2.test.test_evo2_generation --model_name evo2_7b
```

**Scoring (Forward Pass)**:
To score sequences, create a python script (e.g., `score.py`) and run it inside the container:
```bash
docker run --rm --gpus all -v $(pwd):/workspace -w /workspace evo2 python score.py
```

**Persisting Models**:
The models are large. To avoid re-downloading, mount a local cache directory:
```bash
docker run --rm --gpus all \
  -v $(pwd):/workspace \
  -v $HOME/.cache/huggingface:/root/.cache/huggingface \
  -w /workspace \
  evo2 python -m evo2.test.test_evo2_generation --model_name evo2_7b
```

## Checkpoints
- `evo2_7b`: 7B params, 1M context
- `evo2_40b`: 40B params (requires multi-GPU)
- `evo2_7b_base`: 7B params, 8K context
