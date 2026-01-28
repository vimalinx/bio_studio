---
name: rfdiffusion
description: Protein structure generation and design using RFdiffusion. Supports monomer, binder, and motif scaffolding.
allowed-tools: Bash
---

RFdiffusion is a method for structure generation using diffusion models.
This skill uses the `rfdiffusion` Docker container.

## Prerequisites
1. **Docker Image**: Must be built locally.
   `docker build -t rfdiffusion repositories/active/RFdiffusion`
2. **Models**: Must be downloaded to a local directory (e.g., `repositories/active/RFdiffusion/models`).
   `bash repositories/active/RFdiffusion/scripts/download_models.sh repositories/active/RFdiffusion/models`

## Usage Commands

**Basic Monomer Generation**:
```bash
docker run -it --rm --gpus all \
  -v $(pwd)/inputs:/inputs \
  -v $(pwd)/outputs:/outputs \
  -v /run/media/vimalinx/Data/bio_studio/repositories/active/RFdiffusion/models:/models \
  rfdiffusion \
  inference.output_prefix=/outputs/test \
  inference.model_directory_path=/models \
  'contigmap.contigs=[150-150]' \
  inference.num_designs=1
```

**Motif Scaffolding**:
```bash
docker run -it --rm --gpus all \
  -v $(pwd)/inputs:/inputs \
  -v $(pwd)/outputs:/outputs \
  -v /run/media/vimalinx/Data/bio_studio/repositories/active/RFdiffusion/models:/models \
  rfdiffusion \
  inference.output_prefix=/outputs/motif_test \
  inference.model_directory_path=/models \
  inference.input_pdb=/inputs/target.pdb \
  'contigmap.contigs=[10-40/A163-181/10-40]'
```

**Parameters**:
- `inference.num_designs`: Number of designs to generate.
- `contigmap.contigs`: Definition of the protein structure (lengths, motifs).
