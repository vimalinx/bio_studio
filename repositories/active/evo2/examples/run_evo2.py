#!/usr/bin/env python3
"""
Evo2 Examples for DNA Sequence Analysis
Run this inside the Docker container
"""

import torch
from evo2 import Evo2

# Available models:
# - evo2_1b_base: 1B parameters, 8K context (Best for 12GB VRAM)
# - evo2_7b_base: 7B parameters, 8K context (May OOM on 12GB)
# - evo2_7b: 7B parameters, 1M context (Requires multiple GPUs)
# - evo2_40b: 40B parameters (Requires multiple GPUs)

MODEL_NAME = "evo2_1b_base"  # Change based on your VRAM


def example_forward():
    """Score likelihoods across a DNA sequence."""
    print("\n=== Example 1: Forward Pass (Scoring) ===")
    evo2_model = Evo2(MODEL_NAME)

    sequence = 'ACGTACGTACGTACGT'
    input_ids = torch.tensor(
        evo2_model.tokenizer.tokenize(sequence),
        dtype=torch.int,
    ).unsqueeze(0).to('cuda:0')

    outputs, _ = evo2_model(input_ids)
    logits = outputs[0]

    print(f'Sequence: {sequence}')
    print(f'Logits shape: {logits.shape}')
    print(f'Logits: {logits}')


def example_generation():
    """Generate DNA sequences based on prompts."""
    print("\n=== Example 2: DNA Generation ===")
    evo2_model = Evo2(MODEL_NAME)

    # Generate DNA sequence
    output = evo2_model.generate(
        prompt_seqs=["ACGTACGT"],
        n_tokens=100,
        temperature=1.0,
        top_k=4
    )

    print(f'Generated sequence: {output.sequences[0]}')


def example_embeddings():
    """Extract embeddings for downstream tasks."""
    print("\n=== Example 3: Extract Embeddings ===")
    evo2_model = Evo2(MODEL_NAME)

    sequence = 'ACGTACGTACGT'
    input_ids = torch.tensor(
        evo2_model.tokenizer.tokenize(sequence),
        dtype=torch.int,
    ).unsqueeze(0).to('cuda:0')

    layer_name = 'blocks.28.mlp.l3'

    outputs, embeddings = evo2_model(
        input_ids,
        return_embeddings=True,
        layer_names=[layer_name]
    )

    print(f'Embeddings shape: {embeddings[layer_name].shape}')
    print(f'Embeddings (first 10 values): {embeddings[layer_name][0, 0, :10]}')


if __name__ == "__main__":
    print(f"Using model: {MODEL_NAME}")

    example_forward()
    example_generation()
    example_embeddings()

    print("\n=== All examples completed! ===")
