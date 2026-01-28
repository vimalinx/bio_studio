"""
Test batch DNA generation with memory tracking.
Generates 10000 DNA sequences and monitors GPU memory usage.
"""
import torch
from evo2 import Evo2
import time
import psutil
import os

def get_gpu_memory():
    """Get GPU memory usage in MB"""
    device = torch.cuda.current_device()
    allocated = torch.cuda.memory_allocated(device) / 1024**2
    reserved = torch.cuda.memory_reserved(device) / 1024**2
    total = torch.cuda.get_device_properties(device).total_memory / 1024**2
    return {
        'allocated_mb': allocated,
        'reserved_mb': reserved,
        'total_mb': total,
        'free_mb': total - allocated,
        'usage_percent': (allocated / total) * 100
    }

def test_long_sequence_generation(seq_length=10000):
    """
    Generate a single long DNA sequence and track memory usage.

    Args:
        seq_length: Length of the sequence to generate
    """
    print("=" * 60)
    print("Evo2 1B Base - Long Sequence Generation Test")
    print("=" * 60)
    print(f"Sequence length: {seq_length} tokens")
    print()

    # Load model
    print("Loading evo2_1b_base model...")
    model = Evo2("evo2_1b_base")
    print("Model loaded successfully!")

    # Initial memory state
    initial_mem = get_gpu_memory()
    print(f"\nInitial GPU Memory:")
    print(f"  Allocated: {initial_mem['allocated_mb']:.2f} MB")
    print(f"  Reserved: {initial_mem['reserved_mb']:.2f} MB")
    print(f"  Total: {initial_mem['total_mb']:.2f} MB")
    print(f"  Free: {initial_mem['free_mb']:.2f} MB")
    print(f"  Usage: {initial_mem['usage_percent']:.2f}%")
    print()

    # Generation parameters
    temperature = 1.0
    top_k = 4
    prompt = "ATCG"

    print(f"Starting generation ({seq_length} tokens)...")
    print("-" * 60)

    import time
    start_time = time.time()

    result = model.generate(
        prompt_seqs=[prompt],
        n_tokens=seq_length,
        temperature=temperature,
        top_k=top_k,
        verbose=1  # Show progress
    )

    elapsed_time = time.time() - start_time

    # Check memory after generation
    final_mem = get_gpu_memory()

    print("-" * 60)
    print("\nGeneration Complete!")
    print()

    # Statistics
    print("=" * 60)
    print("Generation Statistics")
    print("=" * 60)
    print(f"Sequence length: {len(result.sequences[0])} bp")
    print(f"Total time: {elapsed_time:.2f} seconds")
    print(f"Tokens per second: {seq_length / elapsed_time:.2f}")
    print()

    print("=" * 60)
    print("GPU Memory Statistics")
    print("=" * 60)
    print(f"Initial Allocated: {initial_mem['allocated_mb']:.2f} MB")
    print(f"Final Allocated: {final_mem['allocated_mb']:.2f} MB")
    print(f"Peak Allocated: {final_mem['allocated_mb']:.2f} MB")
    print(f"Memory Increase: {final_mem['allocated_mb'] - initial_mem['allocated_mb']:.2f} MB")
    print(f"Current Usage: {final_mem['usage_percent']:.2f}%")
    print(f"GPU Total Memory: {final_mem['total_mb']:.2f} MB")
    print()

    # Output the sequence
    print("=" * 60)
    print("Generated Sequence (first 500 bp)")
    print("=" * 60)
    print(result.sequences[0][:500])
    print(f"\n... (total length: {len(result.sequences[0])} bp)")

    return result.sequences[0]

def test_batch_generation(num_sequences=10000, batch_size=100):
    print("=" * 60)
    print("Evo2 1B Base - Batch DNA Generation Test")
    print("=" * 60)
    print(f"Total sequences to generate: {num_sequences}")
    print(f"Batch size: {batch_size}")
    print(f"Number of batches: {num_sequences // batch_size}")
    print()

    # Load model
    print("Loading evo2_1b_base model...")
    model = Evo2("evo2_1b_base")
    print("Model loaded successfully!")

    # Initial memory state
    initial_mem = get_gpu_memory()
    print(f"\nInitial GPU Memory:")
    print(f"  Allocated: {initial_mem['allocated_mb']:.2f} MB")
    print(f"  Reserved: {initial_mem['reserved_mb']:.2f} MB")
    print(f"  Total: {initial_mem['total_mb']:.2f} MB")
    print(f"  Free: {initial_mem['free_mb']:.2f} MB")
    print(f"  Usage: {initial_mem['usage_percent']:.2f}%")
    print()

    # Generation parameters
    n_tokens = 100  # Length of each generated sequence
    temperature = 1.0
    top_k = 4

    # Track statistics
    peak_memory = 0
    total_time = 0
    all_generated_sequences = []

    # Generate in batches
    num_batches = num_sequences // batch_size
    prompt_base = "ATCG"  # Base prompt for all sequences

    print(f"Starting generation ({n_tokens} tokens per sequence)...")
    print("-" * 60)

    for batch_idx in range(num_batches):
        batch_start = time.time()

        # Create prompts for this batch
        prompts = [prompt_base] * batch_size

        # Generate
        result = model.generate(
            prompt_seqs=prompts,
            n_tokens=n_tokens,
            temperature=temperature,
            top_k=top_k,
            batched=True,
            cached_generation=True,
            verbose=0
        )

        # The result is a GenerationOutput object with sequences attribute
        generated = result.sequences
        logprobs = result.logprobs_mean

        all_generated_sequences.extend(generated)

        batch_time = time.time() - batch_start
        total_time += batch_time

        # Check memory after this batch
        current_mem = get_gpu_memory()
        peak_memory = max(peak_memory, current_mem['allocated_mb'])

        # Progress update every 100 batches
        if (batch_idx + 1) % 100 == 0 or batch_idx == 0:
            sequences_done = (batch_idx + 1) * batch_size
            throughput = (batch_size / batch_time)
            print(f"Batch {batch_idx + 1}/{num_batches} | "
                  f"Sequences: {sequences_done}/{num_sequences} | "
                  f"Time: {batch_time:.2f}s | "
                  f"Throughput: {throughput:.1f} seq/s | "
                  f"GPU: {current_mem['usage_percent']:.1f}%")

    print("-" * 60)
    print("\nGeneration Complete!")
    print()

    # Final memory state
    final_mem = get_gpu_memory()

    # Statistics
    print("=" * 60)
    print("Generation Statistics")
    print("=" * 60)
    print(f"Total sequences generated: {len(all_generated_sequences)}")
    print(f"Total time: {total_time:.2f} seconds")
    print(f"Average time per batch: {total_time / num_batches:.2f} seconds")
    print(f"Average throughput: {num_sequences / total_time:.2f} sequences/second")
    print()

    print("=" * 60)
    print("GPU Memory Statistics")
    print("=" * 60)
    print(f"Initial Allocated: {initial_mem['allocated_mb']:.2f} MB")
    print(f"Final Allocated: {final_mem['allocated_mb']:.2f} MB")
    print(f"Peak Allocated: {peak_memory:.2f} MB")
    print(f"Memory Increase: {final_mem['allocated_mb'] - initial_mem['allocated_mb']:.2f} MB")
    print(f"Current Usage: {final_mem['usage_percent']:.2f}%")
    print(f"GPU Total Memory: {final_mem['total_mb']:.2f} MB")
    print()

    # Sample outputs
    print("=" * 60)
    print("Sample Generated Sequences (first 5)")
    print("=" * 60)
    for i, seq in enumerate(all_generated_sequences[:5]):
        print(f"\nSequence {i+1}:")
        print(f"  {seq[:150]}...")
        if len(seq) > 150:
            print(f"  (total length: {len(seq)} bp)")

    return all_generated_sequences

if __name__ == "__main__":
    # Generate 100 sequences in batches
    sequences = test_batch_generation(num_sequences=100, batch_size=10)
