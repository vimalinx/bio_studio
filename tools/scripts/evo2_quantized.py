#!/usr/bin/env python3
"""
Evo 2 Quantized Inference (12GB VRAM Optimized)
-----------------------------------------------
This script enables running the Evo 2 7B model on consumer GPUs (e.g., RTX 3080/4070)
by utilizing 8-bit quantization via bitsandbytes.

Usage:
    python tools/scripts/evo2_quantized.py --model_dir /path/to/evo2-7b-hf --prompt "ACGT"

Requirements:
    pip install transformers bitsandbytes accelerate
"""

import torch
import argparse
import os
from transformers import AutoTokenizer, AutoModelForCausalLM, BitsAndBytesConfig, AutoConfig

def setup_model(model_path):
    print(f"🔄 Loading model from: {model_path}")
    print("   Configuring 8-bit quantization...")
    
    # 1. Quantization Config
    bnb_config = BitsAndBytesConfig(
        load_in_8bit=True,
        bnb_8bit_compute_dtype=torch.float16,
        bnb_8bit_use_double_quant=True,
        llm_int8_skip_modules=["lm_head"] # Important for some architectures
    )

    # 2. Register 'evo2' config if needed (AutoConfig fallback)
    # Note: If the HF folder contains modeling_evo2.py, trust_remote_code=True handles it.
    
    # 3. Load Model
    try:
        model = AutoModelForCausalLM.from_pretrained(
            model_path,
            quantization_config=bnb_config,
            device_map="auto",
            torch_dtype=torch.float16,
            trust_remote_code=True
        )
    except Exception as e:
        print(f"❌ Error loading model: {e}")
        print("   Ensure the directory contains a valid HF model (config.json, *.bin/*.safetensors)")
        return None, None

    # 4. Load Tokenizer
    try:
        tokenizer = AutoTokenizer.from_pretrained(model_path, trust_remote_code=True)
    except:
        print("⚠️  Tokenizer loading failed, trying default Evo2 approach...")
        # Fallback logic could go here
        return None, None

    return model, tokenizer

def generate_sequence(model, tokenizer, prompt, max_new_tokens=100):
    print(f"🧬 generating sequence for prompt: {prompt}")
    
    inputs = torch.tensor(
        tokenizer.tokenize(prompt),
        dtype=torch.int32
    ).unsqueeze(0).to("cuda")

    with torch.no_grad():
        outputs = model(inputs)
        # Simple generation logic (forward pass check)
        logits = outputs.logits
        
    # Clean up
    torch.cuda.empty_cache()
    
    return logits

def main():
    parser = argparse.ArgumentParser(description="Run Evo 2 7B with 8-bit quantization")
    parser.add_argument("--model_dir", type=str, required=True, help="Path to HF-converted model directory")
    parser.add_argument("--prompt", type=str, default="ACGT", help="DNA sequence prompt")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.model_dir):
        print(f"❌ Model directory not found: {args.model_dir}")
        print("   Please convert .pt weights using 'convert_evo2_weights.py' first.")
        exit(1)

    model, tokenizer = setup_model(args.model_dir)
    
    if model:
        logits = generate_sequence(model, tokenizer, args.prompt)
        print(f"✅ Inference successful. Logits shape: {logits.shape}")
        print("   (To implement full generation, add sampling loop or use model.generate)")

if __name__ == "__main__":
    main()
