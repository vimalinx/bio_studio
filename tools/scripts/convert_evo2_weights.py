#!/usr/bin/env python3
"""
把单个 evo2_7b.pt 拆成 HuggingFace 目录
"""
import torch, json, os
from collections import OrderedDict

PT_FILE   = "evo2_7b.pt"
HF_FOLDER = "evo2-7b-hf"
os.makedirs(HF_FOLDER, exist_ok=True)

# ---- 加载权重 ----
print(">>> 加载原始权重 ...")
ckpt = torch.load(PT_FILE, map_location="cpu", weights_only=False)
if "model" in ckpt:                      # 有的存档套一层
    ckpt = ckpt["model"]

# ---- 写 config ----
config = {
  "architectures": ["Evo2ForCausalLM"],
  "model_type": "evo2",
  "vocab_size": 512,
  "hidden_size": 4096,
  "num_hidden_layers": 32,
  "num_attention_heads": 32,
  "intermediate_size": 11008,
  "rms_norm_eps": 1e-05,
  "max_position_embeddings": 131072,
  "rope_theta": 10000.0,
  "tie_word_embeddings": False,
  "torch_dtype": "float16"
}
with open(f"{HF_FOLDER}/config.json", "w") as f:
    json.dump(config, f, indent=2)

# ---- 权重改名 ----
new_state = OrderedDict()
for k, v in ckpt.items():
    if not isinstance(v, torch.Tensor):          # 跳过非张量
        continue
    k = k.replace("layers.", "model.layers.")
    k = k.replace("attn.", "self_attn.")
    new_state[k] = v

# ---- 分片保存 ----
max_shard = 5e9
current_size = 0
shard_idx  = 0
shard_dict = OrderedDict()
index = {"metadata": {"total_size": 0}, "weight_map": {}}

for name, tensor in new_state.items():
    tensor_size = tensor.numel() * tensor.element_size()
    if current_size + tensor_size > max_shard and shard_dict:
        shard_name = f"pytorch_model-{shard_idx+1:05d}-of-{(shard_idx+1):05d}.bin"
        torch.save(shard_dict, f"{HF_FOLDER}/{shard_name}")
        print(f"  保存分片 {shard_name}")
        for key in shard_dict.keys():
            index["weight_map"][key] = shard_name
        shard_idx += 1
        current_size = 0
        shard_dict = OrderedDict
