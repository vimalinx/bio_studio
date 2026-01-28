#!/usr/bin/env python3
"""
Evo 2 埃博拉病毒变异分析脚本
使用 evo2_1b_base 模型预测 568 个变异的效应
"""

import sys
import json
import torch
import numpy as np
from Bio import SeqIO
from collections import defaultdict

# 导入 Evo 2
from evo2 import Evo2

print("=" * 80)
print("Evo 2 埃博拉病毒变异效应预测")
print("=" * 80)
print()

# ==========================================
# 1. 初始化 Evo 2 模型
# ==========================================
print("正在初始化 Evo 2 模型...")
print(f"模型: evo2_1b_base")
evo2_model = Evo2("evo2_1b_base")
print("✓ Evo 2 模型加载完成")
print()

# 检查 CUDA 可用性
device = "cuda:0" if torch.cuda.is_available() else "cpu"
print(f"设备: {device}")
print()

# ==========================================
# 2. 读取参考基因组
# ==========================================
print("正在读取参考基因组...")
ref_fasta = (
    sys.argv[1]
    if len(sys.argv) > 1
    else "/workdir/shared_data/references/ebov_zaire.fa"
)

with open(ref_fasta, "r") as f:
    ref_record = SeqIO.read(f, "fasta")
ref_sequence = str(ref_record.seq)

print(f"参考基因组: {ref_record.id}")
print(f"长度: {len(ref_sequence):,} bp")
print()

# ==========================================
# 3. 读取 VCF 文件
# ==========================================
print("正在读取 VCF 文件...")
vcf_file = (
    sys.argv[2]
    if len(sys.argv) > 2
    else "/workdir/projects/test_rnaseq_analysis/data/processed/SRR1972739_variants.vcf"
)

variants = []
with open(vcf_file, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        chrom, pos, id_, ref, alt, qual, filter_, info, format_, sample = fields[:10]
        variants.append(
            {
                "chrom": chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alt,
                "qual": float(qual) if qual != "." else 0,
                "info": info,
                "sample": sample,
            }
        )

print(f"读取到 {len(variants):,} 个变异")
print()

# ==========================================
# 4. 变异效应预测
# ==========================================
print("正在使用 Evo 2 预测变异效应...")
print()

context_size = 100  # 使用 100 bp 上下文
results = []

for i, variant in enumerate(variants):
    if (i + 1) % 50 == 0:
        print(f"进度: {i + 1}/{len(variants)}")

    pos = variant["pos"]
    ref = variant["ref"]
    alt = variant["alt"]

    # 计算上下文范围
    start = max(0, pos - context_size)
    end = min(len(ref_sequence), pos + len(ref) + context_size)

    # 提取参考序列上下文
    ref_context = ref_sequence[start:pos] + ref + ref_sequence[pos + len(ref) : end]

    # 构建变异序列
    if len(alt) == 1:
        alt_context = ref_sequence[start:pos] + alt + ref_sequence[pos + len(ref) : end]
    else:
        # INDEL 处理
        alt_context = ref_sequence[start:pos] + alt + ref_sequence[pos + len(ref) : end]

    # Tokenize 并预测
    ref_tokens = evo2_model.tokenizer.tokenize(ref_context)
    alt_tokens = evo2_model.tokenizer.tokenize(alt_context)

    # 计算参考序列的似然
    ref_input_ids = torch.tensor([ref_tokens], dtype=torch.long).to(device)
    
    # 获取原始输出，不做拆包假设
    model_output = evo2_model(ref_input_ids)
    
    # DEBUG: 打印输出结构 (更深层)
    if i == 0:
        print(f"DEBUG: model_output type: {type(model_output)}")
        if isinstance(model_output, tuple):
             print(f"DEBUG: model_output length: {len(model_output)}")
             for idx, item in enumerate(model_output):
                 print(f"DEBUG: model_output[{idx}] type: {type(item)}")
                 if isinstance(item, tuple):
                     print(f"DEBUG: model_output[{idx}] length: {len(item)}")
                     for sub_idx, sub_item in enumerate(item):
                         print(f"DEBUG: model_output[{idx}][{sub_idx}] type: {type(sub_item)}")
                         if hasattr(sub_item, 'shape'):
                             print(f"DEBUG: model_output[{idx}][{sub_idx}] shape: {sub_item.shape}")

    # 尝试提取 logits
    ref_logits = None
    
    # 递归查找 Tensor 的辅助函数 (深度优先)
    def find_tensor(obj):
        if isinstance(obj, torch.Tensor):
            return obj
        if isinstance(obj, (tuple, list)):
            for item in obj:
                res = find_tensor(item)
                if res is not None: return res
        return None

    # 如果是 HuggingFace 风格对象
    if hasattr(model_output, 'logits'):
        ref_logits = model_output.logits
    # 如果是 tuple，通常 logits 是第一个 tensor
    elif isinstance(model_output, tuple):
        # 针对刚才发现的 tuple in tuple 结构
        # model_output[0] 是 tuple
        if isinstance(model_output[0], tuple):
             # 尝试取 model_output[0][0]
             ref_logits = model_output[0][0]
        elif isinstance(model_output[0], torch.Tensor):
             ref_logits = model_output[0]
        else:
             # 最后的手段：递归搜索第一个 Tensor
             print("DEBUG: Searching for tensor recursively...")
             ref_logits = find_tensor(model_output)

    if ref_logits is None:
        raise ValueError(f"Could not find logits tensor in output of type {type(model_output)}")
        
    # 再次检查 ref_logits 是否还有 batch 维
    if len(ref_logits.shape) == 3:
        ref_logits = ref_logits[0] # remove batch dim -> [seq, vocab]

    # 计算变异序列的似然
    alt_input_ids = torch.tensor([alt_tokens], dtype=torch.long).to(device)
    # 同上处理
    alt_output_raw = evo2_model(alt_input_ids)
    
    # 使用相同的提取逻辑
    if hasattr(alt_output_raw, 'logits'):
        alt_logits = alt_output_raw.logits
    elif isinstance(alt_output_raw, tuple):
        if isinstance(alt_output_raw[0], tuple):
             alt_logits = alt_output_raw[0][0]
        elif isinstance(alt_output_raw[0], torch.Tensor):
             alt_logits = alt_output_raw[0]
        else:
             alt_logits = find_tensor(alt_output_raw)
             
    if len(alt_logits.shape) == 3:
        alt_logits = alt_logits[0]

    # 获取变异位置的 logits
    # 注意：token 长度可能不等于字符长度，需要重新计算位置
    # 这里我们简化处理：假设 token 差不多对应字符，或者取中间位置
    # 更好的做法是找到变异对应的具体 token index
    
    # 对于 Evo2 CharLevelTokenizer，token 基本上是字符级或者是 byte 级
    # 让我们假设输入长度就是输出长度
    
    seq_len = ref_logits.shape[0]
    # 变异位置应该是输入的 ref_context 的中间
    # ref_context = [start...pos] + ref + [pos+len...end]
    # ref 长度为 1 (SNP) 或更多
    # 所以变异起始位置在 context_size
    
    variant_start_idx = context_size
    # 确保索引不越界
    if variant_start_idx >= seq_len:
         variant_start_idx = seq_len - 1
         
    ref_pos_logits = ref_logits[variant_start_idx]
    alt_pos_logits = alt_logits[variant_start_idx]

    # 计算得分
    ref_max_logit = torch.max(ref_pos_logits).item()
    alt_max_logit = torch.max(alt_pos_logits).item()

    # 效应分数 = alt_logit - ref_logit
    effect_score = alt_max_logit - ref_max_logit

    # 计算序列似然（所有位置的平均）
    ref_log_probs = torch.log_softmax(ref_logits, dim=-1)
    alt_log_probs = torch.log_softmax(alt_logits, dim=-1)

    ref_likelihood = torch.sum(ref_log_probs).item()
    alt_likelihood = torch.sum(alt_log_probs).item()

    likelihood_diff = alt_likelihood - ref_likelihood

    results.append(
        {
            "position": pos,
            "ref": ref,
            "alt": alt,
            "qual": variant["qual"],
            "effect_score": effect_score,
            "likelihood_diff": likelihood_diff,
            "ref_max_logit": ref_max_logit,
            "alt_max_logit": alt_max_logit,
        }
    )

print(f"✓ 完成 {len(variants):,} 个变异的效应预测")
print()

# ==========================================
# 5. 统计分析
# ==========================================
print("=" * 80)
print("变异效应统计")
print("=" * 80)
print()

# 按效应分数排序
results_sorted = sorted(results, key=lambda x: abs(x["effect_score"]), reverse=True)

# 前20个最显著变异
print("前 20 个最显著变异（按效应分数绝对值排序）：")
print()
print(
    f"{'位置':>8} | {'REF':>4} | {'ALT':>4} | {'效应分数':>10} | {'似然差异':>10} | {'质量':>8}"
)
print("-" * 80)

for i, r in enumerate(results_sorted[:20]):
    print(
        f"{r['position']:>8} | {r['ref']:>4} | {r['alt']:>4} | {r['effect_score']:>10.3f} | {r['likelihood_diff']:>10.3f} | {r['qual']:>8.1f}"
    )

print()

# 变异类型统计
snp_count = sum(1 for r in results if len(r["ref"]) == 1 and len(r["alt"]) == 1)
indel_count = sum(1 for r in results if len(r["ref"]) != len(r["alt"]))

print(f"变异类型统计：")
print(f"  SNP: {snp_count:,} 个")
print(f"  INDEL: {indel_count:,} 个")
print()

# 效应分数统计
effect_scores = [r["effect_score"] for r in results]
print(f"效应分数统计：")
print(f"  平均: {np.mean(effect_scores):.3f}")
print(f"  标准差: {np.std(effect_scores):.3f}")
print(f"  最小值: {np.min(effect_scores):.3f}")
print(f"  最大值: {np.max(effect_scores):.3f}")
print(f"  中位数: {np.median(effect_scores):.3f}")
print()

# 似然差异统计
likelihood_diffs = [r["likelihood_diff"] for r in results]
print(f"似然差异统计：")
print(f"  平均: {np.mean(likelihood_diffs):.3f}")
print(f"  标准差: {np.std(likelihood_diffs):.3f}")
print(f"  最小值: {np.min(likelihood_diffs):.3f}")
print(f"  最大值: {np.max(likelihood_diffs):.3f}")
print()

# ==========================================
# 6. 保存结果
# ==========================================
print("正在保存结果...")

# 保存 JSON 格式
output_json = (
    sys.argv[3]
    if len(sys.argv) > 3
    else "/workdir/projects/test_rnaseq_analysis/data/results/evo2_variant_effects.json"
)

with open(output_json, "w") as f:
    json.dump(
        {
            "model": "evo2_1b_base",
            "reference": ref_record.id,
            "reference_length": len(ref_sequence),
            "num_variants": len(results),
            "statistics": {
                "snp_count": snp_count,
                "indel_count": indel_count,
                "effect_score": {
                    "mean": float(np.mean(effect_scores)),
                    "std": float(np.std(effect_scores)),
                    "min": float(np.min(effect_scores)),
                    "max": float(np.max(effect_scores)),
                    "median": float(np.median(effect_scores)),
                },
                "likelihood_diff": {
                    "mean": float(np.mean(likelihood_diffs)),
                    "std": float(np.std(likelihood_diffs)),
                    "min": float(np.min(likelihood_diffs)),
                    "max": float(np.max(likelihood_diffs)),
                },
            },
            "variants": results,
        },
        f,
        indent=2,
    )

print(f"✓ 结果保存到: {output_json}")
print()

# ==========================================
# 7. 保存详细表格（CSV格式）
# ==========================================
output_csv = output_json.replace(".json", ".csv")

with open(output_csv, "w") as f:
    f.write(
        "position,ref,alt,qual,effect_score,likelihood_diff,ref_max_logit,alt_max_logit\n"
    )
    for r in results:
        f.write(
            f"{r['position']},{r['ref']},{r['alt']},{r['qual']},{r['effect_score']:.3f},{r['likelihood_diff']:.3f},{r['ref_max_logit']:.3f},{r['alt_max_logit']:.3f}\n"
        )

print(f"✓ CSV 表格保存到: {output_csv}")
print()

print("=" * 80)
print("分析完成！")
print("=" * 80)
