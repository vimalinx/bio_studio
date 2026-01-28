#!/usr/bin/env python3

import sys
import json
import torch
import numpy as np
from Bio import SeqIO
from evo2 import Evo2
import gc

print("="*80)
print("Evo 2 埃博拉病毒 GP 基因饱和突变扫描")
print("="*80)

# 配置
MODEL_NAME = "evo2_1b_base"
DEVICE = "cuda:0" if torch.cuda.is_available() else "cpu"
GP_START = 5900  # 1-based start
GP_END = 8305    # 1-based end
CONTEXT_SIZE = 50 # 减小上下文以节省显存
BATCH_SIZE = 100  # 批处理大小
OUTPUT_FILE = "/workdir/projects/test_rnaseq_analysis/data/results/gp_saturation_scores.csv"

print(f"设备: {DEVICE}")
print(f"GP 基因范围: {GP_START}-{GP_END}")

# 加载模型
print("正在加载模型...")
evo2_model = Evo2(MODEL_NAME)
print("模型加载完成")

# 读取参考基因组
with open("/workdir/shared_data/references/ebov_zaire.fa", "r") as f:
    ref_record = SeqIO.read(f, "fasta")
ref_sequence = str(ref_record.seq)

# 生成突变列表
print("正在生成突变列表...")
mutations = []
bases = ['A', 'C', 'G', 'T']

# 遍历 GP 基因的每个位置
# 注意：GP_START 是 1-based，Python 切片是 0-based
# 0-based 索引: GP_START-1 到 GP_END
for i in range(GP_START-1, GP_END):
    ref_base = ref_sequence[i]
    pos_1based = i + 1
    
    # 对每个位置生成 3 个突变
    for alt_base in bases:
        if alt_base != ref_base:
            mutations.append({
                "pos": i,  # 0-based index in genome
                "ref": ref_base,
                "alt": alt_base,
                "pos_1based": pos_1based
            })

print(f"总突变数: {len(mutations)}")

# 初始化结果文件
with open(OUTPUT_FILE, "w") as f:
    f.write("position,ref,alt,effect_score\n")

# 批处理预测
print("开始预测...")
total_processed = 0

# 为了效率，我们只计算一次参考序列的 logits
# 但由于每个突变的上下文可能略有不同（边缘效应），最准确的方法是每批重新计算
# 这里为了准确性，我们对每个突变计算差异

def process_batch(batch_mutations):
    results = []
    
    for mut in batch_mutations:
        pos = mut["pos"]
        ref = mut["ref"]
        alt = mut["alt"]
        
        start = max(0, pos - CONTEXT_SIZE)
        end = min(len(ref_sequence), pos + 1 + CONTEXT_SIZE)
        
        # 提取上下文
        ref_context = ref_sequence[start:end]
        alt_context = ref_sequence[start:pos] + alt + ref_sequence[pos+1:end]
        
        try:
            # Tokenize
            ref_tokens = evo2_model.tokenizer.tokenize(ref_context)
            alt_tokens = evo2_model.tokenizer.tokenize(alt_context)
            
            ref_input = torch.tensor([ref_tokens], dtype=torch.long).to(DEVICE)
            alt_input = torch.tensor([alt_tokens], dtype=torch.long).to(DEVICE)
            
            # 预测
            with torch.no_grad():
                # Evo 2 returns (logits_tuple, embeddings)
                # logits_tuple is (sequence_logits, None)
                ref_out = evo2_model(ref_input)
                alt_out = evo2_model(alt_input)
                
                ref_logits = ref_out[0][0][0] # [seq_len, vocab_size]
                alt_logits = alt_out[0][0][0]
            
            # 计算效应
            # 找到突变位点在上下文中的索引
            # ref_context = ... + ref + ...
            # 突变位点在 tokens 中的索引大致为 CONTEXT_SIZE
            # 但由于 tokenizer 可能合并字符，我们需要精确计算
            
            # 简化策略：使用序列中间位置
            # 对于字符级 tokenizer (Evo 2 default)，每个字符一个 token
            # 突变位点索引 = pos - start
            
            idx = pos - start
            
            if idx < len(ref_logits) and idx < len(alt_logits):
                ref_score = torch.max(ref_logits[idx]).item()
                alt_score = torch.max(alt_logits[idx]).item()
                effect = alt_score - ref_score
                
                results.append(f"{mut['pos_1based']},{ref},{alt},{effect:.4f}\n")
            
        except Exception as e:
            # print(f"Error at {pos}: {e}")
            pass
            
    return results

# 分批执行
for i in range(0, len(mutations), BATCH_SIZE):
    batch = mutations[i:min(i+BATCH_SIZE, len(mutations))]
    
    batch_results = process_batch(batch)
    
    # 写入结果
    with open(OUTPUT_FILE, "a") as f:
        for line in batch_results:
            f.write(line)
            
    total_processed += len(batch)
    if total_processed % 500 == 0:
        print(f"进度: {total_processed}/{len(mutations)} ({(total_processed/len(mutations))*100:.1f}%)")
        sys.stdout.flush()
        
    # 清理内存
    torch.cuda.empty_cache()
    if total_processed % 1000 == 0:
        gc.collect()

print("扫描完成！")
