# 埃博拉病毒 RNA-seq 分析项目 (SRR1972739)

本项目演示了一个完整的、由 AI 驱动的生物信息学分析流程，从原始测序数据处理到深度学习辅助的功能预测。

## 📊 项目概览

- **样本**: SRR1972739 (Zaire ebolavirus, RNA-seq)
- **分析工具**: FastQC, Bowtie2, SAMtools, bcftools, Evo 2 (1B Base)
- **主要成果**: 
  - 修正了物种不匹配问题
  - 优化了序列比对率 (+10%)
  - 识别了 568 个高质量变异
  - 利用 Evo 2 深度学习模型预测了 GP 基因的功能热点

## 📂 目录结构

```
test_rnaseq_analysis/
├── data/
│   ├── raw/                # 原始测序数据 (FASTQ)
│   ├── processed/          # 中间文件 (BAM, VCF, QC)
│   └── results/            # 最终分析结果 (JSON, CSV)
│       ├── variant_statistics.json    # 变异基础统计
│       ├── evo2_variant_effects.json  # 变异效应预测
│       ├── gp_saturation_scores.csv   # GP 基因饱和突变数据
│       └── gp_analysis.json           # GP 基因功能热点分析
├── scripts/                # 分析脚本
│   └── evo2_gp_saturation.py  # Evo 2 饱和突变扫描脚本
├── ANALYSIS_REPORT.md      # 初始分析报告 (含错误诊断)
├── EBOV_ANALYSIS_REPORT.md # 修正后的比对报告
├── FINAL_ANALYSIS_REPORT.md # 最终综合分析报告
├── EVO2_ENHANCED_REPORT.md # Evo 2 增强分析报告
├── GP_EVO2_ANALYSIS_REPORT.md # GP 基因深度分析报告
└── GP_VALIDATION_REPORT.md # 实验数据验证报告
```

## 📝 主要报告

1. **[Evo 2 增强分析报告](./EVO2_ENHANCED_REPORT.md)**
   - 项目的核心总结，包含所有分析步骤和结果概览。

2. **[GP 基因功能热点报告](./GP_EVO2_ANALYSIS_REPORT.md)**
   - 详细展示了 Evo 2 对埃博拉病毒糖蛋白 (GP) 的饱和突变分析结果。
   - 识别了药物靶点（保守区）和免疫逃逸位点（高容忍区）。

3. **[GP 功能验证报告](./GP_VALIDATION_REPORT.md)**
   - 将 Evo 2 的预测结果与 Uniprot 已知功能数据进行比对验证。
   - 证明了 AI 模型准确"重现"了生物学知识。

## 🚀 复现指南

### 环境准备
项目依赖 Docker 环境运行深度学习模型。请确保安装了 Docker 和 NVIDIA Container Toolkit。

```bash
# 启动 Evo 2 环境
cd ../../
docker compose -f evo2/docker-compose.override.yml run --rm -d --name evo2-analysis evo2 sleep infinity
```

### 运行分析
脚本位于 `scripts/` 目录。

```bash
# 运行 GP 基因饱和突变扫描
docker exec evo2-analysis python3 /workdir/projects/test_rnaseq_analysis/scripts/evo2_gp_saturation.py
```

## 📈 关键发现摘要

- **高病毒载量**: 平均覆盖深度 >3700x。
- **G→A 过渡偏好**: 占所有变异的 76%，提示 RNA 编辑活性。
- **药物靶点**: GP 基因的 Furin 切割位点 (R501) 和关键糖基化位点 (N454) 高度保守。
- **逃逸风险**: GP 基因的 Mucin-like Domain (MLD) 具有极高的突变容忍度。

---
**最后更新**: 2026-01-22
