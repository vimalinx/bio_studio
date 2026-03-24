# yeast_rnaseq_demo

## 项目类型
rnaseq

## 描述
Yeast RNA-seq demonstration pipeline using STAR and featureCounts.

## 项目结构

```
yeast_rnaseq_demo/
├── data/                  # 项目数据（完全独立）
│   ├── raw/               # 原始输入数据
│   ├── processed/         # 中间产物
│   ├── results/           # 最终结果
│   └── references/        # 项目内参考数据
├── scripts/               # 项目脚本
│   ├── config.py          # 项目配置与路径定义
│   ├── validate_project.py # 项目级自检入口
│   └── pipeline.py        # 主分析流程入口
├── notebooks/             # 探索性分析
├── logs/                  # 运行日志与验证报告
└── README.md              # 本文件
```

## 快速开始

### 1. 自检与查看步骤
```bash
cd scripts
python validate_project.py
python pipeline.py --validate
python pipeline.py --steps
```

### 2. 准备数据
将原始数据放入 `data/raw/`，文件名默认匹配 `*_R1.fastq.gz` / `*_R2.fastq.gz`。

### 3. 运行分析
```bash
cd scripts
python pipeline.py
```

## 说明

- `validate_project.py` 负责项目级结构与配置自检
- `pipeline.py --steps` 会列出索引构建、比对和定量三个阶段
- 这个 demo 现在优先通过工作区共享 runtime 调用 `lib/modules/rnaseq.py` / `samtools_wrapper.py`，同时保留本地 `counts_matrix.csv` 整理逻辑
- 如果项目被单独拷出工作区、找不到 `lib/`，会自动回退到脚本内置的 STAR + featureCounts 命令实现

## 项目状态

- [ ] 数据准备
- [ ] 质量控制
- [ ] 主要分析
- [ ] 结果验证
