# example_rnaseq

## 项目类型
rnaseq

## 描述
示例RNA-seq分析项目

## 项目结构

```
example_rnaseq/
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

### 1. 自检模板
```bash
cd scripts
python validate_project.py
python pipeline.py --validate
python pipeline.py --steps
```

### 2. 准备数据
将原始数据放入 `data/raw/`。

### 3. 编辑配置
按项目需要更新 `scripts/config.py` 中的样本、参考和线程数。

### 4. 运行流程
```bash
cd scripts
python pipeline.py
```

## 说明

- `validate_project.py` 负责项目级结构与配置自检
- `logs/validation_report.json` 和 `logs/validation_report.md` 会记录验证结果
- 缺少真实数据或参考文件时默认给出警告，不阻断模板自检

## 类型骨架

这个模板预置了 **RNA-seq** 最小骨架，默认围绕 `FastQC/Fastp`、`STAR`、`featureCounts` 和 `MultiQC` 组织。

- 推荐把参考基因组放在 `data/references/genome.fa`
- 推荐把注释放在 `data/references/annotation.gtf`
- 双端测序默认使用 `*_R1.fastq.gz` / `*_R2.fastq.gz`


## 项目状态

- [ ] 模板自检通过
- [ ] 数据准备
- [ ] 质量控制
- [ ] 主要分析
- [ ] 结果验证
