# ai_design_playground

## 项目类型
generic

## 描述
In-silico playground: generate toy sequences, analyze, and validate with local tools/AI models.

## 项目结构

```
ai_design_playground/
├── data/                  # 项目数据（完全独立）
│   ├── raw/              # 原始数据
│   ├── processed/         # 中间结果
│   ├── results/          # 最终结果
│   └── references/       # 项目特定参考序列
├── scripts/              # 项目脚本（调用lib模块）
│   ├── pipeline.py       # 主要分析流程
│   ├── config.py         # 项目配置
│   └── analysis.py      # 项目特定分析
├── notebooks/           # Jupyter notebooks
├── logs/               # 日志文件
└── README.md           # 本文件
```

## 使用方法

这个项目是“纯干实验/纯仿真”的玩具工作流：生成一批合成 toy DNA（不对应真实生物体），然后用本地工具做简单分析并产出报告。

### 1. 准备数据
不需要手动准备数据。

### 2. 编辑配置
编辑 `scripts/config.py` 设置项目参数

### 3. 运行分析
```bash
cd scripts
python pipeline.py
```

运行后主要输出：
- `data/raw/toy_dna.fa`
- `data/processed/prodigal_proteins.faa` / `data/processed/prodigal_genes.fna`
- `data/results/report.md` / `data/results/report.json`

可选：加一个轻量的 AI 校验步骤（ESM contact 预测，会按需下载模型）：
```bash
python pipeline.py --use-esm-contacts --esm-model esm2_t6_8M_UR50D
```

## 项目状态

- [ ] 数据准备
- [ ] 质量控制
- [ ] 主要分析
- [ ] 结果验证

## 笔记
记录你的分析过程和发现
