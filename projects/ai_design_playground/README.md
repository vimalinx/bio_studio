# ai_design_playground

## 项目类型
generic

## 描述
In-silico playground: generate toy sequences, analyze, and validate with local tools/AI models.

## 项目结构

```
ai_design_playground/
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

这个项目是“纯干实验/纯仿真”的玩具工作流：生成一批合成 toy DNA（不对应真实生物体），然后用本地工具做简单分析并产出报告。

### 1. 自检与查看步骤
```bash
cd scripts
python validate_project.py
python pipeline.py --validate
python pipeline.py --steps
```

### 2. 运行分析
不需要手动准备数据；pipeline 会自行生成 toy 输入。

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

## 说明

- `validate_project.py` 负责项目级结构与配置自检
- `pipeline.py --steps` 会列出 toy 数据生成、本地分析和报告输出三个阶段
- 这个项目保留自定义分析逻辑，不强制收敛到通用模板的四步骨架
- 在工作区内运行时，底层 `seqkit` / `Prodigal` / `RNAfold` 调用会优先走共享 runtime 和 `lib/modules`
- 如果项目被单独拷出工作区、找不到 `lib/`，会自动回退到项目脚本内置实现

## 项目状态

- [ ] 数据准备
- [ ] 质量控制
- [ ] 主要分析
- [ ] 结果验证
