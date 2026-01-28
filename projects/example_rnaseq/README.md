# example_rnaseq

## 项目类型
rnaseq

## 描述
示例RNA-seq分析项目

## 项目结构

```
example_rnaseq/
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

### 1. 准备数据
将原始数据放入 `data/raw/`

### 2. 编辑配置
编辑 `scripts/config.py` 设置项目参数

### 3. 运行分析
```bash
cd scripts
python pipeline.py
```

## 项目状态

- [ ] 数据准备
- [ ] 质量控制
- [ ] 主要分析
- [ ] 结果验证

## 笔记
记录你的分析过程和发现
