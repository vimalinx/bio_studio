# cross_system_benchmark

## 项目类型
special benchmark project

## 描述
Cross-system benchmark panel for Ebola, SARS-CoV-2, and yeast using curated public GEO datasets and lightweight processed-data baselines.

## 项目结构

```text
cross_system_benchmark/
├── data/
│   ├── raw/               # 下载到本地的 GEO processed inputs
│   ├── processed/         # 归一化后的表达长表、样本表和摘要
│   ├── results/           # 单数据集与综合报告
│   └── references/        # 小体量参考与参考 manifest
├── metadata/              # benchmark catalog、样本表、对比定义
├── scripts/
│   ├── config.py
│   ├── validate_project.py
│   ├── pipeline.py
│   ├── common.py
│   ├── fetch_benchmarks.py
│   ├── prepare_references.py
│   ├── run_baseline.py
│   └── summarize_cross_systems.py
├── notebooks/
├── logs/
└── README.md
```

## 快速开始

### 1. 自检与查看步骤
```bash
cd scripts
python validate_project.py
python pipeline.py --validate
python pipeline.py --steps
```

### 2. 生成 benchmark catalog
```bash
python pipeline.py fetch
```

### 3. 生成 reference manifest 并准备默认小体量参考
```bash
python pipeline.py prepare
```

### 4. 下载 processed benchmark 文件并跑 baseline
```bash
python pipeline.py baseline
```

### 5. 生成单数据集与跨系统报告
```bash
python pipeline.py report
```

## 当前固定 benchmark

- Ebola: `GSE114905`
- SARS-CoV-2: `GSE147507`
- Yeast: `GSE67149`

## 说明

- 这是一个 **专项项目**，不是单一 RNA-seq 模板项目
- 第一版主目标是 benchmark catalog 和 processed-data baseline
- 不默认下载全量 raw FASTQ
- 不默认构建大型人类宿主参考
- `metadata/` 是 accession、样本和比较定义的真源

## 结果产物

- `metadata/datasets.yaml`
- `metadata/samples.tsv`
- `metadata/comparisons.tsv`
- `metadata/source_links.md`
- `data/processed/<dataset>/expression.tsv`
- `data/results/per_dataset/*`
- `data/results/integrated/cross_system_summary.md`

## 项目状态

- [ ] 模板自检通过
- [ ] benchmark catalog
- [ ] references prepared
- [ ] baseline outputs
- [ ] integrated summary
