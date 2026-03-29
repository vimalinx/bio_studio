# Cross-System Benchmark Design

> 本文档记录 2026-03-27 这一轮“Ebola + SARS-CoV-2 + 酵母”跨系统 benchmark 项目的设计决策。

## 目标

在 `projects/` 下新增一个可复用的专项项目，用统一目录和统一脚本管理三套公开 benchmark：

- Ebola virus: `GSE114905`
- SARS-CoV-2: `GSE147507`
- *Saccharomyces cerevisiae*: `GSE67149`（隶属 `GSE67151`）

第一版主目标只有两个：

1. 形成可复用 benchmark 资产层
2. 形成可直接运行的 baseline 分析层

附加目标：

- 输出一份把三套系统放在一起的跨系统生物学摘要报告

## 项目定位

这个项目不是单一物种的一条 RNA-seq 流水线，也不是全量 raw FASTQ 重分析仓库。  
它更接近一个“带分析能力的 benchmark panel”。

因此推荐把它作为 **专项项目** 落在：

- `projects/cross_system_benchmark/`

而不是强行塞进完全通用的模板语义。

## 为什么选这三套数据

### Ebola

- GEO: `GSE114905`
- 公开标题: `Transcriptome in Huh7 cells infected with 4 Ebolaviruses`
- 优点:
  - 数据集小
  - 有明确时间点 `1-3 dpi`
  - GEO 页面直接提供 processed Excel 文件
  - 非常适合第一版先走“processed-data baseline”

第一版只固定使用其中的 `EBOV 1/2/3 dpi` 条目，不把四种 ebolavirus 一起混成主比较。

### SARS-CoV-2

- GEO: `GSE147507`
- 公开标题: `SARS-CoV-2 launches a unique transcriptional signature from in vitro, ex vivo, and in vivo systems`
- 优点:
  - 经典 benchmark，复用广
  - GEO 页面直接提供 `RawReadCounts_Human.tsv.gz`
  - 可先从人类样本处理后矩阵切出一个干净子集做 baseline

第一版推荐优先使用 `NHBE mock vs SARS-CoV-2` 这组比较。

### Yeast

- GEO subseries: `GSE67149`
- super series: `GSE67151`
- 公开标题: `Rpd3 drives transcriptional quiescence`
- 优点:
  - 条件定义清楚
  - 有 `log`, `Q`, `DS` 等状态
  - GEO 页面直接提供 `Normalized_FPKM_Values.xlsx`

第一版固定使用 `wild type log vs wild type Q`。

## 数据与参考来源

推荐在项目元数据里固定以下官方来源：

- Ebola benchmark:
  - <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114905>
- SARS-CoV-2 benchmark:
  - <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507>
- Yeast benchmark:
  - <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67149>
- Ebola reference accession:
  - `NC_002549.1`
- SARS-CoV-2 reference accession:
  - `NC_045512.2`
- Yeast reference source:
  - SGD / S288C reference assets

第一版不默认下载大型宿主参考。  
如果后续需要 raw-level 子集重分析，再单独扩展人类参考准备逻辑。

## 目录设计

```text
projects/cross_system_benchmark/
├── data/
│   ├── raw/
│   │   ├── ebola/
│   │   ├── sars_cov_2/
│   │   └── yeast/
│   ├── processed/
│   │   ├── ebola/
│   │   ├── sars_cov_2/
│   │   └── yeast/
│   ├── results/
│   │   ├── per_dataset/
│   │   └── integrated/
│   └── references/
│       ├── ebola/
│       ├── sars_cov_2/
│       └── yeast/
├── metadata/
│   ├── datasets.yaml
│   ├── samples.tsv
│   ├── comparisons.tsv
│   └── source_links.md
├── scripts/
│   ├── config.py
│   ├── validate_project.py
│   ├── pipeline.py
│   ├── fetch_benchmarks.py
│   ├── prepare_references.py
│   ├── run_baseline.py
│   └── summarize_cross_systems.py
├── logs/
├── notebooks/
└── README.md
```

## 工作流设计

### 阶段 1: benchmark catalog

输出：

- `metadata/datasets.yaml`
- `metadata/samples.tsv`
- `metadata/comparisons.tsv`
- `metadata/source_links.md`

这一层的任务是把 accession、条件、推荐对比、参考来源和 processed 文件链接固定下来。  
这一层不依赖大规模下载。

### 阶段 2: reference preparation

输出：

- `data/references/ebola/*`
- `data/references/sars_cov_2/*`
- `data/references/yeast/*`

第一版优先准备：

- Ebola viral reference
- SARS-CoV-2 viral reference
- Yeast reference

对于大型宿主参考，只在元数据里记录推荐来源，不默认拉取。

### 阶段 3: baseline analysis

这一层默认基于 GEO 可直接下载的 processed matrices / counts 做分析，而不是默认全量拉 raw FASTQ。

每个系统至少输出：

- 样本子集说明
- 基础表达矩阵摘要
- 一个二维投影或层次聚类图
- 一个 top feature 表
- 一份 dataset summary

### 阶段 4: integrated summary

输出：

- `data/results/integrated/cross_system_summary.md`

这里不做伪严谨的基因一对一跨物种比较。  
比较维度只放在更稳的生物学主题层面，例如：

- 宿主应激与抗病毒反应
- 转录与翻译负担
- 代谢压制与状态切换
- 病毒感染与静息状态的系统性差异

## 第一版明确做什么

- 新增专项项目 `cross_system_benchmark`
- 固定三套 benchmark 的官方元数据与下载入口
- 提供可重复生成的 benchmark catalog
- 下载或准备小体量参考资产
- 基于 processed data 跑三套 baseline
- 生成单数据集报告与综合报告

## 第一版明确不做什么

- 不默认下载全量 SRA raw FASTQ
- 不做完整的人类宿主参考构建
- 不做复杂宿主-病毒联合比对
- 不做单细胞、多组学或结构联动
- 不做深度学习模型训练 benchmark

## 验证标准

项目第一版完成至少满足：

- `python scripts/project.py validate cross_system_benchmark` 通过
- `python scripts/project.py steps cross_system_benchmark` 能列出阶段
- `pipeline.py fetch` 能生成 benchmark catalog
- `pipeline.py prepare` 能生成参考准备结果或明确的 reference manifest
- `pipeline.py baseline` 能为三个系统各生成一份结果
- `pipeline.py report` 能生成 integrated 附加报告

## 风险与缓解

### 风险 1: GEO 补充文件格式不完全统一

缓解：

- 每个 dataset 用自己独立的 parser
- 在项目级脚本层统一输出到共同的中间格式

### 风险 2: 样本命名不统一

缓解：

- `samples.tsv` 作为规范化真源
- 下游脚本只依赖标准化字段，不直接依赖原始列名

### 风险 3: 跨系统比较容易做成“看起来有道理但方法很虚”

缓解：

- integrated 报告只做主题级比较
- 不做虚假的逐基因直接映射结论

## 设计结论

这不是“一个更大的 RNA-seq 模板项目”，而是一个 **跨系统 benchmark panel**。  
因此第一版的关键不是吞下所有 raw data，而是先把：

- benchmark 资产层
- processed-data baseline
- 跨系统摘要

三层稳稳接起来。
