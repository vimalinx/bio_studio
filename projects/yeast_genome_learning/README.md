# 酵母菌基因组学习项目

## 项目形态

这是一个 **学习项目**。

- 主入口是分阶段 Bash 教学脚本，而不是统一分析模板 pipeline
- 它不强制并入通用模板 shared runtime
- 当前只补齐最小统一入口：`validate_project.py` 与 `pipeline.py --validate/--steps`

## 📚 项目概述
从酿酒酵母 (*Saccharomyces cerevisiae*) 开始学习基因组数据库分析。

## 🎯 学习目标
- ✅ 理解基因组数据库结构
- ✅ 掌握序列分析基础
- ✅ 学习基因注释解读
- ✅ 熟练使用 BLAST 工具

## 📁 项目结构
```
yeast_genome_learning/
├── data/              # 原始数据
├── scripts/           # 分析脚本
├── results/           # 分析结果
└── logs/              # 运行日志
```

## 🚀 快速开始

### 1. 安装酵母菌数据库
```bash
cd scripts
bash 01_setup_database.sh
```

### 2. 验证安装
```bash
bash 02_verify_install.sh
```

### 3. 运行示例分析
```bash
bash 03_extract_gene.sh ACT1
```

### 4. 统一入口
```bash
python pipeline.py --validate
python pipeline.py --steps
```

## 📊 酵母菌基因组信息
- 物种: *Saccharomyces cerevisiae* (酿酒酵母)
- 菌株: S288C
- 版本: R64-1-1
- 基因组大小: ~12.1 Mb
- 染色体数量: 16 + 线粒体
- 基因数量: ~6,000

## 💡 推荐学习顺序
1. 了解数据库结构 (scripts/01_setup_database.sh)
2. 验证安装 (scripts/02_verify_install.sh)
3. 提取特定基因 (scripts/03_extract_gene.sh)
4. 序列分析 (scripts/04_sequence_analysis.sh)
5. BLAST 搜索 (scripts/05_blast_search.sh)

## 📖 参考资源
- SGD: https://www.yeastgenome.org
- NCBI Yeast: https://www.ncbi.nlm.nih.gov/genome/?term=Saccharomyces+cerevisiae

## 📝 日志
- 2026-01-25: 项目创建
