# Bio Studio - AI协同生物工作区 v2.1

> 一个简洁高效的AI驱动的生物信息学研究环境

**核心理念**: 项目隔离，原子化模块，AI负责执行

---

## 📚 关键文档 (必读)

- **[文档索引](docs/README.md)** 📌: 所有文档入口与导航。
- **[避坑指南](docs/BEST_PRACTICES.md)** ⚠️: 包含 PATH 设置、Evo2 模型输出处理等关键经验。
- **[环境清单](docs/ENVIRONMENT.md)** 🛠️: 当前已安装的所有 Python 库和生信工具版本快照。
- **[工具部署协议](docs/AI_TOOL_PROTOCOL.md)** 🤖: AI 部署新工具的完整流程。
- **[分析执行协议](docs/AI_ANALYSIS_PROTOCOL.md)** 🧪: 分析流程规范与归档要求。

---

## 🎯 这是什么？

Bio Studio 是重构后的生物信息学工作区，核心改进：

- ✅ **AI 深度集成** - 内置 Evo 2 (基因组大模型)、Biomni 等 AI Agent。
- ✅ **全栈环境** - 预装 50+ 核心生信工具 (BWA, STAR, MultiQC, etc.)。
- ✅ **项目隔离** - 每个项目完全独立，数据和脚本互不干扰。
- ✅ **原子化模块** - 通用功能封装在 `lib/`，可复用。

**你只需要用自然语言告诉AI要做什么，它会自动选择工具、执行任务、生成报告。**

---

## 🚀 5分钟快速开始

```bash
# 1. 进入工作区
cd ~/bio_studio

# 2. 激活环境 (如果尚未激活)
source ~/miniforge3/bin/activate bio

# 3. 创建新项目
python3 lib/create_project.py my_analysis --type rnaseq

# 4. 运行分析
cd projects/my_analysis/scripts
python pipeline.py
```

---

## 🔧 已安装工具 (v2.1)

> 版本以 `docs/ENVIRONMENT.md` 为准。

### 🧬 核心生信工具 (CLI)
| 类别 | 工具 | 用途 |
|---|---|---|
| **质控** | `fastp`, `FastQC`, `MultiQC` | 极速质控与报告汇总 |
| **比对** | `BWA-MEM`, `Bowtie2`, `STAR`, `HISAT2` | DNA/RNA 序列比对 |
| **处理** | `Samtools`, `Bcftools`, `Bedtools`, `SeqKit` | BAM/VCF/FASTA 操作 |
| **分析** | `FeatureCounts`, `IQ-TREE 2`, `MAFFT` | 定量、进化树、多序列比对 |

### 🤖 AI 与 深度学习
| 模型/库 | 用途 | 状态 |
|---|---|---|
| **Evo 2 (1B)** | 基因组基础模型，预测变异效应 | ✅ 已部署 (Docker) |
| **ESM** | 蛋白质语言模型 | ✅ 已安装 (Python) |
| **Biomni** | 生物医学 AI Agent | ✅ 已部署 |
| **PyTorch 2.5** | 深度学习后端 | ✅ 支持 CUDA (RTX 5070) |

---

## 📁 目录结构

```
bio_studio/
│
├── lib/                           # ⭐ 通用库（原子化模块）
│   ├── modules/                   # 可复用的Python模块 (qc, alignment, variant...)
│   └── create_project.py          # 项目模板生成器
│
├── projects/                      # ⭐ 项目目录（每个项目独立）
│   └── <project_name>/
│       ├── data/                  # 该项目的数据
│       ├── scripts/               # 项目脚本
│       └── README.md              # 项目说明
│
├── scripts/                       # ⭐ 系统脚本
│   └── maintenance/               # 维护与环境脚本 (update_env, generate_report)
│
├── repositories/                  # 外部工具源码
│   └── active/                    # 当前使用的仓库
│       ├── Biomni/                # 生物医学 AI Agent
│       ├── evo2/                  # Evo 2 基因组大模型
│       └── RFdiffusion/           # 蛋白质设计工具
├── shared_data/                   # 共享参考基因组/数据库
└── docs/                          # 系统文档
```

---

## 📝 开发日志

**v2.1 (2026-01-28)**:
- 🚀 **实战验证**: 完成酵母菌全流程 RNA-seq 分析 (STAR + featureCounts + Variant Calling)。
- 🧪 **自动化测试**: 建立 `test_env_validation` 环境，一键验证核心工具链状态。
- 🧹 **架构治理**: 实施严格的根目录洁癖策略，归档旧数据，确立 "Federal" 联邦制项目架构。
- 🐛 **脚本修复**: 修正 Evo 2 显存管理、BCFtools 变异检测流程及多个可视化脚本 Bug。

**v2.0 (2025-01-22)**:
- ✨ 重构：项目独立架构
- ✨ 新增：lib/通用模块库

---

**现在就开始你的生物研究之旅吧！** 🚀

有问题？直接问AI。
