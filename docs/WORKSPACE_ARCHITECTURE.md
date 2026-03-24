# Bio Studio 工作区架构

> 这份文档定义 Bio Studio 当前工作区的结构边界。目标不是描述每个文件，而是说明什么是核心工作面、什么是第三方引擎、什么是共享资源、什么只是遗留或保留位。

## 一句话定位

Bio Studio 是一个面向 AI agent 的生物计算工作台，而不是单一应用仓库。

它由五类东西组成：
- 行为协议：告诉 AI 怎么工作
- 本地骨架：提供项目模板、脚手架和维护脚本
- AI 接口：skills 和 MCP servers
- 项目工作区：真实分析任务与产物
- 外部引擎：Evo2、RFdiffusion、Biomni、LiteFold 等第三方能力

## 顶层目录角色

| 路径 | 角色 | 应该如何使用 |
|---|---|---|
| `README.md` / `CLAUDE.md` | 总入口与操作规则 | 新 agent / 新协作者优先阅读 |
| `docs/` | 协议、环境、架构与记录 | 放工作区级知识，不放项目专属结果 |
| `lib/` | 本地可复用骨架 | 放通用模块、项目模板、CLI 包装层 |
| `scripts/` | 系统维护脚本 | 放环境检查、安装、维护类脚本 |
| `mcp-servers/` | AI 工具接口层 | 放可独立配置的 MCP servers |
| `.claude/skills/` | agent 技能层 | 放工作流与能力提示，不放业务结果 |
| `projects/` | 任务工作区 | 所有真实分析、报告、日志、结果优先落这里 |
| `shared_data/` | 共享参考资源 | 放跨项目可复用的数据和索引 |
| `repositories/active/` | 第三方引擎源码 | 放当前启用的外部仓库，不作为主开发面 |
| `repositories/archived/` | 归档引擎区 | 放退役或备用的外部仓库 |

## 核心工作面

下面这些目录是 Bio Studio 的主工作面：

1. `docs/`
原因：这里定义协议、环境状态和架构边界，是 AI 的规则层。

2. `projects/`
原因：这里是实际分析的现场。数据、脚本、日志、报告都应该优先在这里沉淀。

3. `lib/`、`scripts/`、`mcp-servers/`
原因：这是工作台自己的基础设施层，负责把环境和外部能力编排成可调用工作流。

## 第三方引擎边界

`repositories/active/` 下的内容不属于 Bio Studio 的主产品代码，它们是挂进工作区的外部能力层。

当前主要包括：
- `RFdiffusion`
- `evo2`
- `Biomni`
- `litefold/source`

治理规则：
- 默认不要在这些目录里做本地架构重构
- 如果只是想“使用能力”，优先在 `projects/`、`tools/scripts/`、`mcp-servers/` 或 `lib/` 层编排
- 如果必须修改第三方仓库，先确认这是本地 patch 还是准备长期维护的 fork

## 项目工作区规则

`projects/` 是分析工作的首选落点。

推荐结构：

```text
projects/<name>/
├── data/
│   ├── raw/
│   ├── processed/
│   ├── results/
│   └── references/
├── scripts/
├── logs/
├── notebooks/
└── README.md
```

规则：
- 原始数据放 `data/raw/`
- 中间结果放 `data/processed/`
- 最终结果放 `data/results/`
- 可复现逻辑优先写进 `scripts/`
- 根目录不要新增散落的 `data/`、`notebooks/`、`test/` 等临时目录

并不是 `projects/` 下所有目录都必须收敛成同一类模板项目。当前至少有两类“特殊项目”：
- 专项项目：例如 `projects/test_rnaseq_analysis/`，保留面向具体研究问题的定制流程与结果现场
- 学习项目：例如 `projects/yeast_genome_learning/`，保留教程式脚本与教学材料

这些项目可以补齐最小统一入口（如 `validate_project.py`、`pipeline.py --validate/--steps`），但不应为追求统一而强行重写成模板 runtime 项目。

对应到工作区 CLI：
- `python scripts/project.py validate <项目名>` / `steps <项目名>` 应尽量都可用
- `python scripts/project.py run <学习项目>` 可以只返回教学引导，而不是直接批量执行所有 lesson 脚本
- `python scripts/project.py run <专项项目>` 仍然可以转发到项目自己的定制 pipeline

## 共享资源规则

`shared_data/` 只放跨项目共享的数据，不放某个单独项目的临时结果。

适合放：
- 参考基因组
- 通用索引
- 可复用数据库
- 工具元数据

不适合放：
- 单个项目独有的 FASTQ / BAM / VCF
- 一次性实验结果
- 仍在频繁改写的临时文件

## 保留位与遗留区

### `docs/site/`
正式的工作区介绍页所在位置。它是文档入口的一部分，不是一次性演示文件。

### `mcp-servers/bio-design-mcp/` 与 `mcp-servers/bio-lab-mcp/`
这两个目录已经是已落地的 MCP 接口层。
- `bio-design-mcp` 提供 DNA 设计、mRNA 优化与靶点评估入口
- `bio-lab-mcp` 提供项目发现、项目级自检和工作区级 smoke test 入口
- 两者都通过 `mcp-servers/render_claude_config.py` 进入统一配置渲染链路

### `_legacy_workspace_v1/`
旧工作区结构的保留区，用于保存历史布局痕迹。默认不要把新工作继续写回这里。

### `tmp/`
临时缓冲区。只适合短生命周期文件，不应成为长期存储目录。

### `evo2-7b-hf/`
为 Evo2 HuggingFace 权重转换预留的位置。当前为空时，表示该工作面尚未启用。

## 当前状态说明

截至本次整理：
- `docs/site/index.html` 已成为正式介绍入口
- `scripts/project.py` 已成为工作区级统一项目入口
- `mcp-servers/bio-design-mcp/` 与 `mcp-servers/bio-lab-mcp/` 已可作为真实服务接入 Claude / agent 配置
- `.gitmodules` 已登记 `RFdiffusion` 与 `litefold/source` 两个 gitlink 路径
- 外部仓库仍可能存在本地未提交改动，这属于工作树状态，不属于本次治理变更范围

## 推荐阅读顺序

1. `README.md`
2. `CLAUDE.md`
3. `docs/README.md`
4. `docs/WORKSPACE_ARCHITECTURE.md`
5. `projects/test_env_validation/README.md`
6. `projects/ai_design_playground/README.md`
7. `projects/test_rnaseq_analysis/README.md`

## 设计原则

- AI 优先，不等于人类不可读
- 项目隔离优先于“根目录方便”
- 第三方引擎与本地骨架必须分层
- 工作区知识应沉淀到 `docs/`，项目知识应沉淀到 `projects/`
- 先把边界说清，再做更深的功能扩展
