# Biology Skill Factory Design

> 本文档记录 2026-03-27 这一轮“把生物信息学工具系统化 skill 化”的设计决策。

## 目标

在当前工作区内建立一套可复用的 **biology skill factory**，把现有生物学 skill、工作区已安装的生信工具，以及后续新增工具，统一纳入同一套规范与自动生成流程中。

核心目标有三层：

1. 给当前项目写一份权威的 skill 规范文档
2. 生成一个当前工作区的 biology skill catalog
3. 提供脚本，自动从本地工具帮助信息与文档中提取描述，并批量生成单工具 skill

## 当前现状

当前工作区内已经存在 9 个生物学相关 skill：

- `bioinformatics-toolkit`
- `biomni`
- `evo2`
- `phage-design`
- `protein-structure`
- `rfdiffusion`
- `rnaseq-pipeline`
- `sequence-analysis`
- `yeast_database`

这些 skill 主要是“领域级 skill”，并不是“一工具一 skill”的细粒度布局。  
例如 `bioinformatics-toolkit` 已经在一个 skill 里揉进了 BLAST、HMMER、bowtie2、bwa、samtools、prodigal、ViennaRNA、Biopython 等多类工具。

因此当前问题不是“完全没有 skill”，而是：

- 缺少统一规范
- 缺少统一 catalog
- 缺少自动发现与批量生成
- 缺少把现有领域 skill 与单工具 skill 协同管理的结构

## 设计判断

### 1. 保留现有 skill，不推倒重来

现有 skill 已经承载领域工作流，直接删除或拆碎会造成工作流倒退。  
本轮只在其上增加新的 **单工具 skill 层**，而不是重写已有 skill。

### 2. 单工具 skill 默认偏手动触发

如果一下子生成几十到上百个可自动触发的工具 skill，Claude 的技能描述上下文会被挤爆。  
因此自动生成的单工具 skill 第一版统一采用：

- `disable-model-invocation: true`
- `user-invocable: true`

也就是默认用户手动 `/tool-name` 调用，避免模型自动乱触发。

### 3. “工具来源”先以工作区可发现资产为准

本轮不从全互联网穷举所有生物信息工具。  
第一版只收敛到三个来源：

1. `scripts/maintenance/install_bio_tools.sh` 里列出的工具
2. 当前 `bio` 环境 `PATH` 上可发现的生信工具
3. 现有 `.claude/skills/bioinformatics-toolkit/TOOLS.md` 与相关 skill 里已经文档化的工具

这样可控、可重复、可验证。

### 4. LLM 负责生成草稿，脚本负责结构化

自动生成 skill 的质量不能完全依赖 `man` 或 `-h` 原文。  
因此采用双层策略：

- 脚本负责收集帮助信息、路径、版本、来源、模板骨架
- LLM 负责把这些原始材料压成面向 Claude Skill 的 `SKILL.md`

但输出必须走统一模板，不能每个工具风格漂移。

### 5. 密钥与服务地址不写进仓库

用户提供了自定义 API key、base URL、模型名和思考强度。  
这些值会被工具脚本使用，但不应硬编码进 git 跟踪文件。  
仓库中只保留：

- `config.example`
- 环境变量读取逻辑
- 可选本地未跟踪配置入口

## 目标目录结构

```text
bio_studio/
├── docs/
│   └── skills/
│       ├── biology-skill-spec.md
│       └── biology-skill-catalog.md
├── scripts/
│   └── skills/
│       ├── discover_bio_tools.py
│       ├── generate_bio_skills.py
│       ├── render_skill_from_help.py
│       └── config.example.json
└── .claude/
    └── skills/
        ├── bioinformatics-toolkit/
        ├── ...
        └── <tool-name>/
            ├── SKILL.md
            └── references/
                └── help.md
```

## 自动生成流程

### 阶段 1: discover

扫描并输出候选工具 catalog。

来源包括：

- 安装脚本里的工具包名
- PATH 上可执行命令
- 现有 skill 文档里明确点名的工具

输出：

- 结构化 JSON/Markdown catalog
- 每个工具的命令名、来源、路径、版本、是否有 `man`

### 阶段 2: collect help

为每个工具采集原始材料：

- `tool --help`
- `tool -h`
- `man tool`
- `tool --version`

并写入每个工具自己的 `references/help.md`。

### 阶段 3: render skill

把结构化元数据和帮助文本喂给 LLM，生成统一模板的 `SKILL.md`。

统一约束：

- `name` 使用工具名或规范化后的短名
- `description` 要写明“何时用”
- 默认 `disable-model-invocation: true`
- `allowed-tools` 尽量只开放与该工具直接相关的 Bash 权限
- 正文不超过约定长度

### 阶段 4: catalog refresh

刷新总览文档，列出：

- 已有领域 skill
- 自动生成的单工具 skill
- 每个工具 skill 的来源与状态

## 第一版范围

第一版要完成：

- 写 skill 规范文档
- 写 skill catalog 文档
- 提供 discover / render / generate 三个脚本
- 用当前工作区可发现的一批工具生成第一波 skill

第一版不做：

- 全互联网级工具爬取
- 每个 skill 的深度人工润色
- 复杂多语言 docs parser
- 自动发布到外部 skill registry

## 验证标准

第一版完成至少满足：

- `docs/skills/biology-skill-spec.md` 存在并可作为以后新增 skill 的标准
- `docs/skills/biology-skill-catalog.md` 能列出当前 biologically relevant skills
- `discover_bio_tools.py` 能稳定输出候选工具列表
- `generate_bio_skills.py` 能为至少一批本地工具生成 `.claude/skills/<tool>/SKILL.md`
- 生成结果目录结构符合 Claude Code skills 约定
- 对同一工具重复运行时，结果可重复更新，不会失控追加垃圾文件

## 设计结论

本轮不是“再手写几个 skill”，而是建立一条长期可复用的 **biology skill production line**。  
以后每次新增生信工具，默认都走这条线：

1. discover
2. collect help
3. render skill
4. refresh catalog

这样当前工作区里的生物学能力才会逐步从“零散知识”变成“结构化技能层”。
