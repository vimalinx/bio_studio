# Bio Studio

> 面向 AI agent 的生物计算工作台，不是单一分析脚本仓库，也不是单一 Web 应用。

Bio Studio 把几个层面放在同一个工作区里：

- 项目工作区：真实分析任务、数据、日志和报告都优先落在 `projects/`
- 本地骨架：`lib/`、`scripts/`、项目模板和统一 CLI
- AI 接口：`mcp-servers/`、skills、介绍页和协议文档
- 外部引擎：Evo2、Biomni、RFdiffusion、LiteFold 等第三方能力

如果你是第一次进这个仓库，先看这三个入口：

- [docs/README.md](docs/README.md)
- [docs/WORKSPACE_ARCHITECTURE.md](docs/WORKSPACE_ARCHITECTURE.md)
- [docs/site/index.html](docs/site/index.html)

## 这套工作区现在能做什么

- 用统一入口创建、验证和运行项目：`python scripts/project.py ...`
- 跑工作区级环境 smoke test：`python scripts/project.py workspace-validate`
- 复用共享 runtime 和模块执行常见流程，而不是在每个项目里重复写命令
- 提供 MCP servers 给 AI 直接调用：
  - `bio-design-mcp`：DNA / mRNA 设计、靶点评估
  - `bio-lab-mcp`：项目发现、项目级验证、工作区级验证
  - `bio-sequence-mcp`：序列分析
  - `bio-structure-mcp`：结构分析
  - `bio-database-mcp`：NCBI / UniProt / PubMed 查询
- 保留外部研究引擎源码，但把“如何调用它们”尽量收敛到工作区骨架层

## 快速开始

先进入仓库根目录并激活 `bio` 环境：

```bash
cd /path/to/bio_studio
source <conda_install>/bin/activate bio
```

然后按这个顺序走：

```bash
# 查看环境和目录边界
python scripts/maintenance/generate_env_report.py

# 创建新项目
python scripts/project.py create my_analysis --type rnaseq --description "My analysis"

# 先自检，再看步骤
python scripts/project.py validate my_analysis
python scripts/project.py steps my_analysis

# 再运行
python scripts/project.py run my_analysis

# 需要时做工作区级 smoke test
python scripts/project.py workspace-validate
```

如果你想先看一个已经接好的项目，而不是从空项目开始，优先看：

- [projects/ai_design_playground/README.md](projects/ai_design_playground/README.md)
- [projects/test_env_validation/README.md](projects/test_env_validation/README.md)
- [projects/test_rnaseq_analysis/README.md](projects/test_rnaseq_analysis/README.md)

## 统一入口怎么理解

工作区主入口是 [scripts/project.py](scripts/project.py)。

常用命令：

```bash
python scripts/project.py create <项目名> --type rnaseq
python scripts/project.py validate <项目名>
python scripts/project.py steps <项目名>
python scripts/project.py run <项目名>
python scripts/project.py workspace-validate
```

`run` 的语义不是所有项目都完全一样：

- 模板项目和 demo 项目：通常会直接执行项目 `pipeline.py`
- 专项项目：会转发到项目自己的定制 pipeline
- 学习项目：可能只输出下一步教学命令，而不是替你一键跑完整教程

这不是不一致，而是刻意保留项目边界。详细解释见 [docs/WORKSPACE_ARCHITECTURE.md](docs/WORKSPACE_ARCHITECTURE.md)。

## MCP 和 AI 接口

当前 MCP 配置在 [mcp-servers/](mcp-servers/)。

如果你要把这套 MCP servers 接到 Claude Code 或其他 agent 运行时，直接生成配置：

```bash
cd mcp-servers
python render_claude_config.py --write claude-config.json
cat claude-config.json
```

现在的活跃 MCP servers：

- [mcp-servers/bio-design-mcp/README.md](mcp-servers/bio-design-mcp/README.md)
- [mcp-servers/bio-lab-mcp/README.md](mcp-servers/bio-lab-mcp/README.md)
- [mcp-servers/bio-sequence-mcp/README.md](mcp-servers/bio-sequence-mcp/README.md)
- [mcp-servers/bio-structure-mcp/README.md](mcp-servers/bio-structure-mcp/README.md)
- [mcp-servers/bio-database-mcp/README.md](mcp-servers/bio-database-mcp/README.md)

## 在线 API 控制面

如果你准备把这套工作区升级成在线服务，不建议重写现有骨架，而应该在现有 CLI / MCP / `projects/` 物化模型之上加一层控制面。

已经补好的设计与实施入口：

- [docs/superpowers/specs/2026-03-29-online-api-control-plane-design.md](docs/superpowers/specs/2026-03-29-online-api-control-plane-design.md)
- [docs/superpowers/plans/2026-03-29-online-api-control-plane-phase-1.md](docs/superpowers/plans/2026-03-29-online-api-control-plane-phase-1.md)

这条路线的核心原则是：

- 保留 `scripts/project.py` 作为执行脊柱
- 保留 `mcp-servers/` 作为工具接口层
- 保留 `projects/` 作为可复现实验与产物边界
- 新增 HTTP/API 层负责请求拆解、能力匹配、运行状态和工件索引

当前已经落下一个最小可用入口 [scripts/api_server.py](scripts/api_server.py)，本地可直接启动：

```bash
python scripts/api_server.py
```

第一批可用端点：

- `GET /healthz`
- `GET /v1/capabilities`
- `POST /v1/plans/preview`
- `POST /v1/runs`
- `GET /v1/runs/{run_id}`
- `GET /v1/runs/{run_id}/events`
- `GET /v1/runs/{run_id}/artifacts`

## 安全 CI 和发布护栏

这个仓库现在已经接上最小、只读的 GitHub Actions 稳定测试链：

- workflow 只跑 `Stable Tests`
- 权限是只读 `contents: read`
- 不做自动 commit、自动 push、自动开 PR
- 不使用 `schedule`、`workflow_run`、`pull_request_target`

本地想复用同一条稳定链，直接运行：

```bash
bash scripts/ci/run_stable_tests.sh
```

这条链是“仓库治理护栏”，不是重型生信计算入口。它的目标是稳，不是把整套研究环境搬进 Actions。

## LiteFold 桥接入口

LiteFold 当前通过工作区桥接脚本接入，而不是直接改写 vendored 上游源码。

可以先看工作区识别到的状态：

```bash
python scripts/litefold.py status
python scripts/litefold.py status --json
```

如果你要看它现在到底能不能被当前 Python 环境启动，先跑预检：

```bash
python scripts/litefold.py preflight
python scripts/litefold.py preflight --json
```

预检会把两类东西一起告诉你：

- 当前 Python 环境还缺哪些 LiteFold 依赖
- 主仓库维护的 LiteFold 安装资产在哪里

如果缺包，先看 dry-run 安装计划：

```bash
bash scripts/setup/setup_litefold_env.sh --dry-run
```

真正要补全环境时，再执行：

```bash
bash scripts/setup/setup_litefold_env.sh
```

这套安装脚本默认使用仓库里的 [requirements-litefold-selfhosted.txt](requirements-litefold-selfhosted.txt)。

如果已经有 LiteFold 服务在跑，可以直接探活：

```bash
python scripts/litefold.py health --url http://localhost:7114
```

如果你只是想知道当前工作区建议怎么从 vendored LiteFold 启动 self-hosted 入口：

```bash
python scripts/litefold.py print-selfhosted-command
```

如果你想先看桥接层准备怎么起服务，而不真的执行：

```bash
python scripts/litefold.py start-selfhosted --dry-run
python scripts/litefold.py start-selfhosted --dry-run --json
```

这层桥接目前负责发现、预检、体检和启动包装，仍然遵守一个原则：

- `repositories/active/litefold/source` 是上游引擎工作树
- 工作区优先在外层做封装和治理，不默认修改上游源码

## 当前工作区结构

```text
bio_studio/
├── docs/                    # 协议、架构、环境与站点入口
├── lib/                     # 共享模块、模板 runtime、项目验证辅助
├── mcp-servers/             # MCP server 接口层
├── projects/                # 真实项目工作区
├── repositories/active/     # 第三方引擎源码
├── scripts/                 # 统一 CLI 和维护脚本
├── shared_data/             # 共享参考数据
└── tools/                   # 独立工具脚本
```

几个最重要的约束：

- 真实分析结果优先放 `projects/`，不要把根目录当临时实验场
- `repositories/active/` 是外部能力层，不是主产品代码层
- `docs/` 放工作区知识，项目专属说明放项目自己的 `README.md`
- 想让 AI 可靠调用能力，优先加到 `lib/` 或 `mcp-servers/`，不要只留在一次性脚本里

## 已经沉淀下来的能力面

共享能力不是只在 README 里写写，而是已经有代码和测试支撑：

- 项目模板与统一验证入口
- 工作区 Python / PATH 环境解析
- 项目级 CLI 工具需求校验
- demo 项目与特殊项目的统一入口兼容
- `bio-design-mcp` / `bio-lab-mcp` 的最小可用服务骨架

变更细节可以直接看 [docs/CHANGELOG.md](docs/CHANGELOG.md)。

## 推荐阅读顺序

1. [README.md](README.md)
2. [docs/README.md](docs/README.md)
3. [docs/WORKSPACE_ARCHITECTURE.md](docs/WORKSPACE_ARCHITECTURE.md)
4. [projects/test_env_validation/README.md](projects/test_env_validation/README.md)
5. [projects/ai_design_playground/README.md](projects/ai_design_playground/README.md)

## 注意事项

- 这个仓库包含很多“工作区资产”，不只是源代码，还包括项目脚本、测试、文档和第三方引擎入口
- 一些目录可能存在本地工作树修改，这是研究型工作区的常态，不代表都适合直接公开发布
- 如果要对外发布，先确认哪些结果、日志、站点文件和第三方子模块状态应该一起带出去
