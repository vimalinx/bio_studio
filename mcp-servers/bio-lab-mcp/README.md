# Bio Lab MCP Server

工作区实验编排 MCP server，主要面向 AI agent 调用当前仓库已经存在的项目发现、自检和环境验证入口。

## 当前能力

1. `list_workspace_projects`
   - 枚举 `projects/` 下带统一 `scripts/` 入口的项目
   - 返回 `config.py` / `validate_project.py` / `pipeline.py` 是否存在

2. `get_project_summary`
   - 读取项目配置
   - 返回项目类型、描述、线程数、样本数、所需 CLI 工具
   - 如果已有 `logs/validation_report.json`，会带上最近一次验证摘要

3. `get_project_steps`
   - 调用 `scripts/project.py steps <project_name>`
   - 返回项目 pipeline 当前暴露的步骤说明

4. `run_project_validation`
   - 调用 `scripts/project.py validate <project_name>`
   - 运行项目级自检并返回 stdout / stderr / exit code

5. `run_workspace_validation`
   - 调用 `scripts/project.py workspace-validate`
   - 运行工作区级 smoke test
   - 这个流程会真的执行 toy 数据、比对和 QC，耗时相对更长

6. `preview_workspace_plan`
   - 复用在线控制面的同一套 planner
   - 把自然语言请求拆成 task brief、能力匹配和计划步骤
   - 适合在真正执行前先看系统准备怎么路由

7. `get_server_capabilities`
   - 返回当前 server 的工作区路径、统一入口路径和工具列表

## 对应代码入口

- `scripts/project.py`
- `projects/test_env_validation/scripts/run_validation.py`
- `lib/project_validation.py`
- `lib/workspace_env.py`

## 最小配置

```bash
python mcp-servers/render_claude_config.py --write mcp-servers/claude-config.json
```

```json
{
  "mcpServers": {
    "bio-lab": {
      "command": "/home/vimalinx/miniforge3/envs/bio/bin/python",
      "args": ["/home/vimalinx/Projects/bio_studio/mcp-servers/bio-lab-mcp/lab_server.py"]
    }
  }
}
```

## 示例

### 查看工作区项目

```
用户: 这个工作区有哪些项目已经接上统一入口了？

AI: [调用 list_workspace_projects]
    - ai_design_playground
    - example_rnaseq
    - test_env_validation
    - test_rnaseq_analysis
    - yeast_genome_learning
    - yeast_rnaseq_demo
```

### 跑项目级自检

```
用户: 帮我验证 ai_design_playground

AI: [调用 run_project_validation]
    - returncode: 0
    - stdout: Validation summary: PASS
```

### 跑工作区级 smoke test

```
用户: 跑一下整个工作区的环境验证

AI: [调用 run_workspace_validation]
    - 会触发 scripts/project.py workspace-validate
```

### 预览任务会怎么被拆解

```
用户: 先别跑，预览一下 “analyze SARS-CoV-2 literature and sequence evidence” 会怎么被系统拆解

AI: [调用 preview_workspace_plan]
    - brief.task_type: analysis
    - selected_capabilities:
      - database.search.pubmed
      - sequence.analysis.basic
    - steps:
      - Create or select a workspace project shell
      - Search literature and biomedical databases
      - Basic sequence analysis
```

## 边界

- 这个 server 目前不是仪器控制层，也不是 Opentrons/Benchling 之类的外部 lab 平台接口。
- 它封装的是当前仓库已经存在的项目和验证链路，重点是让 AI 能直接发现和调用这些入口。
