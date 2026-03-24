# Bio Studio Workspace Governance Design

> 本文档记录 2026-03-23 这一轮“架构治理优先”的第一阶段设计，用于后续继续清理与扩展工作区。

## 目标

把 Bio Studio 从“能用但需要口头解释”推进到“目录边界清楚、入口统一、AI 与人都能快速定位”的状态。

## 第一阶段范围

本轮只处理治理与说明层，不改分析逻辑，不触碰第三方仓库内部实现。

包括：
- 把 `docs/site/index.html` 升级为正式入口
- 新增工作区架构文档
- 给几个语义不清目录补 README 或定位说明
- 对齐 `.gitmodules` 与当前 gitlink 记录
- 更新主 README、文档索引和变更日志

## 不在本轮范围

- 实现 `bio-design-mcp` / `bio-lab-mcp`
- 重写 `projects/` 内 pipeline
- 清理第三方仓库的 dirty working tree
- 删除历史项目、大数据或旧结果

## 设计判断

### 1. 文档优先
这个仓库的核心问题不是“功能缺少一个脚本”，而是工作区边界没有完全显式化。  
第一阶段优先修正导航、定位和分层，而不是立刻继续加功能。

### 2. 不破坏现有工作树
`repositories/active/` 下存在真实第三方仓库和本地改动。  
本轮只做元数据补齐，不改仓库内容，不 reset，不迁移。

### 3. 让入口可见
`docs/site/index.html` 已经具备展示价值，但如果 README 和 docs index 不指向它，它仍然只是孤文件。  
因此入口接线是第一优先级。

### 4. 保留位先命名，再决定删不删
`tmp/`、`evo2-7b-hf/`、`_legacy_workspace_v1/` 这类目录先通过 README 明确身份。  
是否进一步清理，留到下一阶段。

## 下一阶段候选

第一阶段完成后，第二阶段可从以下三个方向中任选一个推进：

1. AI 接口增强
- 补 `bio-design-mcp`
- 补 `bio-lab-mcp`
- 整理 `.claude/skills`

2. 工作流标准化
- 统一项目模板
- 规范报告产出
- 强化验证项目与脚本接口

3. 外部引擎治理
- 明确 fork / vendored / submodule 策略
- 为 active repositories 增加统一说明与维护规则
