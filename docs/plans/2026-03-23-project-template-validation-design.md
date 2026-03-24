# Project Template Validation Design

> 第二阶段设计记录：把项目模板从空骨架升级为带项目级自检能力的最小工作单元。

## 目标

- 新建项目后可立即运行项目级验证，而不是先手工补脚手架
- 明确区分“工作区级环境验收”和“项目级模板自检”
- 让模板、样例和协议保持一致

## 范围

本轮聚焦：
- `lib/create_project.py`
- 模板生成的 `README.md` / `scripts/config.py` / `scripts/pipeline.py`
- 新增 `scripts/validate_project.py`
- 同步 `projects/example_rnaseq`
- 同步 `docs/AI_ANALYSIS_PROTOCOL.md`

不包含：
- 重型生信流程实现
- MCP 扩展
- 第三方仓库接入改造

## 设计决策

### 1. 项目模板默认自包含

模板生成的脚本先以项目自身为中心运行，不强依赖向上回溯导入工作区根目录。只有需要复用 `lib/` 时，再显式引入。

### 2. 建立项目级验证脚本

新增 `scripts/validate_project.py`，用于检查：
- 目录结构完整性
- 配置模块可导入
- 配置字段是否齐全
- 路径定义是否落在项目内
- 是否生成验证日志与 JSON 摘要

### 3. 统一入口

`scripts/pipeline.py` 提供：
- `--steps`
- `--validate`
- 正常步骤执行

项目级验证通过 `pipeline.py --validate` 和 `validate_project.py` 两条入口都能跑通。

### 4. 保留双层验证模型

- `projects/test_env_validation`：工作区级环境 smoke test
- 模板内置 `validate_project.py`：项目级自检

二者不合并，避免把环境验收和项目模板绑死。
