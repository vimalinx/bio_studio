# LiteFold Launcher Design

## Goal

把 LiteFold 外层桥接从“只能看状态”补到“能做启动前预检，并能用正确上下文启动 self-hosted 服务”。

## Design

外层包装新增两类能力：

- `preflight`：基于 LiteFold 现有源码结构，检查推荐工作目录、`PYTHONPATH`、候选依赖文件和缺失 Python 模块
- `start-selfhosted`：用桥接层推导出的工作目录与环境变量启动 `selfhosted.py`，并提供 dry-run 模式

这次仍然明确不做：

- 不修改 `repositories/active/litefold/source` 内任何文件
- 不自动安装依赖
- 不把 Docker 包装成默认可用路径

## Why

当前上游 `selfhosted.py` 同时依赖三类导入路径：

- `litefold.selfhosted.*`
- `fold_models`
- `selfhosted.*`

因此直接 `cd selfhosted && python selfhosted.py` 并不稳。外层包装需要把运行目录与 `PYTHONPATH` 固定下来，否则这个入口只是“看起来有脚本”，不等于“能稳定起服务”。

## Interface

在现有 LiteFold 桥接基础上补这些输出：

- 推荐启动工作目录
- 推荐 `PYTHONPATH`
- 候选 requirements 文件
- 从 self-hosted 相关源码中提取的外部模块根名
- 当前 Python 环境中缺失的模块列表

CLI 增加：

- `python scripts/litefold.py preflight`
- `python scripts/litefold.py start-selfhosted --dry-run`

## Validation Strategy

- 用单元测试验证预检结果、推荐命令和 dry-run CLI 输出
- 不启动真实 LiteFold 服务
- 不依赖真实 GPU、Docker 或联网状态
