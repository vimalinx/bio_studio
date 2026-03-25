# LiteFold Env Assets Design

## Goal

给 LiteFold self-hosted 增加一套主仓库维护的环境补全资产，让外层桥接不只知道“缺包”，还知道“该装什么、用什么脚本装”。

## Design

这次补三样东西：

- 仓库自带的 LiteFold self-hosted requirements 文件
- 仓库自带的 LiteFold setup 脚本
- 更完整的 `preflight`，把 `selfhosted` 和 `fold_models` 两边的外部 Python 依赖都纳入扫描

明确不做：

- 不修改 `repositories/active/litefold/source`
- 不让 `scripts/litefold.py` 默认自动改环境
- 不在测试里真的联网安装依赖

## Interface

新增资产：

- `requirements-litefold-selfhosted.txt`
- `scripts/setup/setup_litefold_env.sh`

桥接层输出新增：

- 工作区 requirements 文件路径
- 工作区 setup 脚本路径
- 推荐安装命令
- 来自 `fold_models` 的额外依赖模块

## Validation Strategy

- 用 pytest 验证 `preflight` 是否暴露安装资产与更完整依赖列表
- 用 dry-run 测试 shell 脚本输出，不做真实安装
- 继续保持稳定测试链纯本地、无网络
