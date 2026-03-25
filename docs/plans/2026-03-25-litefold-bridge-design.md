# LiteFold Bridge Design

## Goal

为工作区增加一个最小、可测试、不会改动上游 LiteFold 子模块源码的桥接层。

## Design

这层桥接只负责三件事：

- 发现 `repositories/active/litefold/source` 下的关键文件与目录
- 提供一个轻量健康检查入口，必要时探测远端 `/health`
- 生成建议的 self-hosted 启动命令，帮助用户复用 vendored LiteFold 源码

明确不做：

- 修改 `repositories/active/litefold/source` 内任何源码
- 在仓库测试里真的启动 Docker、下载模型或访问 GPU
- 把 LiteFold 深度整合进统一项目运行入口

## Interface

新增两个入口：

- `lib/litefold_bridge.py`
  - 负责路径发现、状态汇总、健康检查、推荐命令生成
- `scripts/litefold.py`
  - 提供 `status`、`health`、`print-selfhosted-command` 三个命令

## Validation Strategy

- 用纯 Python 单元测试覆盖路径发现、命令生成、健康检查和 CLI 输出
- 用 monkeypatch 隔离 Docker 与 HTTP 依赖
- 把新测试加入 `scripts/ci/run_stable_tests.sh`，保持 CI 可稳定运行

## README Update

根 README 额外补两块信息：

- 现在的 GitHub Actions 是只读、安全、不自触发的稳定测试链
- LiteFold 目前通过桥接脚本暴露最小治理入口，而不是直接改写上游引擎
