# GitHub Actions Safe CI Design

## Goal

为 `bio_studio` 增加一套“默认安全、不自触发、不写回仓库”的 GitHub Actions。

## Design

这套 CI 只承担两件事：
- 在 `push` / `pull_request` 上运行一组稳定、纯 Python 的回归测试
- 在 `workflow_dispatch` 下允许人工手动触发同一条检查链

明确不做：
- 自动提交
- 自动 push
- 自动改版本号
- 自动创建 PR
- `schedule` 定时任务
- `workflow_run` 链式触发
- `pull_request_target`

## Stability Rules

1. Workflow 默认只读权限：`contents: read`
2. `actions/checkout` 禁止保留写凭证：`persist-credentials: false`
3. 使用 `concurrency.cancel-in-progress: true`，同分支重复提交时自动取消旧任务
4. 只监听代码相关路径，不因普通文档改动反复运行
5. 安装最小 CI 依赖，而不是整套研究环境

## Validation Strategy

- 新增 `requirements-ci.txt`，只包含 CI 必需依赖
- 新增 `scripts/ci/run_stable_tests.sh`，固定稳定测试集合
- 新增 `tests/test_github_actions_ci.py`，防回归检查 workflow 是否重新变成危险形态
