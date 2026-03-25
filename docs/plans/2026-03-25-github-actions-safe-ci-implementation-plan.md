# GitHub Actions Safe CI Implementation Plan

## Scope

- In: GitHub Actions workflow、安全约束、最小依赖、回归测试
- Out: 自动部署、自动发布、自动同步子模块、重型生信工具链云端执行

## Files

- Create: `.github/workflows/ci.yml`
- Create: `requirements-ci.txt`
- Create: `scripts/ci/run_stable_tests.sh`
- Create: `tests/test_github_actions_ci.py`
- Update: `docs/CHANGELOG.md`

## Steps

1. 添加只读 `CI` workflow，并限制触发面
2. 抽出最小 CI 依赖，避免安装整套本地研究环境
3. 固定稳定测试集合，避免 Actions 跑到高成本或高波动任务
4. 添加 workflow 安全回归测试，防止后续误改出 `schedule` / `workflow_run` / 自动 push
5. 运行本地 pytest 验证新增测试与稳定测试集合
6. 在隔离虚拟环境中验证 `requirements-ci.txt` 足以支撑 CI
