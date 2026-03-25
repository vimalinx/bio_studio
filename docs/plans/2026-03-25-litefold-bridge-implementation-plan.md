# LiteFold Bridge Implementation Plan

## Scope

- In: LiteFold 桥接模块、CLI、README 补充、稳定测试
- Out: LiteFold 模型下载、Docker 编排、GPU 推理、上游源码改造

## Files

- Create: `lib/litefold_bridge.py`
- Create: `scripts/litefold.py`
- Create: `tests/test_litefold_bridge.py`
- Update: `scripts/ci/run_stable_tests.sh`
- Update: `README.md`

## Steps

1. 先写 LiteFold 桥接层测试，覆盖路径发现、命令生成和健康检查
2. 运行新增测试，确认在实现前按预期失败
3. 实现 `lib/litefold_bridge.py` 的最小行为
4. 实现 `scripts/litefold.py` CLI，并补一个轻量 CLI 测试
5. 更新 README，写明安全 CI 与 LiteFold 桥接入口
6. 把 LiteFold 测试加入稳定测试链路
7. 运行新增测试和稳定测试集合，确认全部通过
