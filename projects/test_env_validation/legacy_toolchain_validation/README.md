# legacy_toolchain_validation

这个目录是从工作区根目录的旧版 `test_env_validation/` 迁移过来的历史验证产物与脚本备份，保留仅用于回溯与对照。

推荐使用当前标准化项目目录：`projects/test_env_validation/`。

说明：
- `data/`：旧版 toy 数据与参考文件
- `results/`：旧版工具链跑出的示例输出（可用于快速比对）
- `scripts/run_validation.py`：旧版验证脚本（已改成不依赖运行时 cwd）
