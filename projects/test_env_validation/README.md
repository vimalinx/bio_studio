# test_env_validation

## 角色
这是 **工作区级环境 smoke test**，不是普通分析项目模板示例。

它的职责是验证：
- 当前 `bio_studio` 环境中的核心生信 CLI 是否可运行
- Python 关键依赖是否可导入
- 一条最小化 RNA-seq 风格链路能否在 toy data 上跑通

如果你只是新建了一个分析项目，要先跑项目级自检：
```bash
python scripts/project.py validate <项目名>
```

如果你要验整个工作区工具链，运行：
```bash
python scripts/project.py workspace-validate
```

工作区 CLI 现在会优先尝试检测本机 `bio` 环境解释器来运行这个 smoke test。
这能避免当前 shell 落在系统 Python 时，误报 `Bio` 一类依赖缺失。

## 目录说明

- `scripts/run_validation.py`: 工作区级环境 smoke test 主入口
- `scripts/validate_project.py`: 项目级结构与配置自检，方便和普通模板保持一致
- `scripts/pipeline.py`: 轻量模板入口，支持 `--steps` / `--validate`
- `legacy_toolchain_validation/`: 早期验证产物归档

如果你需要只检查这个项目本身的结构与路径定义，也可以运行：
```bash
python scripts/project.py validate test_env_validation
```

## 注意

- 这里的通过结果不等于任意项目配置都正确
- 任意项目模板通过 `validate_project.py`，也不等于整个工作区工具链已验收
- `scanpy` 当前属于可选检查项；缺失时会给出警告，但不应让当前 smoke test 直接失败
- smoke test 当前生成的 toy FASTQ 是成对 reads，`featureCounts` 会按 paired-end 输入处理
- `multiqc_report/` 在每次运行前会被重建，并从搜索路径中排除旧报告，避免重复扫描历史 `multiqc_data`
