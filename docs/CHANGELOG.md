# Bio Studio 变更日志

记录工具部署、分析流程和环境变更。每次操作完成后追加一条记录。

## [2026-03-25] - 主仓整洁度：忽略外部引擎子模块内部脏状态
- **操作**: 更新 `.gitmodules`，为 `repositories/active/RFdiffusion` 与 `repositories/active/litefold/source` 增加 `ignore = dirty`；同步更新 `repositories/README.md` 说明主仓与外部引擎仓库的边界
- **结果**: 成功
- **影响**: 主仓 `git status` 不再因外部引擎子模块内部的本地未提交改动而持续显示脏状态；真正的子模块指针变化仍会保留为主仓可见变更
- **验证**: `git status --short --ignore-submodules=dirty`

## [2026-03-24] - 统一 CLI：学习项目 run 返回工作区级引导
- **操作**: 更新 `scripts/project.py`，在 `run` 子命令下向项目 pipeline 透传工作区 CLI 上下文；升级 `projects/yeast_genome_learning/scripts/pipeline.py`，在经由统一入口调用时输出工作区级 `steps` / `validate` 后续命令；补充 `tests/test_workspace_project_cli.py` 对专项项目 `validate` / `steps` 与学习项目 `run` 引导的覆盖；同步更新 `README.md`、`docs/README.md`、`docs/AI_ANALYSIS_PROTOCOL.md` 与 `docs/WORKSPACE_ARCHITECTURE.md`
- **结果**: 成功
- **影响**: `python scripts/project.py run yeast_genome_learning` 不再只提示项目内局部命令，而是明确给出工作区级下一步操作；统一入口对特殊项目的行为语义更清晰
- **验证**: `python -m pytest tests/test_workspace_project_cli.py tests/test_special_project_minimal_integration.py -q`

## [2026-03-24] - MCP 扩展：落地设计/实验编排服务与配置渲染
- **操作**: 新增 `mcp-servers/bio-design-mcp/` 与 `mcp-servers/bio-lab-mcp/` 的服务实现；新增 `mcp-servers/render_claude_config.py` 统一渲染 `claude-config.json`；更新 `mcp-servers/install-all.sh`、`mcp-servers/mcp-requirements.txt` 与数据库/序列/结构 server 的解释器与配置行为；同步更新 `mcp-servers/README.md`、`docs/WORKSPACE_ARCHITECTURE.md` 与相关 README
- **结果**: 成功
- **影响**: 工作区现在不只有序列/结构/数据库查询 MCP，还具备设计入口和项目发现/验证入口；Claude 配置不再依赖手写绝对路径
- **验证**: `python -m pytest tests/test_design_mcp_server.py tests/test_lab_mcp_server.py tests/test_mcp_config_render.py tests/test_mcp_server_entrypoints.py tests/test_mcp_readme_paths.py tests/test_database_mcp_server_config.py tests/test_design_tool_script_compat.py -q`

## [2026-03-24] - 工作区入口：修复 workspace-validate 解释器漂移
- **操作**: 更新 `scripts/project.py`，让项目入口优先解析本机 `bio` 环境 Python；更新 `projects/test_env_validation/scripts/run_validation.py`，将 `scanpy` 调整为可选依赖并在缺少必需库时输出明确解释器提示；新增 `tests/test_workspace_runtime_consistency.py`
- **结果**: 成功
- **影响**: `python scripts/project.py workspace-validate` 不再因为当前 shell 落在系统 Python 而误报 `Bio` 缺失；当前环境未安装 `scanpy` 时，smoke test 也能按预期继续执行
- **验证**: `python -m pytest tests/test_project_template_validation.py tests/test_workspace_project_cli.py tests/test_workspace_validation_project.py tests/test_project_type_templates.py tests/test_example_project_sync.py tests/test_workspace_runtime_consistency.py -q`；`python scripts/project.py workspace-validate`

## [2026-03-24] - 工作区验证：修复 toy paired-end 与 MultiQC 重扫噪音
- **操作**: 更新 `projects/test_env_validation/scripts/run_validation.py`，新增 paired-end mock read 生成辅助函数，使 R2 使用 fragment 末端反向互补；将 `featureCounts` 调整为 `-p`；在执行 `multiqc` 前重建 `multiqc_report/` 并通过 `--ignore '*/multiqc_report/*'` 排除旧报告；新增 `tests/test_workspace_validation_mock_data.py`
- **结果**: 成功
- **影响**: `workspace-validate` 现在会产生可正常双端比对、可计数、可汇总的 toy 数据，MultiQC 不再因旧 `multiqc_data` 被重复扫描而输出 featureCounts / fastp 噪音
- **验证**: `python -m pytest tests/test_workspace_validation_mock_data.py -q`；`python scripts/project.py workspace-validate`

## [2026-03-24] - Demo 项目：对齐模板级验证入口
- **操作**: 为 `projects/ai_design_playground/` 与 `projects/yeast_rnaseq_demo/` 新增 `scripts/validate_project.py`；升级两者的 `scripts/config.py` 与 `scripts/pipeline.py`，统一支持 `--validate` / `--steps`；同步更新对应 README；新增 `tests/test_demo_project_validation_entrypoints.py`
- **结果**: 成功
- **影响**: 两个已有 demo 项目现在也具备和模板项目一致的项目级自检入口，能输出 `logs/validation_report.json` 与 `logs/validation_report.md`，同时保留各自原有的自定义分析流程
- **验证**: `python -m pytest tests/test_demo_project_validation_entrypoints.py -q`；`python -m pytest tests/test_project_template_validation.py tests/test_workspace_project_cli.py tests/test_workspace_validation_project.py tests/test_project_type_templates.py tests/test_example_project_sync.py tests/test_workspace_runtime_consistency.py tests/test_workspace_validation_mock_data.py tests/test_demo_project_validation_entrypoints.py -q`

## [2026-03-24] - 模板运行时：把类型模板接到共享 modules
- **操作**: 修复 `lib/modules/utils.py` 的 shell 命令处理；为 `lib/modules/qc.py` 增加 `fastp` wrapper；新增 `lib/modules/rnaseq.py` 与 `lib/template_runtime.py`；升级 `lib/create_project.py` 生成的模板 `pipeline.py`，在工作区内优先通过共享 runtime 调用 `lib/modules`；同步 `projects/example_rnaseq/scripts/pipeline.py`；新增 `tests/test_project_template_module_integration.py`
- **结果**: 成功
- **影响**: `rnaseq` / `variant` / `phylogeny` 模板不再只是打印步骤提示，而是在具备工作区上下文与输入前提时可直接复用共享 wrapper；项目脱离工作区时仍保留原有轻量模板行为
- **验证**: `python -m pytest tests/test_project_template_module_integration.py -q`

## [2026-03-24] - Demo 运行时：让 yeast RNA-seq demo 复用共享模块
- **操作**: 扩展 `lib/template_runtime.py`，为 RNA-seq 增加 `build_index` / `alignment` / `quantification` 共享步骤；更新 `projects/yeast_rnaseq_demo/scripts/config.py` 暴露 STAR / featureCounts 参数；升级 `projects/yeast_rnaseq_demo/scripts/pipeline.py` 优先走共享 runtime，同时保留 `counts_matrix.csv` 后处理；新增 `tests/test_demo_project_shared_runtime.py`
- **结果**: 成功
- **影响**: `yeast_rnaseq_demo` 不再手写完整的 STAR / featureCounts 调用链作为唯一实现，而是优先复用共享 RNA-seq wrapper；离开工作区时仍可回退到项目内置命令路径
- **验证**: `python -m pytest tests/test_demo_project_shared_runtime.py -q`

## [2026-03-24] - Demo 运行时：让 AI design playground 复用共享模块
- **操作**: 新增 `lib/modules/sequence_analysis.py`，提供 `seqkit stats` 与 `RNAfold` wrapper；扩展 `lib/template_runtime.py`，增加 `run_ai_design_playground_analysis(...)`；升级 `projects/ai_design_playground/scripts/pipeline.py` 优先通过共享 runtime 执行 `seqkit` / `Prodigal` / `RNAfold`，同时保留 toy 数据生成、蛋白摘要与报告拼装逻辑；新增 `tests/test_ai_design_playground_shared_runtime.py`
- **结果**: 成功
- **影响**: `ai_design_playground` 的底层工具调用不再只存在于项目脚本内部，dry-lab 分析链开始沉淀为可复用模块；离开工作区时仍保留原有本地实现作为回退
- **验证**: `python -m pytest tests/test_ai_design_playground_shared_runtime.py -q`

## [2026-03-24] - 特殊项目：补齐最小统一入口并明确分类
- **操作**: 为 `projects/test_rnaseq_analysis/` 与 `projects/yeast_genome_learning/` 新增轻量 `scripts/validate_project.py`；为 `yeast_genome_learning` 新增兼容型 `scripts/pipeline.py`；为 `test_rnaseq_analysis/scripts/pipeline.py` 补充 `--validate` / `--steps`；更新两个项目 README 与 `QUICKSTART.md` 明确“专项项目”/“学习项目”分类；在 `docs/WORKSPACE_ARCHITECTURE.md` 记录特殊项目边界；新增 `tests/test_special_project_minimal_integration.py`
- **结果**: 成功
- **影响**: 这两个项目现在可以被工作区统一入口最小识别和自检，但仍保留各自原有的研究型/教学型执行方式，不被强行改造成模板 runtime 项目
- **验证**: `python -m pytest tests/test_special_project_minimal_integration.py -q`

## [2026-03-23] - 项目模板：引入类型化骨架并锁定示例同步
- **操作**: 扩展 `lib/create_project.py`，为 `rnaseq` / `variant` / `phylogeny` / `genome` 注入类型化 README、配置字段和步骤提示；新增 `tests/test_project_type_templates.py`；同步 `projects/example_rnaseq/` 到新版 RNA-seq skeleton；新增 `tests/test_example_project_sync.py` 防止示例模板漂移
- **结果**: 成功
- **影响**: 新建项目不再只是通用空壳，而是带有最小可运行的分析类型语义；示例项目与生成器保持一致，后续模板升级更容易被测试捕获
- **验证**: `python -m pytest tests/test_project_template_validation.py tests/test_workspace_project_cli.py tests/test_workspace_validation_project.py tests/test_project_type_templates.py tests/test_example_project_sync.py -q`

## [2026-03-23] - 工作区验证项目：对齐新版模板验证接口
- **操作**: 升级 `projects/test_env_validation/scripts/config.py` 为 `Path` 风格配置；新增 `projects/test_env_validation/scripts/validate_project.py`；更新 `projects/test_env_validation/scripts/pipeline.py` 支持 `--validate`；新增 `tests/test_workspace_validation_project.py`
- **结果**: 成功
- **影响**: 工作区环境 smoke test 项目自身也具备了和普通模板一致的项目级自检入口，但真实工具链验收仍以 `run_validation.py` / `python scripts/project.py workspace-validate` 为主
- **验证**: `python -m pytest tests/test_workspace_validation_project.py -q`

## [2026-03-23] - 工作区入口：新增统一项目 CLI
- **操作**: 新增 `scripts/project.py`，提供工作区级 `create` / `validate` / `steps` / `run` / `workspace-validate` 子命令；新增 `tests/test_workspace_project_cli.py`；更新 `README.md`、`docs/AI_ANALYSIS_PROTOCOL.md` 与 `projects/test_env_validation/README.md`
- **结果**: 成功
- **影响**: 常规项目创建、自检和运行不再要求手工 `cd` 到项目脚本目录；工作区级 smoke test 有了统一入口，且与项目模板级验证边界更清晰
- **验证**: `python -m pytest tests/test_workspace_project_cli.py -q`

## [2026-03-23] - 项目模板：补齐项目级自检与统一验证入口
- **操作**: 升级 `lib/create_project.py`，新建项目时同时生成 `scripts/config.py`、`scripts/validate_project.py` 与支持 `--validate` / `--steps` 的 `scripts/pipeline.py`；同步更新 `projects/example_rnaseq/` 模板示例；更新 `docs/AI_ANALYSIS_PROTOCOL.md`
- **结果**: 成功
- **影响**: 新创建项目可在无真实数据前提下完成模板级自检，并输出 `logs/validation_report.json` 与 `logs/validation_report.md`；工作区级环境验收仍保留在 `projects/test_env_validation`
- **验证**: `python -m pytest tests/test_project_template_validation.py -q`

## [2026-03-23] - 架构治理：补齐工作区导航与边界说明
- **操作**: 新增 `docs/WORKSPACE_ARCHITECTURE.md` 说明核心工作面、第三方引擎层、共享资源与遗留区；新增 `docs/site/README.md`；为 `repositories/`、`_legacy_workspace_v1/`、`tmp/`、`evo2-7b-hf/` 补充目录说明；更新 `README.md` 与 `docs/README.md` 接入介绍页和架构文档；新增 `docs/plans/2026-03-23-workspace-governance-design.md`
- **结果**: 成功
- **影响**: 工作区入口更统一，目录边界更明确；未改动分析脚本逻辑与项目结果
- **验证**: `docs/site/index.html` 可直接打开；文档入口已在 `README.md` 与 `docs/README.md` 可见

## [2026-03-23] - Git 元数据：对齐 litefold 子模块映射
- **操作**: 更新 `.gitmodules`，补登记 `repositories/active/litefold/source` 的上游地址 `https://github.com/Anindyadeep/litefold.git`
- **结果**: 成功
- **影响**: 主仓库的 gitlink 记录与 `.gitmodules` 更一致；未修改 `litefold/source` 仓库内部内容
- **验证**: `.gitmodules` 已同时登记 `RFdiffusion` 与 `litefold/source` 两个路径

## [2026-03-23] - 结构标识：补齐占位 MCP 与 LiteFold 外层说明
- **操作**: 为 `mcp-servers/bio-design-mcp/`、`mcp-servers/bio-lab-mcp/` 新增 README；为 `repositories/active/litefold/` 新增外层说明；同步更新 `mcp-servers/README.md` 与 `docs/WORKSPACE_ARCHITECTURE.md`
- **结果**: 成功
- **影响**: 已落地接口与规划槽位边界更清晰；LiteFold 外层包装目录职责更明确
- **验证**: 相关目录均可直接通过 README 理解当前用途与状态

## [2026-02-28] - 环境：创建 conda `bio`
- **操作**: 安装 Miniforge 到 `~/miniforge3`；新增镜像配置 `tools/configs/condarc_cn`；更新 `setup_env.sh` 优先使用 mamba 且将 PyTorch 固定为 CPU build；修正 `requirements.txt` 中不可用的 PyPI 依赖并移除非 PyPI 的 `seqtk`；安装 pip 依赖；生成环境快照
- **结果**: 成功
- **影响**: 新增 conda base 与 `bio` 环境（Python 3.10）；`docs/ENVIRONMENT.md` 更新；未影响现有项目目录
- **验证**: `conda activate bio` 后 `python -V` / `fastqc --version` / `bwa` / `samtools --version` 可用；`python scripts/maintenance/generate_env_report.py` 运行成功

## [2026-02-28] - 工具部署：补齐常用生信 CLI
- **操作**: 在 conda env `bio` 中安装 `multiqc`, `fastp`, `seqkit`, `STAR`, `hisat2`, `subread`(featureCounts), `mafft`, `viennarna`(RNAfold), `prodigal`；更新环境快照
- **结果**: 成功
- **影响**: `bio` 环境新增常用工具；`samtools`/`bcftools`/`htslib` 因依赖被 bioconda 回退到 1.22.x（功能不受影响，属常见依赖约束）
- **验证**: `multiqc --version` / `fastp --version` / `seqkit version` / `STAR --version` / `hisat2 --version` / `featureCounts -v` / `mafft --version` / `RNAfold --version` / `prodigal -v` 均可运行；`docs/ENVIRONMENT.md` 更新

## [2026-02-28] - 项目：新增 AI 设计玩具工作流（纯干实验）
- **操作**: 创建示例项目 `projects/ai_design_playground/`，提供端到端脚本 `projects/ai_design_playground/scripts/pipeline.py`（toy DNA 生成 → Prodigal ORF → RNAfold MFE → 报告输出），支持可选 ESM contacts 校验
- **结果**: 成功
- **影响**: 新增一个可直接运行的 dry-lab playground，不依赖真实实验材料
- **验证**: 在 `conda activate bio` 后运行 `cd projects/ai_design_playground/scripts && python pipeline.py` 生成 `data/results/report.md`

## [2026-02-28] - 环境：补齐 `protein_predict.py` 依赖
- **操作**: 在 conda env `bio` 中安装 `scikit-learn`（用于 `tools/scripts/protein_predict.py` 内部 PCA 降维步骤）；更新环境快照
- **结果**: 成功
- **影响**: `docs/ENVIRONMENT.md` 的 `sklearn` 状态从 Not Installed 变为 Installed
- **验证**: `conda activate bio` 后 `python -c "import sklearn"` 成功；`python scripts/maintenance/generate_env_report.py` 运行成功

## [2026-02-28] - 维护：仓库一致性与目录洁癖
- **操作**: 添加 `.gitmodules` 修复 `repositories/active/RFdiffusion` 子模块映射；将根目录 `test_env_validation/` 迁移到 `projects/test_env_validation/legacy_toolchain_validation/`；移除根目录 `notebooks/`；修复 `projects/yeast_rnaseq_demo/data/references/` 中坏的绝对软链接为相对链接；生成环境快照
- **结果**: 成功
- **影响**: 仅工作区结构/文档变更；未安装/升级工具；`docs/ENVIRONMENT.md` 更新为当前 shell 环境检测结果（未激活 conda 时会显示工具缺失）
- **验证**: `git submodule status` 不再报错；工作区无 broken symlink；`python scripts/maintenance/generate_env_report.py` 运行成功

## [2026-01-26] - 文档与结构同步
- **操作**: 修正文档指向、补齐分析协议与索引、合并工具部署文档、更新共享数据说明
- **结果**: 成功
- **影响**: 仅文档与脚本，无运行环境变更
- **验证**: `docs/ENVIRONMENT.md` 已更新

## [2026-01-26] - 文档清理
- **操作**: 移除 `docs/AI_TOOL_DEPLOYMENT.md`
- **结果**: 成功
- **影响**: 仅文档调整
- **验证**: N/A

## [2026-01-25] - 初始化
- **操作**: 创建变更日志
- **结果**: 成功
- **影响**: 无
- **验证**: N/A
