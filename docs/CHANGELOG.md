# Bio Studio 变更日志

记录工具部署、分析流程和环境变更。每次操作完成后追加一条记录。

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
