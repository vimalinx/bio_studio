# AI 分析执行协议 (AI_ANALYSIS_PROTOCOL.md)

> **版本**: 1.1  
> **最后更新**: 2026-03-23  
> **用途**: 规范分析流程，保证可复现与项目隔离

---

## 🎯 目标

- 分析可复现、可追踪
- 项目数据与脚本隔离
- 结果有明确归档位置
- 新建项目可立即完成模板级自检

---

## ✅ 执行流程

### 1. 明确任务与输入输出
确认：
- 输入数据类型与路径
- 预期输出（报告/图表/表格/模型）
- 关键参数与参考版本

### 2. 创建/确认项目结构
推荐使用工作区入口：
```bash
python scripts/project.py create <项目名> --type <类型> --description "<描述>"
```

创建后先运行项目级自检：
```bash
python scripts/project.py validate <项目名>
python scripts/project.py steps <项目名>
```

如需直接使用底层模板生成器，也可以运行：
```bash
python lib/create_project.py <项目名> --type <类型>
```

已有项目则检查 `projects/<项目名>/README.md`、目录完整性与最近一次验证报告。

### 3. 数据整理
- 原始数据放 `projects/<项目名>/data/raw/`
- 中间结果放 `data/processed/`
- 最终结果放 `data/results/`
- 项目特定参考放 `data/references/`

#### 🚫 严禁操作
- **禁止** 在工作区根目录创建 `data/` 文件夹
- **禁止** 跨项目直接引用数据（使用 `shared_data` 或复制）

### 4. 配置与脚本
在 `projects/<项目名>/scripts/` 下：
- 更新 `config.py`（路径、样本、参考）
- 先确认 `validate_project.py` 通过
- 编写/更新 `pipeline.py`

模板型项目在工作区内运行时，`pipeline.py` 会优先尝试调用 `lib/template_runtime.py`，再由它转发到 `lib/modules/` 中的可复用 wrapper。
如果项目被单独拷出工作区、找不到 `lib/`，则会自动回退到模板内置的提示型步骤，不影响 `--validate` 与 `--steps`。

推荐优先从工作区根目录执行：
```bash
python scripts/project.py run <项目名>
python scripts/project.py run <项目名> --step quality_control
```

注意：
- 对模板项目 / demo 项目，`run` 通常会直接执行 `pipeline.py`
- 对专项项目，`run` 会转发到项目自定义 pipeline
- 对学习项目，`run` 可以只输出推荐的下一步命令与文档入口，不一定直接执行整条教学链

建议：Notebook 仅用于临时探索；可复现逻辑应落在 `scripts/`。
（仓库默认忽略 `*.ipynb`，避免把大量中间探索文件提交进版本库。）

### 5. 运行与记录
- 运行分析脚本
- 项目级验证会写出：
  - `projects/<项目名>/logs/validation_report.json`
  - `projects/<项目名>/logs/validation_report.md`
- 关键输出写入 `projects/<项目名>/logs/`

### 6. 结果归档
建议输出：
- `projects/<项目名>/README.md` 更新结论
- 必要时生成 `ANALYSIS_REPORT.md`（项目根目录）

### 7. 文档同步
- 更新 `docs/CHANGELOG.md`
- 如环境变动，运行 `python scripts/maintenance/generate_env_report.py`
- 遇到坑，记录到 `docs/BEST_PRACTICES.md`

---

## 🧭 双层验证模型

- `projects/<项目名>/scripts/validate_project.py`：项目模板级自检，检查结构、配置和路径边界
- `projects/test_env_validation/`：工作区级环境 smoke test，检查工具链与依赖环境，可通过 `python scripts/project.py workspace-validate` 触发

不要把这两层验证混为一谈。模板验证通过，不代表整个工作区工具链已验收；环境 smoke test 通过，也不代表某个项目配置正确。

---

## 🧭 约定与建议

- 通用参考优先放 `shared_data/references/`
- 路径操作使用 `pathlib.Path`
- 可复用工具调用优先沉淀到 `lib/modules/`，不要在每个模板项目里重复手搓 wrapper
- 旧结果放 `archive/` 避免混淆
