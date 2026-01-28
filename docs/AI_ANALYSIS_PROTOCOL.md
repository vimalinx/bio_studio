# AI 分析执行协议 (AI_ANALYSIS_PROTOCOL.md)

> **版本**: 1.0  
> **最后更新**: 2026-01-25  
> **用途**: 规范分析流程，保证可复现与项目隔离

---

## 🎯 目标

- 分析可复现、可追踪
- 项目数据与脚本隔离
- 结果有明确归档位置

---

## ✅ 执行流程

### 1. 明确任务与输入输出
确认：
- 输入数据类型与路径
- 预期输出（报告/图表/表格/模型）
- 关键参数与参考版本

### 2. 创建/确认项目结构
新项目使用模板：
```bash
python lib/create_project.py <项目名> --type <类型>
```
已有项目则检查 `projects/<项目名>/README.md` 与目录完整性。

### 3. 数据整理
- 原始数据放 `projects/<项目名>/data/raw/`
- 中间结果放 `data/processed/`
- 最终结果放 `data/results/`

#### 🚫 严禁操作
- **禁止** 在工作区根目录创建 `data/` 文件夹
- **禁止** 跨项目直接引用数据（使用 `shared_data` 或复制）

### 4. 配置与脚本
在 `projects/<项目名>/scripts/` 下：
- 更新 `config.py`（路径、样本、参考）
- 编写/更新 `pipeline.py`

### 5. 运行与记录
- 运行分析脚本
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

## 🧭 约定与建议

- 通用参考优先放 `shared_data/references/`
- 路径操作使用 `pathlib.Path`
- 旧结果放 `archive/` 避免混淆
