# AI 工具部署协议 (AI_TOOL_PROTOCOL.md)

> **版本**: 1.0
> **最后更新**: 2026-01-26
> **用途**: AI 部署新工具时必须遵循的完整流程
> **模式**: AI 全托管（用户给出来源，AI完成部署与文档）
>
> **说明**: 当前工作区不包含 `auto-deploy/`。涉及自动测试/回滚的步骤需手动执行，并用以下脚本替代：
> - 环境快照: `python scripts/maintenance/generate_env_report.py`
> - Skill 模板: `python tools/scripts/deploy_tool.py --scan` 或 `--repo <path>`

---

## 🔔 触发条件

当用户指令包含以下意图时，进入本协议流程：
- 明确提出“添加工具 / 部署工具”
- 提供 GitHub 仓库链接并希望使用该工具

---

## 🎯 协议目标

确保每次工具部署：
1. ✅ 工具可正常使用
2. ✅ 不破坏现有环境
3. ✅ 有完整的文档和 Skill
4. ✅ 可回滚

---

## 📋 部署前检查清单

在开始部署前，必须完成以下检查：

```
□ 1. 确认工具名称和来源 (Git URL / Conda 包名 / PyPI 包名)
□ 2. 确认工具用途和类别
□ 3. 生成环境快照: python scripts/maintenance/generate_env_report.py
□ 4. 记录变更计划: docs/CHANGELOG.md (新增待办条目)
□ 5. 阅读工具的官方文档/README
```

---

## 🔄 完整部署流程

### 阶段 1: 信息收集 (5 分钟)

**步骤 1.1: 确认工具信息**
```
工具名称: _______________
来源类型: [ ] Git仓库  [ ] Conda包  [ ] PyPI包  [ ] 二进制文件
来源地址: _______________
工具类别: [ ] 质控  [ ] 比对  [ ] 变异  [ ] 定量  [ ] 可视化  [ ] AI模型  [ ] 其他
依赖要求: _______________
```

**步骤 1.2: 阅读官方文档**
- 找到工具的 README 或文档
- 记录安装命令
- 记录依赖要求
- 记录使用示例

**步骤 1.3: 智能分析 (Context Gathering)**
AI 必须主动阅读以下内容（无需用户手动提供）：
1. **README**: 理解核心功能、安装依赖、关键参数
2. **Examples**: `examples/` 目录下的脚本，提取真实用法模式
3. **Scripts**: 核心 Python/Shell 脚本，理解输入输出格式
4. **Environment**: `requirements.txt` 或 `environment.yml`

**步骤 1.4: 检查兼容性**
```bash
# 检查是否已安装
which <工具名> || conda list | grep <工具名>

# 检查 Python 版本要求
# 检查是否需要特定的库版本（可能冲突）
```

---

### 阶段 2: 环境准备 (2 分钟)

**步骤 2.1: 记录环境快照**
```bash
python scripts/maintenance/generate_env_report.py
```

**步骤 2.2: 激活环境**
```bash
source ~/miniforge3/bin/activate bio
```

---

### 阶段 3: 安装工具 (时间不定)

根据工具来源选择安装方式：

#### 方式 A: Conda/Mamba 安装 (推荐)
```bash
# 优先用 mamba（更快）
mamba install -n bio <包名> -c conda-forge -c bioconda -y

# 或用 conda
conda install -n bio <包名> -c conda-forge -c bioconda -y
```

#### 方式 B: PyPI 安装
```bash
# 激活环境后
pip install <包名>
```

#### 方式 C: Git 仓库安装
```bash
# 1. 克隆到 repositories/active/
cd ~/bio_studio/repositories/active
git clone <git_url> <工具名>

# 2. 安装依赖
cd <工具名>
pip install -r requirements.txt  # 如果有
pip install -e .  # 如果是 Python 包

# 3. 或者按官方文档安装
```

#### 方式 D: 二进制文件
```bash
# 1. 下载到 tools/bin/
cd ~/bio_studio/tools/bin
wget <下载链接>
chmod +x <文件名>

# 2. 添加到 PATH（在 bio 环境中）
# 编辑 ~/miniforge3/envs/bio/etc/conda/activate.d/bio_studio.sh
```

---

### 阶段 4: 验证安装 (5 分钟)

**步骤 4.1: 基础验证**
```bash
# 检查命令是否可用
which <工具名>

# 查看版本/帮助
<工具名> --version
<工具名> --help
```

**步骤 4.2: 功能验证**
```bash
# 用官方示例或最小测试数据运行
# 记录运行命令和预期输出
```

**步骤 4.3: 环境完整性验证**
```bash
# 确认没有破坏其他工具（生成最新快照并对比）
python scripts/maintenance/generate_env_report.py
```

如果任何步骤失败：
```bash
# 手动回滚（卸载包/删除仓库/恢复配置）
conda remove -n bio <包名>
pip uninstall <包名>
rm -rf repositories/active/<工具名>
```

---

### 阶段 5: 创建 Skill 文件 (10 分钟)

推荐先生成模板：
```bash
python tools/scripts/deploy_tool.py --repo ~/bio_studio/repositories/active/<工具名>
```

**步骤 5.1: 创建 Skill 目录**
```bash
mkdir -p ~/bio_studio/.claude/skills/<工具名>
```

**步骤 5.2: 编写 SKILL.md**

使用以下模板创建 `~/.claude/skills/<工具名>/SKILL.md`:

```markdown
# <工具名> Skill

## 概述
<一句话描述工具用途>

## 安装信息
- **安装方式**: <conda/pip/git/binary>
- **安装命令**: `<安装命令>`
- **版本**: <版本号>
- **安装日期**: <日期>

## 用途
<详细描述工具能做什么>

## 基本用法

### 命令格式
```bash
<工具名> [选项] <输入> <输出>
```

### 常用参数
| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-i/--input` | 输入文件 | 必需 |
| `-o/--output` | 输出文件 | 必需 |
| `-t/--threads` | 线程数 | 1 |

### 示例

#### 示例 1: 基本用法
```bash
<具体命令示例>
```

#### 示例 2: 高级用法
```bash
<具体命令示例>
```

## 输入输出

### 输入格式
- <格式说明>

### 输出格式
- <格式说明>

## 常见问题

### Q: <常见问题 1>
A: <解答>

### Q: <常见问题 2>
A: <解答>

## 相关工具
- <相关工具 1>
- <相关工具 2>

## 参考链接
- 官方文档: <链接>
- GitHub: <链接>
```

---

### 阶段 6: 创建/更新 MCP 服务器 (可选, 15 分钟)

如果工具需要 MCP 集成：

**步骤 6.1: 判断是否需要 MCP**
- 工具是否经常被调用？
- 是否需要程序化访问？
- 是否有复杂的参数？

如果是，继续创建 MCP 服务器。

**步骤 6.2: 创建 MCP 服务器目录**
```bash
mkdir -p ~/bio_studio/mcp-servers/<工具名>-mcp
```

**步骤 6.3: 使用模板创建服务器**
复制已有服务器模板（如 `mcp-servers/bio-sequence-mcp/sequence_server.py`）并修改。

**步骤 6.4: 更新 MCP 配置**
编辑 `mcp-servers/claude-config.json`，添加新服务器。

---

### 阶段 7: 更新文档 (5 分钟)

**步骤 7.1: 更新 CHANGELOG**
```bash
# 编辑 docs/CHANGELOG.md，添加:
## [2026-01-25] - 工具部署
- **操作**: 安装 <工具名> <版本>
- **方式**: <conda/pip/git>
- **结果**: 成功
- **验证**: <工具名> --help 正常，ENVIRONMENT.md 已更新
```

**步骤 7.2: 更新环境报告**
```bash
python scripts/maintenance/generate_env_report.py
```

**步骤 7.3: 更新 BEST_PRACTICES (如果踩坑)**
如果安装过程中遇到问题并解决了，记录到 `docs/BEST_PRACTICES.md`。

---

### 阶段 8: 最终验证 (2 分钟)

**完成检查清单:**
```
□ 工具可以运行 --help
□ 工具可以处理测试数据
□ 环境快照已更新 (ENVIRONMENT.md)
□ SKILL.md 已创建
□ CHANGELOG.md 已更新
```

全部完成后，向用户报告：
```
✅ 工具 <名称> 部署完成

安装信息:
- 版本: x.x.x
- 安装方式: conda/pip/git
- Skill 位置: .claude/skills/<工具名>/SKILL.md

验证结果:
- 基础测试: 通过
- 环境完整性: 通过

使用示例:
<一个简单的使用示例>
```

---

## 🚨 失败处理流程

如果任何阶段失败：

### 步骤 1: 停止操作
不要继续尝试，先分析问题。

### 步骤 2: 收集信息
```bash
# 保存错误日志
# 记录失败的命令
# 记录错误信息
```

### 步骤 3: 尝试回滚
```bash
conda remove -n bio <包名>
pip uninstall <包名>
rm -rf repositories/active/<工具名>
```

### 步骤 4: 验证回滚
```bash
python scripts/maintenance/generate_env_report.py
```

### 步骤 5: 报告用户
```
❌ 工具 <名称> 部署失败

失败阶段: <阶段名>
错误信息: <错误信息>
已尝试: <尝试过的解决方法>
当前状态: 已回滚到部署前状态

建议:
- <可能的解决方案 1>
- <可能的解决方案 2>

需要用户协助:
- <需要什么帮助>
```

---

## 📝 常见工具安装速查

### Bioconda 工具 (大多数生信工具)
```bash
mamba install -n bio <工具名> -c bioconda -c conda-forge -y
```

### Python 库
```bash
pip install <库名>
```

### R 包 (在 conda 环境中)
```bash
mamba install -n bio r-<包名> -c conda-forge -y
# 或
R -e "install.packages('<包名>')"
```

### 需要编译的工具
```bash
cd ~/bio_studio/repositories
git clone <url>
cd <工具名>
mkdir build && cd build
cmake ..
make -j$(nproc)
make install  # 或复制到 tools/bin/
```

### Docker 容器化工具
```bash
docker pull <镜像名>
# 在 tools/scripts/ 创建包装脚本
```

---

## ✅ 协议完成确认

---

## 🧪 示例交互

**User**: "把 esm 添加进来"

**AI**:
1. `git clone ...`
2. 阅读 README 和 examples
3. 写入 `.claude/skills/esm/SKILL.md`
4. 报告: "已部署 ESM。它主要用于... 你可以使用 `/esm` 或直接问我 '分析蛋白质序列'。"

部署完成后，AI 必须能回答以下问题：

1. 工具安装在哪里？
2. 如何调用这个工具？
3. 输入输出格式是什么？
4. 有哪些常用参数？
5. 如果出问题怎么办？

如果无法回答任何一个问题，部署**未完成**。
