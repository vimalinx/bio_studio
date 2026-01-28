# 🚀 Biomni 快速上手指南

## ⚡ 5分钟快速开始

### 步骤1: 安装环境

```bash
# 选择一个环境安装方式
conda env create -f biomni_env/environment.yml        # 推荐：基础环境
conda activate biomni_e1

# 安装Biomni
pip install biomni --upgrade
```

### 步骤2: 配置API

```bash
# 运行配置工具
python switch_profile.py

# 选择一个配置：
# 1. anthropic   - Claude模型
# 2. openai      - GPT模型
# 3. custom      - 本地模型
```

### 步骤3: 编辑API密钥

```bash
nano .env
# 填入你的API密钥
```

### 步骤4: 开始使用

```python
from biomni.agent import A1

agent = A1(path='./data', llm='claude-sonnet-4-5')
agent.go("规划一个CRISPR筛选实验")
```

## 📁 项目结构

```
Biomni/
├── profiles/                    # API配置文件目录
│   ├── anthropic.env           # Claude模型配置
│   ├── openai.env              # GPT模型配置
│   ├── azure.env               # Azure OpenAI配置
│   ├── custom.env              # 自定义模型配置
│   ├── biomni-r0.env           # Biomni-R0配置
│   └── default.env             # 默认配置
├── switch_profile.py            # 配置切换工具 ⭐
├── test_config.py               # 配置测试工具
├── quick_start.sh               # 快速启动脚本
├── API_CONFIG_GUIDE.md          # 详细配置指南
└── SETUP_API_CONFIG.md          # 完整部署文档
```

## 🎯 常用命令

### 配置管理

```bash
# 交互式菜单
python switch_profile.py

# 列出所有配置
python switch_profile.py list

# 切换配置
python switch_profile.py switch anthropic
python switch_profile.py switch openai
python switch_profile.py switch custom

# 查看配置详情
python switch_profile.py info anthropic

# 测试配置
python test_config.py
```

### 使用不同模型

```python
# Claude (Anthropic)
from biomni.agent import A1
agent = A1(path='./data', llm='claude-sonnet-4-5')

# GPT-4o (OpenAI)
agent = A1(path='./data', llm='gpt-4o')

# Azure OpenAI
agent = A1(path='./data', llm='azure-gpt-4o')

# 本地模型 (Ollama/vLLM/SGLang)
agent = A1(
    path='./data',
    llm='llama3.2',
    source='Custom',
    base_url='http://localhost:8000/v1'
)

# Biomni-R0 (需要先启动SGLang服务器)
from biomni.config import default_config
default_config.llm = "claude-3-5-sonnet-20241022"
default_config.source = "Anthropic"

agent = A1(
    llm="biomni/Biomni-R0-32B-Preview",
    source="Custom",
    base_url="http://localhost:30000/v1",
    api_key="EMPTY"
)
```

## 🔧 快速切换场景

### 场景1: 开发 → 生产

```bash
# 开发：使用便宜的模型
python switch_profile.py switch openai
# 编辑 .env: BIOMNI_LLM=gpt-4o-mini

# 生产：使用最佳模型
python switch_profile.py switch anthropic
# 编辑 .env: BIOMNI_LLM=claude-sonnet-4-5
```

### 场景2: 云端 → 本地

```bash
# 云端API
python switch_profile.py switch anthropic

# 本地模型
python switch_profile.py switch custom
# 确保本地服务已启动：ollama serve
```

### 场景3: 不同项目

```bash
# 项目A: 使用OpenAI
python switch_profile.py switch openai

# 项目B: 使用Anthropic
python switch_profile.py switch anthropic

# 项目C: 使用本地模型
python switch_profile.py switch custom
```

## 📖 文档索引

| 文档 | 用途 |
|------|------|
| **API_CONFIG_GUIDE.md** | 配置切换系统详细指南 |
| **SETUP_API_CONFIG.md** | 完整部署和使用文档 |
| **README.md** | Biomni官方文档 |
| **profiles/README.md** | 配置文件说明 |

## 🎓 示例任务

```python
from biomni.agent import A1

# 初始化
agent = A1(path='./data', llm='claude-sonnet-4-5')

# CRISPR筛选设计
agent.go("设计CRISPR筛选来识别调节T细胞耗竭的基因")

# 单细胞分析
agent.go("分析 [路径] 的scRNA-seq数据")

# 药物性质预测
agent.go("预测化合物的ADMET性质：CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")

# 基因功能分析
agent.go("分析BRCA1基因的功能和疾病关联")
```

## 🌐 Web界面

```python
from biomni.agent import A1

agent = A1(path='./data', llm='claude-sonnet-4-5')
agent.launch_gradio_demo()
# 访问 http://localhost:7860
```

## ⚙️ 高级配置

### 全局配置

```python
from biomni.config import default_config
from biomni.agent import A1

# 设置全局默认
default_config.llm = "gpt-4o"
default_config.timeout_seconds = 1200

# 所有代理使用这个配置
agent = A1()
```

### 环境变量

```bash
# 在 .bashrc 或 .zshrc 中设置
export ANTHROPIC_API_KEY="your_key"
export OPENAI_API_KEY="your_key"
export LLM_SOURCE="Anthropic"
export BIOMNI_LLM="claude-sonnet-4-5"
```

## 🆘 获取帮助

```bash
# 配置工具帮助
python switch_profile.py --help

# 测试配置
python test_config.py

# 查看文档
cat API_CONFIG_GUIDE.md
cat SETUP_API_CONFIG.md
```

## 🔗 有用的链接

- 🌐 [Biomni官网](https://biomni.stanford.edu)
- 📚 [GitHub仓库](https://github.com/snap-stanford/Biomni)
- 💬 [Slack社区](https://join.slack.com/t/biomnigroup/shared_invite/zt-3avks4913-dotMBt8D_apQnJ3mG~ak6Q)
- 📖 [论文](https://www.biorxiv.org/content/10.1101/2025.05.30.656746v1)

---

**需要更多帮助？**
- 查看 [API_CONFIG_GUIDE.md](API_CONFIG_GUIDE.md) 了解配置详情
- 查看 [SETUP_API_CONFIG.md](SETUP_API_CONFIG.md) 了解完整部署流程
- 运行 `python switch_profile.py` 开始使用配置工具
