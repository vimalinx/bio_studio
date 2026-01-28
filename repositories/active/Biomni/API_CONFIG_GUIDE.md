# Biomni API配置切换系统 🔄

一个简单易用的API配置管理系统，让在多个API提供商之间切换变得轻松。

## 📋 目录

- [快速开始](#快速开始)
- [可用配置](#可用配置)
- [使用方法](#使用方法)
- [创建自定义配置](#创建自定义配置)
- [故障排除](#故障排除)

## 🚀 快速开始

### 1. 安装Biomni

```bash
# 设置环境（选择一个）
conda env create -f biomni_env/environment.yml        # 基础环境
# 或
bash biomni_env/setup.sh                               # 完整环境
# 或
conda env create -f biomni_env/fixed_env.yml           # 精简环境

# 激活环境
conda activate biomni_e1

# 安装Biomni
pip install biomni --upgrade
```

### 2. 切换配置

```bash
# 交互式菜单（推荐）
python switch_profile.py

# 或命令行直接切换
python switch_profile.py switch anthropic

# 或使用快速启动脚本
bash quick_start.sh
```

### 3. 配置API密钥

编辑 `.env` 文件，填入你的实际API密钥：

```bash
nano .env
# 或
vim .env
```

### 4. 开始使用

```python
from biomni.agent import A1

# 初始化代理
agent = A1(path='./data', llm='claude-sonnet-4-5')

# 执行任务
agent.go("你的生物医学任务")
```

## 📦 可用配置

| 配置名称 | 描述 | 适用场景 |
|---------|------|---------|
| **anthropic.env** | Anthropic Claude模型 | 只使用Claude模型 |
| **openai.env** | OpenAI GPT模型 | 只使用GPT-4o等模型 |
| **azure.env** | Azure OpenAI | 使用Azure托管的OpenAI模型 |
| **custom.env** | 自定义模型服务 | 使用Ollama、vLLM、SGLang等 |
| **biomni-r0.env** | Biomni-R0专用配置 | 使用Biomni-R0推理模型 |
| **default.env** | 默认配置 | 支持多个提供商 |

## 💡 使用方法

### 交互式菜单（推荐）

```bash
python switch_profile.py
```

功能：
- 查看所有可用配置
- 切换配置
- 查看配置详情
- 创建新配置
- 直接编辑.env文件

### 命令行模式

```bash
# 列出所有配置
python switch_profile.py list

# 切换到特定配置
python switch_profile.py switch anthropic
python switch_profile.py switch openai
python switch_profile.py switch custom

# 查看配置信息
python switch_profile.py info anthropic

# 查看帮助
python switch_profile.py --help
```

### 手动切换

```bash
# 复制配置文件
cp profiles/anthropic.env .env

# 编辑API密钥
nano .env
```

## 🎯 创建自定义配置

### 方法1：使用工具创建

```bash
python switch_profile.py
# 选择 'c' -> 创建新配置
# 输入配置名称
# 选择模板或创建空白配置
# 编辑配置内容
```

### 方法2：手动创建

```bash
# 1. 复制现有配置作为模板
cp profiles/anthropic.env profiles/myconfig.env

# 2. 编辑新配置
nano profiles/myconfig.env

# 3. 切换到新配置
python switch_profile.py switch myconfig
```

### 配置文件模板

```bash
# 配置文件模板示例

# API密钥（根据你的提供商选择）
ANTHROPIC_API_KEY=your_key_here
OPENAI_API_KEY=your_key_here
GEMINI_API_KEY=your_key_here

# 模型提供商（可选）
LLM_SOURCE=Anthropic  # 或 OpenAI, AzureOpenAI, Custom等

# 默认模型（可选）
BIOMNI_LLM=claude-sonnet-4-5

# 其他配置（可选）
BIOMNI_TEMPERATURE=0.7
BIOMNI_TIMEOUT_SECONDS=600
BIOMNI_DATA_PATH=./data

# 自定义模型配置（如果使用Custom）
# CUSTOM_MODEL_BASE_URL=http://localhost:8000/v1
# CUSTOM_MODEL_API_KEY=your_key_here
```

## 🔧 高级用法

### 使用Biomni-R0

```bash
# 1. 首先启动SGLang服务器
python -m sglang.launch_server \
  --model-path RyanLi0802/Biomni-R0-Preview \
  --port 30000 \
  --host 0.0.0.0 \
  --mem-fraction-static 0.8 \
  --tp 2 \
  --trust-remote-code

# 2. 切换到Biomni-R0配置
python switch_profile.py switch biomni-r0

# 3. 编辑.env文件，填入Anthropic API密钥（用于数据库查询）

# 4. 使用
from biomni.config import default_config
from biomni.agent import A1

default_config.llm = "claude-3-5-sonnet-20241022"
default_config.source = "Anthropic"

agent = A1(
    llm="biomni/Biomni-R0-32B-Preview",
    source="Custom",
    base_url="http://localhost:30000/v1",
    api_key="EMPTY"
)
```

### 使用Ollama本地模型

```bash
# 1. 启动Ollama
ollama serve

# 2. 切换到custom配置
python switch_profile.py switch custom

# 3. 编辑.env文件：
#    CUSTOM_MODEL_BASE_URL=http://localhost:11434/v1
#    LLM_SOURCE=Custom
#    BIOMNI_LLM=llama3.2

# 4. 使用
from biomni.agent import A1
agent = A1(
    llm="llama3.2",
    source="Custom",
    base_url="http://localhost:11434/v1"
)
```

### 环境特定配置

```bash
# 开发环境 - 使用快速/便宜的模型
python switch_profile.py switch openai
# 编辑 .env: BIOMNI_LLM=gpt-4o-mini

# 生产环境 - 使用最佳模型
python switch_profile.py switch anthropic
# 编辑 .env: BIOMNI_LLM=claude-sonnet-4-5

# 本地测试 - 使用自定义模型
python switch_profile.py switch custom
```

## 🧪 测试配置

运行测试脚本验证配置：

```bash
python test_config.py
```

测试内容：
- ✓ 环境文件检查
- ✓ API密钥配置
- ✓ 包导入测试
- ✓ 配置类测试
- ✓ LLM自动检测测试

## 🔍 故障排除

### 问题1: 找不到biomni包

**解决方案:**
```bash
pip install biomni --upgrade
# 或
pip install git+https://github.com/snap-stanford/Biomni.git@main
```

### 问题2: API密钥未生效

**解决方案:**
1. 确认.env文件存在且在项目根目录
2. 检查API密钥是否正确（没有多余的空格）
3. 重启Python解释器
4. 使用 `test_config.py` 验证配置

### 问题3: 无法切换配置

**解决方案:**
```bash
# 检查profiles目录是否存在
ls profiles/

# 确认配置文件存在
ls profiles/*.env

# 手动切换
cp profiles/anthropic.env .env
```

### 问题4: 模型名称错误

**常见模型名称:**
- Anthropic: `claude-sonnet-4-5`, `claude-opus-4-20250514`
- OpenAI: `gpt-4o`, `gpt-4o-mini`, `gpt-4-turbo`
- Azure: `azure-gpt-4o` (注意前缀)
- Custom: 取决于你的模型服务

## 📖 相关文档

- [完整部署指南](SETUP_API_CONFIG.md)
- [Biomni官方文档](README.md)
- [配置选项详解](docs/configuration.md)
- [贡献指南](CONTRIBUTION.md)

## 🤝 贡献

欢迎贡献！查看[贡献指南](CONTRIBUTION.md)了解如何参与。

## 📄 许可证

Biomni使用Apache 2.0许可证。详见[LICENSE](LICENSE)文件。

---

**快速链接:**
- 🔬 [Biomni官网](https://biomni.stanford.edu)
- 📚 [文档](README.md)
- 💬 [Slack社区](https://join.slack.com/t/biomnigroup/shared_invite/zt-3avks4913-dotMBt8D_apQnJ3mG~ak6Q)
- 🐦 [Twitter](https://x.com/ProjectBiomni)
