# Biomni 部署与API切换配置指南

本文档说明如何部署Biomni并方便地切换不同的API配置。

## 快速开始

### 1. 安装环境

首先需要设置Biomni环境：

```bash
# 基础环境（推荐开始使用这个）
conda env create -f biomni_env/environment.yml

# 或者完整环境E1（需要>10小时和>30GB空间）
bash biomni_env/setup.sh

# 或者精简环境（不含R和CLI工具）
conda env create -f biomni_env/fixed_env.yml
```

激活环境：

```bash
conda activate biomni_e1
```

安装Biomni包：

```bash
# 从PyPI安装
pip install biomni --upgrade

# 或从GitHub源安装最新版本
pip install git+https://github.com/snap-stanford/Biomni.git@main
```

### 2. 配置API密钥

使用我们提供的配置切换工具：

```bash
# 交互式菜单模式
python switch_profile.py

# 或者直接列出所有可用配置
python switch_profile.py list

# 切换到特定配置（例如：Anthropic）
python switch_profile.py switch anthropic

# 查看配置信息
python switch_profile.py info openai
```

### 3. 编辑API密钥

切换到配置后，需要编辑 `.env` 文件填入实际的API密钥：

```bash
# 使用文本编辑器
nano .env

# 或使用vim
vim .env
```

或使用配置工具：

```bash
python switch_profile.py
# 然后选择 'e' 选项直接编辑 .env 文件
```

## 可用配置说明

### 1. **anthropic.env** - Anthropic Claude模型
适用场景：只想使用Claude模型
```python
from biomni.agent import A1
agent = A1(path='./data', llm='claude-sonnet-4-5')
```

### 2. **openai.env** - OpenAI模型
适用场景：只想使用OpenAI模型（GPT-4o等）
```python
from biomni.agent import A1
agent = A1(path='./data', llm='gpt-4o')
```

### 3. **azure.env** - Azure OpenAI
适用场景：使用Azure OpenAI服务
```python
from biomni.agent import A1
agent = A1(path='./data', llm='azure-gpt-4o')
```

### 4. **custom.env** - 自定义模型服务
适用场景：使用Ollama、vLLM、SGLang等本地或自定义模型服务
```python
from biomni.agent import A1
agent = A1(
    path='./data',
    llm='llama3.2',  # 或其他自定义模型名称
    source='Custom',
    base_url='http://localhost:8000/v1',
    api_key='EMPTY'
)
```

### 5. **biomni-r0.env** - Biomni-R0模型
适用场景：使用Biomni-R0推理模型（需要先启动SGLang服务器）

首先启动SGLang服务器：
```bash
python -m sglang.launch_server \
  --model-path RyanLi0802/Biomni-R0-Preview \
  --port 30000 \
  --host 0.0.0.0 \
  --mem-fraction-static 0.8 \
  --tp 2 \
  --trust-remote-code \
  --json-model-override-args '{"rope_scaling":{"rope_type":"yarn","factor":1.0,"original_max_position_embeddings":32768}, "max_position_embeddings": 131072}'
```

然后使用：
```python
from biomni.config import default_config
from biomni.agent import A1

# 数据库查询使用Claude
default_config.llm = "claude-3-5-sonnet-20241022"
default_config.source = "Anthropic"

# 代理推理使用Biomni-R0
agent = A1(
    llm="biomni/Biomni-R0-32B-Preview",
    source="Custom",
    base_url="http://localhost:30000/v1",
    api_key="EMPTY"
)

agent.go("Your biomedical task here")
```

## 使用示例

### 基础使用

```python
from biomni.agent import A1

# 初始化代理（首次运行会自动下载~11GB的数据湖）
agent = A1(path='./data', llm='claude-sonnet-4-5')

# 执行生物医学任务
agent.go("规划一个CRISPR筛选实验来识别调节T细胞耗竭的基因")
agent.go("分析[路径]的单细胞RNA测序数据并生成假设")
agent.go("预测这个化合物的ADMET属性：CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")
```

### 启动Gradio界面

```python
from biomni.agent import A1

agent = A1(path='./data', llm='claude-sonnet-4-5')
agent.launch_gradio_demo()
# 访问 http://localhost:7860
```

### 配置管理

```python
from biomni.config import default_config
from biomni.agent import A1

# 全局配置（推荐）
default_config.llm = "gpt-4o"
default_config.timeout_seconds = 1200

# 所有代理和数据库查询都使用这个配置
agent = A1()
```

## 快速切换API配置

### 方法1：命令行切换

```bash
# 列出所有配置
python switch_profile.py

# 切换到OpenAI配置
python switch_profile.py switch openai

# 切换到Anthropic配置
python switch_profile.py switch anthropic

# 切换到自定义模型配置
python switch_profile.py switch custom
```

### 方法2：交互式菜单

```bash
python switch_profile.py
```

然后按照菜单提示选择配置。

### 方法3：手动切换

```bash
# 复制配置文件
cp profiles/openai.env .env

# 编辑API密钥
nano .env
```

## 创建自定义配置

### 使用工具创建

```bash
python switch_profile.py
# 选择 'c' 选项
# 输入配置名称
# 选择模板或创建空配置
# 填入API密钥
```

### 手动创建

1. 复制现有配置：
```bash
cp profiles/anthropic.env profiles/myconfig.env
```

2. 编辑新配置：
```bash
nano profiles/myconfig.env
```

3. 切换到新配置：
```bash
python switch_profile.py switch myconfig
```

## 常见问题

### 1. API密钥安全

- 不要将包含真实API密钥的配置文件提交到Git
- `.env` 文件已在 `.gitignore` 中，不会被提交
- 可以使用 `profiles/` 目录存储不同环境的配置

### 2. 模型名称格式

- **Anthropic**: `claude-sonnet-4-5`, `claude-opus-4-20250514`
- **OpenAI**: `gpt-4o`, `gpt-4o-mini`, `gpt-4-turbo`
- **Azure**: 必须加 `azure-` 前缀，如 `azure-gpt-4o`
- **自定义**: 取决于你的模型服务（如Ollama中的 `llama3.2`）

### 3. 环境变量优先级

1. 代码中的直接参数（优先级最高）
2. `.env` 文件
3. 系统环境变量
4. `default_config` 全局配置

### 4. 多配置管理

可以同时保存多个配置文件，根据需要切换：

```bash
# 开发环境使用便宜的模型
python switch_profile.py switch openai  # gpt-4o-mini

# 生产环境使用更好的模型
python switch_profile.py switch anthropic  # claude-sonnet-4-5

# 本地测试使用自定义模型
python switch_profile.py switch custom
```

## 下一步

1. 阅读 [README.md](README.md) 了解更多功能
2. 查看 [docs/configuration.md](docs/configuration.md) 了解详细配置选项
3. 尝试 [tutorials/biomni_101.ipynb](tutorials/biomni_101.ipynb) 教程
4. 参考 [CONTRIBUTION.md](CONTRIBUTION.md) 贡献工具和数据集

## 支持的模型提供商

- **Anthropic**: Claude Sonnet, Opus, Haiku
- **OpenAI**: GPT-4o, GPT-4 Turbo, GPT-3.5
- **Azure OpenAI**: Azure托管的OpenAI模型
- **Google Gemini**: Gemini Pro, Ultra
- **Groq**: 快速推理服务
- **AWS Bedrock**: 亚马逊云模型服务
- **Ollama**: 本地开源模型
- **自定义**: 任何OpenAI兼容的API（vLLM, SGLang, Text Generation WebUI等）
