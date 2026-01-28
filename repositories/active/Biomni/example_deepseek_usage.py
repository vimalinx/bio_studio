#!/usr/bin/env python3
"""
Biomni + DeepSeek 完整使用示例

演示如何使用DeepSeek API运行Biomni生物医学AI代理
"""

import os
from dotenv import load_dotenv

# 加载环境变量
load_dotenv()

print("=" * 70)
print("Biomni + DeepSeek 使用示例")
print("=" * 70)

# ============================================
# 示例 1: 基础配置
# ============================================
print("\n【示例 1】基础配置\n")
print("""
from biomni.agent import A1

# 方法1: 直接传入参数
agent = A1(
    path='./data',
    llm='deepseek-chat',
    source='Custom',
    base_url='https://api.deepseek.com',
    api_key='sk-6f73c67f11d5469e846aba019b0f3530'
)

# 方法2: 使用配置文件（推荐）
# 先运行: python switch_profile.py switch deepseek
# 然后编辑 .env 文件填入API密钥
# 之后就可以简单使用:
agent = A1(path='./data')

# 执行任务
agent.go("设计一个CRISPR筛选实验来研究T细胞耗竭机制")
""")

# ============================================
# 示例 2: 使用配置管理
# ============================================
print("\n【示例 2】使用全局配置\n")
print("""
from biomni.config import default_config
from biomni.agent import A1

# 设置全局默认配置
default_config.llm = "deepseek-chat"
default_config.source = "Custom"
default_config.base_url = "https://api.deepseek.com"
default_config.api_key = "sk-6f73c67f11d5469e846aba019b0f3530"
default_config.temperature = 0.7

# 所有代理和数据库查询都使用这个配置
agent = A1()

# 现在可以执行任务
agent.go("分析BRCA1基因的功能和疾病关联")
""")

# ============================================
# 示例 3: 生物医学任务
# ============================================
print("\n【示例 3】实际生物医学任务\n")
print("""
from biomni.agent import A1

agent = A1(
    path='./data',
    llm='deepseek-chat',
    source='Custom',
    base_url='https://api.deepseek.com',
    api_key='your_api_key_here'
)

# 任务1: CRISPR筛选设计
agent.go("设计一个CRISPR-Cas9筛选实验，识别调节乳腺癌细胞对药物耐药性的基因")

# 任务2: 单细胞RNA测序分析
agent.go("分析单细胞RNA测序数据，识别肿瘤微环境中的免疫细胞亚群")

# 任务3: 药物性质预测
agent.go("预测这个SMILES字符串化合物的ADMET性质: CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")

# 任务4: 基因功能分析
agent.go("查找TP53基因的功能、相关疾病和最新研究进展")

# 任务5: 生物信息学流程设计
agent.go("设计一个NGS数据分析流程，用于检测体细胞突变")
""")

# ============================================
# 示例 4: Gradio界面
# ============================================
print("\n【示例 4】启动Web界面\n")
print("""
from biomni.agent import A1

agent = A1(
    path='./data',
    llm='deepseek-chat',
    source='Custom',
    base_url='https://api.deepseek.com',
    api_key='your_api_key_here'
)

# 启动Gradio界面
agent.launch_gradio_demo()

# 在浏览器中访问 http://localhost:7860
# 可以通过图形界面与Biomni交互
""")

# ============================================
# 示例 5: 成本优化
# ============================================
print("\n【示例 5】DeepSeek成本优势\n")
print("""
DeepSeek vs 其他模型 (价格对比，每百万tokens):

模型                    | 输入价格   | 输出价格
------------------------|------------|------------
DeepSeek-Chat          | ¥1.0       | ¥2.0
Claude Sonnet 4.5      | $3.0       | $15.0
GPT-4o                 | $2.5       | $10.0

优势:
✓ 成本仅为Claude/GPT的 1/10 到 1/20
✓ 中文支持优秀
✓ 32K上下文窗口
✓ 适合大规模生物医学数据分析

典型任务成本:
- CRISPR筛选设计 (~5000 tokens): ¥0.01
- 单细胞分析 (~10000 tokens): ¥0.02
- 复杂多步推理 (~50000 tokens): ¥0.1
""")

# ============================================
# 配置步骤
# ============================================
print("\n" + "=" * 70)
print("快速配置步骤")
print("=" * 70)
print("""
步骤1: 切换到DeepSeek配置
    python switch_profile.py switch deepseek

步骤2: 编辑API密钥
    nano .env
    # 确认API密钥已正确设置

步骤3: 测试连接
    python test_deepseek.py

步骤4: 运行Biomni
    python
    >>> from biomni.agent import A1
    >>> agent = A1(path='./data')
    >>> agent.go("你的生物医学任务")

注意事项:
1. 确保网络可以访问 api.deepseek.com
2. 如使用代理，请临时取消代理设置:
   unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY
3. 首次运行会下载约11GB的数据湖文件
4. DeepSeek响应速度很快，适合交互式使用
""")

# ============================================
# 实际测试
# ============================================
print("\n" + "=" * 70)
print("实际测试 (可选)")
print("=" * 70)

test_input = input("\n是否要运行实际的Biomni测试？这可能需要几分钟。(y/N): ").strip().lower()

if test_input == 'y':
    try:
        from biomni.llm import get_llm
        from langchain_core.messages import HumanMessage

        print("\n正在测试Biomni + DeepSeek...")

        llm = get_llm(
            model="deepseek-chat",
            source="Custom",
            base_url="https://api.deepseek.com",
            api_key="sk-6f73c67f11d5469e846aba019b0f3530",
            temperature=0.7
        )

        # 测试简单任务
        messages = [HumanMessage(content="简要介绍CRISPR-Cas9技术在基因编辑中的应用（不超过100字）")]
        response = llm.invoke(messages)

        print("\n响应:")
        print("-" * 70)
        print(response.content)
        print("-" * 70)
        print("\n✓ 测试成功！Biomni + DeepSeek 已准备就绪。")

    except Exception as e:
        print(f"\n✗ 测试失败: {e}")
        print("请确保已安装所需依赖: pip install biomni langchain-openai")
else:
    print("\n跳过实际测试。")

print("\n" + "=" * 70)
print("更多信息:")
print("  - DeepSeek API文档: https://platform.deepseek.com/api-docs/")
print("  - Biomni文档: README.md")
print("  - 配置指南: API_CONFIG_GUIDE.md")
print("=" * 70)
