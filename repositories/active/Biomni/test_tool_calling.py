#!/usr/bin/env python3
"""
测试 DeepSeek-Reasoner 工具调用能力
"""
import os

# 取消代理设置
for var in ['http_proxy', 'https_proxy', 'HTTP_PROXY', 'HTTPS_PROXY', 
            'all_proxy', 'ALL_PROXY', 'no_proxy', 'NO_PROXY']:
    os.environ.pop(var, None)

from dotenv import load_dotenv
load_dotenv("/media/vimalinx/Data/bio_studio/Biomni/.env", override=False)

print("=" * 70)
print("🧪 测试 DeepSeek-Reasoner 工具调用能力")
print("=" * 70)

# 配置默认 LLM（用于工具内的数据库查询）
from biomni.config import default_config
default_config.llm = "deepseek-chat"
default_config.source = "Custom"
default_config.base_url = "https://api.deepseek.com"
default_config.api_key = "sk-6f73c67f11d5469e846aba019b0f3530"

from biomni.agent import A1

# 初始化 Agent
print("\n🔧 初始化 Agent...")
agent = A1(
    path='/media/vimalinx/Data/bio_studio/Biomni/data',
    llm='deepseek-reasoner',
    source='Custom',
    base_url='https://api.deepseek.com',
    api_key='sk-6f73c67f11d5469e846aba019b0f3530',
    expected_data_lake_files=[],
    timeout_seconds=1200,
    use_tool_retriever=False
)

print("\n✅ Agent 已就绪")

# 查看可用工具
print("\n📋 可用工具模块:")
for module_name in list(agent.module2api.keys())[:10]:
    tools = agent.module2api[module_name]
    print(f"    - {module_name}: {len(tools)} 个工具")

print("\n" + "=" * 70)
print("🧪 测试: 调用 UniProt 工具查询 BRCA1 蛋白信息")
print("-" * 70)
try:
    # 使用直接 URL 查询（不需要 LLM）
    agent.go("使用 UniProt API 查询 BRCA1 人类蛋白的信息，直接使用这个 URL: https://rest.uniprot.org/uniprotkb/P38398")
except Exception as e:
    print(f"错误: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 70)
print("✅ 测试完成")
print("=" * 70)
