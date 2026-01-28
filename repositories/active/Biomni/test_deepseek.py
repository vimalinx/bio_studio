#!/usr/bin/env python3
"""
测试 Biomni + DeepSeek-Reasoner 配置
"""
import os
import sys

# 取消代理设置
for var in ['http_proxy', 'https_proxy', 'HTTP_PROXY', 'HTTPS_PROXY', 
            'all_proxy', 'ALL_PROXY', 'no_proxy', 'NO_PROXY']:
    os.environ.pop(var, None)

from dotenv import load_dotenv
load_dotenv("/media/vimalinx/Data/bio_studio/Biomni/.env", override=False)

print("=" * 70)
print("🧪 测试 DeepSeek-Reasoner + Biomni 配置")
print("=" * 70)

# 检查环境变量
print("\n📋 环境变量:")
print(f"  LLM_SOURCE: {os.getenv('LLM_SOURCE')}")
print(f"  CUSTOM_MODEL_BASE_URL: {os.getenv('CUSTOM_MODEL_BASE_URL')}")
print(f"  BIOMNI_LLM: {os.getenv('BIOMNI_LLM')}")
print(f"  BIOMNI_TIMEOUT_SECONDS: {os.getenv('BIOMNI_TIMEOUT_SECONDS')}")

try:
    from biomni.llm import get_llm
    from langchain_core.messages import HumanMessage

    print("\n✅ biomni 模块加载成功")

    # 创建 LLM 实例
    print("\n🔧 创建 DeepSeek-Reasoner LLM 实例...")
    llm = get_llm(
        model="deepseek-reasoner",
        source="Custom",
        base_url="https://api.deepseek.com",
        api_key="sk-6f73c67f11d5469e846aba019b0f3530",
        temperature=0.7
    )
    print(f"  模型: deepseek-reasoner")
    print(f"  API: https://api.deepseek.com")

    # 测试简单调用
    print("\n💬 测试 API 调用（可能需要一些时间）...")
    import time
    start = time.time()
    
    test_message = [HumanMessage(content="用通俗易懂的语言解释CRISPR基因编辑的原理")]
    response = llm.invoke(test_message)
    
    elapsed = time.time() - start
    
    print(f"\n⏱️  响应时间: {elapsed:.1f}秒")
    print("\n📤 DeepSeek-Reasoner 回复:")
    print("-" * 70)
    print(response.content)
    print("-" * 70)

    print("\n" + "=" * 70)
    print("✅ 配置测试成功！")
    print("=" * 70)

    print("\n📖 使用示例:")
    print("""
from biomni.agent import A1

# 使用 DeepSeek-Reasoner 初始化 Biomni Agent
agent = A1(
    path='./data',
    llm='deepseek-reasoner',    # 深度推理模型
    source='Custom',
    base_url='https://api.deepseek.com',
    api_key='sk-6f73c67f11d5469e846aba019b0f3530',
    expected_data_lake_files=[],  # 跳过数据湖下载
    timeout_seconds=1200          # 20分钟超时
)

# 执行复杂生物医学任务
agent.go("设计一个CRISPR筛选实验来识别调节T细胞耗竭的关键基因")
""")

except Exception as e:
    print(f"\n❌ 错误: {e}")
    import traceback
    traceback.print_exc()
