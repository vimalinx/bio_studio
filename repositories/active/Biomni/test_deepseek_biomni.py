#!/usr/bin/env python3
"""
测试Biomni与DeepSeek集成
"""

import os
import sys

# 设置环境变量
os.environ['OPENAI_API_KEY'] = 'sk-6f73c67f11d5469e846aba019b0f3530'
os.environ['LLM_SOURCE'] = 'Custom'
os.environ['BIOMNI_LLM'] = 'deepseek-chat'
os.environ['CUSTOM_MODEL_BASE_URL'] = 'https://api.deepseek.com'
os.environ['CUSTOM_MODEL_API_KEY'] = 'sk-6f73c67f11d5469e846aba019b0f3530'

print("=" * 60)
print("测试 Biomni + DeepSeek 集成")
print("=" * 60)

try:
    # 测试导入
    print("\n1. 导入Biomni模块...")
    from biomni.llm import get_llm

    print("   ✓ 导入成功")

    # 创建LLM实例
    print("\n2. 创建DeepSeek LLM实例...")
    llm = get_llm(
        model="deepseek-chat",
        source="Custom",
        base_url="https://api.deepseek.com",
        api_key="sk-6f73c67f11d5469e846aba019b0f3530",
        temperature=0.7
    )

    print("   ✓ LLM实例创建成功")

    # 测试简单的调用
    print("\n3. 测试LLM调用...")
    from langchain_core.messages import HumanMessage

    messages = [HumanMessage(content="请用中文回答：DeepSeek API与Biomni集成测试成功了吗？请简短回答。")]

    response = llm.invoke(messages)

    print("   ✓ LLM调用成功")
    print(f"\n响应内容:")
    print("-" * 60)
    print(response.content)
    print("-" * 60)

    print("\n" + "=" * 60)
    print("✓ Biomni + DeepSeek 集成测试成功！")
    print("=" * 60)

    print("\n下一步:")
    print("1. 切换到DeepSeek配置:")
    print("   python switch_profile.py switch deepseek")
    print("\n2. 使用Biomni:")
    print("   from biomni.agent import A1")
    print("   agent = A1(")
    print("       path='./data',")
    print("       llm='deepseek-chat',")
    print("       source='Custom',")
    print("       base_url='https://api.deepseek.com',")
    print("       api_key='sk-6f73c67f11d5469e846aba019b0f3530'")
    print("   )")
    print("\n3. 执行任务:")
    print("   agent.go('你的生物医学任务')")

except ImportError as e:
    print(f"\n✗ 导入错误: {e}")
    print("\n需要先安装Biomni:")
    print("   pip install biomni --upgrade")
    sys.exit(1)

except Exception as e:
    print(f"\n✗ 锋试失败: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
