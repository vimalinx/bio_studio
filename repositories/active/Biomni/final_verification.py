#!/usr/bin/env python3
"""最终验证所有功能"""

import os
import sys

# 取消代理设置
for var in ['http_proxy', 'https_proxy', 'HTTP_PROXY', 'HTTPS_PROXY', 'all_proxy', 'ALL_PROXY']:
    os.environ.pop(var, None)

print("=" * 70)
print("🎉 DeepSeek + Biomni 最终验证")
print("=" * 70)

# 验证1: 检查文件
print("\n[1/5] 检查文件...")
files_to_check = [
    "profiles/deepseek.env",
    "test_deepseek.py",
    "test_deepseek_biomni.py",
    "example_deepseek_usage.py",
    "DEEPSEEK_GUIDE.md",
    "DEEPSEEK_SUCCESS.md",
    "switch_profile.py"
]

all_exist = True
for file in files_to_check:
    if os.path.exists(file):
        print(f"  ✓ {file}")
    else:
        print(f"  ✗ {file} 不存在")
        all_exist = False

if not all_exist:
    print("\n✗ 部分文件缺失")
    sys.exit(1)

# 验证2: 测试DeepSeek API
print("\n[2/5] 测试DeepSeek API连接...")
try:
    from openai import OpenAI

    client = OpenAI(
        api_key="sk-6f73c67f11d5469e846aba019b0f3530",
        base_url="https://api.deepseek.com"
    )

    response = client.chat.completions.create(
        model="deepseek-chat",
        messages=[{"role": "user", "content": "Hi"}],
        max_tokens=10
    )

    print("  ✓ DeepSeek API连接正常")
    print(f"  响应: {response.choices[0].message.content[:50]}...")
except Exception as e:
    print(f"  ✗ DeepSeek API测试失败: {e}")
    sys.exit(1)

# 验证3: 测试Biomni集成
print("\n[3/5] 测试Biomni集成...")
try:
    from biomni.llm import get_llm
    from langchain_core.messages import HumanMessage

    llm = get_llm(
        model="deepseek-chat",
        source="Custom",
        base_url="https://api.deepseek.com",
        api_key="sk-6f73c67f11d5469e846aba019b0f3530"
    )

    messages = [HumanMessage(content="Say 'OK' in Chinese.")]
    response = llm.invoke(messages)

    print("  ✓ Biomni + DeepSeek集成正常")
    print(f"  响应: {response.content}")
except Exception as e:
    print(f"  ✗ Biomni集成测试失败: {e}")
    sys.exit(1)

# 验证4: 测试配置系统
print("\n[4/5] 测试配置系统...")
try:
    # 检查配置文件
    if os.path.exists("profiles/deepseek.env"):
        with open("profiles/deepseek.env") as f:
            content = f.read()
            if "deepseek-chat" in content and "api.deepseek.com" in content:
                print("  ✓ DeepSeek配置文件正确")
            else:
                print("  ✗ DeepSeek配置文件不完整")
                sys.exit(1)
    else:
        print("  ✗ DeepSeek配置文件不存在")
        sys.exit(1)
except Exception as e:
    print(f"  ✗ 配置系统测试失败: {e}")
    sys.exit(1)

# 验证5: 列出所有配置
print("\n[5/5] 验证配置切换系统...")
try:
    import subprocess
    result = subprocess.run(
        ["python", "switch_profile.py", "list"],
        capture_output=True,
        text=True,
        timeout=10
    )

    if "deepseek" in result.stdout:
        print("  ✓ 配置系统包含DeepSeek")
        print("\n可用配置:")
        for line in result.stdout.split("\n"):
            if "  -" in line:
                print(f"    {line.strip()}")
    else:
        print("  ✗ 配置系统未找到DeepSeek")
        sys.exit(1)
except Exception as e:
    print(f"  ✗ 配置系统测试失败: {e}")
    sys.exit(1)

# 成功！
print("\n" + "=" * 70)
print("🎉 所有验证通过！")
print("=" * 70)

print("\n✅ DeepSeek + Biomni 已完全配置并可以使用！")

print("\n📋 快速开始:")
print("  1. 切换到DeepSeek配置:")
print("     python switch_profile.py switch deepseek")
print("\n  2. 使用Biomni:")
print("     from biomni.agent import A1")
print("     agent = A1(path='./data', llm='deepseek-chat', source='Custom',")
print("                base_url='https://api.deepseek.com',")
print("                api_key='sk-6f73c67f11d5469e846aba019b0f3530')")
print("     agent.go('你的生物医学任务')")

print("\n📚 文档:")
print("  - DEEPSEEK_GUIDE.md     - 完整使用指南")
print("  - DEEPSEEK_SUCCESS.md   - 总结文档")
print("  - API_CONFIG_GUIDE.md   - 配置系统指南")

print("\n🧪 测试脚本:")
print("  - python test_deepseek.py          - 测试API")
print("  - python test_deepseek_biomni.py   - 测试集成")
print("  - python example_deepseek_usage.py - 完整示例")

print("\n💡 提示:")
print("  - DeepSeek成本仅为其他模型的 1/10 到 1/20")
print("  - 适合开发、测试和大规模使用")
print("  - 中文支持优秀")

print("\n" + "=" * 70)
print("开始使用吧！ 🚀")
print("=" * 70)
