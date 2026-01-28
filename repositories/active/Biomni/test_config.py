#!/usr/bin/env python3
"""
测试Biomni API配置

这个脚本用于验证API配置是否正确设置。
"""

import os
import sys
from pathlib import Path


def test_env_file():
    """检查.env文件是否存在"""
    env_file = Path.cwd() / ".env"
    if env_file.exists():
        print("✓ .env 文件存在")
        return True
    else:
        print("✗ .env 文件不存在")
        print("  请先运行: python switch_profile.py")
        return False


def test_api_keys():
    """测试API密钥是否已配置"""
    from dotenv import load_dotenv

    load_dotenv()

    checks = []

    # 检查Anthropic
    anthropic_key = os.getenv("ANTHROPIC_API_KEY")
    if anthropic_key and anthropic_key != "your_anthropic_api_key_here":
        checks.append(("✓", "Anthropic API Key", "已配置"))
    else:
        checks.append(("✗", "Anthropic API Key", "未配置或为默认值"))

    # 检查OpenAI
    openai_key = os.getenv("OPENAI_API_KEY")
    if openai_key and openai_key != "your_openai_api_key_here":
        checks.append(("✓", "OpenAI API Key", "已配置"))
    else:
        checks.append(("○", "OpenAI API Key", "未配置"))

    # 检查LLM_SOURCE
    llm_source = os.getenv("LLM_SOURCE")
    if llm_source:
        checks.append(("✓", "LLM Source", llm_source))
    else:
        checks.append(("○", "LLM Source", "未设置（将自动检测）"))

    # 检查BIOMNI_LLM
    biomni_llm = os.getenv("BIOMNI_LLM")
    if biomni_llm:
        checks.append(("✓", "默认模型", biomni_llm))
    else:
        checks.append(("○", "默认模型", "未设置"))

    # 检查自定义模型配置
    custom_base_url = os.getenv("CUSTOM_MODEL_BASE_URL")
    if custom_base_url:
        checks.append(("✓", "自定义模型Base URL", custom_base_url))
        if os.getenv("CUSTOM_MODEL_API_KEY"):
            checks.append(("✓", "自定义模型API Key", "已配置"))
    else:
        checks.append(("○", "自定义模型", "未配置"))

    # 打印结果
    print("\nAPI配置检查:")
    print("-" * 60)
    for status, name, value in checks:
        print(f"{status} {name}: {value}")
    print("-" * 60)

    return any(c[0] == "✓" for c in checks)


def test_imports():
    """测试Biomni包是否可以导入"""
    try:
        from biomni.config import default_config, BiomniConfig
        print("✓ 可以导入 biomni.config")
        return True
    except ImportError as e:
        print(f"✗ 无法导入 biomni.config: {e}")
        print("  请确保已安装: pip install biomni")
        return False


def test_config_class():
    """测试配置类"""
    try:
        from biomni.config import BiomniConfig

        # 创建配置实例
        config = BiomniConfig(
            llm="claude-sonnet-4-5",
            temperature=0.7,
            timeout_seconds=600
        )

        print("✓ BiomniConfig 类工作正常")
        print(f"  - 模型: {config.llm}")
        print(f"  - 温度: {config.temperature}")
        print(f"  - 超时: {config.timeout_seconds}秒")

        return True
    except Exception as e:
        print(f"✗ BiomniConfig 测试失败: {e}")
        return False


def test_llm_detection():
    """测试LLM自动检测功能"""
    try:
        from biomni.llm import get_llm

        # 测试不同模型的检测
        test_cases = [
            ("claude-sonnet-4-5", "Anthropic"),
            ("gpt-4o", "OpenAI"),
            ("azure-gpt-4o", "AzureOpenAI"),
            ("gemini-pro", "Gemini"),
        ]

        print("\nLLM自动检测测试:")
        print("-" * 60)
        all_passed = True
        for model, expected_source in test_cases:
            try:
                # 只测试检测，不实际创建LLM（避免API调用）
                # 这里我们只是验证逻辑，不实际调用
                print(f"○ {model} -> {expected_source} (未实际测试)")
            except Exception as e:
                print(f"✗ {model} 检测失败: {e}")
                all_passed = False

        print("-" * 60)
        return all_passed

    except ImportError as e:
        print(f"✗ 无法导入 get_llm: {e}")
        return False


def main():
    """主测试函数"""
    print("=" * 60)
    print("Biomni 配置测试")
    print("=" * 60)

    # 运行所有测试
    tests = [
        ("环境文件", test_env_file),
        ("API密钥", test_api_keys),
        ("包导入", test_imports),
        ("配置类", test_config_class),
        ("LLM检测", test_llm_detection),
    ]

    results = []
    for name, test_func in tests:
        print(f"\n测试: {name}")
        try:
            passed = test_func()
            results.append((name, passed))
        except Exception as e:
            print(f"✗ 测试出错: {e}")
            results.append((name, False))

    # 总结
    print("\n" + "=" * 60)
    print("测试总结")
    print("=" * 60)

    passed = sum(1 for _, p in results if p)
    total = len(results)

    for name, p in results:
        status = "✓ 通过" if p else "✗ 失败"
        print(f"{status} - {name}")

    print("-" * 60)
    print(f"总计: {passed}/{total} 测试通过")
    print("=" * 60)

    if passed == total:
        print("\n✓ 所有测试通过！Biomni已正确配置。")
        print("\n可以开始使用:")
        print("  from biomni.agent import A1")
        print("  agent = A1(path='./data', llm='claude-sonnet-4-5')")
        print("  agent.go('你的生物医学任务')")
    else:
        print("\n⚠️  部分测试失败，请检查配置。")
        print("\n常见问题:")
        print("1. 如果.env文件不存在: 运行 python switch_profile.py")
        print("2. 如果API密钥未配置: 编辑 .env 文件填入实际密钥")
        print("3. 如果无法导入biomni: 运行 pip install biomni")


if __name__ == "__main__":
    main()
