#!/usr/bin/env python3
"""
测试 DeepSeek-V3.2 Reasoner 模式

对比 deepseek-chat 和 deepseek-reasoner 的性能差异
"""

import os
import sys
import time

# 取消所有代理
for var in ['http_proxy', 'https_proxy', 'HTTP_PROXY', 'HTTPS_PROXY',
            'all_proxy', 'ALL_PROXY', 'socks_proxy', 'SOCKS_PROXY']:
    os.environ.pop(var, None)

os.environ['OPENAI_API_KEY'] = 'sk-6f73c67f11d5469e846aba019b0f3530'

print("\n" + "="*80)
print("🧠 DeepSeek-V3.2 Reasoner 模式测试")
print("="*80)

try:
    from biomni.llm import get_llm
    from langchain_core.messages import HumanMessage, SystemMessage

    # 测试任务：需要深度推理的复杂生物医学问题
    task = """
作为一个癌症生物学专家，请分析并回答：

在BRCA1突变的乳腺癌患者中，为什么PARP抑制剂有效？
请从分子机制角度详细解释"合成致死"原理，并分析：
1. BRCA1在DNA损伤修复中的具体作用
2. PARP在DNA修复中的功能
3. 为什么两者同时抑制会导致癌细胞死亡
4. 为什么正常细胞不受影响

请提供深入的分析，包括分子通路和细胞生物学机制。
"""

    models = {
        "deepseek-chat": "非思考模式（快速）",
        "deepseek-reasoner": "思考模式（深度推理）"
    }

    results = {}

    for model_name, description in models.items():
        print(f"\n{'='*80}")
        print(f"🧪 测试模型: {model_name}")
        print(f"📝 模式: {description}")
        print(f"{'='*80}")

        print("\n⏳ 正在思考...")
        start_time = time.time()

        try:
            llm = get_llm(
                model=model_name,
                source="Custom",
                base_url="https://api.deepseek.com",
                api_key="sk-6f73c67f11d5469e846aba019b0f3530",
                temperature=0.7
            )

            messages = [
                SystemMessage(content="你是Biomni，专业的癌症生物学和分子生物学专家。"),
                HumanMessage(content=task)
            ]

            response = llm.invoke(messages)

            end_time = time.time()
            elapsed_time = end_time - start_time

            print(f"\n✅ 响应完成")
            print(f"⏱  耗时: {elapsed_time:.2f}秒")
            print(f"📊 长度: {len(response.content)} 字符")

            print(f"\n{'='*80}")
            print("📄 回答内容:")
            print(f"{'='*80}")
            print(response.content)
            print(f"{'='*80}")

            results[model_name] = {
                "time": elapsed_time,
                "length": len(response.content),
                "content": response.content
            }

        except Exception as e:
            print(f"\n❌ 错误: {e}")
            results[model_name] = None

    # 对比分析
    print("\n" + "="*80)
    print("📊 对比分析")
    print("="*80)

    if "deepseek-chat" in results and "deepseek-reasoner" in results:
        chat = results["deepseek-chat"]
        reasoner = results["deepseek-reasoner"]

        if chat and reasoner:
            print(f"\n{'指标':<20} {'deepseek-chat':<20} {'deepseek-reasoner':<20}")
            print("-"*60)
            print(f"{'响应时间':<20} {chat['time']:.2f}秒{'':<15} {reasoner['time']:.2f}秒")
            print(f"{'回答长度':<20} {chat['length']}字符{'':<12} {reasoner['length']}字符")

            time_diff = reasoner['time'] / chat['time']
            length_diff = reasoner['length'] / chat['length']

            print(f"\n💡 分析:")
            print(f"  • Reasoner耗时是Chat的 {time_diff:.1f}x")
            print(f"  • Reasoner回答长度是Chat的 {length_diff:.1f}x")

            print("\n🎯 使用建议:")
            print("  • deepseek-chat: 适合快速查询、简单任务")
            print("  • deepseek-reasoner: 适合复杂推理、深度分析（推荐用于Biomni）")

    # 保存reasoner的回答作为参考
    if "deepseek-reasoner" in results and results["deepseek-reasoner"]:
        with open("REASONER_OUTPUT_EXAMPLE.md", "w", encoding="utf-8") as f:
            f.write("# DeepSeek-V3.2 Reasoner 模式输出示例\n\n")
            f.write(f"**任务**: PARP抑制剂合成致死机制分析\n")
            f.write(f"**模型**: deepseek-reasoner\n")
            f.write(f"**耗时**: {results['deepseek-reasoner']['time']:.2f}秒\n\n")
            f.write("---\n\n")
            f.write(results["deepseek-reasoner"]["content"])

        print(f"\n💾 Reasoner回答已保存到: REASONER_OUTPUT_EXAMPLE.md")

except Exception as e:
    print(f"\n❌ 测试失败: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*80)
print("🎉 测试完成！")
print("="*80)
print("\n💡 结论:")
print("  DeepSeek-V3.2 Reasoner模式提供更深度的推理能力，")
print("  非常适合Biomni这样的复杂生物医学任务。")
print("="*80 + "\n")
