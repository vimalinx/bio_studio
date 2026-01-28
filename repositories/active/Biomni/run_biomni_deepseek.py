#!/usr/bin/env python3
"""
Biomni + DeepSeek-Reasoner 快速启动脚本
"""
import os

# 取消代理设置（避免代理冲突）
for var in ['http_proxy', 'https_proxy', 'HTTP_PROXY', 'HTTPS_PROXY', 
            'all_proxy', 'ALL_PROXY', 'no_proxy', 'NO_PROXY']:
    os.environ.pop(var, None)

# 加载配置
from dotenv import load_dotenv
load_dotenv("/media/vimalinx/Data/bio_studio/Biomni/.env", override=False)

from biomni.agent import A1

print("=" * 70)
print("🧬 Biomni + DeepSeek-Reasoner 生物医学 AI 助手")
print("=" * 70)
print("\n⚠️  使用 deepseek-reasoner 模型（推理更深，响应较慢）")

# 初始化 Agent
print("\n🔧 初始化 Agent...")
agent = A1(
    path='/media/vimalinx/Data/bio_studio/Biomni/data',
    llm='deepseek-reasoner',
    source='Custom',
    base_url='https://api.deepseek.com',
    api_key='sk-6f73c67f11d5469e846aba019b0f3530',
    expected_data_lake_files=[],  # 跳过数据湖下载，加快启动
    timeout_seconds=1200  # 20分钟超时
)

print("\n✅ Agent 已就绪！")
print("\n💡 输入你的生物医学问题，或输入 'quit' 退出")
print("-" * 70)

# 交互循环
while True:
    try:
        user_input = input("\n📝 请输入问题: ").strip()
        
        if not user_input:
            continue
            
        if user_input.lower() in ['quit', 'exit', 'q', '退出']:
            print("\n👋 再见！")
            break
        
        print(f"\n🤔 DeepSeek-Reasoner 正在深度思考...")
        print("   (这可能需要一些时间...)")
        print("-" * 70)
        
        agent.go(user_input)
        
        print("-" * 70)
        
    except KeyboardInterrupt:
        print("\n\n👋 再见！")
        break
    except Exception as e:
        print(f"\n❌ 错误: {e}")
        import traceback
        traceback.print_exc()
