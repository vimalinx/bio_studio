#!/usr/bin/env python3
"""
Test script for Biomni Agent (A1)
Task: Retrieve SARS-CoV-2 Spike protein info
"""

import sys
import os
import subprocess
from pathlib import Path

# 1. 动态定位 Bio Studio 根目录
try:
    root = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip()
    BIO_STUDIO_ROOT = Path(root)
except:
    # 回退到当前脚本的上上上级目录 (projects/test.../scripts/ -> root)
    BIO_STUDIO_ROOT = Path(__file__).resolve().parent.parent.parent.parent

print(f"📂 Bio Studio Root: {BIO_STUDIO_ROOT}")

# 2. 添加 Biomni 到路径
BIOMNI_PATH = BIO_STUDIO_ROOT / "repositories/active/Biomni"
if not BIOMNI_PATH.exists():
    print(f"❌ Biomni not found at {BIOMNI_PATH}")
    sys.exit(1)

sys.path.append(str(BIOMNI_PATH))
print(f"✅ Biomni path added: {BIOMNI_PATH}")

# 3. 导入 Biomni
try:
    from biomni.agent import A1
except ImportError as e:
    print(f"❌ Failed to import Biomni: {e}")
    # 尝试打印 sys.path 调试
    print("Python Path:", sys.path)
    sys.exit(1)

# 4. 初始化 Agent
print("🤖 Initializing Biomni Agent...")
# 注意：我们传入 expected_data_lake_files=[] 以跳过巨大的数据下载
# 假设我们只用它的推理和搜索能力
try:
    agent = A1(
        path=str(BIOMNI_PATH / 'data'), 
        llm='gpt-4o',  # 尝试使用 GPT-4o，或者根据 .env 配置自动选择
        expected_data_lake_files=[] # Skip large download
    )
except Exception as e:
    print(f"❌ Agent init failed: {e}")
    sys.exit(1)

# 5. 执行任务
task = "Find the amino acid sequence length of the SARS-CoV-2 Spike protein (Wuhan-Hu-1) and the range of its Receptor Binding Domain (RBD)."
print(f"\n🚀 Executing Task: {task}\n")

try:
    result = agent.go(task)
    print("\n✨ Result:")
    print(result)
except Exception as e:
    print(f"\n❌ Execution failed: {e}")
