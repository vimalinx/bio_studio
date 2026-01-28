#!/usr/bin/env python3
"""
生物设计自动化演示 - 实际可运行版本
"""
import os

# 取消代理设置
for var in ['http_proxy', 'https_proxy', 'HTTP_PROXY', 'HTTPS_PROXY', 
            'all_proxy', 'ALL_PROXY', 'no_proxy', 'NO_PROXY']:
    os.environ.pop(var, None)

from dotenv import load_dotenv
load_dotenv("/media/vimalinx/Data/bio_studio/Biomni/.env", override=False)

from biomni.config import default_config
default_config.llm = "deepseek-chat"
default_config.source = "Custom"
default_config.base_url = "https://api.deepseek.com"
default_config.api_key = "sk-6f73c67f11d5469e846aba019b0f3530"

from biomni.llm import get_llm
from langchain_core.messages import HumanMessage, SystemMessage

print("=" * 80)
print("🧬 生物设计自动化演示")
print("=" * 80)

# 创建 LLM
llm = get_llm(
    model="deepseek-reasoner",
    source="Custom",
    base_url="https://api.deepseek.com",
    api_key="sk-6f73c67f11d5469e846aba019b0f3530",
    temperature=0.7
)

# 示例需求
user_request = "设计一个能特异性识别并杀伤肺癌细胞的 CAR-T 细胞治疗方案"

print(f"\n📝 用户需求: {user_request}")
print("\n" + "=" * 80)
print("🔍 DeepSeek-Reasoner 正在分析...")
print("=" * 80)

# 系统提示
system_prompt = """你是一个专业的生物医学设计专家。你的任务是将用户的自然语言需求转化为结构化的生物设计方案。

你需要考虑：
1. 蛋白质设计和选择
2. 分子作用机制
3. mRNA 和载体设计
4. 安全性和法规考虑

请提供详细、专业、可操作的设计方案。"""

# 阶段 1: 需求分析
print("\n📋 阶段 1: 需求分析")
print("-" * 80)

stage1_prompt = f"""
请分析以下生物医学需求，提取关键信息：

用户需求: {user_request}

请提取并输出：
1. 主要目标
2. 目标疾病/组织
3. 作用机制
4. 关键技术挑战
5. 建议的实施方案（至少2个）
"""

msg1 = [SystemMessage(content=system_prompt), HumanMessage(content=stage1_prompt)]
response1 = llm.invoke(msg1)
print(response1.content)

# 阶段 2: 蛋白质设计
print("\n" + "=" * 80)
print("🔗 阶段 2: 蛋白质组件设计")
print("-" * 80)

stage2_prompt = f"""
基于上述需求分析，设计 CAR-T 细胞的蛋白质组件。

请详细说明：
1. scFv（单链抗体）的选择 - 针对肺癌细胞表面抗原
2. 铰链区选择
3. 跨膜区选择（CD8α 或 CD28）
4. 共刺激结构域（4-1BB 或 CD28）
5. 激活结构域（CD3ζ）

对于每个组件，请提供：
- 推荐的蛋白名称
- UniProt ID
- 功能理由
- 备选方案
"""

msg2 = [SystemMessage(content=system_prompt), HumanMessage(content=stage2_prompt)]
response2 = llm.invoke(msg2)
print(response2.content)

# 阶段 3: mRNA 设计考虑
print("\n" + "=" * 80)
print("🧬 阶段 3: mRNA 设计考虑")
print("-" * 80)

stage3_prompt = """
针对上述 CAR-T 设计，说明 mRNA 序列设计的关键考虑：

1. 密码子优化策略
2. GC 含量控制
3. 二级结构避免
4. 5' 和 3' UTR 选择
5. 修饰核苷酸建议（如 N1-methyl-pseudouridine）
"""

msg3 = [SystemMessage(content=system_prompt), HumanMessage(content=stage3_prompt)]
response3 = llm.invoke(msg3)
print(response3.content)

# 阶段 4: 载体选择
print("\n" + "=" * 80)
print("📦 阶段 4: 载体系统选择")
print("-" * 80)

stage4_prompt = """
比较 CAR-T 细胞生产的不同载体系统：

1. 慢病毒载体
2. 逆转录病毒载体
3. 转座子系统（PiggyBac, Sleeping Beauty）
4. mRNA 电转

请从以下方面评估：
- 转导效率
- 整合风险
- 生产成本
- 临床应用成熟度
- 监管考虑

给出最终推荐和理由。
"""

msg4 = [SystemMessage(content=system_prompt), HumanMessage(content=stage4_prompt)]
response4 = llm.invoke(msg4)
print(response4.content)

# 阶段 5: 验证策略
print("\n" + "=" * 80)
print("✅ 阶段 5: 验证策略")
print("-" * 80)

stage5_prompt = """
制定完整的验证策略，包括：

1. 计算验证（预测）
   - 蛋白质结构预测
   - 免疫原性预测
   - 脱靶分析

2. 体外验证
   - 细胞毒性实验
   - 特异性检测
   - 细胞因子释放

3. 体内验证
   - 动物模型选择
   - 安全性评估
   - 药效学评估

请提供详细的实验设计。
"""

msg5 = [SystemMessage(content=system_prompt), HumanMessage(content=stage5_prompt)]
response5 = llm.invoke(msg5)
print(response5.content)

print("\n" + "=" * 80)
print("📊 演示完成！")
print("=" * 80)
print("\n💡 这只是一个演示框架。完整的系统需要：")
print("  - 集成专业数据库（UniProt, IEDB, Addgene）")
print("  - 添加序列分析工具（Biopython, ViennaRNA）")
print("  - 连接结构预测 API（ESMFold, AlphaFold）")
print("  - 实现完整的工作流自动化")
