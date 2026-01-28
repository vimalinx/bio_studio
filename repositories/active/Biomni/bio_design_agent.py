#!/usr/bin/env python3
"""
生物设计自动化 Agent - 概念验证原型

从自然语言功能描述 → 蛋白质链 → mRNA 设计 → 载体筛选 → DNA 设计 → 验证
"""
import os
import json

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

from biomni.agent import A1

print("=" * 80)
print("🧬 生物设计自动化 Agent")
print("=" * 80)

class BioDesignAgent:
    """生物设计自动化 Agent"""
    
    def __init__(self):
        print("\n🔧 初始化 BioDesign Agent...")
        
        # 初始化 Biomni Agent with DeepSeek-Reasoner
        self.agent = A1(
            path='/media/vimalinx/Data/bio_studio/Biomni/data',
            llm='deepseek-reasoner',
            source='Custom',
            base_url='https://api.deepseek.com',
            api_key='sk-6f73c67f11d5469e846aba019b0f3530',
            expected_data_lake_files=[],
            timeout_seconds=1200
        )
        
        # 存储设计状态
        self.design_state = {
            'user_description': None,
            'functional_spec': None,
            'protein_chain': None,
            'mrna_design': None,
            'vector_selection': None,
            'dna_design': None,
            'validation_report': None
        }
        
        print("✅ BioDesign Agent 已就绪")
    
    def parse_requirements(self, description: str) -> dict:
        """阶段 1: 解析用户需求"""
        print("\n" + "=" * 80)
        print("📋 阶段 1: 需求分析")
        print("=" * 80)
        print(f"\n用户输入: {description}")
        
        prompt = f"""
请分析以下生物医学需求，提取关键信息：

用户需求: {description}

请提取并输出以下信息（JSON格式）:
{{
    "target": "主要目标（疾病/应用）",
    "mechanism": "作用机制",
    "target_tissue": "目标组织/细胞",
    "constraints": "约束条件",
    "delivery_method": "递送方式要求（如果有）",
    "suggested_approaches": ["建议方案1", "建议方案2", "建议方案3"]
}}
"""
        
        result = self._execute_with_prompt(prompt)
        self.design_state['user_description'] = description
        self.design_state['functional_spec'] = result
        
        return result
    
    def design_protein_chain(self, spec: dict) -> dict:
        """阶段 2: 设计蛋白质作用链"""
        print("\n" + "=" * 80)
        print("🔗 阶段 2: 蛋白质作用链设计")
        print("=" * 80)
        
        prompt = f"""
基于以下功能规范，设计蛋白质作用链：

功能规范:
- 目标: {spec.get('target', 'N/A')}
- 机制: {spec.get('mechanism', 'N/A')}
- 目标组织: {spec.get('target_tissue', 'N/A')}

请设计蛋白质作用链，包括：
1. 核心效应蛋白（实现主要功能）
2. 辅助蛋白（增强/调节功能）
3. 报告/标记蛋白（可选，用于追踪）

对于每个蛋白，提供：
- 蛋白名称和UniProt ID
- 功能描述
- 相互作用关系
- 推荐的优化方向

输出格式:
{{
    "approach": "选择的方案",
    "protein_chain": [
        {{
            "name": "蛋白名称",
            "uniprot_id": "UniProt ID",
            "function": "功能描述",
            "role": "核心/辅助/标记"
        }}
    ],
    "interactions": "蛋白间相互作用描述"
}}
"""
        
        result = self._execute_with_prompt(prompt)
        self.design_state['protein_chain'] = result
        return result
    
    def design_mrna(self, protein_chain: dict) -> dict:
        """阶段 3: mRNA 序列设计"""
        print("\n" + "=" * 80)
        print("🧬 阶段 3: mRNA 序列设计")
        print("=" * 80)
        
        prompt = f"""
基于以下蛋白质链，设计 mRNA 序列：

蛋白质信息: {json.dumps(protein_chain, ensure_ascii=False)}

请设计 mRNA 序列，考虑：
1. 密码子优化（人类细胞表达）
2. 5' UTR 和 3' UTR 选择
3. GC 含量优化
4. 二级结构避免
5. 修饰核苷酸建议

输出设计规范，包括：
- 每个蛋白的优化策略
- 推荐的 UTR
- 密码子适应性目标
- 特殊注意事项

（注意：不需要输出完整序列，只需要设计规范和策略）
"""
        
        result = self._execute_with_prompt(prompt)
        self.design_state['mrna_design'] = result
        return result
    
    def select_vector(self, mrna_design: dict, target_tissue: str) -> dict:
        """阶段 4: 载体筛选"""
        print("\n" + "=" * 80)
        print("📦 阶段 4: 载体筛选与评估")
        print("=" * 80)
        
        prompt = f"""
为以下设计筛选合适的递送载体：

目标组织: {target_tissue}
mRNA 设计信息: {json.dumps(mrna_design, ensure_ascii=False)[:500]}

请评估并推荐载体系统：
1. LNP（脂质纳米颗粒）
2. AAV（腺相关病毒）
3. 腺病毒
4. 其他新兴载体

对每个载体评估：
- 适用性评分（1-10）
- 优点
- 缺点
- 容量限制
- 免疫原性风险
- 临床应用状态

输出最终推荐及理由：
{{
    "recommended_vector": "推荐的载体",
    "alternatives": ["备选方案1", "备选方案2"],
    "rationale": "推荐理由",
    "specifications": "载体规格要求"
}}
"""
        
        result = self._execute_with_prompt(prompt)
        self.design_state['vector_selection'] = result
        return result
    
    def design_dna(self, vector_selection: dict) -> dict:
        """阶段 5: DNA 设计"""
        print("\n" + "=" * 80)
        print("🧪 阶段 5: DNA 载体设计")
        print("=" * 80)
        
        prompt = f"""
基于载体选择，设计完整的 DNA 载体：

载体信息: {json.dumps(vector_selection, ensure_ascii=False)}

请设计 DNA 载体，包括：
1. 载体骨架结构
2. 启动子选择（考虑组织特异性）
3. 多克隆位点
4. 选择标记
5. 增强子/调控元件
6. 终止信号
7. 克隆策略

输出设计概要和克隆方案：
{{
    "vector_map": "载体图谱描述",
    "key_components": ["关键元件1", "关键元件2"],
    "cloning_strategy": "克隆策略",
    "restriction_sites": "推荐酶切位点",
    "synthesis_plan": "基因合成计划"
}}
"""
        
        result = self._execute_with_prompt(prompt)
        self.design_state['dna_design'] = result
        return result
    
    def validate_design(self) -> dict:
        """阶段 6: 计算验证"""
        print("\n" + "=" * 80)
        print("✅ 阶段 6: 设计验证")
        print("=" * 80)
        
        prompt = f"""
对以下完整设计进行验证和风险评估：

设计摘要:
{json.dumps(self.design_state, ensure_ascii=False, indent=2)}

请进行以下验证：
1. 蛋白质结构和稳定性
2. mRNA 表达效率预测
3. 载体安全性评估
4. 免疫原性风险
5. 脱靶效应分析
6. 法规符合性

输出验证报告：
{{
    "overall_score": "综合评分 (1-10)",
    "risk_assessment": {{"high": [], "medium": [], "low": []}},
    "validation_steps": ["验证步骤1", "验证步骤2"],
    "experimental_validation": "建议的湿实验验证方案",
    "next_steps": "后续行动建议"
}}
"""
        
        result = self._execute_with_prompt(prompt)
        self.design_state['validation_report'] = result
        return result
    
    def _execute_with_prompt(self, prompt: str) -> str:
        """执行单步任务并返回结果"""
        # 这里简化处理，实际需要解析 agent.go() 的输出
        # 对于 PoC，我们直接返回模拟的结构化数据
        return {"status": "designed", "prompt": prompt}
    
    def design_from_description(self, description: str) -> dict:
        """完整工作流：从描述到设计"""
        print(f"\n{'=' * 80}")
        print(f"🚀 开始生物设计工作流")
        print(f"{'=' * 80}")
        print(f"\n原始需求: {description}")
        
        # 执行各阶段
        spec = self.parse_requirements(description)
        protein = self.design_protein_chain(spec)
        mrna = self.design_mrna(protein)
        vector = self.select_vector(mrna, spec.get('target_tissue', 'general'))
        dna = self.design_dna(vector)
        validation = self.validate_design()
        
        # 生成最终报告
        print("\n" + "=" * 80)
        print("📊 设计完成！生成最终报告")
        print("=" * 80)
        
        return self.design_state


# 演示用例
if __name__ == "__main__":
    designer = BioDesignAgent()
    
    # 示例需求
    example_descriptions = [
        "我想要一个能特异性识别并杀伤肺癌细胞的 CAR-T 细胞治疗方案",
        "设计一个 mRNA 疫苗，用于预防 COVID-19 变异株",
        "开发一个肝脏特异性基因治疗方案，用于治疗血友病 B"
    ]
    
    print("\n📝 示例需求:")
    for i, desc in enumerate(example_descriptions, 1):
        print(f"\n{i}. {desc}")
    
    print("\n💡 要运行完整工作流，请使用:")
    print("  designer.design_from_description('你的需求描述')")
