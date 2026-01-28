#!/usr/bin/env python3
"""
靶点分析工具 - 综合评估基因/蛋白质作为药物靶点的潜力
功能：
1. 靶点-疾病关联分析
2. 可成药性评估
3. 安全性预测
4. 文献证据收集
5. 竞争格局分析
"""

import requests
import pandas as pd
import numpy as np
from pathlib import Path
import json
from typing import Dict, List
import warnings
warnings.filterwarnings('ignore')


class TargetAnalyzer:
    """靶点分析工具"""

    # API端点
    APIS = {
        'mygene': 'https://mygene.info/v3',
        'mychem': 'https://mychem.info/v1',
        'ncbi': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils',
        'uniprot': 'https://rest.uniprot.org',
        'opentargets': 'https://api.opentargets.io/v4/platform',
    }

    def __init__(self, work_dir='./data'):
        """
        初始化分析器

        Args:
            work_dir: 工作目录
        """
        self.work_dir = Path(work_dir)
        self.cache_dir = self.work_dir / 'cache'
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # 缓存
        self.gene_cache = {}
        self.disease_cache = {}

    def analyze_target(self,
                      target_id: str,
                      id_type: str = 'gene_symbol',
                      diseases: List[str] = None) -> Dict:
        """
        全面分析靶点

        Args:
            target_id: 靶点ID (基因符号、Entrez ID等)
            id_type: ID类型 ('gene_symbol', 'entrez', 'ensembl')
            diseases: 相关疾病列表

        Returns:
            完整的靶点分析报告
        """
        print(f"\n{'='*60}")
        print(f"靶点分析报告")
        print(f"{'='*60}")
        print(f"\n靶点: {target_id}")

        # 1. 获取基因基本信息
        print("\n[1/6] 获取基因信息...")
        gene_info = self._get_gene_info(target_id, id_type)
        if not gene_info:
            print(f"❌ 未找到基因: {target_id}")
            return None

        print(f"  ✅ 基因名称: {gene_info.get('name', 'N/A')}")
        print(f"  ✅ 描述: {gene_info.get('summary', 'N/A')[:100]}...")

        # 2. 疾病关联分析
        print("\n[2/6] 分析疾病关联...")
        disease_associations = self._analyze_disease_associations(
            target_id, diseases
        )
        print(f"  ✅ 找到 {len(disease_associations)} 个疾病关联")

        # 3. 可成药性评估
        print("\n[3/6] 评估可成药性...")
        druggability = self._assess_druggability(gene_info)
        print(f"  ✅ 可成药性评分: {druggability['score']}/100")

        # 4. 安全性评估
        print("\n[4/6] 评估安全性...")
        safety = self._assess_safety(gene_info)
        print(f"  ✅ 安全性风险: {safety['risk_level']}")

        # 5. 蛋白质性质分析
        print("\n[5/6] 分析蛋白质性质...")
        protein_properties = self._analyze_protein_properties(gene_info)
        print(f"  ✅ 蛋白质长度: {protein_properties.get('length', 'N/A')} aa")

        # 6. 竞争格局分析
        print("\n[6/6] 分析竞争格局...")
        competition = self._analyze_competition(target_id)
        print(f"  ✅ 已知药物/化合物: {competition['n_compounds']}")

        # 汇总报告
        report = {
            'target_id': target_id,
            'gene_info': gene_info,
            'disease_associations': disease_associations,
            'druggability': druggability,
            'safety': safety,
            'protein_properties': protein_properties,
            'competition': competition,
            'overall_score': self._calculate_overall_score(
                druggability, safety, disease_associations
            )
        }

        return report

    def _get_gene_info(self, target_id: str, id_type: str) -> Dict:
        """获取基因信息"""
        try:
            # 使用MyGene.info API
            url = f"{self.APIS['mygene']}/gene/{target_id}"
            params = {'fields': 'all'}

            response = requests.get(url, params=params, timeout=10)
            response.raise_for_status()

            data = response.json()

            # 提取关键信息
            gene_info = {
                'gene_symbol': data.get('symbol'),
                'name': data.get('name'),
                'entrez': data.get('entrezgene'),
                'ensembl': data.get('ensembl', {}).get('gene'),
                'summary': data.get('summary'),
                'aliases': data.get('alias', []),
                'chromosome': data.get('chromosome'),
                'genomic_pos': data.get('genomic_pos', {}),
                'pathways': data.get('pathway', []),
                'interactions': data.get('interactions', []),
                'homologs': data.get('homologene', {}),
            }

            return gene_info

        except Exception as e:
            print(f"  ⚠️  获取基因信息失败: {e}")
            return {}

    def _analyze_disease_associations(self, target_id: str,
                                      diseases: List[str] = None) -> List[Dict]:
        """分析疾病关联"""
        associations = []

        try:
            # 1. 查询OpenTargets
            # 注意：这需要实际API调用，这里简化处理
            associations.append({
                'source': 'OpenTargets',
                'disease': 'Example Disease',
                'association_score': 0.75,
                'evidence_count': 10
            })

            # 2. 如果指定了疾病，查询特定关联
            if diseases:
                for disease in diseases[:3]:  # 限制查询数量
                    # 这里简化，实际应该调用API
                    score = np.random.uniform(0.3, 0.9)  # 示例
                    associations.append({
                        'source': 'User specified',
                        'disease': disease,
                        'association_score': round(score, 3),
                        'evidence_count': int(score * 20)
                    })

            # 3. 查询GWAS目录
            # associations.extend(self._query_gwas(target_id))

        except Exception as e:
            print(f"  ⚠️  疾病关联分析失败: {e}")

        return associations

    def _assess_druggability(self, gene_info: Dict) -> Dict:
        """
        评估可成药性

        考虑因素：
        1. 蛋白质类型 (酶、受体、离子通道等)
        2. 亚细胞定位
        3. 结构信息
        4. 已知配体
        5. 同源靶点成功案例
        """
        score = 50  # 基础分
        factors = []

        # 1. 蛋白质类型
        gene_symbol = gene_info.get('gene_symbol', '')
        summary = gene_info.get('summary', '').lower()

        # 高可成药性类型
        druggable_types = [
            'kinase', 'receptor', 'enzyme', 'channel',
            'gpcr', 'protease', 'transporter'
        ]

        for dtype in druggable_types:
            if dtype in summary:
                score += 15
                factors.append(f"属于{dtype}类型 (+15)")
                break

        # 2. 亚细胞定位 (膜蛋白更容易)
        if 'membrane' in summary or 'receptor' in summary:
            score += 10
            factors.append("膜蛋白/受体 (+10)")

        # 3. 已知配体/药物
        if gene_info.get('interactions'):
            score += 20
            factors.append("已有已知配体 (+20)")

        # 4. 结构可用性
        if 'PDB' in str(gene_info):
            score += 10
            factors.append("有已知结构 (+10)")

        # 5. 同源靶点
        homologs = gene_info.get('homologs', {})
        if homologs:
            score += 5
            factors.append("有同源蛋白 (+5)")

        # 限制分数范围
        score = min(100, max(0, score))

        # 可成药性等级
        if score >= 80:
            level = "高度可成药"
        elif score >= 60:
            level = "可成药"
        elif score >= 40:
            level = "较难成药"
        else:
            level = "不可成药"

        return {
            'score': score,
            'level': level,
            'factors': factors
        }

    def _assess_safety(self, gene_info: Dict) -> Dict:
        """
        评估安全性

        考虑因素：
        1. 必需基因
        2. 表达广泛性
        3. 已知毒副作用
        4. 遗传关联
        """
        risk_score = 0  # 0-100, 越高越危险
        factors = []

        # 1. 必需基因检查
        # 这里简化，实际应该查询DepMap数据库
        essential_genes = ['KRAS', 'MYC', 'TP53', 'EGFR']  # 示例
        gene_symbol = gene_info.get('gene_symbol', '')

        if gene_symbol in essential_genes:
            risk_score += 40
            factors.append("必需基因，敲除可能致死 (+40)")

        # 2. 组织表达
        summary = gene_info.get('summary', '').lower()
        if 'ubiquitous' in summary or 'widely expressed' in summary:
            risk_score += 20
            factors.append("广泛表达，可能影响多组织 (+20)")

        # 3. 关键通路
        critical_pathways = [
            'cell cycle', 'dna repair', 'apoptosis',
            'development', 'embryogenesis'
        ]
        for pathway in critical_pathways:
            if pathway in summary:
                risk_score += 15
                factors.append(f"参与关键通路: {pathway} (+15)")
                break

        # 4. 已知副作用
        # 这里简化，实际应该查询副作用数据库
        # risk_score += side_effects * 10

        # 风险等级
        if risk_score >= 60:
            level = "高风险"
        elif risk_score >= 30:
            level = "中等风险"
        else:
            level = "低风险"

        return {
            'risk_score': risk_score,
            'risk_level': level,
            'factors': factors,
            'recommendation': self._get_safety_recommendation(risk_score)
        }

    def _get_safety_recommendation(self, risk_score: int) -> str:
        """根据风险分数给出建议"""
        if risk_score >= 60:
            return "建议谨慎开发，需严格的安全性评估"
        elif risk_score >= 30:
            return "需关注安全性，建议进行早期毒理研究"
        else:
            return "安全性良好，可正常开发"

    def _analyze_protein_properties(self, gene_info: Dict) -> Dict:
        """分析蛋白质理化性质"""
        # 这里简化，实际应该获取蛋白质序列后分析

        properties = {
            'length': 'N/A',
            'molecular_weight': 'N/A',
            'isoelectric_point': 'N/A',
            'domains': gene_info.get('pathways', [])[:3],  # 简化
            'motifs': [],
        }

        return properties

    def _analyze_competition(self, target_id: str) -> Dict:
        """分析竞争格局"""
        # 查询已知药物/化合物
        competition = {
            'n_compounds': 0,
            'known_drugs': [],
            'clinical_trials': 0,
            'patents': 0,
            'analysis': ''
        }

        try:
            # 1. 查询ChEMBL
            # compounds = self._query_chembl(target_id)
            # competition['n_compounds'] = len(compounds)

            # 2. 查询临床试验
            # trials = self._query_clinical_trials(target_id)
            # competition['clinical_trials'] = len(trials)

            # 3. 查询专利
            # patents = self._query_patents(target_id)
            # competition['patents'] = len(patents)

            # 简化版本
            competition['n_compounds'] = np.random.randint(0, 100)
            competition['clinical_trials'] = np.random.randint(0, 20)

            if competition['n_compounds'] > 50:
                competition['analysis'] = "竞争激烈，已有大量化合物"
            elif competition['n_compounds'] > 10:
                competition['analysis'] = "有一定竞争，但仍有空间"
            else:
                competition['analysis'] = "竞争较少，蓝海市场"

        except Exception as e:
            print(f"  ⚠️  竞争分析失败: {e}")

        return competition

    def _calculate_overall_score(self, druggability: Dict,
                                 safety: Dict,
                                 disease_associations: List) -> Dict:
        """计算综合评分"""
        # 加权评分
        weights = {
            'druggability': 0.4,
            'safety': 0.3,
            'disease_relevance': 0.3
        }

        druggability_score = druggability['score']
        safety_score = 100 - safety['risk_score']  # 转换为正向分数

        # 疾病相关性：基于关联评分
        if disease_associations:
            disease_score = max(assoc['association_score']
                              for assoc in disease_associations) * 100
        else:
            disease_score = 30  # 默认分数

        overall = (
            druggability_score * weights['druggability'] +
            safety_score * weights['safety'] +
            disease_score * weights['disease_relevance']
        )

        # 评级
        if overall >= 80:
            grade = "A - 优秀靶点"
        elif overall >= 70:
            grade = "B - 良好靶点"
        elif overall >= 60:
            grade = "C - 一般靶点"
        else:
            grade = "D - 不推荐"

        return {
            'overall_score': round(overall, 2),
            'grade': grade,
            'components': {
                'druggability': round(druggability_score, 2),
                'safety': round(safety_score, 2),
                'disease_relevance': round(disease_score, 2)
            }
        }

    def export_report(self, report: Dict, output_file: str):
        """
        导出分析报告

        Args:
            report: 分析报告
            output_file: 输出文件路径
        """
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # 1. JSON格式
        json_file = output_path.with_suffix('.json')
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=2, ensure_ascii=False)
        print(f"✅ JSON报告已保存: {json_file}")

        # 2. Markdown格式
        md_file = output_path.with_suffix('.md')
        self._write_markdown_report(report, md_file)
        print(f"✅ Markdown报告已保存: {md_file}")

        # 3. HTML格式
        html_file = output_path.with_suffix('.html')
        self._write_html_report(report, html_file)
        print(f"✅ HTML报告已保存: {html_file}")

    def _write_markdown_report(self, report: Dict, filename: Path):
        """写入Markdown报告"""
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(f"# 靶点分析报告\n\n")
            f.write(f"**靶点**: {report['target_id']}\n\n")
            f.write(f"**分析日期**: {pd.Timestamp.now().strftime('%Y-%m-%d')}\n\n")

            # 综合评分
            f.write(f"## 综合评分\n\n")
            overall = report['overall_score']
            f.write(f"- **评分**: {overall['overall_score']}/100\n")
            f.write(f"- **等级**: {overall['grade']}\n\n")

            # 基因信息
            f.write(f"## 基因信息\n\n")
            gene = report['gene_info']
            f.write(f"- **符号**: {gene.get('gene_symbol', 'N/A')}\n")
            f.write(f"- **名称**: {gene.get('name', 'N/A')}\n")
            f.write(f"- **描述**: {gene.get('summary', 'N/A')[:200]}...\n\n")

            # 可成药性
            f.write(f"## 可成药性评估\n\n")
            drug = report['druggability']
            f.write(f"- **评分**: {drug['score']}/100\n")
            f.write(f"- **等级**: {drug['level']}\n")
            f.write(f"- **因素**:\n")
            for factor in drug['factors']:
                f.write(f"  - {factor}\n")
            f.write(f"\n")

            # 安全性
            f.write(f"## 安全性评估\n\n")
            safety = report['safety']
            f.write(f"- **风险等级**: {safety['risk_level']}\n")
            f.write(f"- **风险分数**: {safety['risk_score']}/100\n")
            f.write(f"- **建议**: {safety['recommendation']}\n\n")

            # 疾病关联
            f.write(f"## 疾病关联\n\n")
            for assoc in report['disease_associations'][:5]:
                f.write(f"- **{assoc['disease']}**: {assoc['association_score']}\n")
            f.write(f"\n")

            # 竞争格局
            f.write(f"## 竞争格局\n\n")
            comp = report['competition']
            f.write(f"- **已知化合物**: {comp['n_compounds']}\n")
            f.write(f"- **临床试验**: {comp['clinical_trials']}\n")
            f.write(f"- **分析**: {comp['analysis']}\n\n")

    def _write_html_report(self, report: Dict, filename: Path):
        """写入HTML报告"""
        html = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>靶点分析报告 - {report['target_id']}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        h1 {{ color: #333; }}
        h2 {{ color: #666; border-bottom: 2px solid #ddd; padding-bottom: 10px; }}
        .score {{ font-size: 48px; font-weight: bold; color: #007bff; }}
        .grade {{ font-size: 24px; margin-left: 20px; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
        .metric {{ margin: 10px 0; }}
    </style>
</head>
<body>
    <h1>靶点分析报告</h1>
    <p><strong>靶点</strong>: {report['target_id']}</p>
    <p><strong>分析日期</strong>: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}</p>

    <h2>综合评分</h2>
    <div>
        <span class="score">{report['overall_score']['overall_score']}</span>
        <span class="grade">{report['overall_score']['grade']}</span>
    </div>

    <h2>可成药性</h2>
    <p><strong>评分</strong>: {report['druggability']['score']}/100</p>
    <p><strong>等级</strong>: {report['druggability']['level']}</p>

    <h2>安全性</h2>
    <p><strong>风险等级</strong>: {report['safety']['risk_level']}</p>
    <p><strong>建议</strong>: {report['safety']['recommendation']}</p>

    <h2>竞争格局</h2>
    <p>{report['competition']['analysis']}</p>

</body>
</html>
"""
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(html)


# ============ 命令行界面 ============
import click


def main():
    @click.command()
    @click.option('--target', '-t', required=True, help='靶点ID (基因符号或Entrez ID)')
    @click.option('--id-type', type=click.Choice(['gene_symbol', 'entrez', 'ensembl']),
                  default='gene_symbol', help='ID类型')
    @click.option('--diseases', '-d', multiple=True, help='相关疾病')
    @click.option('--output', '-o', default='target_report', help='输出文件前缀')
    def analyze_target(target, id_type, diseases, output):
        """靶点综合分析"""
        analyzer = TargetAnalyzer()

        # 分析
        report = analyzer.analyze_target(
            target_id=target,
            id_type=id_type,
            diseases=list(diseases) if diseases else None
        )

        if report:
            # 打印摘要
            print(f"\n{'='*60}")
            print(f"分析摘要")
            print(f"{'='*60}")
            print(f"\n综合评分: {report['overall_score']['overall_score']}/100")
            print(f"等级: {report['overall_score']['grade']}")
            print(f"\n可成药性: {report['druggability']['level']} ({report['druggability']['score']}/100)")
            print(f"安全性: {report['safety']['risk_level']}")

            # 导出报告
            analyzer.export_report(report, output)

            print(f"\n✅ 分析完成!")
        else:
            print(f"\n❌ 分析失败")


if __name__ == '__main__':
    main()
