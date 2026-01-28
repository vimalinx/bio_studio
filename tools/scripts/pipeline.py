#!/usr/bin/env python3
"""
完整工作流：从DNA设计到蛋白质验证和靶点分析
整合所有工具的端到端流程
"""

import sys
from pathlib import Path
import click
import json

# 添加scripts目录到路径
sys.path.insert(0, str(Path(__file__).parent))

from dna_design import DNADesigner
from mrna_optimize import MRNAOptimizer
from protein_predict import ProteinPredictor
from target_analysis import TargetAnalyzer


class BioDesignPipeline:
    """生物设计完整工作流"""

    def __init__(self, work_dir='./data'):
        """
        初始化工作流

        Args:
            work_dir: 工作目录
        """
        self.work_dir = Path(work_dir)
        self.results_dir = self.work_dir / 'results'
        self.results_dir.mkdir(parents=True, exist_ok=True)

        # 初始化各个工具
        print("初始化生物设计工具...")
        self.dna_designer = DNADesigner()
        self.mrna_optimizer = MRNAOptimizer()
        # 注意：蛋白质预测器会下载模型，比较耗时
        # self.protein_predictor = ProteinPredictor()
        self.protein_predictor = None
        self.target_analyzer = TargetAnalyzer(work_dir)
        print("✅ 工具初始化完成\n")

    def run_full_pipeline(self,
                         protein_sequence: str,
                         target_name: str = 'unknown_target',
                         species: str = 'E.coli',
                         optimize_for: str = 'high'):
        """
        运行完整流程

        Args:
            protein_sequence: 蛋白质序列
            target_name: 靶点名称
            species: 表达物种
            optimize_for: 优化目标 ('high', 'medium', 'low')

        Returns:
            完整的结果字典
        """
        print(f"\n{'='*70}")
        print(f"开始完整生物设计流程")
        print(f"{'='*70}")
        print(f"\n靶点名称: {target_name}")
        print(f"序列长度: {len(protein_sequence)} aa")
        print(f"表达宿主: {species}")
        print(f"优化水平: {optimize_for}\n")

        results = {
            'target_name': target_name,
            'input_sequence': protein_sequence,
            'steps': {}
        }

        # ============ 步骤1: 靶点分析 ============
        print("\n" + "="*70)
        print("步骤 1/5: 靶点分析")
        print("="*70)

        try:
            target_report = self.target_analyzer.analyze_target(
                target_id=target_name,
                id_type='gene_symbol'
            )

            if target_report:
                results['steps']['target_analysis'] = target_report
                print(f"\n✅ 靶点分析完成")
                print(f"   综合评分: {target_report['overall_score']['overall_score']}/100")
                print(f"   可成药性: {target_report['druggability']['level']}")
                print(f"   安全性: {target_report['safety']['risk_level']}")
            else:
                print(f"\n⚠️  靶点分析失败，继续其他步骤...")

        except Exception as e:
            print(f"\n⚠️  靶点分析出错: {e}")

        # ============ 步骤2: DNA设计 ============
        print("\n" + "="*70)
        print("步骤 2/5: DNA序列设计")
        print("="*70)

        try:
            # 蛋白质转DNA
            dna_sequence = self.dna_designer.protein_to_dna(
                protein_sequence,
                optimize=True
            )

            # 分析DNA序列
            dna_analysis = self.dna_designer.analyze_sequence(dna_sequence)

            # 设计引物
            primers = self.dna_designer.design_primers(dna_sequence)

            # 保存结果
            dna_output = self.results_dir / f"{target_name}_dna.fa"
            self.dna_designer.export_sequence(
                dna_sequence,
                dna_output,
                description=f'Designed for {species}'
            )

            results['steps']['dna_design'] = {
                'sequence': dna_sequence,
                'analysis': dna_analysis,
                'primers': primers,
                'output_file': str(dna_output)
            }

            print(f"\n✅ DNA设计完成")
            print(f"   序列长度: {len(dna_sequence)} bp")
            print(f"   GC含量: {dna_analysis['gc_content']}%")
            print(f"   CAI: {dna_analysis['codon_adaptation_index']}")
            if 'error' not in primers:
                print(f"   引物对: 已设计")
                print(f"   正向引物: {primers['forward_primer']['sequence']}")
                print(f"   反向引物: {primers['reverse_primer']['sequence']}")
                print(f"   产物大小: {primers['product_size']} bp")

        except Exception as e:
            print(f"\n❌ DNA设计失败: {e}")
            import traceback
            traceback.print_exc()
            return None

        # ============ 步骤3: mRNA优化 ============
        print("\n" + "="*70)
        print("步骤 3/5: mRNA序列优化")
        print("="*70)

        try:
            # mRNA设计
            mrna_result = self.mrna_optimizer.optimize_for_expression(
                protein_sequence,
                expression_level=optimize_for
            )

            # 保存结果
            mrna_output = self.results_dir / f"{target_name}_mrna"
            self.mrna_optimizer.export_mrna(mrna_result, mrna_output)

            results['steps']['mrna_design'] = mrna_result

            print(f"\n✅ mRNA设计完成")
            print(f"   序列长度: {mrna_result['length']} nt")
            print(f"   GC含量: {mrna_result['analysis']['gc_content_total']}%")
            print(f"   稳定性评分: {mrna_result['analysis']['stability_score']}")
            print(f"   折叠能: {mrna_result['analysis']['estimated_folding_energy']:.2f} kcal/mol")

        except Exception as e:
            print(f"\n❌ mRNA设计失败: {e}")
            import traceback
            traceback.print_exc()
            return None

        # ============ 步骤4: 蛋白质结构预测 ============
        print("\n" + "="*70)
        print("步骤 4/5: 蛋白质结构预测")
        print("="*70)

        if self.protein_predictor is None:
            print("\n⚠️  蛋白质预测器未初始化")
            print("   提示：首次使用需要下载ESM模型 (~3GB)")
            if click.confirm("是否现在初始化蛋白质预测器？"):
                try:
                    self.protein_predictor = ProteinPredictor()
                except Exception as e:
                    print(f"❌ 初始化失败: {e}")
                    print("   跳过蛋白质预测步骤...")
        else:
            try:
                # 结构预测
                protein_output = self.results_dir / f"{target_name}_protein"
                protein_result = self.protein_predictor.predict_structure(
                    protein_sequence,
                    str(protein_output)
                )

                # 溶解度预测
                solubility = self.protein_predictor.predict_solubility(protein_sequence)

                # 稳定性分析
                stability = self.protein_predictor.analyze_stability(protein_sequence)

                results['steps']['protein_analysis'] = {
                    'structure_file': str(protein_output / "predicted_structure.pdb"),
                    'solubility': solubility,
                    'stability': stability
                }

                print(f"\n✅ 蛋白质分析完成")
                print(f"   溶解度: {solubility['prediction']}")
                print(f"   稳定性: {stability['stability']}")
                print(f"   分子量: {stability['molecular_weight']:.2f} Da")
                print(f"   等电点: {stability['isoelectric_point']:.2f}")

            except Exception as e:
                print(f"\n❌ 蛋白质分析失败: {e}")
                import traceback
                traceback.print_exc()

        # ============ 步骤5: 生成完整报告 ============
        print("\n" + "="*70)
        print("步骤 5/5: 生成完整报告")
        print("="*70)

        try:
            # 保存JSON报告
            report_file = self.results_dir / f"{target_name}_full_report.json"
            with open(report_file, 'w', encoding='utf-8') as f:
                json.dump(results, f, indent=2, ensure_ascii=False)

            # 生成Markdown报告
            self._generate_markdown_report(results, target_name)

            # 生成HTML报告
            self._generate_html_report(results, target_name)

            print(f"\n✅ 报告生成完成")
            print(f"   JSON: {report_file}")
            print(f"   Markdown: {self.results_dir / f'{target_name}_report.md'}")
            print(f"   HTML: {self.results_dir / f'{target_name}_report.html'}")

        except Exception as e:
            print(f"\n❌ 报告生成失败: {e}")

        # 打印总结
        print("\n" + "="*70)
        print("流程完成!")
        print("="*70)
        print(f"\n所有结果已保存到: {self.results_dir}")
        print(f"\n主要输出文件:")
        print(f"  - DNA序列: {self.results_dir / f'{target_name}_dna.fa'}")
        print(f"  - mRNA序列: {self.results_dir / f'{target_name}_mrna.fa'}")
        print(f"  - 完整报告: {self.results_dir / f'{target_name}_report.html'}")

        return results

    def _generate_markdown_report(self, results: dict, target_name: str):
        """生成Markdown报告"""
        md_file = self.results_dir / f"{target_name}_report.md"

        with open(md_file, 'w', encoding='utf-8') as f:
            f.write(f"# 生物设计完整报告\n\n")
            f.write(f"**靶点**: {results['target_name']}\n")
            f.write(f"**序列长度**: {len(results['input_sequence'])} aa\n\n")

            # 靶点分析
            if 'target_analysis' in results['steps']:
                f.write(f"## 1. 靶点分析\n\n")
                ta = results['steps']['target_analysis']
                f.write(f"- **综合评分**: {ta['overall_score']['overall_score']}/100\n")
                f.write(f"- **等级**: {ta['overall_score']['grade']}\n")
                f.write(f"- **可成药性**: {ta['druggability']['level']}\n")
                f.write(f"- **安全性**: {ta['safety']['risk_level']}\n\n")

            # DNA设计
            if 'dna_design' in results['steps']:
                f.write(f"## 2. DNA设计\n\n")
                dna = results['steps']['dna_design']
                f.write(f"- **序列长度**: {len(dna['sequence'])} bp\n")
                f.write(f"- **GC含量**: {dna['analysis']['gc_content']}%\n")
                f.write(f"- **CAI**: {dna['analysis']['codon_adaptation_index']}\n")
                if 'error' not in dna['primers']:
                    f.write(f"\n### 引物序列\n\n")
                    f.write(f"- 正向: `{dna['primers']['forward_primer']['sequence']}`\n")
                    f.write(f"- 反向: `{dna['primers']['reverse_primer']['sequence']}`\n")
                f.write(f"\n")

            # mRNA设计
            if 'mrna_design' in results['steps']:
                f.write(f"## 3. mRNA设计\n\n")
                mrna = results['steps']['mrna_design']
                f.write(f"- **序列长度**: {mrna['length']} nt\n")
                f.write(f"- **GC含量**: {mrna['analysis']['gc_content_total']}%\n")
                f.write(f"- **稳定性评分**: {mrna['analysis']['stability_score']}\n\n")

            # 蛋白质分析
            if 'protein_analysis' in results['steps']:
                f.write(f"## 4. 蛋白质分析\n\n")
                prot = results['steps']['protein_analysis']
                if 'solubility' in prot:
                    f.write(f"- **溶解度**: {prot['solubility']['prediction']}\n")
                if 'stability' in prot:
                    f.write(f"- **稳定性**: {prot['stability']['stability']}\n")
                    f.write(f"- **分子量**: {prot['stability']['molecular_weight']:.2f} Da\n")
                    f.write(f"- **等电点**: {prot['stability']['isoelectric_point']:.2f}\n\n")

    def _generate_html_report(self, results: dict, target_name: str):
        """生成HTML报告"""
        html_file = self.results_dir / f"{target_name}_report.html"

        html = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>生物设计报告 - {target_name}</title>
    <style>
        body {{ font-family: 'Segoe UI', Arial, sans-serif; margin: 40px; background: #f5f5f5; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; padding: 30px; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; border-left: 4px solid #3498db; padding-left: 15px; }}
        .header {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; border-radius: 5px; margin-bottom: 30px; }}
        .metric {{ display: inline-block; margin: 15px 30px 15px 0; }}
        .metric-label {{ font-size: 14px; color: #7f8c8d; }}
        .metric-value {{ font-size: 32px; font-weight: bold; color: #2c3e50; }}
        .success {{ color: #27ae60; }}
        .warning {{ color: #f39c12; }}
        .error {{ color: #e74c3c; }}
        .sequence {{ background: #f8f9fa; padding: 15px; border-radius: 5px; font-family: monospace; font-size: 12px; overflow-x: auto; }}
        table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
        th {{ background: #34495e; color: white; padding: 12px; text-align: left; }}
        td {{ border-bottom: 1px solid #ddd; padding: 12px; }}
        tr:hover {{ background: #f5f5f5; }}
        .grade {{ font-size: 48px; font-weight: bold; }}
        .grade-A {{ color: #27ae60; }}
        .grade-B {{ color: #3498db; }}
        .grade-C {{ color: #f39c12; }}
        .grade-D {{ color: #e74c3c; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1 style="border: none; color: white;">生物设计完整报告</h1>
            <p style="font-size: 18px;">靶点: {results['target_name']} | 序列长度: {len(results['input_sequence'])} aa</p>
        </div>

"""

        # 靶点分析部分
        if 'target_analysis' in results['steps']:
            ta = results['steps']['target_analysis']
            grade_class = f"grade-{ta['overall_score']['grade'][0]}"
            html += f"""
        <h2>1. 靶点分析</h2>
        <div>
            <span class="metric">
                <div class="metric-value {grade_class}">{ta['overall_score']['overall_score']}</div>
                <div class="metric-label">综合评分</div>
            </span>
            <span class="metric">
                <div class="metric-value">{ta['overall_score']['grade']}</div>
                <div class="metric-label">等级</div>
            </span>
            <span class="metric">
                <div class="metric-value">{ta['druggability']['score']}</div>
                <div class="metric-label">可成药性</div>
            </span>
            <span class="metric">
                <div class="metric-value">{ta['safety']['risk_level']}</div>
                <div class="metric-label">安全性</div>
            </span>
        </div>
"""

        # DNA设计部分
        if 'dna_design' in results['steps']:
            dna = results['steps']['dna_design']
            html += f"""
        <h2>2. DNA设计</h2>
        <table>
            <tr><th>参数</th><th>值</th></tr>
            <tr><td>序列长度</td><td>{len(dna['sequence'])} bp</td></tr>
            <tr><td>GC含量</td><td>{dna['analysis']['gc_content']}%</td></tr>
            <tr><td>CAI</td><td>{dna['analysis']['codon_adaptation_index']}</td></tr>
            <tr><td>限制性位点</td><td>{len(dna['analysis']['restriction_sites'])} 个</td></tr>
        </table>
"""

        # mRNA设计部分
        if 'mrna_design' in results['steps']:
            mrna = results['steps']['mrna_design']
            html += f"""
        <h2>3. mRNA设计</h2>
        <table>
            <tr><th>参数</th><th>值</th></tr>
            <tr><td>序列长度</td><td>{mrna['length']} nt</td></tr>
            <tr><td>GC含量</td><td>{mrna['analysis']['gc_content_total']}%</td></tr>
            <tr><td>稳定性评分</td><td>{mrna['analysis']['stability_score']}</td></tr>
            <tr><td>折叠能</td><td>{mrna['analysis']['estimated_folding_energy']:.2f} kcal/mol</td></tr>
        </table>
"""

        # 蛋白质分析部分
        if 'protein_analysis' in results['steps']:
            prot = results['steps']['protein_analysis']
            if 'stability' in prot:
                html += f"""
        <h2>4. 蛋白质分析</h2>
        <table>
            <tr><th>参数</th><th>值</th></tr>
            <tr><td>分子量</td><td>{prot['stability']['molecular_weight']:.2f} Da</td></tr>
            <tr><td>等电点</td><td>{prot['stability']['isoelectric_point']:.2f}</td></tr>
            <tr><td>稳定性</td><td>{prot['stability']['stability']}</td></tr>
"""
                if 'solubility' in prot:
                    html += f"<tr><td>溶解度</td><td>{prot['solubility']['prediction']}</td></tr>"
                html += "</table>"

        html += """
    </div>
</body>
</html>
"""
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html)


# ============ 命令行界面 ============
@click.command()
@click.option('--sequence', '-s', required=True, help='蛋白质序列或FASTA文件')
@click.option('--target', '-t', default='unknown_target', help='靶点名称')
@click.option('--species', default='E.coli', help='表达宿主')
@click.option('--optimize', type=click.Choice(['low', 'medium', 'high']),
              default='high', help='表达优化水平')
@click.option('--work-dir', default='./data', help='工作目录')
def run_pipeline(sequence, target, species, optimize, work_dir):
    """
    运行完整生物设计流程

    示例:
        python pipeline.py -s "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
        -t "hemoglobin" --optimize high
    """
    pipeline = BioDesignPipeline(work_dir=work_dir)

    # 读取序列
    if Path(sequence).exists():
        from Bio import SeqIO
        record = SeqIO.read(sequence, 'fasta')
        protein_seq = str(record.seq)
        if target == 'unknown_target':
            target = record.id
    else:
        protein_seq = sequence

    # 运行流程
    results = pipeline.run_full_pipeline(
        protein_sequence=protein_seq,
        target_name=target,
        species=species,
        optimize_for=optimize
    )

    if results:
        print("\n✅ 流程执行成功!")
    else:
        print("\n❌ 流程执行失败")
        sys.exit(1)


if __name__ == '__main__':
    run_pipeline()
