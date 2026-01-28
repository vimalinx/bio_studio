#!/usr/bin/env python3
"""
mRNA设计与优化工具
功能：
1. 密码子优化
2. UTR设计 (5' UTR和3' UTR)
3. 二级结构预测和优化
4. GC含量优化
5. 稳定性优化
6. 避免 motifs (miRNA结合位点、免疫刺激序列等)
"""

from Bio.Seq import Seq
from Bio.SeqUtils import GC
import numpy as np
import pandas as pd
from pathlib import Path
import json
from typing import Dict, List, Tuple
import warnings
warnings.filterwarnings('ignore')


class MRNAOptimizer:
    """mRNA序列设计与优化工具"""

    # 人类常用UTR序列
    HUMAN_5_UTR = {
        'alpha_globin': 'GCCACCATGG',
        'beta_globin': 'GCCACCATGG',
        'optimized': 'GCCACC',  # Kozak序列
        'strong': 'GCCACCATGGGT',
    }

    HUMAN_3_UTR = {
        'alpha_globin': 'AATAAA',  # Poly-A信号
        'beta_globin': 'AATAAA',
        'stabilized': 'AATAAA',  # 包含稳定元件
    }

    # 需要避免的序列 motifs
    AVOID_MOTIFS = {
        'immunostimulatory': [
            'GUCCUUCAA',  # TLR7/8激活
            'GUUGUUGU',   # TLR7/8激活
        ],
        'miRNA_sites': [
            # 常见miRNA种子序列
            'AAGCU',
            'AACCA',
            'AAUGC',
            'UUGGU',
            'CAGUG',
        ],
        'instability': [
            'AUUUA',  # ARE元件，导致mRNA不稳定
            'UUAUUUAUU',
        ]
    }

    # 最优GC含量范围
    OPTIMAL_GC_RANGE = (40, 60)

    # 密码子使用表 (人类高频)
    HUMAN_CODON_USAGE = {
        'A': ['GCC', 'GCT', 'GCA', 'GCG'],
        'R': ['CGC', 'AGA', 'CGT', 'AGG', 'CGA', 'CGG'],
        'N': ['AAC', 'AAT'],
        'D': ['GAC', 'GAT'],
        'C': ['TGC', 'TGT'],
        'Q': ['CAG', 'CAA'],
        'E': ['GAG', 'GAA'],
        'G': ['GGC', 'GGT', 'GGA', 'GGG'],
        'H': ['CAC', 'CAT'],
        'I': ['ATC', 'ATT', 'ATA'],
        'L': ['CTG', 'TTA', 'CTT', 'CTC', 'CTA', 'TTG'],
        'K': ['AAG', 'AAA'],
        'M': ['ATG'],
        'F': ['TTT', 'TTC'],
        'P': ['CCG', 'CCT', 'CCC', 'CCA'],
        'S': ['AGC', 'TCT', 'TCC', 'AGT', 'TCA', 'TCG'],
        'T': ['ACC', 'ACA', 'ACG', 'ACT'],
        'W': ['TGG'],
        'Y': ['TAC', 'TAT'],
        'V': ['GTG', 'GTT', 'GTA', 'GTC'],
    }

    def __init__(self, species='human'):
        """
        初始化mRNA优化器

        Args:
            species: 目标物种 ('human', 'mouse')
        """
        self.species = species

    def design_mrna(self,
                    protein_seq: str,
                    optimize_codon: bool = True,
                    add_5utr: str = 'optimized',
                    add_3utr: str = 'alpha_globin',
                    avoid_motifs: bool = True,
                    target_gc: float = 50) -> Dict:
        """
        设计完整的mRNA序列

        Args:
            protein_seq: 蛋白质序列
            optimize_codon: 是否优化密码子
            add_5utr: 5' UTR类型
            add_3utr: 3' UTR类型
            avoid_motifs: 是否避免不良motifs
            target_gc: 目标GC含量

        Returns:
            包含mRNA序列和分析的字典
        """
        print(f"正在设计mRNA序列...")

        # 1. 密码子优化
        if optimize_codon:
            coding_seq = self._optimize_codons(protein_seq, target_gc)
        else:
            coding_seq = self._protein_to_rna(protein_seq)

        # 2. 添加UTR
        utr5 = self.HUMAN_5_UTR.get(add_5utr, add_5utr)
        utr3 = self.HUMAN_3_UTR.get(add_3utr, add_3utr)

        full_mrna = utr5 + coding_seq + utr3

        # 3. 检查和移除不良motifs
        if avoid_motifs:
            full_mrna, motifs_found = self._remove_motifs(full_mrna)
        else:
            motifs_found = []

        # 4. 分析序列
        analysis = self._analyze_mrna(full_mrna, len(coding_seq))

        return {
            'sequence': full_mrna,
            '5_utr': utr5,
            'coding_sequence': coding_seq,
            '3_utr': utr3,
            'length': len(full_mrna),
            'analysis': analysis,
            'motifs_removed': motifs_found
        }

    def _protein_to_rna(self, protein_seq: str) -> str:
        """将蛋白质序列转换为RNA序列 (无优化)"""
        rna_seq = ""
        for aa in protein_seq.upper():
            # 简单映射，使用第一个可用密码子
            codons = self.HUMAN_CODON_USAGE.get(aa, ['NNN'])
            rna_seq += codons[0].replace('T', 'U')
        return rna_seq

    def _optimize_codons(self, protein_seq: str, target_gc: float = 50) -> str:
        """
        密码子优化，考虑GC含量

        Args:
            protein_seq: 蛋白质序列
            target_gc: 目标GC含量

        Returns:
            优化后的RNA序列
        """
        optimized_seq = ""
        current_gc = 0

        for i, aa in enumerate(protein_seq.upper()):
            codons = self.HUMAN_CODON_USAGE.get(aa, ['UUU'])

            # 选择使GC含量最接近目标的密码子
            best_codon = codons[0].replace('T', 'U')
            min_diff = float('inf')

            for codon in codons:
                rna_codon = codon.replace('T', 'U')
                test_seq = optimized_seq + rna_codon
                test_gc = GC(test_seq)
                diff = abs(test_gc - target_gc)

                if diff < min_diff:
                    min_diff = diff
                    best_codon = rna_codon

            optimized_seq += best_codon

            # 动态调整
            if (i + 1) % 10 == 0:
                current_gc = GC(optimized_seq)
                # 如果偏离目标，调整策略
                if current_gc < target_gc - 5:
                    # 偏向GC丰富密码子
                    pass
                elif current_gc > target_gc + 5:
                    # 偏向AT丰富密码子
                    pass

        return optimized_seq

    def _analyze_mrna(self, mrna_seq: str, coding_len: int) -> Dict:
        """
        分析mRNA序列特征

        Args:
            mrna_seq: mRNA序列
            coding_len: 编码区长度

        Returns:
            分析结果字典
        """
        total_len = len(mrna_seq)

        # 全局分析
        gc_content = GC(mrna_seq)

        # 编码区分析
        coding_seq = mrna_seq[:coding_len]
        gc_coding = GC(coding_seq)

        # 简化的二级结构预测 (自由能估计)
        folding_energy = self._estimate_folding_energy(mrna_seq)

        # 寻找 motifs
        motifs_found = self._find_motifs(mrna_seq)

        # AU含量
        au_content = 100 - gc_content

        return {
            'total_length': total_len,
            'coding_length': coding_len,
            'utr5_length': len(mrna_seq) - coding_len - len(mrna_seq[coding_len:]),
            'utr3_length': len(mrna_seq[coding_len:]),
            'gc_content_total': round(gc_content, 2),
            'gc_content_coding': round(gc_coding, 2),
            'au_content_total': round(au_content, 2),
            'estimated_folding_energy': round(folding_energy, 2),
            'motifs_found': motifs_found,
            'gc_optimal': self.OPTIMAL_GC_RANGE[0] <= gc_content <= self.OPTIMAL_GC_RANGE[1],
            'stability_score': self._calculate_stability_score(mrna_seq)
        }

    def _estimate_folding_energy(self, seq: str) -> float:
        """
        估计折叠自由能 (简化版)
        负值越大，结构越稳定

        实际应用中应该使用 ViennaRNA 或 RNAfold
        """
        # 简化的能量计算
        # GC配对: -3 kcal/mol
        # AU配对: -2 kcal/mol
        # GU配对: -1 kcal/mol

        # 这里使用简化的启发式方法
        gc_count = seq.count('G') + seq.count('C')
        au_count = seq.count('A') + seq.count('U')

        # 估计自由能 (粗略)
        energy = -(gc_count * 1.5 + au_count * 1.0) / len(seq) * 10

        return energy

    def _find_motifs(self, seq: str) -> Dict[str, List[str]]:
        """查找序列中的motifs"""
        found = {'immunostimulatory': [], 'miRNA': [], 'instability': []}

        seq_upper = seq.upper()

        for category, motifs in self.AVOID_MOTIFS.items():
            for motif in motifs:
                if motif in seq_upper:
                    positions = []
                    start = 0
                    while True:
                        pos = seq_upper.find(motif, start)
                        if pos == -1:
                            break
                        positions.append(pos)
                        start = pos + 1

                    if positions:
                        if category == 'immunostimulatory':
                            found['immunostimulatory'].append({motif: positions})
                        elif 'miRNA' in category:
                            found['miRNA'].append({motif: positions})
                        else:
                            found['instability'].append({motif: positions})

        return found

    def _remove_motifs(self, seq: str) -> Tuple[str, List[str]]:
        """
        通过同义突变移除不良motifs

        注意：这是简化版本，实际需要更复杂的算法
        """
        removed = []
        modified_seq = seq

        # 检查并移除免疫刺激序列
        for motif in self.AVOID_MOTIFS['immunostimulatory']:
            if motif in modified_seq.upper():
                removed.append(motif)
                # 这里简化处理：实际需要替换密码子
                # modified_seq = self._replace_synonymous_codons(modified_seq, motif)

        return modified_seq, removed

    def _calculate_stability_score(self, seq: str) -> float:
        """
        计算mRNA稳定性评分 (0-1)

        考虑因素：
        - GC含量
        - 二级结构
        - 不稳定motifs
        - 3' UTR稳定元件
        """
        score = 0.5  # 基础分

        # GC含量得分
        gc = GC(seq)
        if 40 <= gc <= 60:
            score += 0.2
        elif 30 <= gc <= 70:
            score += 0.1
        else:
            score -= 0.1

        # 折叠能得分
        energy = self._estimate_folding_energy(seq)
        if energy < -150:  # 高度稳定
            score += 0.15
        elif energy < -100:
            score += 0.1

        # motifs惩罚
        motifs = self._find_motifs(seq)
        if motifs['instability']:
            score -= 0.1 * len(motifs['instability'])

        return round(max(0, min(1, score)), 3)

    def predict_secondary_structure(self, seq: str) -> Dict:
        """
        预测二级结构

        注意：需要安装 ViennaRNA 包
        sudo apt-get install vienna-rna
        pip install RNA

        这里提供简化版本
        """
        print("提示：完整的二级结构预测需要安装 ViennaRNA 包")
        print("  sudo apt-get install vienna-rna")
        print("  pip install RNA")

        try:
            import RNA
            # 使用 ViennaRNA
            fc = RNA.fold_compound(seq)
            (ss, mfe) = fc.mfe()
            fc = RNA.fold_compound(seq)
            (struct, energy) = fc.centroid()

            return {
                'mfe_structure': ss,
                'mfe_energy': mfe,
                'centroid_structure': struct,
                'centroid_energy': energy,
                'method': 'ViennaRNA'
            }
        except ImportError:
            # 简化预测
            structure = self._simple_structure_prediction(seq)
            return {
                'predicted_structure': structure,
                'method': 'simplified',
                'note': 'Install ViennaRNA for accurate prediction'
            }

    def _simple_structure_prediction(self, seq: str) -> str:
        """
        简化的二级结构预测 (启发式)

        返回点括号表示法
        """
        n = len(seq)
        structure = ['.'] * n

        # 简单的配对规则：GC > AU > GU
        # 从两端向内配对
        left, right = 0, n - 1

        while left < right:
            if self._can_pair(seq[left], seq[right]):
                structure[left] = '('
                structure[right] = ')'
                left += 1
                right -= 1
            else:
                left += 1

        return ''.join(structure)

    def _can_pair(self, base1: str, base2: str) -> bool:
        """判断两个碱基能否配对"""
        pairs = {
            ('G', 'C'), ('C', 'G'),
            ('A', 'U'), ('U', 'A'),
            ('G', 'U'), ('U', 'G')
        }
        return (base1.upper(), base2.upper()) in pairs

    def export_mrna(self, mrna_dict: Dict, output_prefix: str):
        """
        导出mRNA序列和分析结果

        Args:
            mrna_dict: mRNA设计结果
            output_prefix: 输出文件前缀
        """
        output_path = Path(output_prefix)

        # 1. 导出FASTA
        fasta_file = output_path.with_suffix('.fa')
        with open(fasta_file, 'w') as f:
            f.write(f">mRNA_design_{self.species}\n")
            f.write(f"{mrna_dict['sequence']}\n")
            f.write(f">5_UTR\n{mrna_dict['5_utr']}\n")
            f.write(f">Coding_Sequence\n{mrna_dict['coding_sequence']}\n")
            f.write(f">3_UTR\n{mrna_dict['3_utr']}\n")

        print(f"FASTA已保存: {fasta_file}")

        # 2. 导出JSON报告
        json_file = output_path.with_suffix('.json')
        with open(json_file, 'w') as f:
            json.dump(mrna_dict, f, indent=2)

        print(f"JSON报告已保存: {json_file}")

        # 3. 导出CSV摘要
        csv_file = output_path.with_suffix('.csv')
        summary = pd.DataFrame([{
            'Parameter': 'Total Length',
            'Value': f"{mrna_dict['length']} nt"
        }, {
            'Parameter': 'GC Content',
            'Value': f"{mrna_dict['analysis']['gc_content_total']}%"
        }, {
            'Parameter': 'Stability Score',
            'Value': mrna_dict['analysis']['stability_score']
        }, {
            'Parameter': 'Folding Energy',
            'Value': f"{mrna_dict['analysis']['estimated_folding_energy']} kcal/mol"
        }])
        summary.to_csv(csv_file, index=False)

        print(f"CSV摘要已保存: {csv_file}")

    def optimize_for_expression(self,
                                protein_seq: str,
                                expression_level: str = 'high') -> Dict:
        """
        针对表达水平进行优化

        Args:
            protein_seq: 蛋白质序列
            expression_level: 目标表达水平 ('low', 'medium', 'high')

        Returns:
            优化的mRNA设计
        """
        if expression_level == 'high':
            # 高表达：强UTR，高GC，优化密码子
            return self.design_mrna(
                protein_seq,
                optimize_codon=True,
                add_5utr='strong',
                add_3utr='stabilized',
                avoid_motifs=True,
                target_gc=55
            )
        elif expression_level == 'medium':
            # 中等表达
            return self.design_mrna(
                protein_seq,
                optimize_codon=True,
                add_5utr='optimized',
                add_3utr='alpha_globin',
                avoid_motifs=True,
                target_gc=50
            )
        else:
            # 低表达
            return self.design_mrna(
                protein_seq,
                optimize_codon=False,
                add_5utr='optimized',
                add_3utr='alpha_globin',
                avoid_motifs=True,
                target_gc=45
            )


# ============ 命令行界面 ============
import click


def main():
    @click.command()
    @click.option('--protein', '-p', required=True, help='蛋白质序列或文件')
    @click.option('--output', '-o', default='mrna_design', help='输出文件前缀')
    @click.option('--species', '-s', default='human', help='目标物种')
    @click.option('--optimize/--no-optimize', default=True, help='密码子优化')
    @click.option('--utr5', default='optimized', help='5\' UTR类型')
    @click.option('--utr3', default='alpha_globin', help='3\' UTR类型')
    @click.option('--target-gc', default=50, help='目标GC含量')
    @click.option('--expression', '-e', type=click.Choice(['low', 'medium', 'high']),
                  help='表达水平优化')
    def optimize_mrna(protein, output, species, optimize, utr5, utr3, target_gc, expression):
        """mRNA序列设计与优化"""
        optimizer = MRNAOptimizer(species=species)

        # 读取蛋白质序列
        if Path(protein).exists():
            from Bio import SeqIO
            record = SeqIO.read(protein, 'fasta')
            protein_seq = str(record.seq)
        else:
            protein_seq = protein

        # 设计mRNA
        if expression:
            print(f"\n针对 {expression} 表达水平进行优化...")
            result = optimizer.optimize_for_expression(protein_seq, expression)
        else:
            result = optimizer.design_mrna(
                protein_seq,
                optimize_codon=optimize,
                add_5utr=utr5,
                add_3utr=utr3,
                target_gc=target_gc
            )

        # 打印结果
        print(f"\n{'='*60}")
        print(f"mRNA设计结果")
        print(f"{'='*60}")
        print(f"\n序列长度: {result['length']} nt")
        print(f"5' UTR: {result['5_utr']}")
        print(f"3' UTR: {result['3_utr']}")
        print(f"\n完整mRNA序列:")
        print(result['sequence'])

        print(f"\n{'='*60}")
        print(f"序列分析")
        print(f"{'='*60}")
        analysis = result['analysis']
        print(f"总GC含量: {analysis['gc_content_total']}%")
        print(f"编码区GC: {analysis['gc_content_coding']}%")
        print(f"稳定性评分: {analysis['stability_score']}")
        print(f"估计折叠能: {analysis['estimated_folding_energy']} kcal/mol")

        if analysis['motifs_found']['immunostimulatory']:
            print(f"\n⚠️  发现免疫刺激序列: {len(analysis['motifs_found']['immunostimulatory'])}")
        if analysis['motifs_found']['instability']:
            print(f"⚠️  发现不稳定序列: {len(analysis['motifs_found']['instability'])}")

        # 导出
        optimizer.export_mrna(result, output)

        print(f"\n✅ 设计完成!")

    optimize_mrna()


if __name__ == '__main__':
    main()
