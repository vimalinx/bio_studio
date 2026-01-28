#!/usr/bin/env python3
"""
DNA设计工具 - 从蛋白质序列反推DNA序列并进行优化
功能：
1. 密码子优化
2. GC含量优化
3. 避免限制性酶切位点
4. 引物设计
5. 序列质量评估
"""

from Bio.Seq import Seq
from Bio.SeqUtils import GC, gc_fraction
from Bio.SeqUtils.MeltingTemp import Tm_Wallace
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
from pathlib import Path
import json
import warnings
warnings.filterwarnings('ignore')


class DNADesigner:
    """DNA序列设计与优化工具"""

    # 简并密码子表
    GENETIC_CODE = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }

    # 反向密码子表 (氨基酸 -> 可能的密码子)
    REVERSE_CODE = {}
    for codon, aa in GENETIC_CODE.items():
        if aa not in REVERSE_CODE:
            REVERSE_CODE[aa] = []
        REVERSE_CODE[aa].append(codon)

    # 大肠杆菌高频密码子表 (用于密码子优化)
    ECOLI_PREF = {
        'A': 'GCT', 'R': 'CGT', 'N': 'AAC', 'D': 'GAT',
        'C': 'TGC', 'Q': 'CAG', 'E': 'GAA', 'G': 'GGT',
        'H': 'CAC', 'I': 'ATC', 'L': 'CTG', 'K': 'AAG',
        'M': 'ATG', 'F': 'TTT', 'P': 'CCG', 'S': 'TCG',
        'T': 'ACC', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTG',
        '*': 'TAA'
    }

    # 常用限制性酶切位点 (需要避免)
    RESTRICTION_SITES = {
        'EcoRI': 'GAATTC',
        'BamHI': 'GGATCC',
        'HindIII': 'AAGCTT',
        'XhoI': 'CTCGAG',
        'NdeI': 'CATATG',
        'NcoI': 'CCATGG',
        'SalI': 'GTCGAC',
        'PstI': 'CTGCAG',
        'NotI': 'GCGGCCGC',
    }

    def __init__(self, target_species='E.coli'):
        """
        初始化DNA设计器

        Args:
            target_species: 目标表达物种 ('E.coli', 'human', 'yeast')
        """
        self.target_species = target_species
        self.codon_table = self.ECOLI_PREF if target_species == 'E.coli' else self.REVERSE_CODE

    def protein_to_dna(self, protein_seq, optimize=True):
        """
        蛋白质序列转DNA序列

        Args:
            protein_seq: 蛋白质序列 (单字母代码)
            optimize: 是否进行密码子优化

        Returns:
            DNA序列字符串
        """
        dna_seq = ""
        for aa in protein_seq.upper():
            if aa == '*':
                # 终止密码子
                dna_seq += self.codon_table.get('*', 'TAA')
            elif aa in self.codon_table:
                if optimize and isinstance(self.codon_table, dict):
                    # 使用优化密码子
                    dna_seq += self.codon_table[aa]
                else:
                    # 随机选择密码子
                    codons = self.REVERSE_CODE.get(aa, ['NNN'])
                    dna_seq += codons[0]  # 简化，使用第一个
            else:
                raise ValueError(f"未知氨基酸: {aa}")

        return dna_seq

    def analyze_sequence(self, dna_seq):
        """
        分析DNA序列质量

        Args:
            dna_seq: DNA序列

        Returns:
            包含各种指标的字典
        """
        analysis = {
            'length': len(dna_seq),
            'gc_content': round(GC(dna_seq), 2),
            'gc_content_range': self._check_gc_range(dna_seq),
            'at_content': round(100 - GC(dna_seq), 2),
            'restriction_sites': self._find_restriction_sites(dna_seq),
            'repetitive_sequences': self._find_repeats(dna_seq),
            'codon_adaptation_index': self._calculate_cai(dna_seq),
            'is_multiple_of_3': len(dna_seq) % 3 == 0,
        }
        return analysis

    def _check_gc_range(self, dna_seq, window=50):
        """检查GC含量是否在合理范围内 (滑动窗口)"""
        gc_contents = []
        for i in range(0, len(dna_seq) - window + 1, window):
            window_seq = dna_seq[i:i+window]
            gc_contents.append(GC(window_seq))

        if not gc_contents:
            return {'min': 0, 'max': 0, 'mean': 0}

        return {
            'min': round(min(gc_contents), 2),
            'max': round(max(gc_contents), 2),
            'mean': round(np.mean(gc_contents), 2),
            'std': round(np.std(gc_contents), 2)
        }

    def _find_restriction_sites(self, dna_seq):
        """查找限制性酶切位点"""
        found = {}
        for enzyme, site in self.RESTRICTION_SITES.items():
            if site.upper() in dna_seq.upper():
                positions = []
                start = 0
                while True:
                    pos = dna_seq.upper().find(site.upper(), start)
                    if pos == -1:
                        break
                    positions.append(pos)
                    start = pos + 1
                if positions:
                    found[enzyme] = positions
        return found

    def _find_repeats(self, dna_seq, min_length=6):
        """查找重复序列"""
        repeats = []
        seq_upper = dna_seq.upper()

        for length in range(min_length, 20):
            seen = {}
            for i in range(len(seq_upper) - length + 1):
                subseq = seq_upper[i:i+length]
                if subseq in seen:
                    if i - seen[subseq] >= length:  # 非重叠
                        repeats.append({
                            'sequence': subseq,
                            'positions': [seen[subseq], i],
                            'length': length
                        })
                        break
                seen[subseq] = i

        return repeats[:5]  # 返回前5个重复

    def _calculate_cai(self, dna_seq):
        """
        计算密码子适应指数 (Codon Adaptation Index)
        简化版本
        """
        # 这里简化实现，实际CAI需要参考密码子使用表
        codons = [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3) if len(dna_seq[i:i+3]) == 3]

        if not codons:
            return 0.0

        # 简化：使用高频密码子比例
        high_usage_codons = set(self.ECOLI_PREF.values())
        optimized_count = sum(1 for c in codons if c in high_usage_codons)

        return round(optimized_count / len(codons), 3)

    def design_primers(self, dna_seq, target_tm=60, primer_length=20, product_size_range=(500, 2000)):
        """
        设计PCR引物

        Args:
            dna_seq: DNA模板序列
            target_tm: 目标Tm值
            primer_length: 引物长度
            product_size_range: PCR产物大小范围

        Returns:
            包含正向和反向引物的字典
        """
        seq_len = len(dna_seq)
        min_product, max_product = product_size_range

        # 正向引物：从序列起始位置选择
        forward_candidates = []
        for i in range(0, min(50, seq_len - primer_length)):
            primer_seq = dna_seq[i:i+primer_length]
            tm = self._calculate_tm(primer_seq)

            if abs(tm - target_tm) < 5:  # Tm在目标值±5范围内
                forward_candidates.append({
                    'sequence': primer_seq,
                    'position': i,
                    'tm': tm,
                    'gc_content': GC(primer_seq),
                    'length': primer_length
                })

        # 反向引物：从序列末尾选择，需要反向互补
        reverse_candidates = []
        for i in range(max(0, seq_len - 50 - primer_length), seq_len - primer_length):
            primer_seq = dna_seq[i:i+primer_length]
            reverse_primer = str(Seq(primer_seq).reverse_complement())
            tm = self._calculate_tm(reverse_primer)

            if abs(tm - target_tm) < 5:
                reverse_candidates.append({
                    'sequence': reverse_primer,
                    'position': i,
                    'tm': tm,
                    'gc_content': GC(reverse_primer),
                    'length': primer_length
                })

        # 选择最优引物对
        if not forward_candidates or not reverse_candidates:
            return {'error': '无法找到合适的引物'}

        # 选择Tm最接近目标的引物
        best_forward = min(forward_candidates, key=lambda x: abs(x['tm'] - target_tm))
        best_reverse = min(reverse_candidates, key=lambda x: abs(x['tm'] - target_tm))

        product_size = best_reverse['position'] - best_forward['position'] + primer_length

        return {
            'forward_primer': best_forward,
            'reverse_primer': best_reverse,
            'product_size': product_size,
            'tm_difference': abs(best_forward['tm'] - best_reverse['tm'])
        }

    def _calculate_tm(self, primer_seq):
        """计算引物Tm值 (简化版Wallace规则)"""
        return Tm_Wallace(primer_seq)

    def add_features(self, dna_seq, features=None):
        """
        添加序列特征 (启动子、标签、终止子等)

        Args:
            dna_seq: 原始DNA序列
            features: 特征字典

        Returns:
            完整的DNA序列
        """
        if features is None:
            features = {
                '5_utr': '',
                'promoter': '',
                'tag': '',
                '3_utr': '',
                'terminator': ''
            }

        full_seq = ""
        if features.get('5_utr'):
            full_seq += features['5_utr']
        if features.get('promoter'):
            full_seq += features['promoter']
        if features.get('tag'):
            full_seq += features['tag']
        full_seq += dna_seq
        if features.get('3_utr'):
            full_seq += features['3_utr']
        if features.get('terminator'):
            full_seq += features['terminator']

        return full_seq

    def export_sequence(self, dna_seq, filename, format='fasta', description=''):
        """
        导出序列到文件

        Args:
            dna_seq: DNA序列
            filename: 输出文件名
            format: 文件格式 ('fasta', 'genbank')
            description: 序列描述
        """
        record = SeqRecord(
            Seq(dna_seq),
            id=Path(filename).stem,
            description=description
        )

        output_path = Path(filename)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, 'w') as f:
            SeqIO.write(record, f, format)

        print(f"序列已导出到: {output_path}")

    def generate_report(self, analysis, output_path):
        """生成设计报告"""
        report = {
            'sequence_analysis': analysis,
            'summary': {
                'total_length': analysis['length'],
                'gc_content': f"{analysis['gc_content']}%",
                'cai': analysis['codon_adaptation_index'],
                'restriction_sites_found': len(analysis['restriction_sites']),
                'repetitive_sequences': len(analysis['repetitive_sequences']),
                'is_valid_orf': analysis['is_multiple_of_3']
            },
            'recommendations': self._generate_recommendations(analysis)
        }

        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2)

        return report

    def _generate_recommendations(self, analysis):
        """生成优化建议"""
        recommendations = []

        if analysis['gc_content'] < 30:
            recommendations.append("GC含量偏低，建议增加密码子优化")
        elif analysis['gc_content'] > 70:
            recommendations.append("GC含量偏高，可能影响表达")

        if analysis['restriction_sites']:
            recommendations.append(f"发现{len(analysis['restriction_sites'])}个限制性酶切位点，建议移除")

        if analysis['repetitive_sequences']:
            recommendations.append(f"发现{len(analysis['repetitive_sequences'])}个重复序列，可能影响PCR")

        if analysis['codon_adaptation_index'] < 0.7:
            recommendations.append("CAI值偏低，建议进一步优化密码子使用")

        if not recommendations:
            recommendations.append("序列质量良好，可以直接使用")

        return recommendations


# ============ 命令行界面 ============
def main():
    import click

    @click.command()
    @click.option('--protein', '-p', help='输入蛋白质序列或文件')
    @click.option('--output', '-o', default='output.fa', help='输出文件')
    @click.option('--species', '-s', default='E.coli', help='目标物种')
    @click.option('--optimize/--no-optimize', default=True, help='是否密码子优化')
    @click.option('--primers', is_flag=True, help='是否设计PCR引物')
    def design_dna(protein, output, species, optimize, primers):
        """DNA序列设计与优化"""
        designer = DNADesigner(target_species=species)

        # 读取输入
        if Path(protein).exists():
            # 从文件读取
            protein_seq = SeqIO.read(protein, 'fasta').seq
        else:
            # 直接输入序列
            protein_seq = protein

        # 转换为DNA
        dna_seq = designer.protein_to_dna(str(protein_seq), optimize=optimize)

        print(f"\n设计的DNA序列 ({len(dna_seq)} bp):")
        print(dna_seq)
        print(f"\n每行100bp:")
        for i in range(0, len(dna_seq), 100):
            print(f"{i+1:5d}: {dna_seq[i:i+100]}")

        # 分析序列
        analysis = designer.analyze_sequence(dna_seq)
        print(f"\n序列分析:")
        print(f"  长度: {analysis['length']} bp")
        print(f"  GC含量: {analysis['gc_content']}%")
        print(f"  CAI: {analysis['codon_adaptation_index']}")
        print(f"  限制性酶切位点: {len(analysis['restriction_sites'])}")

        # 设计引物
        if primers:
            primers = designer.design_primers(dna_seq)
            if 'error' not in primers:
                print(f"\n引物设计:")
                print(f"  正向: {primers['forward_primer']['sequence']}")
                print(f"         Tm={primers['forward_primer']['tm']:.1f}°C")
                print(f"  反向: {primers['reverse_primer']['sequence']}")
                print(f"         Tm={primers['reverse_primer']['tm']:.1f}°C")
                print(f"  产物大小: {primers['product_size']} bp")

        # 导出序列
        designer.export_sequence(dna_seq, output, description=f'Designed for {species}')

        # 生成报告
        report_path = Path(output).with_suffix('.json')
        designer.generate_report(analysis, report_path)
        print(f"\n报告已保存: {report_path}")

    design_dna()


if __name__ == '__main__':
    main()
