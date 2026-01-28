#!/usr/bin/env python3
"""
Bio Sequence MCP Server
提供序列分析功能的MCP服务器
"""

import asyncio
import logging
from typing import Any
from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import Tool, TextContent
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import io

# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# 创建MCP服务器实例
server = Server("bio-sequence-mcp")

# ============================================================
# 工具函数定义
# ============================================================

def analyze_dna(sequence: str) -> dict:
    """分析DNA序列"""
    try:
        seq = Seq(sequence)
        return {
            "sequence": str(seq),
            "length": len(seq),
            "gc_content": round(gc_fraction(seq) * 100, 2),
            "reverse_complement": str(seq.reverse_complement()),
            "translation": str(seq.translate()),
            "complement": str(seq.complement()),
            "transcribe": str(seq.transcribe())
        }
    except Exception as e:
        return {"error": str(e)}

def analyze_rna(sequence: str) -> dict:
    """分析RNA序列"""
    try:
        seq = Seq(sequence)
        return {
            "sequence": str(seq),
            "length": len(seq),
            "gc_content": round(gc_fraction(seq) * 100, 2),
            "translation": str(seq.translate()),
            "back_transcribe": str(seq.back_transcribe())
        }
    except Exception as e:
        return {"error": str(e)}

def analyze_protein(sequence: str) -> dict:
    """分析蛋白质序列"""
    try:
        seq = Seq(sequence)
        analysed_seq = ProteinAnalysis(str(seq))
        return {
            "sequence": str(seq),
            "length": len(seq),
            "molecular_weight": round(analysed_seq.molecular_weight(), 2),
            "isoelectric_point": round(analysed_seq.isoelectric_point(), 2),
            "instability_index": round(analysed_seq.instability_index(), 2),
            "amino_acids_composition": analysed_seq.get_amino_acids_percent(),
            "gravy": round(analysed_seq.gravy(), 2)  # 疏水性
        }
    except Exception as e:
        return {"error": str(e)}

def find_orfs(sequence: str, min_length: int = 100) -> list:
    """查找开放阅读框(ORF)"""
    try:
        seq = Seq(sequence)
        orfs = []

        # 正向链
        for frame in range(3):
            for i in range(frame, len(seq) - 2, 3):
                codon = str(seq[i:i+3])
                if codon == "ATG":  # 起始密码子
                    for j in range(i+3, len(seq) - 2, 3):
                        codon = str(seq[j:j+3])
                        if codon in ["TAA", "TAG", "TGA"]:  # 终止密码子
                            orf_length = j - i
                            if orf_length >= min_length:
                                orfs.append({
                                    "frame": f"+{frame+1}",
                                    "start": i,
                                    "end": j+3,
                                    "length": orf_length,
                                    "sequence": str(seq[i:j+3]),
                                    "protein": str(seq[i:j+3].translate())
                                })
                            break

        return {"orfs_found": len(orfs), "orfs": orfs[:10]}  # 限制返回前10个
    except Exception as e:
        return {"error": str(e)}

def translate_sequence(sequence: str, genetic_code: int = 1) -> dict:
    """翻译序列为蛋白质"""
    try:
        seq = Seq(sequence)
        return {
            "original": str(seq),
            "translated": str(seq.translate(table=genetic_code, to_stop=True)),
            "genetic_code": genetic_code,
            "length": len(seq.translate(table=genetic_code, to_stop=True))
        }
    except Exception as e:
        return {"error": str(e)}

def reverse_complement(sequence: str) -> dict:
    """生成反向互补序列"""
    try:
        seq = Seq(sequence)
        return {
            "original": str(seq),
            "reverse_complement": str(seq.reverse_complement()),
            "complement": str(seq.complement()),
            "reverse": str(seq[::-1])
        }
    except Exception as e:
        return {"error": str(e)}

# ============================================================
# MCP工具定义
# ============================================================

@server.list_tools()
async def list_tools() -> list[Tool]:
    """列出所有可用的工具"""
    return [
        Tool(
            name="analyze_dna",
            description="分析DNA序列 - GC含量、翻译、互补序列等",
            inputSchema={
                "type": "object",
                "properties": {
                    "sequence": {
                        "type": "string",
                        "description": "DNA序列 (A, T, G, C)"
                    }
                },
                "required": ["sequence"]
            }
        ),
        Tool(
            name="analyze_rna",
            description="分析RNA序列 - GC含量、翻译等",
            inputSchema={
                "type": "object",
                "properties": {
                    "sequence": {
                        "type": "string",
                        "description": "RNA序列 (A, U, G, C)"
                    }
                },
                "required": ["sequence"]
            }
        ),
        Tool(
            name="analyze_protein",
            description="分析蛋白质序列 - 分子量、等电点、疏水性等",
            inputSchema={
                "type": "object",
                "properties": {
                    "sequence": {
                        "type": "string",
                        "description": "蛋白质序列 (单字母氨基酸代码)"
                    }
                },
                "required": ["sequence"]
            }
        ),
        Tool(
            name="find_orfs",
            description="查找DNA序列中的开放阅读框(ORF)",
            inputSchema={
                "type": "object",
                "properties": {
                    "sequence": {
                        "type": "string",
                        "description": "DNA序列"
                    },
                    "min_length": {
                        "type": "integer",
                        "description": "最小ORF长度 (默认100)",
                        "default": 100
                    }
                },
                "required": ["sequence"]
            }
        ),
        Tool(
            name="translate_sequence",
            description="将DNA/RNA序列翻译为蛋白质",
            inputSchema={
                "type": "object",
                "properties": {
                    "sequence": {
                        "type": "string",
                        "description": "DNA或RNA序列"
                    },
                    "genetic_code": {
                        "type": "integer",
                        "description": "遗传密码表编号 (默认1 - 标准密码表)",
                        "default": 1
                    }
                },
                "required": ["sequence"]
            }
        ),
        Tool(
            name="reverse_complement",
            description="生成DNA序列的反向互补序列",
            inputSchema={
                "type": "object",
                "properties": {
                    "sequence": {
                        "type": "string",
                        "description": "DNA序列"
                    }
                },
                "required": ["sequence"]
            }
        )
    ]

@server.call_tool()
async def call_tool(name: str, arguments: Any) -> list[TextContent]:
    """执行工具调用"""
    logger.info(f"调用工具: {name} 参数: {arguments}")

    try:
        if name == "analyze_dna":
            result = analyze_dna(arguments["sequence"])
        elif name == "analyze_rna":
            result = analyze_rna(arguments["sequence"])
        elif name == "analyze_protein":
            result = analyze_protein(arguments["sequence"])
        elif name == "find_orfs":
            result = find_orfs(
                arguments["sequence"],
                arguments.get("min_length", 100)
            )
        elif name == "translate_sequence":
            result = translate_sequence(
                arguments["sequence"],
                arguments.get("genetic_code", 1)
            )
        elif name == "reverse_complement":
            result = reverse_complement(arguments["sequence"])
        else:
            result = {"error": f"未知工具: {name}"}

        import json
        return [TextContent(
            type="text",
            text=json.dumps(result, indent=2, ensure_ascii=False)
        )]

    except Exception as e:
        logger.error(f"工具执行错误: {e}")
        return [TextContent(
            type="text",
            text=str({"error": str(e)})
        )]

# ============================================================
# 主程序
# ============================================================

async def main():
    """启动MCP服务器"""
    async with stdio_server() as (read_stream, write_stream):
        await server.run(
            read_stream,
            write_stream,
            server.create_initialization_options()
        )

if __name__ == "__main__":
    asyncio.run(main())
