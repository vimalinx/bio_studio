#!/usr/bin/env python3
"""
Bio Structure MCP Server
提供蛋白质结构分析功能的MCP服务器
"""

import asyncio
import logging
import json
from typing import Any, Literal
from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import Tool, TextContent
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.Polypeptide import PPBuilder
import os
import tempfile

# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# 创建MCP服务器实例
server = Server("bio-structure-mcp")

# ============================================================
# 工具函数定义
# ============================================================

def parse_pdb(pdb_content: str) -> dict:
    """解析PDB文件内容"""
    try:
        # 写入临时文件
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            temp_file = f.name

        parser = PDBParser()
        structure = parser.get_structure("protein", temp_file)

        # 清理临时文件
        os.unlink(temp_file)

        # 提取基本信息
        chains = []
        residues = 0
        atoms = 0

        for model in structure:
            for chain in model:
                chain_info = {"id": chain.id, "residues": []}
                for residue in chain:
                    if residue.has_id('CA'):
                        residues += 1
                        chain_info["residues"].append({
                            "name": residue.resname,
                            "number": residue.id[1],
                            "atoms": len(residue)
                        })
                    atoms += len(residue)
                chains.append(chain_info)

        return {
            "title": structure.header.get("name", "Unknown"),
            "chains": len(chains),
            "total_residues": residues,
            "total_atoms": atoms,
            "chains_detail": chains,
            "resolution": structure.header.get("resolution", "N/A"),
            "experiment_type": structure.header.get("structure_method", "N/A")
        }

    except Exception as e:
        return {"error": str(e)}

def analyze_secondary_structure(pdb_content: str) -> dict:
    """分析蛋白质二级结构"""
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            temp_file = f.name

        parser = PDBParser()
        structure = parser.get_structure("protein", temp_file)
        os.unlink(temp_file)

        ppb = PPBuilder()
        secondary_structure = {"helices": [], "sheets": [], "coils": []}

        for pp in ppb.build_peptides(structure):
            # 简化的二级结构分析
            seq = str(pp.get_sequence())
            secondary_structure["sequence"] = seq
            secondary_structure["length"] = len(seq)

        # DSSP分析（如果可用）
        try:
            from Bio.PDB.DSSP import DSSP
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
                f.write(pdb_content)
                temp_file = f.name

            dssp = DSSP(structure, temp_file)
            os.unlink(temp_file)

            ss_counts = {"H": 0, "E": 0, "C": 0}  # H: helix, E: sheet, C: coil
            for key in dssp:
                ss = dssp[key][2]
                if ss in ss_counts:
                    ss_counts[ss] += 1

            secondary_structure["dssp_analysis"] = {
                "helix_residues": ss_counts["H"],
                "sheet_residues": ss_counts["E"],
                "coil_residues": ss_counts["C"],
                "helix_percent": round(ss_counts["H"] / sum(ss_counts.values()) * 100, 2) if ss_counts else 0,
                "sheet_percent": round(ss_counts["E"] / sum(ss_counts.values()) * 100, 2) if ss_counts else 0,
                "coil_percent": round(ss_counts["C"] / sum(ss_counts.values()) * 100, 2) if ss_counts else 0
            }
        except:
            secondary_structure["dssp_analysis"] = "DSSP not available"

        return secondary_structure

    except Exception as e:
        return {"error": str(e)}

def calculate_structure_metrics(pdb_content: str) -> dict:
    """计算结构指标"""
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            temp_file = f.name

        parser = PDBParser()
        structure = parser.get_structure("protein", temp_file)
        os.unlink(temp_file)

        # 计算基本几何信息
        atoms = list(structure.get_atoms())
        coordinates = [atom.get_coord() for atom in atoms]

        if len(coordinates) > 0:
            import numpy as np
            coords_array = np.array(coordinates)

            # 质心
            centroid = np.mean(coords_array, axis=0)

            # 尺寸
            max_coords = np.max(coords_array, axis=0)
            min_coords = np.min(coords_array, axis=0)
            dimensions = max_coords - min_coords

            return {
                "num_atoms": len(atoms),
                "centroid": centroid.tolist(),
                "dimensions": {
                    "x": round(float(dimensions[0]), 2),
                    "y": round(float(dimensions[1]), 2),
                    "z": round(float(dimensions[2]), 2)
                },
                "max_coordinates": max_coords.tolist(),
                "min_coordinates": min_coords.tolist()
            }
        else:
            return {"error": "No atoms found in structure"}

    except Exception as e:
        return {"error": str(e)}

def extract_sequence(pdb_content: str) -> dict:
    """从PDB文件提取蛋白质序列"""
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            temp_file = f.name

        parser = PDBParser()
        structure = parser.get_structure("protein", temp_file)
        os.unlink(temp_file)

        ppb = PPBuilder()
        sequences = []

        for pp in ppb.build_peptides(structure):
            sequences.append({
                "chain_id": pp[0].get_parent().id,
                "sequence": str(pp.get_sequence()),
                "length": len(pp.get_sequence())
            })

        return {
            "total_sequences": len(sequences),
            "sequences": sequences
        }

    except Exception as e:
        return {"error": str(e)}

def predict_structure_esm(sequence: str) -> dict:
    """使用ESM预测蛋白质结构（需要ESM库）"""
    try:
        import esm
        import torch
        import subprocess
        import sys
 
        # 检查是否安装了依赖
        try:
            import openfold
        except ImportError:
            return {
                "error": "Missing dependency: openfold. Install with: pip install openfold",
                "alternative": "Use AlphaFold or upload PDB file directly",
                "debug_info": "openfold is required for ESM-Fold v1"
            }
 
        # 加载ESM模型
        try:
            model = esm.pretrained.esmfold_v1()
            model = model.eval().cuda() if torch.cuda.is_available() else model.eval()
        except Exception as e:
            return {
                "error": f"Failed to load ESM-Fold model: {str(e)}",
                "alternative": "Use AlphaFold or upload PDB file directly",
                "debug_info": f"Model loading error: {type(e).__name__}"
            }
 
        # 检查序列长度
        if len(sequence) > 700:
            return {
                "error": f"Sequence too long. Maximum length is 700, got {len(sequence)}",
                "alternative": "Use shorter sequence (< 700 aa) or AlphaFold",
                "debug_info": "ESM-Fold v1 has a maximum sequence length limit"
            }
 
        # 预测结构
        try:
            with torch.no_grad():
                output = model.infer_pdb(sequence)
 
            # 返回PDB格式
            return {
                "status": "success",
                "sequence": sequence,
                "sequence_length": len(sequence),
                "pdb_content": output,
                "pdb_length": len(output),
                "model": "ESM-Fold v1",
                "note": "Structure predicted successfully. Use the pdb_content to save or analyze."
            }
 
        except RuntimeError as e:
            error_msg = str(e)
            if "out of memory" in error_msg.lower():
                return {
                    "error": f"GPU memory error: {error_msg}",
                    "alternative": "Use CPU mode or try a shorter sequence",
                    "debug_info": "Reduce sequence length or use CPU (slower but requires less memory)"
                }
            else:
                return {
                    "error": f"Runtime error: {error_msg}",
                    "alternative": "Use AlphaFold or upload PDB file directly",
                    "debug_info": f"Error type: {type(e).__name__}"
                }
        except Exception as e:
            return {
                "error": f"Prediction error: {str(e)}",
                "alternative": "Use AlphaFold or upload PDB file directly",
                "debug_info": f"Error type: {type(e).__name__}"
            }
 
    except ImportError as e:
        return {
            "error": f"ESM library not available: {str(e)}",
            "alternative": "Use AlphaFold or upload PDB file directly",
            "note": "Install with: pip install fair-esm openfold omegaconf"
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
            name="parse_pdb",
            description="解析PDB文件，提取蛋白质结构基本信息",
            inputSchema={
                "type": "object",
                "properties": {
                    "pdb_content": {
                        "type": "string",
                        "description": "PDB文件内容（纯文本）"
                    }
                },
                "required": ["pdb_content"]
            }
        ),
        Tool(
            name="analyze_secondary_structure",
            description="分析蛋白质二级结构（α螺旋、β折叠等）",
            inputSchema={
                "type": "object",
                "properties": {
                    "pdb_content": {
                        "type": "string",
                        "description": "PDB文件内容"
                    }
                },
                "required": ["pdb_content"]
            }
        ),
        Tool(
            name="calculate_structure_metrics",
            description="计算结构几何指标（质心、尺寸、原子数等）",
            inputSchema={
                "type": "object",
                "properties": {
                    "pdb_content": {
                        "type": "string",
                        "description": "PDB文件内容"
                    }
                },
                "required": ["pdb_content"]
            }
        ),
        Tool(
            name="extract_sequence",
            description="从PDB文件提取蛋白质序列",
            inputSchema={
                "type": "object",
                "properties": {
                    "pdb_content": {
                        "type": "string",
                        "description": "PDB文件内容"
                    }
                },
                "required": ["pdb_content"]
            }
        ),
        Tool(
            name="predict_structure_esm",
            description="使用ESM-Fold预测蛋白质3D结构",
            inputSchema={
                "type": "object",
                "properties": {
                    "sequence": {
                        "type": "string",
                        "description": "蛋白质序列（单字母代码）"
                    }
                },
                "required": ["sequence"]
            }
        )
    ]

@server.call_tool()
async def call_tool(name: str, arguments: Any) -> list[TextContent]:
    """执行工具调用"""
    logger.info(f"调用工具: {name}")

    try:
        if name == "parse_pdb":
            result = parse_pdb(arguments["pdb_content"])
        elif name == "analyze_secondary_structure":
            result = analyze_secondary_structure(arguments["pdb_content"])
        elif name == "calculate_structure_metrics":
            result = calculate_structure_metrics(arguments["pdb_content"])
        elif name == "extract_sequence":
            result = extract_sequence(arguments["pdb_content"])
        elif name == "predict_structure_esm":
            result = predict_structure_esm(arguments["sequence"])
        else:
            result = {"error": f"未知工具: {name}"}

        return [TextContent(
            type="text",
            text=json.dumps(result, indent=2, ensure_ascii=False)
        )]

    except Exception as e:
        logger.error(f"工具执行错误: {e}")
        return [TextContent(
            type="text",
            text=json.dumps({"error": str(e)}, indent=2)
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
