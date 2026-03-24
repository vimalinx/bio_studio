#!/usr/bin/env python3
"""
Bio Design MCP Server
提供 DNA/mRNA 设计与靶点评估能力的 MCP 服务器
"""

from __future__ import annotations

import asyncio
import importlib.util
import json
import logging
from pathlib import Path
from typing import Any

from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import TextContent, Tool


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

server = Server("bio-design-mcp")

ROOT = Path(__file__).resolve().parents[2]
SERVER_DIR = Path(__file__).resolve().parent
TOOLS_DIR = ROOT / "tools" / "scripts"
RUNTIME_DIR = SERVER_DIR / ".runtime"

SCRIPT_PATHS = {
    "DNADesigner": TOOLS_DIR / "dna_design.py",
    "MRNAOptimizer": TOOLS_DIR / "mrna_optimize.py",
    "TargetAnalyzer": TOOLS_DIR / "target_analysis.py",
}

_MODULE_CACHE: dict[str, Any] | None = None


def _load_symbol(path: Path, symbol: str) -> Any:
    spec = importlib.util.spec_from_file_location(f"bio_design_{path.stem}", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module spec for {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    try:
        return getattr(module, symbol)
    except AttributeError as exc:
        raise RuntimeError(f"{path.name} does not define {symbol}") from exc


def load_design_modules() -> dict[str, Any]:
    global _MODULE_CACHE
    if _MODULE_CACHE is None:
        _MODULE_CACHE = {
            symbol: _load_symbol(path, symbol)
            for symbol, path in SCRIPT_PATHS.items()
        }
    return _MODULE_CACHE


def get_server_capabilities() -> dict[str, Any]:
    modules: dict[str, Any] = {}
    for symbol, path in SCRIPT_PATHS.items():
        try:
            _load_symbol(path, symbol)
        except Exception as exc:
            modules[symbol] = {"status": "error", "script": str(path), "error": str(exc)}
        else:
            modules[symbol] = {"status": "ready", "script": str(path)}

    return {
        "server": "bio-design-mcp",
        "tools": ["design_dna", "design_mrna", "analyze_target"],
        "runtime_dir": str(RUNTIME_DIR),
        "modules": modules,
    }


def design_dna(
    protein_sequence: str,
    species: str = "E.coli",
    optimize: bool = True,
    include_primers: bool = False,
) -> dict[str, Any]:
    modules = load_design_modules()
    designer = modules["DNADesigner"](target_species=species)
    dna_sequence = designer.protein_to_dna(protein_sequence, optimize=optimize)
    payload = {
        "tool": "design_dna",
        "protein_sequence": protein_sequence,
        "species": species,
        "optimize": optimize,
        "dna_sequence": dna_sequence,
        "analysis": designer.analyze_sequence(dna_sequence),
    }
    if include_primers:
        payload["primers"] = designer.design_primers(dna_sequence)
    return payload


def design_mrna(
    protein_sequence: str,
    species: str = "human",
    expression_level: str = "high",
) -> dict[str, Any]:
    modules = load_design_modules()
    optimizer = modules["MRNAOptimizer"](species=species)
    result = optimizer.optimize_for_expression(
        protein_sequence, expression_level=expression_level
    )
    return {
        "tool": "design_mrna",
        "protein_sequence": protein_sequence,
        "species": species,
        "expression_level": expression_level,
        "mrna_sequence": result.get("sequence"),
        "length": result.get("length"),
        "analysis": result.get("analysis", {}),
        "design": result,
    }


def analyze_target(
    target_id: str,
    id_type: str = "gene_symbol",
    diseases: list[str] | None = None,
) -> dict[str, Any]:
    modules = load_design_modules()
    analyzer = modules["TargetAnalyzer"](work_dir=str(RUNTIME_DIR))
    report = analyzer.analyze_target(target_id, id_type=id_type, diseases=diseases)
    return {
        "tool": "analyze_target",
        "target_id": target_id,
        "id_type": id_type,
        "report": report,
    }


@server.list_tools()
async def list_tools() -> list[Tool]:
    return [
        Tool(
            name="design_dna",
            description="根据蛋白质序列生成 DNA 设计，并返回基础质量分析",
            inputSchema={
                "type": "object",
                "properties": {
                    "protein_sequence": {
                        "type": "string",
                        "description": "蛋白质序列（单字母氨基酸代码）",
                    },
                    "species": {
                        "type": "string",
                        "description": "目标表达物种，默认 E.coli",
                        "default": "E.coli",
                    },
                    "optimize": {
                        "type": "boolean",
                        "description": "是否进行密码子优化",
                        "default": True,
                    },
                    "include_primers": {
                        "type": "boolean",
                        "description": "是否附带 PCR 引物建议",
                        "default": False,
                    },
                },
                "required": ["protein_sequence"],
            },
        ),
        Tool(
            name="design_mrna",
            description="针对表达水平优化 mRNA 设计，返回序列与稳定性分析",
            inputSchema={
                "type": "object",
                "properties": {
                    "protein_sequence": {
                        "type": "string",
                        "description": "蛋白质序列（单字母氨基酸代码）",
                    },
                    "species": {
                        "type": "string",
                        "description": "目标物种，默认 human",
                        "default": "human",
                    },
                    "expression_level": {
                        "type": "string",
                        "description": "目标表达水平",
                        "default": "high",
                        "enum": ["low", "medium", "high"],
                    },
                },
                "required": ["protein_sequence"],
            },
        ),
        Tool(
            name="analyze_target",
            description="对候选靶点做可成药性、安全性与竞争格局的基础评估",
            inputSchema={
                "type": "object",
                "properties": {
                    "target_id": {
                        "type": "string",
                        "description": "靶点 ID，例如 TP53",
                    },
                    "id_type": {
                        "type": "string",
                        "description": "ID 类型",
                        "default": "gene_symbol",
                        "enum": ["gene_symbol", "entrez", "ensembl"],
                    },
                    "diseases": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "可选疾病上下文列表",
                    },
                },
                "required": ["target_id"],
            },
        ),
        Tool(
            name="get_server_capabilities",
            description="返回 bio-design-mcp 当前可用的设计模块与脚本路径",
            inputSchema={"type": "object", "properties": {}},
        ),
    ]


@server.call_tool()
async def call_tool(name: str, arguments: Any) -> list[TextContent]:
    logger.info("调用工具: %s", name)

    try:
        arguments = arguments or {}

        if name == "design_dna":
            result = design_dna(
                arguments["protein_sequence"],
                arguments.get("species", "E.coli"),
                arguments.get("optimize", True),
                arguments.get("include_primers", False),
            )
        elif name == "design_mrna":
            result = design_mrna(
                arguments["protein_sequence"],
                arguments.get("species", "human"),
                arguments.get("expression_level", "high"),
            )
        elif name == "analyze_target":
            result = analyze_target(
                arguments["target_id"],
                arguments.get("id_type", "gene_symbol"),
                arguments.get("diseases"),
            )
        elif name == "get_server_capabilities":
            result = get_server_capabilities()
        else:
            result = {"error": f"未知工具: {name}"}

        return [
            TextContent(
                type="text",
                text=json.dumps(result, indent=2, ensure_ascii=False),
            )
        ]
    except Exception as exc:
        logger.error("工具执行错误: %s", exc)
        return [
            TextContent(
                type="text",
                text=json.dumps({"error": str(exc)}, indent=2, ensure_ascii=False),
            )
        ]


async def main_async():
    RUNTIME_DIR.mkdir(parents=True, exist_ok=True)
    async with stdio_server() as (read_stream, write_stream):
        await server.run(
            read_stream,
            write_stream,
            server.create_initialization_options(),
        )


def main():
    asyncio.run(main_async())


if __name__ == "__main__":
    main()
