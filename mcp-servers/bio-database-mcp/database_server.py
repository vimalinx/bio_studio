#!/usr/bin/env python3
"""
Bio Database MCP Server
提供生物数据库查询功能的MCP服务器
"""

import asyncio
import logging
import json
import os
from typing import Any, Optional
from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import Tool, TextContent
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW
import time

# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# 创建MCP服务器实例
server = Server("bio-database-mcp")

ENTREZ_EMAIL_ENV_VARS = (
    "BIO_STUDIO_ENTREZ_EMAIL",
    "NCBI_EMAIL",
    "ENTREZ_EMAIL",
)
ENTREZ_API_KEY_ENV_VARS = (
    "BIO_STUDIO_NCBI_API_KEY",
    "NCBI_API_KEY",
    "ENTREZ_API_KEY",
)


def _resolve_env_var(names: tuple[str, ...]) -> tuple[str | None, str | None]:
    for name in names:
        value = os.environ.get(name, "").strip()
        if value:
            return value, name
    return None, None


def configure_entrez() -> dict[str, Any]:
    """从环境变量配置 NCBI Entrez 身份。"""
    email, email_source = _resolve_env_var(ENTREZ_EMAIL_ENV_VARS)
    api_key, api_key_source = _resolve_env_var(ENTREZ_API_KEY_ENV_VARS)

    Entrez.email = email
    Entrez.api_key = api_key

    if email:
        logger.debug("Configured Entrez email from %s", email_source)
    else:
        logger.warning(
            "NCBI Entrez email is not configured; set one of %s",
            ", ".join(ENTREZ_EMAIL_ENV_VARS),
        )

    return {
        "email": email,
        "email_source": email_source,
        "api_key_configured": bool(api_key),
        "api_key_source": api_key_source,
    }


def append_entrez_warning(payload: dict[str, Any], config: dict[str, Any]) -> dict[str, Any]:
    if not config.get("email"):
        payload["config_warning"] = (
            "NCBI Entrez email is not configured. "
            "Set BIO_STUDIO_ENTREZ_EMAIL or NCBI_EMAIL before using this server heavily."
        )
    return payload


def get_server_config_status() -> dict[str, Any]:
    """返回当前数据库 MCP 的配置状态，不暴露敏感值。"""
    config = configure_entrez()
    payload = {
        "server": "bio-database-mcp",
        "entrez_email_configured": bool(config.get("email")),
        "entrez_email_source": config.get("email_source"),
        "ncbi_api_key_configured": bool(config.get("api_key_configured")),
        "ncbi_api_key_source": config.get("api_key_source"),
        "recommended_env_vars": {
            "email": list(ENTREZ_EMAIL_ENV_VARS),
            "api_key": list(ENTREZ_API_KEY_ENV_VARS),
        },
        "notes": [
            "NCBI Entrez requests should include a contact email.",
            "NCBI API key is optional but raises the default rate limit.",
        ],
    }
    return append_entrez_warning(payload, config)


configure_entrez()

# ============================================================
# 工具函数定义
# ============================================================

def search_pubmed(query: str, max_results: int = 10) -> dict:
    """搜索PubMed文献数据库"""
    config = configure_entrez()
    try:
        # 搜索PubMed
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=max_results,
            sort="relevance"
        )
        record = Entrez.read(handle)
        handle.close()

        ids = record.get("IdList", [])
        total_count = int(record.get("Count", 0))

        if not ids:
            return append_entrez_warning({
                "query": query,
                "total_found": total_count,
                "results": [],
                "note": "No results found"
            }, config)

        # 获取详细信息
        handle = Entrez.efetch(
            db="pubmed",
            id=",".join(ids),
            rettype="medline",
            retmode="text"
        )
        records = handle.read()
        handle.close()

        # 简单解析（实际应用中需要更复杂的解析）
        articles = []
        for pubmed_id in ids[:max_results]:
            try:
                handle = Entrez.efetch(db="pubmed", id=pubmed_id, retmode="xml")
                article_data = Entrez.read(handle)[0]
                handle.close()

                # 提取关键信息
                article = {
                    "pubmed_id": pubmed_id,
                    "title": article_data.get("MedlineCitation", {})
                                              .get("Article", {})
                                              .get("ArticleTitle", "N/A"),
                    "journal": article_data.get("MedlineCitation", {})
                                              .get("Article", {})
                                              .get("Journal", {})
                                              .get("Title", "N/A"),
                    "year": article_data.get("MedlineCitation", {})
                                           .get("Article", {})
                                           .get("Journal", {})
                                           .get("JournalIssue", {})
                                           .get("PubDate", {})
                                           .get("Year", "N/A"),
                    "authors": [author.get("LastName", "") + " " + author.get("Initials", "")
                              for author in article_data.get("MedlineCitation", {})
                                                          .get("Article", {})
                                                          .get("AuthorList", [])[:5]]
                }
                articles.append(article)
            except:
                pass

        return append_entrez_warning({
            "query": query,
            "total_found": total_count,
            "returned": len(articles),
            "articles": articles
        }, config)

    except Exception as e:
        return append_entrez_warning({"error": str(e)}, config)

def search_nucleotide(query: str, max_results: int = 5) -> dict:
    """搜索NCBI核酸数据库"""
    config = configure_entrez()
    try:
        handle = Entrez.esearch(
            db="nucleotide",
            term=query,
            retmax=max_results
        )
        record = Entrez.read(handle)
        handle.close()

        ids = record.get("IdList", [])

        if not ids:
            return append_entrez_warning(
                {"query": query, "results": [], "note": "No results found"},
                config,
            )

        # 获取序列摘要
        handle = Entrez.esummary(db="nucleotide", id=",".join(ids))
        summaries = Entrez.read(handle)
        handle.close()

        results = []
        for i, uid in enumerate(ids):
            if i < len(summaries):
                summary = summaries[uid]
                results.append({
                    "accession": summary.get("AccessionVersion", uid),
                    "title": summary.get("Title", "N/A"),
                    "organism": summary.get("Organism", "N/A"),
                    "length": summary.get("Length", 0),
                    "gburl": summary.get("GiToTaxon", {}).get("GbUrl", "")
                })

        return append_entrez_warning({
            "query": query,
            "found": len(results),
            "sequences": results
        }, config)

    except Exception as e:
        return append_entrez_warning({"error": str(e)}, config)

def search_protein(query: str, max_results: int = 5) -> dict:
    """搜索UniProt/NCBI蛋白质数据库"""
    config = configure_entrez()
    try:
        handle = Entrez.esearch(
            db="protein",
            term=query,
            retmax=max_results
        )
        record = Entrez.read(handle)
        handle.close()

        ids = record.get("IdList", [])

        if not ids:
            return append_entrez_warning(
                {"query": query, "results": [], "note": "No results found"},
                config,
            )

        # 获取蛋白质摘要
        handle = Entrez.esummary(db="protein", id=",".join(ids))
        summaries = Entrez.read(handle)
        handle.close()

        results = []
        for i, uid in enumerate(ids):
            if i < len(summaries):
                summary = summaries[uid]
                results.append({
                    "accession": summary.get("AccessionVersion", uid),
                    "title": summary.get("Title", "N/A"),
                    "organism": summary.get("Organism", "N/A"),
                    "length": summary.get("Length", 0),
                })

        return append_entrez_warning({
            "query": query,
            "found": len(results),
            "proteins": results
        }, config)

    except Exception as e:
        return append_entrez_warning({"error": str(e)}, config)

def get_gene_info(gene_symbol: str) -> dict:
    """获取基因信息"""
    config = configure_entrez()
    try:
        # 搜索基因
        handle = Entrez.esearch(
            db="gene",
            term=f"{gene_symbol}[gene] AND human[organism]",
            retmax=1
        )
        record = Entrez.read(handle)
        handle.close()

        ids = record.get("IdList", [])
        if not ids:
            return append_entrez_warning({"error": f"Gene {gene_symbol} not found"}, config)

        # 获取详细信息
        gene_id = ids[0]
        handle = Entrez.esummary(db="gene", id=gene_id)
        summary = Entrez.read(handle)
        handle.close()

        gene_data = summary["DocumentSummarySet"]["DocumentSummary"][0]

        return append_entrez_warning({
            "gene_symbol": gene_data.get("name", "N/A"),
            "description": gene_data.get("description", "N/A"),
            "chromosome": gene_data.get("chromosome", "N/A"),
            "location": gene_data.get("map_location", "N/A"),
            "organism": gene_data.get("organism", {}).get("scientificname", "N/A"),
            "alias": gene_data.get("other_aliases", "").split(",") if gene_data.get("other_aliases") else [],
            "summary": gene_data.get("summary", "N/A")
        }, config)

    except Exception as e:
        return append_entrez_warning({"error": str(e)}, config)

def run_blast(sequence: str, program: str = "blastp", database: str = "nr") -> dict:
    """运行NCBI BLAST搜索（异步，可能需要等待）"""
    try:
        # 限制序列长度避免超时
        if len(sequence) > 1000:
            return {
                "error": "Sequence too long for online BLAST",
                "suggestion": "Use local BLAST or shorter sequence",
                "max_length": 1000
            }

        # 运行BLAST
        result_handle = NCBIWWW.qblast(
            program=program,
            database=database,
            sequence=sequence
        )

        # 解析结果
        blast_records = list(SeqIO.parse(result_handle, "blast-xml"))
        result_handle.close()

        if not blast_records:
            return {"sequence": sequence, "hits": []}

        # 提取前10个hit
        hits = []
        record = blast_records[0]

        for alignment in record.alignments[:10]:
            hit = {
                "accession": alignment.hit_id,
                "description": alignment.hit_def,
                "length": alignment.length,
                "e_value": min(hsp.expect for hsp in alignment.hsps),
                "identity": alignment.hsps[0].identities / alignment.length * 100 if alignment.hsps else 0
            }
            hits.append(hit)

        return {
            "program": program,
            "database": database,
            "query_length": len(sequence),
            "hits_found": len(hits),
            "hits": hits
        }

    except Exception as e:
        return {"error": str(e)}

def get_uniprot_info(uniprot_id: str) -> dict:
    """获取UniProt蛋白质信息"""
    try:
        import urllib.request
        import json

        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        with urllib.request.urlopen(url) as response:
            data = json.loads(response.read())

        # 提取关键信息
        info = {
            "accession": data.get("primaryAccession", uniprot_id),
            "name": data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "N/A"),
            "gene": data.get("genes", [{}])[0].get("geneName", {}).get("value", "N/A") if data.get("genes") else "N/A",
            "organism": data.get("organism", {}).get("scientificName", "N/A"),
            "length": data.get("sequence", {}).get("length", 0),
            "mass": data.get("sequence", {}).get("molWeight", 0),
            "function": data.get("comments", [{}])[0].get("texts", [{}])[0].get("value", "N/A") if data.get("comments") else "N/A",
            "subcellular_location": [loc.get("location", {}).get("value", "N/A")
                                     for loc in data.get("comments", [])
                                     if loc.get("commentType") == "SUBCELLULAR_LOCATION"]
        }

        return info

    except Exception as e:
        return {"error": str(e), "note": "Make sure UniProt ID is valid (e.g., P12345)"}

# ============================================================
# MCP工具定义
# ============================================================

@server.list_tools()
async def list_tools() -> list[Tool]:
    """列出所有可用的工具"""
    return [
        Tool(
            name="search_pubmed",
            description="搜索PubMed科学文献数据库",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "搜索关键词（支持布尔运算符: AND, OR, NOT）"
                    },
                    "max_results": {
                        "type": "integer",
                        "description": "最大返回结果数（默认10）",
                        "default": 10
                    }
                },
                "required": ["query"]
            }
        ),
        Tool(
            name="search_nucleotide",
            description="搜索NCBI核酸序列数据库（GenBank）",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "搜索关键词或序列ID"
                    },
                    "max_results": {
                        "type": "integer",
                        "description": "最大返回结果数（默认5）",
                        "default": 5
                    }
                },
                "required": ["query"]
            }
        ),
        Tool(
            name="search_protein",
            description="搜索蛋白质数据库（UniProt/NCBI）",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "蛋白质名称、ID或关键词"
                    },
                    "max_results": {
                        "type": "integer",
                        "description": "最大返回结果数（默认5）",
                        "default": 5
                    }
                },
                "required": ["query"]
            }
        ),
        Tool(
            name="get_gene_info",
            description="获取基因的详细信息",
            inputSchema={
                "type": "object",
                "properties": {
                    "gene_symbol": {
                        "type": "string",
                        "description": "基因符号（如 TP53, EGFR, BRCA1）"
                    }
                },
                "required": ["gene_symbol"]
            }
        ),
        Tool(
            name="run_blast",
            description="运行NCBI BLAST序列比对搜索",
            inputSchema={
                "type": "object",
                "properties": {
                    "sequence": {
                        "type": "string",
                        "description": "要搜索的序列"
                    },
                    "program": {
                        "type": "string",
                        "description": "BLAST程序: blastp (蛋白), blastn (核酸), blastx (翻译)",
                        "default": "blastp",
                        "enum": ["blastp", "blastn", "blastx", "tblastn"]
                    },
                    "database": {
                        "type": "string",
                        "description": "数据库: nr (非冗余蛋白), nt (核酸), refseq等",
                        "default": "nr"
                    }
                },
                "required": ["sequence"]
            }
        ),
        Tool(
            name="get_uniprot_info",
            description="获取UniProt蛋白质详细信息",
            inputSchema={
                "type": "object",
                "properties": {
                    "uniprot_id": {
                        "type": "string",
                        "description": "UniProt ID (如 P04637, P12345)"
                    }
                },
                "required": ["uniprot_id"]
            }
        ),
        Tool(
            name="get_server_config_status",
            description="查看当前数据库 MCP 的 Entrez / NCBI 配置状态（不返回敏感值）",
            inputSchema={
                "type": "object",
                "properties": {}
            }
        )
    ]

@server.call_tool()
async def call_tool(name: str, arguments: Any) -> list[TextContent]:
    """执行工具调用"""
    logger.info(f"调用工具: {name}")

    try:
        if name == "search_pubmed":
            result = search_pubmed(
                arguments["query"],
                arguments.get("max_results", 10)
            )
        elif name == "search_nucleotide":
            result = search_nucleotide(
                arguments["query"],
                arguments.get("max_results", 5)
            )
        elif name == "search_protein":
            result = search_protein(
                arguments["query"],
                arguments.get("max_results", 5)
            )
        elif name == "get_gene_info":
            result = get_gene_info(arguments["gene_symbol"])
        elif name == "run_blast":
            result = run_blast(
                arguments["sequence"],
                arguments.get("program", "blastp"),
                arguments.get("database", "nr")
            )
        elif name == "get_uniprot_info":
            result = get_uniprot_info(arguments["uniprot_id"])
        elif name == "get_server_config_status":
            result = get_server_config_status()
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

async def main_async():
    """启动MCP服务器"""
    async with stdio_server() as (read_stream, write_stream):
        await server.run(
            read_stream,
            write_stream,
            server.create_initialization_options()
        )


def main():
    """同步入口，供 console script 调用。"""
    asyncio.run(main_async())


if __name__ == "__main__":
    main()
