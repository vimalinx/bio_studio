#!/usr/bin/env python3
"""
根据当前仓库位置和 bio 环境生成 Claude MCP 配置。
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from lib.workspace_env import resolve_workspace_python


MCP_DIR = Path(__file__).resolve().parent


def build_config() -> dict[str, object]:
    python_path = str(resolve_workspace_python())
    return {
        "mcpServers": {
            "bio-design": {
                "command": python_path,
                "args": [str(MCP_DIR / "bio-design-mcp" / "design_server.py")],
                "description": "DNA/mRNA 设计与靶点评估",
            },
            "bio-lab": {
                "command": python_path,
                "args": [str(MCP_DIR / "bio-lab-mcp" / "lab_server.py")],
                "description": "工作区项目发现、验证与实验链路 smoke test",
            },
            "bio-sequence": {
                "command": python_path,
                "args": [str(MCP_DIR / "bio-sequence-mcp" / "sequence_server.py")],
                "description": "DNA/RNA/蛋白质序列分析",
            },
            "bio-structure": {
                "command": python_path,
                "args": [str(MCP_DIR / "bio-structure-mcp" / "structure_server.py")],
                "description": "蛋白质结构分析和预测",
            },
            "bio-database": {
                "command": python_path,
                "args": [str(MCP_DIR / "bio-database-mcp" / "database_server.py")],
                "description": "生物数据库查询（NCBI、UniProt等）",
            },
        }
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Render Bio Studio Claude MCP config")
    parser.add_argument("--write", type=Path, help="Write JSON to file")
    parser.add_argument(
        "--python-path",
        action="store_true",
        help="Only print the resolved bio python path",
    )
    args = parser.parse_args()

    if args.python_path:
        print(resolve_workspace_python())
        return 0

    payload = build_config()
    rendered = json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    if args.write:
        args.write.write_text(rendered, encoding="utf-8")
    else:
        sys.stdout.write(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
