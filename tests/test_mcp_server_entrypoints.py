from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_mcp_server_files_use_sync_main_wrapper() -> None:
    server_files = [
        ROOT / "mcp-servers" / "bio-design-mcp" / "design_server.py",
        ROOT / "mcp-servers" / "bio-lab-mcp" / "lab_server.py",
        ROOT / "mcp-servers" / "bio-sequence-mcp" / "sequence_server.py",
        ROOT / "mcp-servers" / "bio-structure-mcp" / "structure_server.py",
        ROOT / "mcp-servers" / "bio-database-mcp" / "database_server.py",
    ]

    for path in server_files:
        text = path.read_text(encoding="utf-8")
        assert "async def main_async()" in text, path.name
        assert "def main():" in text, path.name
        assert "if __name__ == \"__main__\":" in text, path.name
