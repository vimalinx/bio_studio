from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MCP_README = ROOT / "mcp-servers" / "README.md"
DATABASE_README = ROOT / "mcp-servers" / "bio-database-mcp" / "README.md"
LAB_README = ROOT / "mcp-servers" / "bio-lab-mcp" / "README.md"


def test_mcp_readmes_do_not_reference_legacy_mount_paths() -> None:
    for path in (MCP_README, DATABASE_README, LAB_README):
        text = path.read_text(encoding="utf-8")
        assert "/media/vimalinx/Data/bio_studio" not in text, path.name
        assert "/run/media/vimalinx" not in text, path.name


def test_database_readme_uses_env_var_based_entrez_setup() -> None:
    text = DATABASE_README.read_text(encoding="utf-8")

    assert "BIO_STUDIO_ENTREZ_EMAIL" in text
    assert "get_server_config_status" in text
    assert "your-email@example.com" not in text


def test_mcp_root_readme_mentions_renderer() -> None:
    text = MCP_README.read_text(encoding="utf-8")

    assert "render_claude_config.py" in text
    assert "claude-config.json" in text


def test_lab_readme_mentions_workspace_validation_tools() -> None:
    text = LAB_README.read_text(encoding="utf-8")

    assert "run_workspace_validation" in text
    assert "scripts/project.py" in text
