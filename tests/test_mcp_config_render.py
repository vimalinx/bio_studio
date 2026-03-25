from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RENDER_SCRIPT = ROOT / "mcp-servers" / "render_claude_config.py"
INSTALL_SCRIPT = ROOT / "mcp-servers" / "install-all.sh"
CONFIG_PATH = ROOT / "mcp-servers" / "claude-config.json"


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_render_claude_config_builds_repo_relative_server_paths(monkeypatch) -> None:
    module = _load_module(RENDER_SCRIPT, "bio_studio_render_claude_config_test")
    monkeypatch.setattr(module, "resolve_workspace_python", lambda: Path("/tmp/bio/bin/python"))

    config = module.build_config()
    servers = config["mcpServers"]

    assert servers["bio-design"]["command"] == "/tmp/bio/bin/python"
    assert servers["bio-design"]["args"] == [
        str(ROOT / "mcp-servers" / "bio-design-mcp" / "design_server.py")
    ]
    assert servers["bio-lab"]["command"] == "/tmp/bio/bin/python"
    assert servers["bio-lab"]["args"] == [
        str(ROOT / "mcp-servers" / "bio-lab-mcp" / "lab_server.py")
    ]
    assert servers["bio-sequence"]["command"] == "/tmp/bio/bin/python"
    assert servers["bio-sequence"]["args"] == [
        str(ROOT / "mcp-servers" / "bio-sequence-mcp" / "sequence_server.py")
    ]
    assert servers["bio-structure"]["args"] == [
        str(ROOT / "mcp-servers" / "bio-structure-mcp" / "structure_server.py")
    ]
    assert servers["bio-database"]["args"] == [
        str(ROOT / "mcp-servers" / "bio-database-mcp" / "database_server.py")
    ]


def test_install_script_uses_dynamic_repo_paths() -> None:
    text = INSTALL_SCRIPT.read_text(encoding="utf-8")

    assert "/media/vimalinx/Data/bio_studio" not in text
    assert "render_claude_config.py" in text
    assert 'SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"' in text


def test_checked_in_claude_config_matches_current_repo() -> None:
    payload = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    servers = payload["mcpServers"]

    assert Path(servers["bio-design"]["args"][0]).as_posix().endswith(
        "/mcp-servers/bio-design-mcp/design_server.py"
    )
    assert Path(servers["bio-lab"]["args"][0]).as_posix().endswith(
        "/mcp-servers/bio-lab-mcp/lab_server.py"
    )
    assert Path(servers["bio-sequence"]["args"][0]).as_posix().endswith(
        "/mcp-servers/bio-sequence-mcp/sequence_server.py"
    )
    assert Path(servers["bio-structure"]["args"][0]).as_posix().endswith(
        "/mcp-servers/bio-structure-mcp/structure_server.py"
    )
    assert Path(servers["bio-database"]["args"][0]).as_posix().endswith(
        "/mcp-servers/bio-database-mcp/database_server.py"
    )
