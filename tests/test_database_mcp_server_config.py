from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SERVER_PATH = ROOT / "mcp-servers" / "bio-database-mcp" / "database_server.py"


def _install_fake_mcp_modules() -> None:
    mcp_module = types.ModuleType("mcp")
    mcp_server_module = types.ModuleType("mcp.server")
    mcp_stdio_module = types.ModuleType("mcp.server.stdio")
    mcp_types_module = types.ModuleType("mcp.types")

    class FakeServer:
        def __init__(self, name: str) -> None:
            self.name = name

        def list_tools(self):
            def decorator(func):
                return func
            return decorator

        def call_tool(self):
            def decorator(func):
                return func
            return decorator

        def create_initialization_options(self):
            return {}

        async def run(self, *args, **kwargs):
            return None

    async def fake_stdio_server():
        raise RuntimeError("not used in tests")

    class Tool:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

    class TextContent:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

    mcp_server_module.Server = FakeServer
    mcp_stdio_module.stdio_server = fake_stdio_server
    mcp_types_module.Tool = Tool
    mcp_types_module.TextContent = TextContent

    sys.modules["mcp"] = mcp_module
    sys.modules["mcp.server"] = mcp_server_module
    sys.modules["mcp.server.stdio"] = mcp_stdio_module
    sys.modules["mcp.types"] = mcp_types_module


def _install_fake_bio_modules() -> None:
    bio_module = types.ModuleType("Bio")
    bio_blast_module = types.ModuleType("Bio.Blast")

    class FakeEntrez:
        email = None
        api_key = None

    class FakeNCBIWWW:
        @staticmethod
        def qblast(*args, **kwargs):
            raise RuntimeError("not used in tests")

    bio_module.Entrez = FakeEntrez
    bio_module.SeqIO = types.SimpleNamespace(parse=lambda *args, **kwargs: [])
    bio_blast_module.NCBIWWW = FakeNCBIWWW

    sys.modules["Bio"] = bio_module
    sys.modules["Bio.Blast"] = bio_blast_module


def _load_module(name: str):
    _install_fake_mcp_modules()
    _install_fake_bio_modules()
    spec = importlib.util.spec_from_file_location(name, SERVER_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_configure_entrez_prefers_bio_studio_env_vars(monkeypatch) -> None:
    monkeypatch.setenv("BIO_STUDIO_ENTREZ_EMAIL", "lab@example.org")
    monkeypatch.setenv("NCBI_EMAIL", "fallback@example.org")
    monkeypatch.setenv("BIO_STUDIO_NCBI_API_KEY", "secret-key")

    module = _load_module("bio_database_server_config_test")
    config = module.configure_entrez()

    assert config["email"] == "lab@example.org"
    assert config["email_source"] == "BIO_STUDIO_ENTREZ_EMAIL"
    assert config["api_key_configured"] is True
    assert module.Entrez.email == "lab@example.org"
    assert module.Entrez.api_key == "secret-key"


def test_append_entrez_warning_when_email_missing(monkeypatch) -> None:
    monkeypatch.delenv("BIO_STUDIO_ENTREZ_EMAIL", raising=False)
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    monkeypatch.delenv("ENTREZ_EMAIL", raising=False)

    module = _load_module("bio_database_server_warning_test")
    payload = module.append_entrez_warning({"ok": True}, {"email": None})

    assert payload["ok"] is True
    assert "config_warning" in payload
    assert "BIO_STUDIO_ENTREZ_EMAIL" in payload["config_warning"]


def test_main_wrapper_uses_asyncio_run(monkeypatch) -> None:
    module = _load_module("bio_database_server_main_test")
    calls: list[object] = []

    monkeypatch.setattr(module.asyncio, "run", lambda coro: calls.append(coro))

    module.main()

    assert len(calls) == 1
    calls[0].close()


def test_get_server_config_status_reports_sanitized_state(monkeypatch) -> None:
    monkeypatch.setenv("BIO_STUDIO_ENTREZ_EMAIL", "lab@example.org")
    monkeypatch.setenv("BIO_STUDIO_NCBI_API_KEY", "secret-key")

    module = _load_module("bio_database_server_status_test")
    status = module.get_server_config_status()

    assert status["server"] == "bio-database-mcp"
    assert status["entrez_email_configured"] is True
    assert status["ncbi_api_key_configured"] is True
    assert status["entrez_email_source"] == "BIO_STUDIO_ENTREZ_EMAIL"
    assert status["ncbi_api_key_source"] == "BIO_STUDIO_NCBI_API_KEY"
    assert "secret-key" not in str(status)
