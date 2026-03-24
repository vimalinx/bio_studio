from __future__ import annotations

import asyncio
import importlib.util
import json
import sys
import types
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SERVER_PATH = ROOT / "mcp-servers" / "bio-design-mcp" / "design_server.py"


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


def _load_module(name: str):
    _install_fake_mcp_modules()
    spec = importlib.util.spec_from_file_location(name, SERVER_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_list_tools_exposes_design_workflows() -> None:
    module = _load_module("bio_design_server_tools_test")

    tools = asyncio.run(module.list_tools())
    names = [tool.kwargs["name"] for tool in tools]

    assert "design_dna" in names
    assert "design_mrna" in names
    assert "analyze_target" in names
    assert "get_server_capabilities" in names


def test_call_tool_design_dna_uses_loaded_modules(monkeypatch) -> None:
    module = _load_module("bio_design_server_design_dna_test")

    class FakeDNADesigner:
        def __init__(self, target_species: str = "E.coli") -> None:
            self.target_species = target_species

        def protein_to_dna(self, protein_sequence: str, optimize: bool = True) -> str:
            assert protein_sequence == "MKT"
            assert optimize is True
            return "ATGAAAACC"

        def analyze_sequence(self, dna_sequence: str) -> dict[str, object]:
            assert dna_sequence == "ATGAAAACC"
            return {"length": 9, "gc_content": 33.3}

        def design_primers(self, dna_sequence: str) -> dict[str, object]:
            assert dna_sequence == "ATGAAAACC"
            return {"forward_primer": {"sequence": "ATGAAAACC"}}

    fake_modules = {"DNADesigner": FakeDNADesigner}
    monkeypatch.setattr(module, "load_design_modules", lambda: fake_modules)

    result = asyncio.run(
        module.call_tool(
            "design_dna",
            {
                "protein_sequence": "MKT",
                "species": "E.coli",
                "optimize": True,
                "include_primers": True,
            },
        )
    )

    payload = json.loads(result[0].kwargs["text"])
    assert payload["tool"] == "design_dna"
    assert payload["dna_sequence"] == "ATGAAAACC"
    assert payload["analysis"]["length"] == 9
    assert payload["primers"]["forward_primer"]["sequence"] == "ATGAAAACC"


def test_call_tool_design_mrna_uses_expression_level(monkeypatch) -> None:
    module = _load_module("bio_design_server_design_mrna_test")

    class FakeMRNAOptimizer:
        def __init__(self, species: str = "human") -> None:
            self.species = species

        def optimize_for_expression(
            self, protein_sequence: str, expression_level: str = "high"
        ) -> dict[str, object]:
            assert protein_sequence == "MKT"
            assert expression_level == "medium"
            return {
                "sequence": "GCCACCAUGAAAACC",
                "analysis": {"stability_score": 0.9},
                "length": 15,
            }

    fake_modules = {"MRNAOptimizer": FakeMRNAOptimizer}
    monkeypatch.setattr(module, "load_design_modules", lambda: fake_modules)

    result = asyncio.run(
        module.call_tool(
            "design_mrna",
            {
                "protein_sequence": "MKT",
                "species": "human",
                "expression_level": "medium",
            },
        )
    )

    payload = json.loads(result[0].kwargs["text"])
    assert payload["tool"] == "design_mrna"
    assert payload["mrna_sequence"] == "GCCACCAUGAAAACC"
    assert payload["analysis"]["stability_score"] == 0.9


def test_call_tool_analyze_target_uses_loaded_modules(monkeypatch) -> None:
    module = _load_module("bio_design_server_target_test")

    class FakeTargetAnalyzer:
        def __init__(self, work_dir: str) -> None:
            self.work_dir = work_dir

        def analyze_target(
            self,
            target_id: str,
            id_type: str = "gene_symbol",
            diseases: list[str] | None = None,
        ) -> dict[str, object]:
            assert target_id == "TP53"
            assert id_type == "gene_symbol"
            assert diseases == ["cancer"]
            return {
                "target_id": "TP53",
                "overall_score": {"overall_score": 82.5, "grade": "A - 优秀靶点"},
            }

    fake_modules = {"TargetAnalyzer": FakeTargetAnalyzer}
    monkeypatch.setattr(module, "load_design_modules", lambda: fake_modules)

    result = asyncio.run(
        module.call_tool(
            "analyze_target",
            {"target_id": "TP53", "id_type": "gene_symbol", "diseases": ["cancer"]},
        )
    )

    payload = json.loads(result[0].kwargs["text"])
    assert payload["tool"] == "analyze_target"
    assert payload["report"]["target_id"] == "TP53"
    assert payload["report"]["overall_score"]["overall_score"] == 82.5


def test_main_wrapper_uses_asyncio_run(monkeypatch) -> None:
    module = _load_module("bio_design_server_main_test")
    calls: list[object] = []

    monkeypatch.setattr(module.asyncio, "run", lambda coro: calls.append(coro))

    module.main()

    assert len(calls) == 1
    calls[0].close()
