from __future__ import annotations

from lib.service_api.registry import build_default_registry


def test_default_registry_exposes_workspace_and_bio_capabilities() -> None:
    registry = build_default_registry()
    capability_ids = {item.capability_id for item in registry.list_capabilities()}

    assert "workspace.project.create" in capability_ids
    assert "workspace.project.validate" in capability_ids
    assert "workspace.project.run" in capability_ids
    assert "database.search.pubmed" in capability_ids
    assert "sequence.analysis.basic" in capability_ids
