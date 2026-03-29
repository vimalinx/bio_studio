from __future__ import annotations

from dataclasses import dataclass

from .models import CapabilityRecord


@dataclass(slots=True)
class CapabilityRegistry:
    _capabilities: dict[str, CapabilityRecord]

    def list_capabilities(self) -> list[CapabilityRecord]:
        return list(self._capabilities.values())

    def get(self, capability_id: str) -> CapabilityRecord | None:
        return self._capabilities.get(capability_id)

    def find_by_task_type(self, task_type: str) -> list[CapabilityRecord]:
        return [
            capability
            for capability in self._capabilities.values()
            if task_type in capability.task_types
        ]

    def find_by_keyword(self, keyword: str) -> list[CapabilityRecord]:
        lowered = keyword.lower().strip()
        if not lowered:
            return []
        return [
            capability
            for capability in self._capabilities.values()
            if lowered in capability.capability_id.lower()
            or lowered in capability.display_name.lower()
            or any(lowered in entry.lower() for entry in capability.keywords)
        ]


def build_default_registry() -> CapabilityRegistry:
    capabilities = [
        CapabilityRecord(
            capability_id="workspace.project.create",
            display_name="Create workspace project shell",
            task_types=["project_setup", "analysis", "design"],
            keywords=["create", "project", "setup", "workspace"],
            transport="cli",
            backend_target="scripts/project.py create",
            supports_preview=False,
        ),
        CapabilityRecord(
            capability_id="workspace.project.validate",
            display_name="Validate workspace project",
            task_types=["validation", "analysis"],
            keywords=["validate", "project", "check"],
            transport="cli",
            backend_target="scripts/project.py validate",
        ),
        CapabilityRecord(
            capability_id="workspace.project.run",
            display_name="Run workspace project pipeline",
            task_types=["analysis", "design"],
            keywords=["run", "pipeline", "execute", "project"],
            transport="cli",
            backend_target="scripts/project.py run",
            supports_async=True,
        ),
        CapabilityRecord(
            capability_id="workspace.validation.smoke",
            display_name="Workspace validation smoke test",
            task_types=["workspace_validation"],
            keywords=["workspace", "smoke", "validate", "environment"],
            transport="cli",
            backend_target="scripts/project.py workspace-validate",
            supports_async=True,
        ),
        CapabilityRecord(
            capability_id="database.search.pubmed",
            display_name="Search literature and biomedical databases",
            task_types=["analysis", "research", "literature"],
            keywords=["pubmed", "paper", "literature", "database", "search"],
            transport="mcp",
            backend_target="bio-database-mcp",
            supports_async=True,
        ),
        CapabilityRecord(
            capability_id="sequence.analysis.basic",
            display_name="Basic sequence analysis",
            task_types=["analysis", "sequence"],
            keywords=["sequence", "fasta", "blast", "alignment", "genome", "rna", "dna"],
            transport="mcp",
            backend_target="bio-sequence-mcp",
            supports_async=True,
        ),
        CapabilityRecord(
            capability_id="structure.analysis.fold",
            display_name="Protein structure analysis",
            task_types=["analysis", "structure", "design"],
            keywords=["structure", "fold", "protein", "alphafold"],
            transport="mcp",
            backend_target="bio-structure-mcp",
            supports_async=True,
            resource_class="heavy",
        ),
        CapabilityRecord(
            capability_id="design.mrna.optimize",
            display_name="mRNA and construct optimization",
            task_types=["design"],
            keywords=["design", "mrna", "construct", "optimize", "codon"],
            transport="mcp",
            backend_target="bio-design-mcp",
            supports_async=True,
        ),
    ]
    return CapabilityRegistry({item.capability_id: item for item in capabilities})
