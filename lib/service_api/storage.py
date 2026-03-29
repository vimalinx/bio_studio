from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .models import PlanPreview, RunArtifact, RunRecord, utc_now_iso


class FileRunStore:
    def __init__(self, workspace_root: Path) -> None:
        self.workspace_root = Path(workspace_root).resolve()
        self.control_root = self.workspace_root / ".codex-runs" / "service_api"
        self.index_dir = self.control_root / "index"
        self.index_dir.mkdir(parents=True, exist_ok=True)

    def write_run(self, record: RunRecord) -> RunRecord:
        project_root = self.workspace_root / "projects" / record.project_name
        if project_root.exists():
            run_dir = self._resolve_project_run_dir(record.project_name)
        else:
            run_dir = self.control_root / "runs" / record.run_id
        run_dir.mkdir(parents=True, exist_ok=True)

        updated_record = record.model_copy(
            update={
                "artifact_dir": str(run_dir),
                "updated_at": utc_now_iso(),
            }
        )
        manifest_path = run_dir / "run_manifest.json"
        manifest_path.write_text(
            updated_record.model_dump_json(indent=2),
            encoding="utf-8",
        )
        self._index_path(updated_record.run_id).write_text(
            json.dumps(
                {
                    "run_id": updated_record.run_id,
                    "project_name": updated_record.project_name,
                    "run_dir": str(run_dir),
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        return updated_record

    def read_run(self, run_id: str) -> RunRecord:
        manifest_path = self.run_dir(run_id) / "run_manifest.json"
        if not manifest_path.exists():
            raise FileNotFoundError(f"Run manifest not found for {run_id}")
        return RunRecord.model_validate_json(manifest_path.read_text(encoding="utf-8"))

    def append_event(self, run_id: str, payload: dict[str, Any]) -> Path:
        run_dir = self.run_dir(run_id)
        run_dir.mkdir(parents=True, exist_ok=True)
        event_path = run_dir / "events.jsonl"
        event_payload = {"timestamp": utc_now_iso(), **payload}
        with event_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(event_payload, ensure_ascii=False) + "\n")
        return event_path

    def write_plan_preview(self, run_id: str, preview: PlanPreview) -> Path:
        run_dir = self.run_dir(run_id)
        run_dir.mkdir(parents=True, exist_ok=True)
        preview_path = run_dir / "plan_preview.json"
        preview_path.write_text(preview.model_dump_json(indent=2), encoding="utf-8")
        return preview_path

    def read_events(self, run_id: str) -> list[dict[str, Any]]:
        event_path = self.run_dir(run_id) / "events.jsonl"
        if not event_path.exists():
            return []
        return [
            json.loads(line)
            for line in event_path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]

    def list_artifacts(self, run_id: str) -> list[RunArtifact]:
        run_dir = self.run_dir(run_id)
        if not run_dir.exists():
            raise FileNotFoundError(f"Run artifacts not found for {run_id}")

        artifacts: list[RunArtifact] = []
        for child in sorted(run_dir.iterdir(), key=lambda item: item.name):
            if not child.is_file():
                continue
            artifacts.append(
                RunArtifact(
                    kind=self._artifact_kind(child.name),
                    path=str(child),
                    label=child.name,
                )
            )
        return artifacts

    def run_dir(self, run_id: str) -> Path:
        index_path = self._index_path(run_id)
        if not index_path.exists():
            raise FileNotFoundError(f"Run index not found for {run_id}")
        payload = json.loads(index_path.read_text(encoding="utf-8"))
        return Path(payload["run_dir"])

    def _resolve_run_dir(self, project_name: str) -> Path:
        return self._resolve_project_run_dir(project_name)

    def _resolve_project_run_dir(self, project_name: str) -> Path:
        return self.workspace_root / "projects" / project_name / "api"

    def _index_path(self, run_id: str) -> Path:
        return self.index_dir / f"{run_id}.json"

    @staticmethod
    def _artifact_kind(filename: str) -> str:
        if filename == "run_manifest.json":
            return "run_manifest"
        if filename == "plan_preview.json":
            return "plan_preview"
        if filename == "events.jsonl":
            return "events"
        return "file"
