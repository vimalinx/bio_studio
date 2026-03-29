#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

import uvicorn
from fastapi import FastAPI, HTTPException


REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from lib.service_api.execution import default_workspace_root, get_run, submit_run
from lib.service_api.models import RunCreateRequest
from lib.service_api.registry import build_default_registry
from lib.service_api.storage import FileRunStore


def create_app(workspace_root: Path | None = None) -> FastAPI:
    workspace_root = Path(workspace_root or default_workspace_root()).resolve()
    registry = build_default_registry()
    store = FileRunStore(workspace_root)

    app = FastAPI(
        title="Bio Studio Online Control Plane",
        version="0.1.0",
    )

    @app.get("/healthz")
    def healthz() -> dict[str, str]:
        return {
            "status": "ok",
            "workspace_root": str(workspace_root),
        }

    @app.get("/v1/capabilities")
    def list_capabilities():
        return registry.list_capabilities()

    @app.post("/v1/plans/preview")
    def preview_plan(request: RunCreateRequest):
        preview_request = request.model_copy(update={"execution_mode": "plan"})
        from lib.service_api.planner import plan_request

        return plan_request(preview_request, registry)

    @app.post("/v1/runs")
    def create_run(request: RunCreateRequest):
        record, preview = submit_run(request, workspace_root=workspace_root, registry=registry)
        return {"run": record, "plan_preview": preview}

    @app.get("/v1/runs/{run_id}")
    def read_run(run_id: str):
        try:
            return get_run(run_id, workspace_root=workspace_root)
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc

    @app.get("/v1/runs/{run_id}/events")
    def read_run_events(run_id: str):
        try:
            return store.read_events(run_id)
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc

    @app.get("/v1/runs/{run_id}/artifacts")
    def read_run_artifacts(run_id: str):
        try:
            return store.list_artifacts(run_id)
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc

    return app


app = create_app()


def main() -> None:
    uvicorn.run("scripts.api_server:app", host="127.0.0.1", port=8000, reload=False)


if __name__ == "__main__":
    main()
