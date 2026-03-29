from __future__ import annotations

from fastapi.testclient import TestClient

from scripts.api_server import create_app


def test_healthz_returns_ok(tmp_path) -> None:
    client = TestClient(create_app(tmp_path))

    response = client.get("/healthz")

    assert response.status_code == 200
    assert response.json()["status"] == "ok"


def test_plan_preview_endpoint_returns_preview(tmp_path) -> None:
    client = TestClient(create_app(tmp_path))

    response = client.post(
        "/v1/plans/preview",
        json={"prompt": "analyze ebola benchmarks"},
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["brief"]["task_type"] == "analysis"
    assert payload["selected_capabilities"]


def test_runs_endpoint_creates_and_returns_run(tmp_path) -> None:
    client = TestClient(create_app(tmp_path))

    create_response = client.post(
        "/v1/runs",
        json={"prompt": "analyze ebola benchmarks", "project_name": "api_http_run"},
    )

    assert create_response.status_code == 200
    payload = create_response.json()
    run_id = payload["run"]["run_id"]

    read_response = client.get(f"/v1/runs/{run_id}")
    assert read_response.status_code == 200
    assert read_response.json()["project_name"] == "api_http_run"

    events_response = client.get(f"/v1/runs/{run_id}/events")
    assert events_response.status_code == 200
    assert events_response.json()
