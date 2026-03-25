from __future__ import annotations

from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
WORKFLOW_PATH = ROOT / ".github" / "workflows" / "ci.yml"
SCRIPT_PATH = ROOT / "scripts" / "ci" / "run_stable_tests.sh"
REQUIREMENTS_PATH = ROOT / "requirements-ci.txt"


def _load_workflow() -> tuple[str, dict]:
    text = WORKFLOW_PATH.read_text(encoding="utf-8")
    data = yaml.safe_load(text)
    return text, data


def test_safe_ci_workflow_exists() -> None:
    assert WORKFLOW_PATH.exists()
    assert SCRIPT_PATH.exists()
    assert REQUIREMENTS_PATH.exists()


def test_safe_ci_workflow_uses_read_only_permissions() -> None:
    _text, data = _load_workflow()

    assert data["permissions"] == {"contents": "read"}
    assert data["concurrency"]["cancel-in-progress"] is True


def test_safe_ci_workflow_only_uses_safe_triggers() -> None:
    _text, data = _load_workflow()

    triggers = data["on"]
    assert "workflow_dispatch" in triggers
    assert "push" in triggers
    assert "pull_request" in triggers
    assert "schedule" not in triggers
    assert "workflow_run" not in triggers
    assert "pull_request_target" not in triggers

    assert triggers["push"]["branches"] == ["master"]
    assert triggers["pull_request"]["branches"] == ["master"]
    assert ".github/workflows/**" in triggers["push"]["paths"]
    assert "tests/**" in triggers["push"]["paths"]
    assert "docs/**" not in triggers["push"]["paths"]


def test_safe_ci_workflow_does_not_include_self_mutation_commands() -> None:
    text, data = _load_workflow()
    job = data["jobs"]["stable-tests"]
    steps = job["steps"]

    assert job["if"] == "github.actor != 'github-actions[bot]'"
    assert any(step.get("uses") == "actions/checkout@v4" for step in steps)
    checkout = next(step for step in steps if step.get("uses") == "actions/checkout@v4")
    assert checkout["with"]["persist-credentials"] is False
    assert checkout["with"]["fetch-depth"] == 1

    forbidden_snippets = [
        "git push",
        "git commit",
        "gh pr",
        "gh repo",
        "peter-evans/create-pull-request",
        "workflow_run",
        "pull_request_target",
        "schedule:",
    ]
    for snippet in forbidden_snippets:
        assert snippet not in text


def test_safe_ci_requirements_stay_minimal() -> None:
    requirements = REQUIREMENTS_PATH.read_text(encoding="utf-8").splitlines()
    lines = [line.strip() for line in requirements if line.strip() and not line.startswith("#")]

    assert lines == [
        "pytest>=8.4.0",
        "PyYAML>=6.0.0",
    ]


def test_safe_ci_script_runs_python_pytest_only() -> None:
    script = SCRIPT_PATH.read_text(encoding="utf-8")

    assert "python -m pytest" in script
    assert "git push" not in script
    assert "git commit" not in script
    assert "gh pr" not in script
    assert "test_github_actions_ci.py" in script
