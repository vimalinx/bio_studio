from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Mapping, Sequence

from lib.workspace_env import build_subprocess_env, resolve_workspace_python
from lib.service_api.adapters.base import AdapterResult
from lib.service_api.models import utc_now_iso


class WorkspaceCLIAdapter:
    def __init__(self, workspace_root: Path, cli_path: Path | None = None) -> None:
        self.workspace_root = Path(workspace_root).resolve()
        self.cli_path = (
            Path(cli_path).resolve()
            if cli_path is not None
            else self.workspace_root / "scripts" / "project.py"
        )

    def invoke(
        self,
        arguments: Sequence[str],
        *,
        extra_env: Mapping[str, str] | None = None,
    ) -> AdapterResult:
        command = [
            str(resolve_workspace_python()),
            str(self.cli_path),
            "--workspace-root",
            str(self.workspace_root),
            *[str(item) for item in arguments],
        ]
        started_at = utc_now_iso()
        completed = subprocess.run(
            command,
            cwd=self.workspace_root,
            env=build_subprocess_env(extra_env=extra_env),
            capture_output=True,
            text=True,
        )
        finished_at = utc_now_iso()
        return AdapterResult(
            command=command,
            returncode=completed.returncode,
            stdout=completed.stdout,
            stderr=completed.stderr,
            started_at=started_at,
            finished_at=finished_at,
        )
