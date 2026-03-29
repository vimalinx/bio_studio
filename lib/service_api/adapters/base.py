from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class AdapterResult:
    command: list[str]
    returncode: int
    stdout: str
    stderr: str
    started_at: str
    finished_at: str
