from __future__ import annotations

import io
import importlib.util
import json
import urllib.error
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "scripts" / "skills" / "render_skill_from_help.py"


def _load_module(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_render_skill_from_help_outputs_valid_skill_markdown() -> None:
    module = _load_module(SCRIPT_PATH, "render_skill_from_help")
    text = module.render_skill_markdown(
        {
            "name": "samtools",
            "summary": "Utilities for working with SAM/BAM/CRAM files.",
            "help_text": "samtools view\nsamtools sort\nsamtools index",
            "path": "/usr/bin/samtools",
        },
        mode="offline",
    )

    assert text.startswith("---")
    assert "name: samtools" in text
    assert "description:" in text
    assert "disable-model-invocation: true" in text
    assert "samtools" in text


def test_load_render_config_keeps_integer_api_tuning_fields(tmp_path: Path) -> None:
    module = _load_module(SCRIPT_PATH, "render_skill_from_help")
    config_path = tmp_path / "config.local.json"
    config_path.write_text(
        json.dumps(
            {
                "base_url": "http://example.test/v1",
                "model": "gpt-5.4",
                "api_key": "test-key",
                "max_completion_tokens": 4096,
                "request_retries": 5,
                "request_timeout_seconds": 90,
            }
        ),
        encoding="utf-8",
    )

    config = module.load_render_config(config_path=config_path, env={})

    assert config["max_completion_tokens"] == 4096
    assert config["request_retries"] == 5
    assert config["request_timeout_seconds"] == 90


def test_extract_response_text_supports_sse_chunks() -> None:
    module = _load_module(SCRIPT_PATH, "render_skill_from_help")
    body = "\n".join(
        [
            'data: {"choices":[{"delta":{"content":"---\\n"}}]}',
            'data: {"choices":[{"delta":{"content":"name: samtools\\n"}}]}',
            'data: {"choices":[{"delta":{"content":"description: Use when ..."}}]}',
            "data: [DONE]",
        ]
    )

    text = module._extract_response_text(body)

    assert text.startswith("---")
    assert "name: samtools" in text


def test_render_skill_from_help_online_retries_until_valid_response(tmp_path: Path, monkeypatch) -> None:
    module = _load_module(SCRIPT_PATH, "render_skill_from_help")
    config_path = tmp_path / "config.local.json"
    config_path.write_text(
        json.dumps(
            {
                "base_url": "http://example.test/v1",
                "model": "gpt-5.4",
                "api_key": "test-key",
                "reasoning_effort": "xhigh",
                "max_completion_tokens": 2048,
                "request_retries": 1,
            }
        ),
        encoding="utf-8",
    )

    class _FakeResponse:
        def __init__(self, body: str):
            self._body = body.encode("utf-8")

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def read(self) -> bytes:
            return self._body

    seen_efforts: list[str] = []

    def fake_urlopen(request, timeout=120):
        payload = json.loads(request.data.decode("utf-8"))
        effort = payload.get("reasoning_effort", "")
        seen_efforts.append(effort)
        if effort in {"xhigh", "high"}:
            raise urllib.error.HTTPError(
                request.full_url,
                500,
                "Internal Server Error",
                hdrs=None,
                fp=io.BytesIO(b'{"error":{"message":"upstream failed"}}'),
            )
        if effort in {"medium", "low"}:
            return _FakeResponse('data: {"choices":[]}\n\ndata: [DONE]\n')
        return _FakeResponse(
            json.dumps(
                {
                    "choices": [
                        {
                            "message": {
                                "content": (
                                    "---\n"
                                    "name: samtools\n"
                                    "description: Use when you need samtools for BAM work.\n"
                                    "disable-model-invocation: true\n"
                                    "user-invocable: true\n"
                                    "---\n\n"
                                    "# samtools\n\n"
                                    "## Quick Start\n\n"
                                    "- Command: `samtools`\n"
                                    "- Local executable: `/usr/bin/samtools`\n"
                                    "- Raw local help snapshot: [references/help.md](references/help.md)\n\n"
                                    "## What This Tool Is Good For\n\n"
                                    "- BAM and CRAM inspection\n"
                                    "- Sorting and indexing alignments\n"
                                    "- Region-based extraction\n\n"
                                    "## Recommended Workflow\n\n"
                                    "1. Read `references/help.md`.\n"
                                    "2. Check inputs.\n"
                                    "3. Run the smallest command.\n"
                                    "4. Validate outputs.\n\n"
                                    "## Guardrails\n\n"
                                    "- Use explicit input paths.\n"
                                    "- Record version details.\n"
                                    "- Validate output files.\n"
                                )
                            }
                        }
                    ]
                }
            )
        )

    monkeypatch.setattr(module.urllib.request, "urlopen", fake_urlopen)

    text = module.render_skill_markdown(
        {
            "name": "samtools",
            "summary": "Utilities for SAM/BAM/CRAM files.",
            "help_text": "samtools view\nsamtools sort\nsamtools index",
            "path": "/usr/bin/samtools",
        },
        mode="online",
        config_path=config_path,
        env={},
    )

    assert text.startswith("---")
    assert "## Guardrails" in text
    assert seen_efforts == ["xhigh", "high", "medium", "low", ""]
