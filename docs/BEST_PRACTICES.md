# Bio Studio Best Practices & Lessons Learned

> **Last Updated**: 2026-01-24
> **Context**: Summary of experiences from the environment upgrade and Ebola analysis workflow implementation.

This document records critical technical insights, common pitfalls, and architectural decisions to guide future development in Bio Studio.

---

## 1. Pipeline Robustness (The `$PATH` Trap)

### Problem
When running Python scripts that invoke shell commands (via `subprocess`) in non-interactive environments (e.g., CI/CD, Docker exec, or automated agents), the shell often fails to inherit the user's Conda environment variables. This results in `command not found` errors even if tools are installed.

### Solution
**Always explicitly inject the Conda bin directory into `os.environ["PATH"]` at the script's start.**

```python
import os

# Crucial for non-interactive shells
CONDA_BIN = "/home/vimalinx/miniforge3/envs/bio/bin"
if CONDA_BIN not in os.environ["PATH"]:
    os.environ["PATH"] = f"{CONDA_BIN}:{os.environ.get('PATH', '')}"
```

**Rule**: Do not assume the environment is activated. Enforce it in the code.

---

## 2. AI Model Integration (The "Black Box" Output)

### Problem
Foundation models (like Evo 2, StripedHyena) often have complex, undocumented return types from their `forward()` methods. They may return `CausalLMOutput`, `tuple`, or nested `tuple` structures depending on the library version or configuration.
*Case Study*: Evo 2 returned a `tuple` containing another `tuple` containing the `logits` tensor (`((logits, ...), ...)`), causing `AttributeError` and `IndexError` during standard access.

### Best Practice
Write **defensive, adaptive extraction logic** for model outputs. Do not hardcode indices (e.g., `output[0]`).

```python
def find_tensor(obj):
    """Recursively search for the first Tensor in a nested structure."""
    if isinstance(obj, torch.Tensor):
        return obj
    if isinstance(obj, (tuple, list)):
        for item in obj:
            res = find_tensor(item)
            if res is not None: return res
    return None

# Usage
output = model(input_ids)
logits = find_tensor(output)
```

**Rule**: Always inspect `type(output)` during the first integration attempt.

---

## 3. Tool Location Strategy (Conda vs Mamba)

### Problem
Not all tools reside in the active Conda environment's `bin` directory.
- **Environment Tools** (e.g., `samtools`, `python`): Located in `/envs/<name>/bin/`.
- **Infrastructure Tools** (e.g., `mamba`, `conda`): Often located in the base installation's `bin/` or `condabin/`, not the environment.

### Solution
1.  **Do not hardcode paths** like `/envs/bio/bin/mamba`.
2.  **Dynamic Discovery**: Use `which mamba` or look in standard base locations (`~/miniforge3/bin/mamba`).
3.  **Cross-Environment Calls**: To install into an env from a script, call the base `mamba` with `-n <env_name>`.

```bash
# Correct
/path/to/base/bin/mamba install -n bio ...
```

---

## 4. Leveraging Biomni (Agent vs Library)

### Insight
Biomni is architected as a collection of standalone Python tools wrapped by an AI agent. You can bypass the Agent layer to use the tools directly.

### Usage Modes
1.  **Library Mode (Recommended)**: Import `biomni.tool.*` directly. fast, deterministic, no API key required.
    - *Requirement*: Install `langchain` and `beautifulsoup4` even if not using the Agent (imports are coupled).
2.  **Agent Mode**: Use `biomni.agent.A1`. Good for open-ended planning but requires API keys and downloading data lakes.

---

## 5. Project Structure

### Workspace Organization
- **Root**: Keep clean. Only entry points (`start.sh`, `CLAUDE.md`) and config dirs.
- **Scripts**:
    - `projects/<name>/scripts/`: Project-specific logic.
    - `scripts/maintenance/`: System-wide admin scripts.
- **Archiving**: When upgrading a project, move old reports/results to `archive/` immediately to prevent confusion between "current" and "obsolete" findings.

---

## 5. Development Workflow

1.  **Verification First**: Before running a heavy analysis, create a "Toy Project" (e.g., `test_env_validation`) with synthetic data (kb-sized) to test the entire toolchain in seconds.
2.  **Path Handling**: Always use `pathlib.Path` objects instead of strings for file manipulation to avoid OS-specific separator issues and enable LSP support.
    - *Bad*: `RAW_DIR = "data/raw"`
    - *Good*: `RAW_DIR = Path("data") / "raw"`
