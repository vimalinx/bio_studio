---
name: biomni
description: General-purpose biomedical AI agent for research tasks. execute complex biomedical queries and analysis.
allowed-tools: Bash, Python
---

Biomni is an AI agent that can execute biomedical tasks.

## Usage Patterns

### 1. Library Mode (Recommended for Tools)
Use Biomni's powerful tools directly as a Python library. **No API Key required.**

```python
import sys
import subprocess
from pathlib import Path

# 1. Locate Biomni
try:
    root = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip()
    BIOMNI_PATH = Path(root) / "repositories/active/Biomni"
    sys.path.append(str(BIOMNI_PATH))
except:
    pass

# 2. Import Tools
from biomni.tool.molecular_biology import annotate_open_reading_frames, find_restriction_sites

# 3. Use Tools
seq = "ATG...TGA"
orfs = annotate_open_reading_frames(seq, min_length=30)
sites = find_restriction_sites(seq, ["EcoRI"])
```

**Available Modules**:
- `biomni.tool.molecular_biology`: ORFs, Primers, Restriction Sites, Golden Gate Assembly
- `biomni.tool.genomics`: Gene lookup, Sequence alignment
- `biomni.tool.protocols`: Lab protocol generation

### 2. Agent Mode (Full AI)
Initialize the autonomous agent to plan and execute complex tasks. **Requires API Key.**

```python
from biomni.agent import A1

# Initialize agent (set expected_data_lake_files=[] to skip 11GB download)
agent = A1(path=str(BIOMNI_PATH / 'data'), 
           llm='claude-3-5-sonnet-20240620',
           expected_data_lake_files=[])

# Execute task
agent.go("Plan a CRISPR screen to identify genes that regulate T cell exhaustion")
```

## Setup
1. **Dependencies**: `pip install langchain langchain-core beautifulsoup4`
2. **API Key** (Agent Mode only): Set `ANTHROPIC_API_KEY` in `repositories/active/Biomni/.env`

