# API Configuration Profiles

This directory contains multiple API configuration profiles for easy switching.

## Profile Files

Each profile is a `.env` file containing API keys and configuration for a specific setup.

## Available Profiles

- `default.env` - Default configuration (copy from .env.example)
- `anthropic.env` - Anthropic/Claude models only
- `openai.env` - OpenAI models only
- `azure.env` - Azure OpenAI configuration
- `custom.env` - Custom model serving (Ollama, local models, etc.)
- `biomni-r0.env` - Biomni-R0 with SGLang setup

## Usage

Use the `switch_profile.py` script to switch between profiles:

```bash
python switch_profile.py
```

Or manually:

```bash
# Activate a specific profile
cp profiles/anthropic.env .env
```

## Creating New Profiles

To create a new profile:

1. Copy an existing profile: `cp profiles/default.env profiles/myprofile.env`
2. Edit the new file with your API keys
3. Switch to it: `python switch_profile.py` and select your profile
