# Codex supervisor automation

This folder contains local automation for connecting an OpenAI API supervisor loop to Codex CLI.

## Files

- `supervisor_loop.py`: asks an OpenAI model to supervise a coding task, delegates concrete implementation steps to `codex exec`, then feeds Codex output, Git diff, and test output back into the supervisor.

## Requirements

Run from an activated Python environment with:

```powershell
python -m pip install --upgrade openai
```

The terminal must have:

```powershell
$env:OPENAI_API_KEY="..."
```

Do not commit the real API key.

Codex CLI must also be installed and logged in:

```powershell
codex --version
codex login status
```

## First safe run

From the repository root:

```powershell
python automation\supervisor_loop.py . --sandbox read-only --max-turns 1
```

A good first goal is:

```text
Inspect this repository and summarize the structure. Do not edit files.
```

## Local Python APY migration gate

From the repository root, run:

```powershell
python scripts\check_python_apy_migration.py
```

This is a local Python-only gate for APY migration hygiene. It checks that `PythonApyBackend` can be imported without loading MATLAB modules, and that the Python APY capability and reference-coverage payloads are JSON-serialisable.

It does not require or validate MATLAB, OpenAI/Codex, API keys, secrets, full MATLAB parity, or MATLAB reference workflows. Use the separate MATLAB/reference validation flow when parity or fixture/reference behavior is the question.

When `pytest` is installed, `python -m pytest tests -q` can be run separately. GitHub Actions/CI for this gate is intentionally deferred until the repo has a clear Python version and dependency contract.

## Normal cautious run

From the repository root:

```powershell
python automation\supervisor_loop.py . --sandbox workspace-write --max-turns 4
```

The script stops when the supervisor says the task is done, when it needs human guidance, or when the maximum number of Codex turns is reached.

## Safety notes

- Start from a clean Git working tree.
- Work on a branch, not directly on main.
- Use `read-only` for inspection.
- Use `workspace-write` for ordinary local coding edits.
- Avoid broader sandbox modes unless working in an isolated environment.
- Review `git diff` before committing.
