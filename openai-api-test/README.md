# OpenAI API smoke test

This folder contains a minimal local smoke test confirming that the active Python environment can call the OpenAI Responses API.

Setup used:

```powershell
$env:OPENAI_API_KEY="...put key here..."
python -m pip install --upgrade openai
python test_api.py