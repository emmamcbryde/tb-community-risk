# TB Codex Quick Start Cheat Sheet

**Search tag:** `#TB-CODEX-QUICKSTART`

Use this tag in ChatGPT, VS Code, or file search to find this note again.

---

## 1. Open the right things

Open:

```text
1. VS Code
2. ChatGPT in browser
3. PowerShell terminal inside VS Code
4. Optional: File Explorer to your repo folder
```

Your repo folder is:

```powershell
C:\Users\emmas\Documents\GITHUB\tb-community-risk
```

In VS Code:

```text
File → Open Folder → C:\Users\emmas\Documents\GITHUB\tb-community-risk
```

Then open a terminal:

```text
Terminal → New Terminal
```

Make sure it is PowerShell.

---

## 2. Get into the repo

In PowerShell:

```powershell
cd C:\Users\emmas\Documents\GITHUB\tb-community-risk
```

Then check where you are:

```powershell
git branch --show-current
git status --short
```

Most recently, the branch you have been working on is:

```text
python-apy-abm-port
```

If you need to switch to it:

```powershell
git checkout python-apy-abm-port
git pull origin python-apy-abm-port
```

---

## 3. Activate your Python environment

Usually:

```powershell
conda activate tbmodel
```

Then quick check:

```powershell
python --version
```

Optional but useful:

```powershell
python -m unittest discover -s tests
```

---

## 4. Resume Codex

From the repo root:

```powershell
codex resume --last
```

If that does not find the right session:

```powershell
codex resume
```

If you were in a different folder when you started the old session:

```powershell
codex resume --all
```

If you have a specific session ID:

```powershell
codex resume <session-id>
```

Example:

```powershell
codex resume 019d94fb-b8a6-7ad1-93f8-c546cb18eb88
```

---

## 5. If Codex opens but seems confused

First tell it where it is:

```text
We are in the tb-community-risk repository on branch python-apy-abm-port.

Please first run:
git branch --show-current
git status --short
git log --oneline -5

Do not make changes yet. Report status only.
```

This prevents Codex from charging ahead in the wrong branch or folder.

---

## 6. If Codex is stuck or has been running too long

Inside Codex, ask:

```text
Please stop any long-running commands and report status only:
1. Are you still running a command?
2. What command?
3. git status --short
4. files changed
5. tests run
6. anything committed?
Do not make new changes.
```

If you think Streamlit is stuck, in a separate PowerShell window:

```powershell
Get-Process python,streamlit -ErrorAction SilentlyContinue
Get-NetTCPConnection -LocalPort 8501,8599 -ErrorAction SilentlyContinue
```

If needed:

```powershell
Stop-Process -Id <PID>
```

---

## 7. Usual validation commands

Before trusting any Codex change, run or ask Codex to run:

```powershell
python -m unittest discover -s tests
```

Then:

```powershell
$files = @('streamlit_app.py') + (Get-ChildItem app,adapters,engine,pages,ui,scripts -Recurse -Filter *.py | ForEach-Object { $_.FullName })
python -m py_compile @files
```

And:

```powershell
python -m py_compile engine/apy/*.py
```

---

## 8. Streamlit manual test

From repo root:

```powershell
streamlit cache clear
streamlit run streamlit_app.py
```

If port `8501` is busy:

```powershell
streamlit run streamlit_app.py --server.port 8599
```

Then check in browser:

```text
Scenario page
Run Model page
Results page
Dynamic + ABM Compare page
Economics page
```

For current work, the key check is:

```text
Python APY backend works without MATLAB.
Economics/report outputs are offline and not Streamlit-dependent unless explicitly asked.
```

---

## 9. Where files usually live

### Core repo

```text
C:\Users\emmas\Documents\GITHUB\tb-community-risk
```

### Python APY model

```text
engine/apy/
```

### Scripts

```text
scripts/
```

### Outputs

```text
outputs/
```

### Word reports

Likely:

```text
paper/
```

or repo root if you placed them there.

### Excel economics workbooks

Recommended location:

```text
paper/excel/
```

### Generated figures

Recommended location:

```text
paper/figures/
```

### Validation fixtures

```text
validation/matlab_reference/
```

---

## 10. Files not to commit unless you mean to

Be careful with:

```text
outputs/
*.docx
*.xlsx
*.png
*.pdf
openaikey*.txt
reboot_text.txt
index.html
```

Check before committing:

```powershell
git status --short
```

If you see something like:

```text
?? openaikeyAPI20260522.txt
```

do **not** commit it.

For local-only ignore rules:

```powershell
Add-Content .git\info\exclude "`n# Local-only files"
Add-Content .git\info\exclude "/openaikey*.txt"
Add-Content .git\info\exclude "/outputs/"
Add-Content .git\info\exclude "/paper/excel/*.xlsx"
Add-Content .git\info\exclude "/*.docx"
```

Use `.git/info/exclude` for personal local exclusions because it does not change the repo.

---

## 11. Standard “start Codex safely” prompt

Paste this whenever you resume:

```text
We are working in the tb-community-risk repository.

Before doing anything, please run:
git branch --show-current
git status --short
git log --oneline -5

Report the current branch and whether the working tree is clean.

Do not make changes yet.
Do not commit anything yet.
```

Then, after it reports status, give the actual task.

---

## 12. Current project mental map

Your current main working branch has been:

```text
python-apy-abm-port
```

Current broad state:

```text
Python APY backend: working experimentally
MATLAB: still reference implementation
Streamlit online tool: working, but not the current focus
Economics paper module: active focus
Word report: preferred writing/editing format
Excel workbook: audit/decision-tree layer
```

The division of labour is:

```text
APY Python model:
  epidemiological outputs, stochastic simulation, cases prevented

Excel:
  editable costing assumptions, DALYs/QALYs, ICERs, NMB, break-even, charts

Word:
  report writing and interpretation

Codex:
  code/document/workbook updates inside repo

ChatGPT:
  planning, prompts, interpretation, checking logic
```

---

## 13. If you want to continue the economics workbook task

Your next prompt to Codex should probably start like this:

```text
We are working on branch python-apy-abm-port.

This is offline paper/economics work, not online-tool work.

Do not modify:
- Streamlit pages
- MATLAB model code
- dynamic model equations
- APY epidemiological model equations

Current focus:
Update the Excel economics workbook and Word report so QALYs include mortality consistently with DALYs, while preserving Dale-compatible QALYs and adding GBD-aligned QALY sensitivity.

Before making changes, run:
git branch --show-current
git status --short
git log --oneline -5
```

Then paste the detailed task.

---

## 14. Quick emergency reset

If everything feels confusing:

```powershell
cd C:\Users\emmas\Documents\GITHUB\tb-community-risk
git branch --show-current
git status --short
```

If there are unexpected untracked files, do not delete immediately. Ask:

```powershell
git status --short
```

Then paste the output into ChatGPT or Codex and ask for a status-only diagnosis.

---

## Reminder about uploaded files

Some older files uploaded into ChatGPT may expire. If you want ChatGPT to inspect the latest Word report, Excel workbook, or uploaded papers again, re-upload the current versions.
