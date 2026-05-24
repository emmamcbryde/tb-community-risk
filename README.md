# TB Community Risk Models

This repository contains community-level tuberculosis (TB) risk models for planning, prioritisation, and scenario comparison.

The current `model-split` branch separates the user interface, Python modelling code, the MATLAB APY v9 reference model, adapters, and shared app state. The Streamlit app entry point for this branch is:

```bash
streamlit run streamlit_app.py
```

The model is intended to support public health planning and sequencing. It is not a diagnostic tool, a clinical decision rule, or a tool for denying screening or treatment to any individual.

---

## What is included

### 1. APY v9 agent-based model

The APY v9 workflow supports LTBI screening and preventive treatment scenarios. The normal Streamlit APY run path is Python-first and uses the Python APY backend by default. The MATLAB APY v9 implementation remains the reference engine and fixture source, and is available as an optional backend for local reference checks where a MATLAB-capable environment is configured.

The APY workflow includes:

- scenario setup and validation;
- Python APY model execution by default;
- results display;
- partial Python health economics inputs and summaries;
- partial targeting and baseline-vs-comparator strategy comparison;
- partial attributable-risk reporting;
- partial chart/export helpers and downloadable summaries;
- MATLAB reference fixtures and optional MATLAB-backed checks.

Python APY coverage is still being expanded. Treat MATLAB as the reference backend for parity validation and fixture generation, not as a requirement for the normal Streamlit APY workflow.

The core MATLAB simulation engine is:

```text
abm/tb_screening_mc_model_v9.m
```

The main user-facing MATLAB scenario script remains:

```text
abm/run_tb_screening_user_options_v9.m
```

Generated MATLAB/APY outputs should remain under:

```text
abm/output/
```

### 2. Python dynamic LTBI-to-TB model

The repository also contains a Python dynamic model for LTBI-to-active-TB projection. It supports age-structured LTBI seeding, risk-factor multipliers, transmission feedback, intervention simulation, and annualised outputs.

Relevant files include:

```text
engine/dynamic/dynamic_model.py
engine/dynamic/exec_dynamic.py
ui/dynamic_ui.py
```

### 3. Python static model

The older/static Python model components are retained for rapid scenario exploration and for comparison with the dynamic workflow.

Relevant files include:

```text
engine/static/
ui/static_ui.py
```

---

## Repository layout

```text
.
├── streamlit_app.py              # Current Streamlit entry point for model-split
├── requirements.txt              # Python dependencies
├── AGENTS.md                     # Development rules and modelling guardrails
│
├── app/                          # Shared Streamlit state and display helpers
│   ├── state.py
│   └── display.py
│
├── adapters/                     # Python adapters around external/model backends
│   ├── matlab_backend.py        # Optional MATLAB/reference backend adapter
│   ├── paths.py
│   └── serialization.py
│
├── pages/                        # Multipage Streamlit workflow
│   ├── 1_Scenario.py             # APY scenario setup, defaults, validation, save/load
│   ├── 2_Run_Model.py            # Run Python APY by default; optional MATLAB/reference backend
│   ├── 3_Results.py              # Display portable results bundle
│   ├── 4_Economics.py            # Partial Python economics assumptions and outputs
│   ├── 5_Dynamic_Model.py        # Dynamic model workflow page
│   ├── 6_Compare.py              # Partial APY baseline vs comparator comparison
│   └── 7_Dynamic_ABM_Compare.py  # Dynamic model + APY ABM comparison
│
├── abm/                          # MATLAB APY v9 reference engine and fixture source
│   ├── tb_screening_mc_model_v9.m
│   ├── run_scenario_v9.m
│   ├── run_scenario_bundle_json_v9.m
│   ├── run_health_economics_json_v9.m
│   ├── build_default_config_v9.m
│   ├── build_results_bundle_v9.m
│   ├── scenarios/
│   └── output/
│
├── engine/                       # Python APY, dynamic, static, and integration code
│   ├── apy/
│   ├── dynamic/
│   ├── integration/
│   ├── static/
│   ├── infection_backcast.py
│   ├── intervention.py
│   ├── model.py
│   └── params.py
│
├── ui/                           # Dynamic/static UI modules retained from earlier app
│   ├── dynamic_ui.py
│   ├── static_ui.py
│   └── app.py                    # Legacy entry point; not the model-split entry point
│
└── data/                         # Input data used by Python workflows
```

---

## Running locally

### 1. Clone the repository and select the branch

```bash
git clone https://github.com/emmamcbryde/tb-community-risk.git
cd tb-community-risk
git checkout model-split
```

### 2. Create and activate a Python environment

Python 3.11 is the recommended local development and test version for this branch. Local checks have also passed on Python 3.12.4.

On macOS/Linux:

```bash
python -m venv .venv
source .venv/bin/activate
```

On Windows PowerShell:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
```

### 3. Install Python requirements

For runtime use:

```bash
pip install -r requirements.txt
```

For local development and tests:

```bash
pip install -r requirements-dev.txt
```

The main Python requirements are Streamlit, pandas, numpy, Altair, matplotlib, scipy, and requests.

### 4. Run the Streamlit app

```bash
streamlit run streamlit_app.py
```

The app should open in a browser at:

```text
http://localhost:8501
```

---

## APY backend requirements

The normal APY v9 workflow in the Streamlit app uses the Python APY backend by default and should run with the Python requirements installed above.

MATLAB is only required when selecting and running the optional MATLAB/reference APY backend. For MATLAB-backed reference checks, fixture generation, or workflows that explicitly need the MATLAB reference implementation, the local environment needs:

- MATLAB installed;
- the MATLAB Engine API for Python installed in the same Python environment used by Streamlit;
- access to the repository `abm/` folder;
- no hard-coded local machine paths in MATLAB scripts.

The MATLAB adapter starts MATLAB and adds the `abm/` folder to the MATLAB path when MATLAB backend functions are called.

If MATLAB is not available, the default Python APY workflow and the Python dynamic/static components should remain runnable. Only actions that explicitly select and call the MATLAB backend require MATLAB.

The Python APY reference coverage helper is available as `engine.apy.reference_coverage.get_reference_coverage()` and through `PythonApyBackend.reference_coverage()`.

---

## Streamlit Cloud deployment

For the `model-split` branch, set the Streamlit Cloud entry point to:

```text
streamlit_app.py
```

Do not use the old entry point:

```text
ui/app.py
```

That file is retained from the earlier shared dynamic/static app structure, but it is not the main entry point for the current multipage `model-split` branch.

Important deployment note: Streamlit Cloud should not be assumed to have MATLAB installed. The hosted normal run path should use the Python APY backend. MATLAB-backed APY execution is primarily a local or separately configured MATLAB-capable reference workflow.

The Python-only APY migration gate also runs in GitHub Actions via `.github/workflows/python-apy-migration.yml`. The workflow mirrors the documented local checks by running `python scripts/check_python_apy_migration.py` and `python -m pytest tests -q`. It does not cover full MATLAB parity or MATLAB-only reference workflows. The malformed `runtime.txt` deployment contract is a deployment follow-up and is not fixed here.

---

## APY v9 Streamlit workflow

The APY v9 workflow is organised as a multipage Streamlit app.

### Scenario

Use the Scenario page to:

- inspect selected APY backend status;
- load backend defaults;
- edit core scenario assumptions;
- edit advanced run controls;
- validate the scenario;
- save or load scenario JSON files.

Current advanced run controls include:

```text
N
nReps
seed
screenWindow
followHorizon
```

### Run Model

Use the Run Model page to:

- validate the current scenario;
- run the Python APY v9 backend by default, or the optional MATLAB/reference backend when selected;
- store the portable results bundle in Streamlit session state;
- flag old results as stale when assumptions change.

### Results

Use the Results page to:

- view metadata;
- view headline and summary outputs;
- inspect technical details;
- download available partial CSV/JSON outputs;
- check whether economics outputs are current or stale.

### Economics

Use the Economics page to:

- load default economics assumptions;
- load the KWAB150 economics preset;
- edit test, regimen, program, and active TB cost assumptions;
- run the partial Python economics subset against the current model results;
- download economics summaries and assumptions.

The Python economics implementation is partial. Use MATLAB-backed reference outputs where full APY health-economics parity is required.

### Compare

Use the Compare page to compare a baseline APY scenario with a comparator APY scenario. This is a partial Python strategy-comparison workflow, not a complete port of all MATLAB targeting and strategy add-ons.

The Compare page supports:

- loading a default baseline;
- copying the current scenario to the baseline;
- copying the baseline to the comparator;
- applying named stress-test presets;
- editing comparator assumptions;
- running baseline and comparator scenarios;
- running partial economics comparisons;
- displaying assumption, outcome, partial attributable-risk, and partial economics differences where supported;
- downloading comparison CSV/JSON outputs.

If assumptions change after comparison outputs have been generated, old comparison outputs are cleared and the page asks the user to rerun the comparison before showing current outputs.

---

## Dynamic model workflow

The Python dynamic model is retained for continuous-time LTBI-to-TB modelling with transmission feedback. The dynamic workflow supports:

- population size and age distribution inputs;
- historical incidence patterns or uploaded incidence CSVs;
- LTBI back-calculation by age;
- progression risk-factor multipliers;
- calibration of transmission parameters;
- baseline and intervention projection;
- annual incidence and cases-averted outputs.

---

## Integrated workflow

The integrated Dynamic + APY ABM Compare page compares the latest Python dynamic model projection with the latest APY v9 ABM results bundle.

The comparison layer:

- reads existing Streamlit session-state bundles;
- does not run MATLAB unless the user explicitly runs APY from the APY Run Model page;
- shows side-by-side metrics only where a conservative alignment exists;
- flags missing and structurally non-comparable metrics instead of forcing equivalence;
- exposes debug panels for the dynamic bundle and APY bundle.

Current aligned metrics include horizon, population, cumulative cases averted, and, when APY `technical.dynamicComparison` is available, cumulative baseline/intervention active TB cases and relative reduction. This comparison is intentionally partial and conservative.

---

## Input file formats

### Historical incidence CSV

For dynamic model incidence history uploads:

```csv
year,incidence
2010,40
2011,38
2012,36
```

`incidence` should be annual TB incidence per 100,000 population.

### Custom age distribution CSV

For dynamic model custom age distributions:

```csv
AgeGroup,Proportion
0,0.012
1,0.012
2,0.012
```

`Proportion` values should sum to approximately 1.

### APY scenario JSON

APY scenarios are saved and loaded through the Scenario page. Scenario JSON files should be kept under:

```text
abm/scenarios/
```

unless another relative path is explicitly selected.

---

## Outputs

Generated outputs should not be treated as source code.

Recommended output locations:

```text
abm/output/                 # MATLAB/APY generated outputs
abm/scenarios/              # Saved APY scenario JSON files
```

The app may also expose downloadable CSV or JSON outputs from Streamlit session state, including:

- APY summary CSVs;
- key metrics CSVs;
- partial economics summary CSVs;
- assumptions JSON files;
- partial comparison diff CSVs;
- comparison JSON bundles.

---

## Development rules

When editing this branch:

- treat `abm/tb_screening_mc_model_v9.m` as the reference MATLAB APY engine;
- do not create new MATLAB model versions unless explicitly required;
- prefer minimal edits over rewrites;
- keep generated outputs in `abm/output/`;
- avoid hard-coded local machine paths;
- preserve compatibility with `run_tb_screening_user_options_v9.m`;
- keep the model framed as a planning and prioritisation tool, not a tool for refusing care.

---

## Recommended validation after changes

For the local Python-only APY migration gate:

```bash
python scripts/check_python_apy_migration.py
```

This checks MATLAB-free `PythonApyBackend` import behavior and JSON-serialisable Python APY capability/reference-coverage payloads. It does not require or validate MATLAB, OpenAI/Codex, secrets, full MATLAB parity, or MATLAB-only reference workflows. When `pytest` is installed, run `python -m pytest tests -q` separately.

The Python-only APY migration CI check in `.github/workflows/python-apy-migration.yml` mirrors those local checks by running `python scripts/check_python_apy_migration.py` and `python -m pytest tests -q`. It does not cover full MATLAB parity or MATLAB-only reference workflows.

The local migration check is MATLAB-free, Codex-free, OpenAI-key-free, internet-free, and secret-free.

For Python/Streamlit syntax checks:

```powershell
$files = @('streamlit_app.py') + (Get-ChildItem app,adapters,engine,pages,ui -Recurse -Filter *.py | ForEach-Object { $_.FullName })
python -m py_compile @files
```

For broader Python checks, run pytest from the development/test environment:

```bash
python -m pytest tests -q
```

For MATLAB-backed changes, use lightweight MATLAB validation only when MATLAB is available. At minimum, check that the relevant APY v9 functions are on the MATLAB path and that generated outputs are written to the expected output folder.

Example MATLAB path hygiene commands:

```matlab
clear functions
rehash
which tb_screening_mc_model_v9 -all
which run_scenario_v9 -all
```

---

## Troubleshooting

### Optional MATLAB/reference APY actions fail

Likely causes:

- MATLAB is not installed;
- MATLAB Engine API for Python is not installed in the active Python environment;
- MATLAB cannot find the `abm/` folder;
- the app is running in a hosted environment without MATLAB.

Check the backend status shown in the Streamlit sidebar and on the Scenario page. These issues should only block actions that explicitly use the MATLAB/reference backend; the normal Python APY run path should not require MATLAB.

### The app does not open locally

Try:

```bash
streamlit run streamlit_app.py
```

Then open:

```text
http://localhost:8501
```

If that fails, check:

- the active Python environment;
- `pip install -r requirements.txt` completed successfully;
- your browser allows localhost pages;
- endpoint protection or antivirus software is not blocking local Streamlit apps.

### Results or comparison outputs look stale

The app tracks stale state. If assumptions change after a run, rerun the model or comparison before interpreting outputs.

### Streamlit Cloud uses the wrong file

Set the app entry point to:

```text
streamlit_app.py
```

not:

```text
ui/app.py
```

---

## Current APY integration status

The APY migration status is conservative:

- Normal APY Streamlit workflow: Python-first, using `PythonApyBackend` by default.
- MATLAB APY v9: reference backend and fixture source, with `tb_screening_mc_model_v9.m` treated as the reference engine.
- Python economics: partial subset, not full MATLAB health-economics parity.
- Targeting and strategy comparison: partial baseline/comparator support, not a complete port of all MATLAB targeting add-ons.
- Attributable-risk reporting: partial reporting and pass-through/status handling, not complete reactivation attributable-risk parity.
- Chart/export parity: partial JSON/CSV/chart-series helpers, not full MATLAB chart/export parity.
- MATLAB App Designer helper files: documentation/reference inventory only for the Streamlit migration unless explicitly revived; normal APY Streamlit behavior should come from Python/Streamlit code.

The branch also includes a thin dynamic-ABM comparison layer:

```text
engine/integration/risk_factor_crosswalk_v9.csv
engine/integration/compare_dynamic_abm_v9.py
engine/integration/dynamic_output_contract_v9.py
pages/7_Dynamic_ABM_Compare.py
```

The integration target is side-by-side comparison, not forced equivalence. The dynamic Python model and MATLAB APY v9 ABM use different structures and partly different risk-factor definitions, so the comparison page flags which metrics are directly comparable and which are model-specific.
