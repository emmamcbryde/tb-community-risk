# TB Community Risk Models

This repository contains community-level tuberculosis (TB) risk models for planning, prioritisation, and scenario comparison.

The current `model-split` branch separates the user interface, Python modelling code, MATLAB APY v9 agent-based model, adapters, and shared app state. The Streamlit app entry point for this branch is:

```bash
streamlit run streamlit_app.py
```

The model is intended to support public health planning and sequencing. It is not a diagnostic tool, a clinical decision rule, or a tool for denying screening or treatment to any individual.

---

## What is included

### 1. APY v9 agent-based model

The APY v9 workflow supports LTBI screening and preventive treatment scenarios. The normal Streamlit run path uses the Python APY backend by default; the MATLAB APY v9 implementation remains available as an optional/reference backend for parity checks and economics workflows where applicable.

The APY workflow includes:

- scenario setup and validation;
- Python APY model execution by default;
- results display;
- health economics inputs and summaries;
- baseline vs comparator scenario comparison;
- downloadable summaries and comparison outputs.

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
│   ├── matlab_backend.py
│   ├── paths.py
│   └── serialization.py
│
├── pages/                        # Multipage Streamlit workflow
│   ├── 1_Scenario.py             # APY scenario setup, defaults, validation, save/load
│   ├── 2_Run_Model.py            # Run APY v9 backend
│   ├── 3_Results.py              # Display portable results bundle
│   ├── 4_Economics.py            # Economics assumptions and outputs
│   ├── 5_Dynamic_Model.py        # Dynamic model workflow page
│   ├── 6_Compare.py              # APY baseline vs comparator comparison
│   └── 7_Dynamic_ABM_Compare.py  # Dynamic model + APY ABM comparison
│
├── abm/                          # MATLAB APY v9 ABM and related files
│   ├── tb_screening_mc_model_v9.m
│   ├── run_scenario_v9.m
│   ├── run_scenario_bundle_json_v9.m
│   ├── run_health_economics_json_v9.m
│   ├── build_default_config_v9.m
│   ├── build_results_bundle_v9.m
│   ├── scenarios/
│   └── output/
│
├── engine/                       # Python dynamic/static model code
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

```bash
pip install -r requirements.txt
```

The main Python requirements are Streamlit, pandas, numpy, Altair, matplotlib, and scipy.

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

MATLAB is only required when selecting and running the optional MATLAB/reference APY backend. For MATLAB-backed parity checks or economics workflows where applicable, the local environment needs:

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

Important deployment note: Streamlit Cloud should not be assumed to have MATLAB installed. The hosted normal run path should use the Python APY backend. MATLAB-backed APY execution is primarily a local or separately configured MATLAB-capable backend workflow.

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
- download available CSV outputs;
- check whether economics outputs are current or stale.

### Economics

Use the Economics page to:

- load default economics assumptions;
- load the KWAB150 economics preset;
- edit test, regimen, program, and active TB cost assumptions;
- run economics against the current model results;
- download economics summaries and assumptions.

### Compare

Use the Compare page to compare a baseline APY scenario with a comparator APY scenario.

The Compare page supports:

- loading a default baseline;
- copying the current scenario to the baseline;
- copying the baseline to the comparator;
- applying named stress-test presets;
- editing comparator assumptions;
- running baseline and comparator scenarios;
- running economics comparisons;
- displaying assumption, outcome, and economics differences;
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

Current aligned metrics include horizon, population, cumulative cases averted, and, when APY `technical.dynamicComparison` is available, cumulative baseline/intervention active TB cases and relative reduction.

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
- economics summary CSVs;
- assumptions JSON files;
- comparison diff CSVs;
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

For Python/Streamlit syntax checks:

```powershell
$files = @('streamlit_app.py') + (Get-ChildItem app,adapters,engine,pages,ui -Recurse -Filter *.py | ForEach-Object { $_.FullName })
python -m py_compile @files
```

For broader Python checks, add or run tests when available:

```bash
python -m unittest discover -s tests
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

### The app opens but APY scenario actions fail

Likely causes:

- MATLAB is not installed;
- MATLAB Engine API for Python is not installed in the active Python environment;
- MATLAB cannot find the `abm/` folder;
- the app is running in a hosted environment without MATLAB.

Check the backend status shown in the Streamlit sidebar and on the Scenario page.

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

## Current integration status

The branch includes a thin dynamic-ABM comparison layer:

```text
engine/integration/risk_factor_crosswalk_v9.csv
engine/integration/compare_dynamic_abm_v9.py
engine/integration/dynamic_output_contract_v9.py
pages/7_Dynamic_ABM_Compare.py
```

The integration target is side-by-side comparison, not forced equivalence. The dynamic Python model and MATLAB APY v9 ABM use different structures and partly different risk-factor definitions, so the comparison page flags which metrics are directly comparable and which are model-specific.
