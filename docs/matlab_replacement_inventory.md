# MATLAB Replacement Inventory

This inventory summarizes the current APY v9 MATLAB replacement status. MATLAB
v9, especially `abm/tb_screening_mc_model_v9.m`, remains the reference engine
until parity and validation criteria are explicitly updated.

## Python-Native Components

- Core APY package: `engine/apy/`
- Config, validation, data loading, age distribution, regimen, calibration, and
  cohort primitives: `engine/apy/config.py`, `engine/apy/validation.py`,
  `engine/apy/data.py`, `engine/apy/age_distribution.py`,
  `engine/apy/regimen.py`, `engine/apy/calibration.py`, and
  `engine/apy/cohort.py`
- Simulation and runner: `engine/apy/simulation.py` and
  `engine/apy/runner.py`
- Result bundle, summary, natural history, parity, and distributional
  validation helpers: `engine/apy/results_bundle.py`, `engine/apy/summary.py`,
  `engine/apy/natural_history.py`, `engine/apy/parity.py`, and
  `engine/apy/distributional_validation.py`
- Streamlit/backend integration surface: `adapters/python_apy_backend.py`,
  `adapters/backend.py`, `app/`, `pages/`, and `ui/`
- Dynamic/static Python model areas: `engine/dynamic/`, `engine/static/`, and
  `engine/integration/`
- Python validation entry point:
  `scripts/run_apy_distributional_validation.py`

## Remaining MATLAB-Dependent Areas

- Reference ABM engine and MATLAB scenario wrappers in `abm/`, including
  `tb_screening_mc_model_v9.m`, `run_scenario_v9.m`, and
  `run_scenario_bundle_json_v9.m`
- MATLAB App Designer callback and state helpers in `abm/app_*_v9.m`
- MATLAB health economics functions, including `run_health_economics_v9.m`,
  `run_health_economics_json_v9.m`, and economics config helpers
- MATLAB add-on analyses for targeting, attributable risk, natural-history
  add-ons, charting, and exports
- MATLAB reference fixture generation via
  `abm/export_matlab_reference_scenarios_v9.m`
- Full claim of MATLAB distributional parity across all workflows; current
  validation is fixture-based and diagnostic

## Replacement Status

| Category | Components | Status |
|---|---|---|
| already ported | APY config/defaults, validation, data loading, age distribution, regimen handling, calibration, cohort primitives, single-cohort simulation, replicate runner, summary rows, partial results bundle, natural-history dynamic comparison support, parity helpers, distributional validation runner | Python-native implementation exists under `engine/apy/` with tests and compact MATLAB fixture validation. |
| still MATLAB-dependent | MATLAB reference engine, MATLAB JSON wrappers, health economics backend, App Designer callbacks, targeting and attributable-risk add-ons, MATLAB chart/export helpers, MATLAB fixture exporter | Keep MATLAB as the reference and compatibility source until each area is explicitly ported and validated. |
| fixture-only | `validation/matlab_reference/` scenario suite, compact `scenario_config.json`, `matlab_summary.csv`, `matlab_dynamic_comparison.json`, `manifest.json` files | Retain as validation inputs for Python parity checks; these are not ordinary model outputs. |
| documentation-only | Existing port plan, distributional validation note, manual smoke-test docs, and this inventory | Use for implementation sequencing and validation context; documentation is not an executable replacement. |

## Recommended Next Implementation Order

1. Tighten and expand Python APY fixture validation across the committed
   scenario suite.
2. Implement the Streamlit switch to use `adapters/python_apy_backend.py`
   without removing the MATLAB backend path.
3. Port scenario save/load helpers if Python scenarios need MATLAB-compatible
   round trips.
4. Port health economics defaults, validation, and execution.
5. Port targeting, attributable-risk, and broader natural-history add-ons.
6. Port export/chart helper behavior only where the current Streamlit app still
   needs MATLAB-equivalent file outputs.
7. Retire MATLAB-dependent paths only after documented parity gates are met and
   the user explicitly approves removal.

## Fixture Retention

Retain MATLAB reference fixtures for parity and validation unless they are
explicitly removed later. They provide the stable comparison baseline for the
Python port while preserving the planning and sequencing purpose of the model.
