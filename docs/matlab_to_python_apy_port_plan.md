# MATLAB to Python APY v9 Port Plan

## Scope

This inventory covers the MATLAB files in `abm/` on the `model-split` branch.
The Python port should preserve `abm/tb_screening_mc_model_v9.m` as the reference implementation until validation passes.
This plan does not propose changes to the MATLAB engine, the Python dynamic model equations, or the current Streamlit behavior.

The first Python-native APY target is a backend MVP that can run APY scenarios without MATLAB and produce the same JSON-like results contract used by the existing Streamlit APY Results and Dynamic + APY ABM Compare pages.

## Recommended Python Package Layout

```text
engine/apy/
  __init__.py
  config.py                  # default APY config and config normalization
  validation.py              # structured validation report
  simulation.py              # Python port of tb_screening_mc_model_v9 runtime
  scenario.py                # save/load helpers
  results.py                 # summary rows and results bundle contract
  do_nothing.py              # natural-history/do-nothing comparator
  economics.py               # later health economics backend
  targeting.py               # later targeted/random strategy add-ons
  attributable_risk.py       # later attributable-risk add-on
  data/
    default_data.csv
    default_age_distribution.csv
tests/apy/
  fixtures/matlab_reference_bundles/
  test_config_parity.py
  test_validation_parity.py
  test_simulation_smoke.py
  test_results_bundle_contract.py
  test_matlab_reference_parity.py
```

Keep the existing MATLAB-backed adapter while the Python backend is introduced. The future switch should be an adapter choice, not a Streamlit page rewrite.

## Inventory Table

| MATLAB file | Category | Python target module | Priority | Validation target | Notes |
|---|---|---|---|---|---|
| `app_apply_economics_ui_state_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Existing Streamlit economics page remains source of UI behavior | App Designer state writer only. |
| `app_apply_ui_state_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Existing Streamlit Scenario page remains source of UI behavior | App Designer component writer only. |
| `app_collect_ui_state_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Existing Streamlit state helpers cover this role | App Designer component reader only. |
| `app_economics_input_changed_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Streamlit economics dirty/stale state | UI-state orchestration, not model runtime. |
| `app_export_v9.m` | MATLAB App Designer/UI callback | `engine.apy.results` later, if needed | P3 | Compare output files against MATLAB exports if export parity is needed | App wrapper around export. Python downloads currently handle most Streamlit needs. |
| `app_has_fresh_results_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Streamlit session-state stale flags | UI state helper. |
| `app_input_changed_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Streamlit dirty/stale state | UI orchestration. |
| `app_load_economics_preset_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Streamlit economics preset behavior | UI wrapper around economics preset. |
| `app_load_scenario_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Streamlit scenario load flow | UI wrapper around scenario load. |
| `app_refresh_bundle_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Python bundle builder tests | UI state wrapper around bundle construction. |
| `app_run_attributable_risk_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Add-on parity later through `engine.apy.attributable_risk` | UI wrapper around add-on. |
| `app_run_charts_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Chart parity only if MATLAB chart exports remain required | UI wrapper around chart add-on. |
| `app_run_do_nothing_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Do-nothing parity covered by Python `do_nothing.py` | UI wrapper around do-nothing add-on. |
| `app_run_economics_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Economics parity covered later by `engine.apy.economics` | UI wrapper around economics backend. |
| `app_run_targeted_vs_random_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Targeting parity later through `engine.apy.targeting` | UI wrapper around targeting add-ons. |
| `app_run_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Python adapter-level smoke tests | UI wrapper around validate/run/bundle. |
| `app_save_scenario_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Streamlit scenario save flow | UI wrapper around scenario save. |
| `app_startup_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Streamlit app startup state | App Designer startup only. |
| `app_validate_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Streamlit validation flow | UI wrapper around validation. |
| `build_default_config_v9.m` | config/defaults | `engine.apy.config` | P0 | Python default config equals MATLAB JSON-normalized default except portable paths/table representation | Canonical APY scenario defaults. |
| `build_default_economics_config_v9.m` | health economics | `engine.apy.economics` | P1 | Default economics config parity | Needed for Python economics page after core MVP. |
| `build_economics_preset_kwab150_v9.m` | health economics | `engine.apy.economics` | P1 | KWAB150 preset parity | Preserve visible economics preset values. |
| `build_results_bundle_v9.m` | result bundle/output contract | `engine.apy.results` | P0 | JSON bundle schema and summary rows match MATLAB reference bundles | Includes `technical.dynamicComparison`; essential for existing Results and integrated compare pages. |
| `build_ui_schema_v9.m` | MATLAB App Designer/UI callback | none, or `app/schema.py` later | P3 | Field mapping smoke tests only if schema-driven Streamlit forms are added | Current Streamlit forms are hand-built. |
| `collect_validation_issues_json_v9.m` | validation | none | do-not-port | MATLAB adapter bridge remains until Python backend replaces MATLAB | JSON wrapper for MATLAB Engine boundary, not Python runtime logic. |
| `collect_validation_issues_v9.m` | validation | `engine.apy.validation` | P0 | Validation report parity for valid defaults and representative invalid configs | Structured warnings/errors are required before Python runs. |
| `config_to_ui_state_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Existing Streamlit widget sync tests | UI conversion helper. |
| `configs_substantively_equal_v9.m` | validation | `engine.apy.validation` or `app/state.py` if needed | P2 | Dirty/stale equivalence tests if Python backend needs MATLAB-equivalent comparison | Useful but not required for backend MVP. |
| `economics_config_to_ui_state_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Streamlit economics widget sync | UI conversion helper. |
| `export_results_v9.m` | charts/export helper | `engine.apy.results` or `app/downloads.py` | P2 | CSV export parity for summary/key metrics | Python Streamlit already has simple downloads; full export parity can wait. |
| `get_output_dir_v9.m` | charts/export helper | `engine.apy.paths` | P2 | Path resolves to `abm/output/` without hard-coded local paths | Needed only for file-output parity. |
| `initialise_app_state_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Streamlit `app/state.py` coverage | App Designer state initializer. |
| `just_for_executing_codex_checks.m` | legacy/test/generated/ignore | none | do-not-port | None | Temporary check file; do not port. |
| `load_scenario_v9.m` | scenario save/load | `engine.apy.scenario` | P1 | Round-trip saved scenario JSON with MATLAB examples | Useful immediately after core MVP. |
| `regimen_cost_field_v9.m` | health economics | `engine.apy.economics` | P1 | Regimen cost field mapping parity | Small helper for economics assumptions. |
| `results_bundle_to_display_v9.m` | result bundle/output contract | `app/display.py` or `engine.apy.results` | P2 | Existing Streamlit displays continue to consume bundle directly | Mostly MATLAB display shaping; not needed if bundle contract is Python-native. |
| `run_health_economics_json_v9.m` | health economics | none | do-not-port | MATLAB adapter bridge remains until Python economics is ported | JSON wrapper for MATLAB Engine boundary. |
| `run_health_economics_v9.m` | health economics | `engine.apy.economics` | P1 | Economics summary rows and totals match MATLAB reference bundles | Important after epidemiologic MVP. |
| `run_scenario_bundle_json_v9.m` | core simulation engine | none | do-not-port | MATLAB adapter bridge remains during transition | JSON wrapper for MATLAB Engine boundary. |
| `run_scenario_v9.m` | core simulation engine | `engine.apy.simulation` | P0 | Python run resolves config defaults and raw outputs like MATLAB wrapper | Thin but important config-to-engine adapter. |
| `run_tb_screening_compare_strategies_v9.m` | targeting/strategy add-on | `engine.apy.targeting` | P2 | Strategy-comparison output parity | Not required for first single-scenario MVP. |
| `run_tb_screening_do_nothing_v9.m` | natural-history/do-nothing add-on | `engine.apy.do_nothing` | P0 | `technical.dynamicComparison` source values match MATLAB reference | Needed to expose baseline/intervention cumulative active TB metrics conservatively. |
| `run_tb_screening_igra_addons_v9.m` | charts/export helper | `engine.apy.results` later | P3 | Add-on output parity only if retained | Chart/report add-on; not needed for MVP. |
| `run_tb_screening_igra_charts_v9.m` | charts/export helper | `engine.apy.results` or visualization layer later | P3 | Chart data parity only if retained | MATLAB chart export helper, not backend MVP. |
| `run_tb_screening_natural_history_addons_v9.m` | natural-history/do-nothing add-on | `engine.apy.do_nothing` later | P2 | Natural-history add-on table parity | Broader add-on; do-nothing subset is P0. |
| `run_tb_screening_reactivation_attributable_v9.m` | attributable-risk add-on | `engine.apy.attributable_risk` | P2 | Attributable-risk table parity | Useful analytical output, but not required for first Python run. |
| `run_tb_screening_targeted_gradient_v9.m` | targeting/strategy add-on | `engine.apy.targeting` | P2 | Targeting gradient parity | Later targeted strategy exploration. |
| `run_tb_screening_targeting_optima_v9.m` | targeting/strategy add-on | `engine.apy.targeting` | P2 | Targeting optima parity | Later targeted strategy exploration. |
| `run_tb_screening_targeting_profile_v9.m` | targeting/strategy add-on | `engine.apy.targeting` | P2 | Targeting profile parity | Later targeted strategy exploration. |
| `run_tb_screening_user_options_v9.m` | config/defaults | `engine.apy.config` examples/docs | P1 | Compatibility examples reproduce user-options defaults | Preserve compatibility, but do not make it the Python engine entry point. |
| `save_scenario_v9.m` | scenario save/load | `engine.apy.scenario` | P1 | Scenario JSON round-trip parity | Needed for full Streamlit scenario workflow. |
| `schema_component_name_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | None unless schema-driven UI is revived | App Designer naming helper. |
| `schema_field_key_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | None unless schema-driven UI is revived | App Designer naming helper. |
| `summarise_results_v9.m` | result bundle/output contract | `engine.apy.results` | P0 | Summary/key metric values match MATLAB for stored reference runs | May duplicate local helper in engine; still part of output parity. |
| `tb_screening_mc_model_v9.m` | core simulation engine | `engine.apy.simulation` | P0 | Statistical parity against MATLAB reference outputs across fixed seeds/scenarios | Reference engine; do not edit in this migration pass. |
| `ui_state_to_config_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Streamlit form-to-config behavior | UI conversion helper. |
| `ui_state_to_economics_config_v9.m` | MATLAB App Designer/UI callback | none | do-not-port | Streamlit economics widget sync | UI conversion helper. |
| `v9main.m` | legacy/test/generated/ignore | none | do-not-port | None | Convenience script; useful as manual reference only. |
| `validate_config_v9.m` | validation | `engine.apy.validation` | P0 | Validated/normalized config parity | Applies metadata/text defaults and raises model-blocking validation errors. |
| `validate_economics_config_v9.m` | health economics | `engine.apy.economics` | P1 | Economics validation report parity | Needed when Python economics backend is enabled. |
| `validation_report_to_display_v9.m` | validation | `app/display.py` | P2 | Existing validation display remains readable | Display shaping, not backend validation logic. |

## Non-MATLAB Assets Inspected

The P0 port should also treat these files as runtime fixtures or reference artifacts:

| File or group | Role | Port handling |
|---|---|---|
| `abm/default_data.csv` | Primary APY parameter data | Copy or load from `engine/apy/data/default_data.csv`; avoid hard-coded local paths. |
| `abm/default_age_distribution.csv` | Default APY age distribution | Copy or load from `engine/apy/data/default_age_distribution.csv`. |
| `abm/default_age_distribution.xlsx` | Alternate age distribution source | P2 only; CSV should be preferred for Python MVP. |
| `abm/*scenario*.json`, `abm/scenarios/*.json` | Scenario fixtures | Use as P1 scenario load/save and config compatibility fixtures. |
| `abm/output/*.csv` | Generated/reference outputs | Do not treat as canonical unless explicitly captured as reference fixtures. Keep generated APY outputs in `abm/output/`. |
| `abm/APYV9WebApp.mlapp`, `abm/app_code.txt` | Frozen MATLAB UI artifacts | Do not port; useful only as historical UI reference. |

## Minimum P0 Set for Python ABM MVP

The minimum Python-native APY MVP should port:

- `build_default_config_v9.m`
- `validate_config_v9.m`
- `collect_validation_issues_v9.m`
- `run_scenario_v9.m`
- `tb_screening_mc_model_v9.m`
- `summarise_results_v9.m`
- `build_results_bundle_v9.m`
- `run_tb_screening_do_nothing_v9.m`
- Runtime data assets: `default_data.csv` and `default_age_distribution.csv`

This P0 set is enough to:

- load the canonical default APY config;
- validate and normalize a config;
- run one APY scenario in Python;
- produce a JSON-like results bundle compatible with the existing Streamlit Results page;
- populate `technical.dynamicComparison`;
- compare Python outputs against stored MATLAB reference bundles.

## Suggested Port Order

1. Port `engine.apy.config` and `engine.apy.validation` first.
2. Add MATLAB reference fixtures by running small scenarios through the existing MATLAB JSON wrappers and committing only compact JSON bundles under tests or docs fixtures.
3. Port deterministic helpers from `tb_screening_mc_model_v9.m`: regimen library, safe math, age distribution loading, discrete draws, test performance, and summary helpers.
4. Port `simulate_one_cohort` and the replicate runner in `tb_screening_mc_model_v9.m`.
5. Port `run_scenario_v9.m` as the Python adapter from config to simulation args.
6. Port `summarise_results_v9.m` and `build_results_bundle_v9.m`.
7. Port `run_tb_screening_do_nothing_v9.m` sufficiently to produce `technical.dynamicComparison`.
8. Add a Python backend adapter behind the existing Streamlit backend interface, but keep MATLAB as the default until parity is acceptable.

Current implementation status:

- Config/default/validation parity is implemented in `engine.apy.config` and `engine.apy.validation`.
- Deterministic data loading, age distribution, regimen, summary, calibration, and cohort primitives are implemented.
- Single-cohort simulation is implemented in `engine.apy.simulation`.
- The full replicate runner, result bundle builder, do-nothing simulation, economics, and Streamlit Python-backend switch are still pending.
- The Python ABM does not yet claim MATLAB distributional parity for full scenario summaries.

## Validation Strategy

Use MATLAB v9 as the reference until explicitly retired.

P0 validation should be MATLAB-free by default and use stored reference fixtures:

- Default config parity: compare Python default config to JSON-normalized MATLAB default config, allowing portable path differences.
- Validation parity: compare validation issues for valid defaults and representative invalid configs.
- Contract parity: assert Python bundles expose the same top-level keys and key rows consumed by Streamlit.
- Dynamic-comparison parity: compare `technical.dynamicComparison` fields for stored MATLAB reference runs.
- Statistical output parity: for fixed seeds and small `N`/`nReps`, compare medians and intervals within pre-agreed tolerances rather than exact row-by-row identity unless random streams are intentionally matched.
- Regression fixtures: store compact MATLAB reference bundles for at least default APY, IGRA vs TST, 3HP vs 4R, low vs high coverage, and cascade-improvement scenarios.

Manual MATLAB validation remains useful during development but should not be required for ordinary Python CI.

## Priority Summary

### P0 Files

- `build_default_config_v9.m`
- `build_results_bundle_v9.m`
- `collect_validation_issues_v9.m`
- `run_scenario_v9.m`
- `run_tb_screening_do_nothing_v9.m`
- `summarise_results_v9.m`
- `tb_screening_mc_model_v9.m`
- `validate_config_v9.m`

### P1 Files

- `build_default_economics_config_v9.m`
- `build_economics_preset_kwab150_v9.m`
- `load_scenario_v9.m`
- `regimen_cost_field_v9.m`
- `run_health_economics_v9.m`
- `run_tb_screening_user_options_v9.m`
- `save_scenario_v9.m`
- `validate_economics_config_v9.m`

### P2 Files

- `configs_substantively_equal_v9.m`
- `export_results_v9.m`
- `get_output_dir_v9.m`
- `results_bundle_to_display_v9.m`
- `run_tb_screening_compare_strategies_v9.m`
- `run_tb_screening_natural_history_addons_v9.m`
- `run_tb_screening_reactivation_attributable_v9.m`
- `run_tb_screening_targeted_gradient_v9.m`
- `run_tb_screening_targeting_optima_v9.m`
- `run_tb_screening_targeting_profile_v9.m`
- `validation_report_to_display_v9.m`

### P3 Files

- `app_export_v9.m`
- `build_ui_schema_v9.m`
- `run_tb_screening_igra_addons_v9.m`
- `run_tb_screening_igra_charts_v9.m`

### Do-Not-Port Files

- `app_apply_economics_ui_state_v9.m`
- `app_apply_ui_state_v9.m`
- `app_collect_ui_state_v9.m`
- `app_economics_input_changed_v9.m`
- `app_has_fresh_results_v9.m`
- `app_input_changed_v9.m`
- `app_load_economics_preset_v9.m`
- `app_load_scenario_v9.m`
- `app_refresh_bundle_v9.m`
- `app_run_attributable_risk_v9.m`
- `app_run_charts_v9.m`
- `app_run_do_nothing_v9.m`
- `app_run_economics_v9.m`
- `app_run_targeted_vs_random_v9.m`
- `app_run_v9.m`
- `app_save_scenario_v9.m`
- `app_startup_v9.m`
- `app_validate_v9.m`
- `collect_validation_issues_json_v9.m`
- `config_to_ui_state_v9.m`
- `economics_config_to_ui_state_v9.m`
- `initialise_app_state_v9.m`
- `just_for_executing_codex_checks.m`
- `run_health_economics_json_v9.m`
- `run_scenario_bundle_json_v9.m`
- `schema_component_name_v9.m`
- `schema_field_key_v9.m`
- `ui_state_to_config_v9.m`
- `ui_state_to_economics_config_v9.m`
- `v9main.m`

## Remaining Planning Notes

- The Python port should keep APY outputs JSON-like at the backend boundary. Raw numpy/pandas objects should not be stored in Streamlit session state.
- Randomness parity should be defined explicitly. Exact MATLAB/Python RNG matching may not be worth the complexity; distributional and fixed-fixture parity is the pragmatic target.
- Path handling should remain portable. Python should resolve bundled data relative to the package/repo root and continue writing generated APY outputs under `abm/output/` only when file export is requested.
- The model must remain framed as a planning and prioritisation tool, not a tool for denying care.
