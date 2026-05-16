# MATLAB Reference Outputs

MATLAB APY v9 remains the reference implementation while the Python-native APY ABM port is developed.

Reference outputs in this folder are intended for Python validation. They should be compact, portable, and committed only when they are useful as stable test fixtures.

The Python APY port must match MATLAB reference outputs within documented tolerances.
Exact individual-level equality is not expected because MATLAB and NumPy random streams differ.

Validation should focus on:

- summary statistics;
- medians and uncertainty intervals;
- stable output-contract fields;
- `technical.dynamicComparison` values where available;
- compact deterministic parameter snapshots where available;
- scenario configuration compatibility.

Generated reference outputs should not be treated as general model outputs.
Keep ordinary generated APY files in `abm/output/`, and only commit reference files here when they are intentionally small and useful for tests.

## Scenario Suite

The compact validation suite is defined in:

`validation/matlab_reference/scenario_suite_v1.json`

It is intentionally small and covers random, targeted, TST, alternative
regimen, zero-coverage, and high-coverage scenarios.

## Exporting MATLAB Fixtures

From MATLAB at the repository root or with `abm/` on the path:

```matlab
clear functions
rehash
export_matlab_reference_scenarios_v9
```

The exporter reads the scenario suite and writes each scenario to:

`validation/matlab_reference/<scenario_id>/`

Commit compact files only:

- `scenario_config.json`
- `matlab_dynamic_comparison.json`
- `matlab_summary.csv`
- `manifest.json`

Do not commit large generated files:

- `matlab_results_bundle.json`
- `matlab_raw_replicates.csv`
- large example cohort files

## Running Python Validation

The Python validator does not require MATLAB. It uses already-exported fixtures:

```powershell
python scripts/run_apy_distributional_validation.py `
  --reference-root validation/matlab_reference `
  --suite-file validation/matlab_reference/scenario_suite_v1.json `
  --output-dir validation/output/apy_python_matlab_distributional_validation
```

Use `--quick` to override Python `nReps` for local smoke checks. Missing
fixture directories are reported in `scenario_summary.csv` rather than treated
as Python validation crashes.
