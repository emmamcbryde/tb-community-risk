# APY v9 ABM + Dynamic Model Metric Alignment

This note documents the current output-contract overlap between the APY v9 MATLAB ABM Streamlit bundle and the Python dynamic model bundle. It is intentionally conservative: structural differences should be flagged rather than hidden.

## APY results bundle shape

The APY bundle is assembled by `abm/build_results_bundle_v9.m`.

Top-level sections currently include:

- `metadata`: contract/model/config/scenario labels from `results.interfaceConfig`.
- `validation`: optional validation report rows.
- `headline`: `strategy`, `calibration`, `keyMetricsRows`, and `summaryRows`.
- `technical`: calibration metadata, table metadata, and `interfaceConfig`.
- `targetedVsRandom`, `doNothing`, `attributableRisk`, `charts`, `economics`, `downloads`: optional add-on sections.

Current APY `headline.keyMetricsRows` are selected in `abm/summarise_results_v9.m`:

- `nScreened`
- `nTestPositiveNonActive`
- `nFalsePositiveTreated`
- `nTotalCoursesStarted`
- `nTotalCoursesCompleted`
- `nCuredInfection`
- `nPreventedActiveTB`
- `NNS_cureInfection`
- `NNS_preventActiveTB`
- `NNT_started_cureInfection`
- `NNT_started_preventActiveTB`

Current APY technical fields useful for alignment:

- `technical.interfaceConfig.N`
- `technical.interfaceConfig.followHorizon`
- `metadata.scenarioLabel`
- `headline.strategy.testType`
- `headline.strategy.screeningStrategy`
- `headline.strategy.regimen`

## Dynamic results bundle shape

The dynamic bundle is assembled by `engine/integration/dynamic_output_contract_v9.py`.

Current dynamic headline metrics:

- `projection_horizon`
- `final_year_baseline_incidence_per100k`
- `final_year_intervention_incidence_per100k`
- `cumulative_baseline_active_tb_cases`
- `cumulative_intervention_active_tb_cases`
- `cumulative_cases_averted`
- `relative_reduction_cumulative_active_tb_cases`

Current dynamic projection rows are stored at:

- `projection.annualRows`

Rows use these dataframe columns when available:

- `Year`
- `Baseline_inc_per100k`
- `Intervention_inc_per100k`
- `Baseline_annual_count`
- `Intervention_annual_count`
- `Baseline_cum_count`
- `Intervention_cum_count`
- `Cases_averted_per100k`
- `Cases_averted_annual_count`
- `Cases_averted_cum_count`

## Currently comparable metrics

These are comparable now when both bundles expose numeric values:

- Horizon: dynamic `projection_horizon` vs APY `technical.interfaceConfig.followHorizon`.
- Population: dynamic `metadata.population` vs APY `technical.interfaceConfig.N`.
- Cumulative baseline active TB cases: dynamic `cumulative_baseline_active_tb_cases` vs APY `technical.dynamicComparison.cumulative_baseline_active_tb_cases`, when available.
- Cumulative intervention active TB cases: dynamic `cumulative_intervention_active_tb_cases` vs APY `technical.dynamicComparison.cumulative_intervention_active_tb_cases`, when available.
- Cases averted: dynamic `cumulative_cases_averted` vs APY `technical.dynamicComparison.cumulative_cases_averted`, or fallback APY `nPreventedActiveTB`.
- Relative reduction: dynamic `relative_reduction_cumulative_active_tb_cases` vs APY `technical.dynamicComparison.relative_reduction_cumulative_active_tb_cases`, when available.

The cases-averted fallback through `nPreventedActiveTB` is directionally useful but less transparent than the full `technical.dynamicComparison` block because it does not expose the paired baseline/intervention cumulative active TB counts.

## APY dynamic-comparison block

`abm/build_results_bundle_v9.m` now exposes:

- `cumulative_baseline_active_tb_cases` over `followHorizon`.
- `cumulative_intervention_active_tb_cases` over `followHorizon`.
- `cumulative_cases_averted` over `followHorizon`, explicitly aligned to the baseline/intervention cumulative definitions.
- `relative_reduction_cumulative_active_tb_cases`.
- `population` copied from `results.interfaceConfig.N`.
- `followHorizon` copied from `results.interfaceConfig.followHorizon`.

These are stored at:

```text
bundle.technical.dynamicComparison
```

Preferred source order:

1. `doNothing.derived` when available.
2. `results.raw.nActiveBy20y` and `results.raw.nPreventedActiveTB` fallback.
3. unavailable with explicit `missingFields` and `notes`.

No disease equations or simulation logic are changed.

## Metrics missing from dynamic for APY comparison

The dynamic model does not currently expose APY pathway metrics:

- `nScreened`
- `nTestPositiveNonActive`
- `nFalsePositiveTreated`
- `nTotalCoursesStarted`
- `nTotalCoursesCompleted`
- `nCuredInfection`
- `NNS_cureInfection`
- `NNS_preventActiveTB`
- `NNT_started_cureInfection`
- `NNT_started_preventActiveTB`

These should remain flagged as structurally non-comparable unless the dynamic model later gains equivalent screening/cascade states.

## Recommendation

Use `technical.dynamicComparison` for dynamic/APY comparison where available. Continue to flag APY pathway metrics, economics outputs, and risk-factor attributable outputs as structurally non-comparable unless explicit equivalent dynamic outputs are added later.
