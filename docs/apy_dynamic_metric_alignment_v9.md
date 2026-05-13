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
- Cases averted: dynamic `cumulative_cases_averted` vs APY `nPreventedActiveTB`.

The cases-averted comparison is directionally useful but not structurally identical unless the APY bundle explicitly exposes baseline and intervention cumulative active TB cases over the same horizon.

## Metrics missing from APY for better dynamic comparison

Recommended APY-side additions to `abm/build_results_bundle_v9.m` or a non-engine bundle helper:

- `cumulative_baseline_active_tb_cases` over `followHorizon`.
- `cumulative_intervention_active_tb_cases` over `followHorizon`.
- `cumulative_cases_averted` over `followHorizon`, explicitly aligned to the baseline/intervention cumulative definitions.
- `relative_reduction_cumulative_active_tb_cases`.
- `population` copied from `results.interfaceConfig.N`.
- `followHorizon` copied from `results.interfaceConfig.followHorizon`.
- `scenarioLabel`.
- `screeningStrategy`.
- `testType`.
- `regimen`.

Do not add these by changing `abm/tb_screening_mc_model_v9.m`. If needed, add derived display fields in `abm/build_results_bundle_v9.m` or another app-facing wrapper.

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

The smallest useful next APY output-contract addition is a derived comparison section that exposes population, horizon, cumulative intervention active TB cases, cumulative baseline/no-intervention active TB cases, cumulative cases averted, and relative reduction. This should be added at the bundle/output layer only, without changing APY simulation logic.
