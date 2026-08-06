# Milestone 5 Clinician Decision Analysis

Milestone 5 adds the versioned `apy_clinician_decision_analysis_v1`
workflow over the existing APY path:

```text
scenario/configuration
  -> APY expected-value or individual-based model
  -> event ledger
  -> event-ledger economics
  -> decision analysis
  -> UI and workbook
```

It does not create a second epidemiological or economic engine. Transmission
effects, PSA, early-review stopping rules and clinician recommendations remain
outside scope.

## Scenario Comparison

`engine.apy.decision_analysis.run_scenario_comparison` compares explicitly
defined strategies using validated adapters only. Supported scenario changes
include test, regimen, screening strategy, coverage, programme and follow-up
timing, TPT initiation, reviewed test characteristics, reviewed regimen
characteristics and explicit cost-item edits. Arbitrary dotted-path mutation is
rejected.

Expected-value comparisons preserve fractional events. Agent-based comparisons
reuse fixed replicate seeds and record a cohort fingerprint so strategy-only
changes can be checked for the same baseline population, infection state and
untreated comparator history.

## Sensitivity And Threshold Analysis

`engine.apy.sensitivity` defines `apy_one_way_sensitivity_v1` and
`apy_threshold_analysis_v1`. Sensitivity analysis uses deterministic
expected-value runs and only executes parameters with explicit supplied low and
high values. The repository `data/apy_sensitivity_spec.csv` intentionally keeps
APY reference ranges unresolved rather than fabricating ranges.

Threshold analysis evaluates an explicit grid before any root finding. It
reports no crossing, one crossing or multiple crossings and records the complete
grid. Monotonicity is diagnosed, not assumed. NMB thresholds are unavailable
unless the economic willingness-to-pay threshold is valid.

## Early Review

`engine.apy.early_review.run_early_screening_review` implements
`apy_early_screening_review_v1`. It supports explicit beta and discrete-grid
priors only; no hidden prior strength is supplied. The likelihood is an
aggregate binomial approximation:

```text
positive tests ~ Binomial(screened to date, model-predicted positive probability)
```

The predicted positive probability is model-derived for the screened tranche and
includes test sensitivity, specificity and false positives. Zero prevalence is
handled explicitly for false-positive likelihood checks without forcing the APY
calibration through an invalid zero-prevalence target.

Stop-versus-continue compares completed coverage with planned total coverage
using the same deterministic targeting order. Completed screening is retained
and continuing adds only remaining programme activity. Exact within-window
review-time scheduling is not yet represented; outputs set
`timingApproximation=true` and document this approximation.

## Clinician Interface And Export

`pages/5_Decision_Analysis.py` provides three tabs: strategy comparison,
sensitivity/threshold analysis and early screening review. It displays evidence
readiness and provisional status, reports modelled consequences, and avoids
automatic recommendations such as stop/continue or preferred strategy.

The consolidated workbook can include:

- `Decision_scenarios`
- `Scenario_comparison`
- `Scenario_comparison_replicates`
- `Sensitivity_spec`
- `Sensitivity_results`
- `Threshold_results`
- `Threshold_crossings`
- `Early_review_inputs`
- `Early_review_prior_posterior`
- `Early_review_projection`
- `Early_review_summary`
- `Decision_analysis_validation`

Unavailable economics remain blank/null and unresolved evidence remains
unresolved. Synthetic validation scripts are technical fixtures only and do not
establish APY reference evidence.
