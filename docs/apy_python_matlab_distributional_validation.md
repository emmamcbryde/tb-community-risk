# APY Python vs MATLAB Distributional Validation

This note records the current diagnostic comparison between the Python APY port
and the committed MATLAB APY v9 compact reference fixture:

`validation/matlab_reference/default_random_igra_3hp_N1500_seed1/`

The comparison uses the same scenario controls as the fixture:

- `N = 1500`
- `nReps = 2000`
- `seed = 1`
- `screenWindow = 2`
- `followHorizon = 20`
- `screeningStrategy = random`
- `testType = IGRA`
- `regimen = 3HP`

The helper intentionally uses generous initial tolerances because exact
replicate-level equality is not expected between MATLAB and NumPy random
streams. This is a diagnostic parity check, not yet a strict release gate.

## Current Result

Run date: 2026-05-15

All 22 compared metrics passed the current initial tolerance rule.

| Metric | PythonMedian | MatlabMedian | AbsoluteDifference | RelativeDifference | Pass |
| --- | ---: | ---: | ---: | ---: | --- |
| nScreened | 450.0 | 450.0 | 0.0 | 0.0 | true |
| nInfected | 112.0 | 113.0 | -1.0 | -0.0088495575 | true |
| nLatentAtScreen | 31.0 | 30.0 | 1.0 | 0.0333333333 | true |
| nActiveAtScreen | 3.0 | 3.0 | 0.0 | 0.0 | true |
| nTestPositive | 40.0 | 40.0 | 0.0 | 0.0 | true |
| nTestPositiveNonActive | 37.0 | 37.0 | 0.0 | 0.0 | true |
| nFalsePositiveTests | 8.0 | 8.0 | 0.0 | 0.0 | true |
| nFalsePositiveTreated | 7.0 | 7.0 | 0.0 | 0.0 | true |
| nTotalCoursesStarted | 32.0 | 31.0 | 1.0 | 0.0322580645 | true |
| nTotalCoursesCompleted | 25.0 | 25.0 | 0.0 | 0.0 | true |
| nCuredInfection | 16.0 | 16.0 | 0.0 | 0.0 | true |
| nPreventedActiveTB | 4.0 | 4.0 | 0.0 | 0.0 | true |
| nActiveBy2y | 19.0 | 19.0 | 0.0 | 0.0 | true |
| nActiveBy20y | 37.0 | 37.0 | 0.0 | 0.0 | true |
| NNS_cureInfection | 28.125 | 28.125 | 0.0 | 0.0 | true |
| NNS_preventActiveTB | 112.5 | 112.5 | 0.0 | 0.0 | true |
| NNT_started_cureInfection | 1.9333333333 | 1.9333333333 | 0.0 | 0.0 | true |
| NNT_started_preventActiveTB | 8.1 | 8.3333333333 | -0.2333333333 | -0.028 | true |
| cumulative_baseline_active_tb_cases | 37.0 | 37.0 | 0.0 | 0.0 | true |
| cumulative_intervention_active_tb_cases | 33.0 | 33.0 | 0.0 | 0.0 | true |
| cumulative_cases_averted | 4.0 | 4.0 | 0.0 | 0.0 | true |
| relative_reduction_cumulative_active_tb_cases | 0.1052631579 | 0.1 | 0.0052631579 | 0.0526315789 | true |

## Likely Causes Of Residual Differences

- MATLAB and NumPy use different random-number generators and draw ordering, so
  exact replicate-level equality is not expected.
- Median summaries can still agree closely even when individual replicate rows
  differ.
- Small differences in treatment-start and prevention-derived ratios are
  expected until broader distributional validation is run across multiple
  scenarios and seeds.
- The current fixture is one scenario only; it does not validate every strategy,
  test type, regimen, targeting mode, or edge case.

## Current Limitations

- This is not a full distributional validation suite yet.
- Tolerances are intentionally broad and should be tightened after more MATLAB
  reference fixtures are added.
- Do-nothing and `technical.dynamicComparison` parity are now covered for this
  fixture, but economics, attributable-risk add-ons, and Streamlit Python
  backend switching remain pending.
