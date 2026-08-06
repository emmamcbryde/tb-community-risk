# Overnight Goal M5 Status

## Objective

Build and harden a readiness-aware clinician decision-analysis workflow v1 for
the APY LTBI screening model, using the existing APY configuration, event-ledger,
health-economics and evidence-readiness contracts.

## Starting Commit

`a5cc400`

## Current Checkpoint

Checkpoint 5 - workbook, downloads and reproducibility complete.

## Changes Completed

- Created working branch `codex/apy-clinician-decision-analysis-m5`.
- Added `engine/apy/decision_analysis.py` with the
  `apy_clinician_decision_analysis_v1` scenario-comparison contract.
- Added validated scenario field adapters for test, regimen, screening strategy,
  coverage, timing, TPT initiation, test characteristics and regimen
  characteristics.
- Added deterministic expected-value and stochastic APY scenario-comparison
  orchestration over the existing event ledger and health-economics engine.
- Added paired stochastic cohort fingerprints to detect unintended changes in
  baseline population, infection state and untreated comparator history when
  comparing strategy-only changes.
- Added `data/apy_sensitivity_spec.csv` with unresolved APY sensitivity ranges
  rather than fabricated low/high values.
- Added `engine/apy/sensitivity.py` with one-way deterministic sensitivity,
  explicit economic cost-item adapters, threshold-analysis grids, monotonicity
  diagnostics, crossing detection and bracketed root finding only when a valid
  crossing exists.
- Added `engine/apy/early_review.py` with explicit beta/discrete-grid priors,
  aggregate-binomial likelihood, false-positive-aware posterior updating,
  posterior projection summaries and stop-versus-continue coverage comparison.
- Added zero-prevalence likelihood handling for false positives without forcing
  the APY calibration through an invalid zero-prevalence target.
- Added `scripts/validate_apy_early_review.py` using labelled synthetic
  validation inputs.
- Added `pages/5_Decision_Analysis.py` with tabs for strategy comparison,
  one-way sensitivity/threshold analysis and early screening review. The page
  calls the authoritative decision-analysis bundle functions and displays
  unresolved-readiness warnings rather than conclusions.
- Extended `build_results_workbook` with optional decision-analysis output
  sheets for scenarios, scenario summaries, replicate economics, sensitivity,
  threshold grids/crossings, early-review inputs/posteriors/projections and
  validation.
- Wired the Results page workbook download to include decision-analysis outputs
  stored in Streamlit session state.

## Focused Tests Run

- `conda run --no-capture-output -n tbmodel python -m unittest tests.test_apy_decision_analysis -v`
  - 4 tests passed.
- `conda run --no-capture-output -n tbmodel python -m unittest tests.test_apy_decision_analysis -v`
  - Initial Checkpoint 2 run timed out during an unnecessarily expensive
    threshold unit test.
- `conda run --no-capture-output -n tbmodel python -m unittest tests.test_apy_decision_analysis -v`
  - 9 tests passed.
- `conda run --no-capture-output -n tbmodel python -m unittest tests.test_apy_decision_analysis.ApyDecisionAnalysisEarlyReviewTests -v`
  - Initial run timed out before duplicate identical-coverage projections were
    cached.
- `conda run --no-capture-output -n tbmodel python -m unittest tests.test_apy_decision_analysis.ApyDecisionAnalysisEarlyReviewTests -v`
  - 8 tests passed.
- `conda run --no-capture-output -n tbmodel python scripts/validate_apy_early_review.py`
  - Passed. Low-yield posterior mean `0.03547451899578974`; high-yield
    posterior mean `0.09998081433369738`.
- `conda run --no-capture-output -n tbmodel python -m py_compile pages/5_Decision_Analysis.py engine/apy/decision_analysis.py engine/apy/sensitivity.py engine/apy/early_review.py`
  - Passed.
- `conda run --no-capture-output -n tbmodel python -m unittest tests.test_apy_decision_analysis.ApyDecisionAnalysisPageSmokeTests -v`
  - 1 test passed.
- `conda run --no-capture-output -n tbmodel python -m unittest tests.test_apy_decision_analysis.ApyDecisionAnalysisWorkbookTests -v`
  - 2 tests passed.

## Defects Found During Review

- None in Checkpoint 1 focused verification.
- Checkpoint 2: threshold test initially performed unnecessary bisection work
  and exceeded the focused-test timeout. Corrected the test to verify the audit
  grid and precondition behavior without expensive root finding.
- Checkpoint 3: zero prevalence could not be passed through APY calibration.
  Added a narrow zero-prevalence projection path for no-infection false-positive
  likelihood checks.
- Checkpoint 3: identical completed/planned coverage performed duplicate
  deterministic projections. Added a local projection cache.
- Checkpoint 3: several Bayesian-mechanics unit tests were initially too
  calibration-heavy. Reworked them to use a deterministic projection stub while
  retaining model-facing boundary coverage and a separate validation script.
- Checkpoint 4: replaced deprecated `pd.np` page code with an explicit NumPy
  import before committing.
- Checkpoint 5: no defects in focused workbook verification.

## Checkpoint Commit Hashes

- Checkpoint 1: `860d080`
- Checkpoint 2: `d4c7888`
- Checkpoint 3: `9e5f4a4`
- Checkpoint 4: `13cc727`
- Checkpoint 5: pending commit.

## Remaining Work

- Checkpoint 6: adversarial review, validation script, final hardening.

## Blockers

- External scientific evidence remains unresolved by design; software work can
  continue with unresolved outputs clearly marked.

## Final Verification Results

- Pending.
