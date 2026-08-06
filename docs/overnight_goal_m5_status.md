# Overnight Goal M5 Status

## Objective

Build and harden a readiness-aware clinician decision-analysis workflow v1 for
the APY LTBI screening model, using the existing APY configuration, event-ledger,
health-economics and evidence-readiness contracts.

## Starting Commit

`a5cc400`

## Current Checkpoint

Checkpoint 1 - scenario comparison contract and runner complete.

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

## Focused Tests Run

- `conda run --no-capture-output -n tbmodel python -m unittest tests.test_apy_decision_analysis -v`
  - 4 tests passed.

## Defects Found During Review

- None in Checkpoint 1 focused verification.

## Checkpoint Commit Hashes

- Checkpoint 1: pending commit.

## Remaining Work

- Checkpoint 2: one-way sensitivity and threshold analysis.
- Checkpoint 3: early screening review and prevalence updating.
- Checkpoint 4: clinician decision-analysis interface.
- Checkpoint 5: workbook/download reproducibility.
- Checkpoint 6: adversarial review, validation script, final hardening.

## Blockers

- External scientific evidence remains unresolved by design; software work can
  continue with unresolved outputs clearly marked.

## Final Verification Results

- Pending.
