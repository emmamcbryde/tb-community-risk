# Explicit recent/remote TBI calibration - Milestone 1 specification

This branch starts a new model-development line derived from the frozen SA
Health Streamlit release. It does not modify the submitted SA Health release
and does not expose a new pathway in the standard Streamlit workflow.

The Milestone 1 scope is limited to audit, mathematical specification, pure
calibration functions and tests. The epidemiological runner, event ledger,
health economics, DALYs, MATLAB, dynamic transmission model, frozen artifacts,
deployment, release branch and release tag are intentionally unchanged.

## Selected identifiers

- Analysis basis: `explicit_recent_remote_tbi_foundation_v1`
- Natural-history semantics: `explicit_recent_remote_tbi_history_v1`
- Calibration method: `constant_hazard_age_window_recent_remote_v1`
- Infection-history contract: `recent_remote_tbi_history_contract_v1`
- Calibration contract: `explicit_recent_remote_tbi_calibration_v1`
- Active-TB observation schema: `active_tb_observation_targets_v1`

These are new identifiers. They do not reuse or reactivate the retired
`continuous_markov_recent_remote` pathway and do not replace the frozen
`sa_health_matlab_v9_compatibility_reference` / `matlab_v9_implicit_early_late`
metadata.

## Definitions

At baseline, each person has one effective infection state:

1. `uninfected`: no infection in either exposure window.
2. `remote_only`: at least one infection acquired more than five years before
   baseline and no infection during the most recent five years.
3. `recent`: at least one TB infection acquired during the five years
   immediately before baseline.

If both remote and recent infection occurred, the effective state is `recent`.
The implementation retains the ability to record prior remote exposure, but
does not double-count that person in population prevalence totals.

The user-facing target labels are:

- Recently infected within 5 years
- Remote infection only — infected 5 or more years ago with no infection in
  the last 5 years

## Age windows

Let current age be `a_i` years and `u` be years before baseline. No exposure is
allowed before birth.

Recent cumulative infection hazard:

```text
H_recent,i = integral_0^min(5, a_i) lambda_recent,i(u) du
```

Remote cumulative infection hazard:

```text
H_remote,i = integral_5^min(100, a_i) lambda_remote,i(u) du
```

For people aged five years or younger, `H_remote,i = 0`.

The initial implementation uses constant hazards within each window:

```text
H_recent,i = lambda_recent * min(5, a_i)
H_remote,i = lambda_remote * max(0, min(100, a_i) - 5)
```

The 100-year limit is a practical historical cap. Ages greater than 100 retain
five recent years and 95 remote years; no additional remote exposure is added
beyond the cap.

## Mutually exclusive probabilities

For each person or age stratum:

```text
P(recent) = 1 - exp(-H_recent)

P(remote only) = exp(-H_recent) * [1 - exp(-H_remote)]

P(uninfected) = exp[-(H_recent + H_remote)]
```

The implementation validates that probabilities are finite, lie in `[0, 1]`,
sum to one within numerical tolerance and have no remote exposure at age five
or younger.

Prior remote exposure is also available as an internal diagnostic:

```text
P(prior remote exposure) = 1 - exp(-H_remote)
P(recent with prior remote exposure) = P(recent) * P(prior remote exposure)
```

Those quantities are not added to total TBI prevalence.

## Population targets and calibration

Future UI inputs are two mutually exclusive proportions of the total
population:

1. Proportion recently infected within five years.
2. Proportion with remote-only infection.

They are not proportions among people with TBI.

For an age distribution with weights `w_i`, the pure calibration solves for
`lambda_recent` first, then `lambda_remote` conditional on the fitted recent
hazard:

```text
sum_i w_i * P_i(recent) = target_recent

sum_i w_i * P_i(remote only) = target_remote_only
```

This sequence is mathematically consistent because `P(remote only)` includes
the probability of no recent infection. The solver does not treat the remote
target as an independent raw prevalence.

The root-finding is deterministic bisection with explicit tolerances and no
random draws. Zero recent and zero remote targets return zero hazards. Feasible
near-boundary targets are bracketed by expanding the upper hazard bound.
Impossible targets raise `CalibrationError` with a feasibility assessment;
targets are not clipped.

The structured result records requested targets, fitted hazards, achieved
recent, remote-only, total TBI and uninfected prevalences, absolute residuals,
convergence and feasibility status, diagnostic messages, selected identifiers
and target definitions.

## Feasibility

The age distribution must be non-empty, finite, non-negative and sum to one
within `1e-9`.

Feasibility checks include:

- `target_recent` and `target_remote_only` must be probabilities.
- Their sum must not exceed one.
- `target_recent` cannot exceed the population weight with any recent exposure
  duration.
- After fitting the recent hazard, `target_remote_only` cannot exceed the
  population weight with remote exposure, weighted by the probability of no
  recent infection.

The exact upper boundary where a target equals the limiting prevalence may
require an infinite hazard. This implementation supports robust near-boundary
finite targets and reports unbracketed exact-limit cases rather than silently
inventing a finite value.

## Age interpretation

Two quantities must not be described as equivalent:

- probability of recent infection in the total age group;
- proportion recent among all people with TBI in that age group.

Under a constant recent hazard, people older than five have the same
total-population probability of recent infection, because each has five recent
exposure years. However, the proportion recent among all people with TBI
generally declines with age because remote infection accumulates.

## Infection timing specification

Milestone 1 specifies timing but does not wire it into the runner.

Conditional on the effective state:

- `recent`: draw the most recent infection event within `[0, min(5, age)]`
  years before baseline from the constant-hazard distribution truncated to the
  available recent window.
- `remote_only`: draw the most recent remote-window infection event within
  `[5, min(100, age)]` years before baseline from the constant-hazard
  distribution truncated to the remote window, conditional on no recent event.
- both remote and recent: classify effective state as `recent`, retain
  `priorRemoteExposure=true`, and treat recent reinfection as resetting the
  higher-progression-risk state in future integration unless a later scientific
  decision chooses a different mechanism.

Audit finding: the inherited runner currently uses a baseline recent flag plus
a Markov recent-to-remote progression compartment. It does not distinguish
first infection, most recent infection and reinfection reset mechanisms. That
mechanism remains unchanged in this milestone.

## Risk-factor separation

Milestone 1 does not use non-age risk factors in the new acquisition
calibration. Age affects exposure only through time alive.

Existing code contains several different semantics:

- `riskPrev`: prevalence of risk factors in the generated population.
- `infOR`: used in the inherited compatibility infection-prevalence model for
  marijuana, contact and renal flags.
- `diseaseOR` / `disOR`: named as odds ratios but applied as multiplicative
  progression hazard multipliers in calibration, cohort simulation and
  expected-value pathways.
- `targetAgeOR`: age odds ratio target for inherited LTBI prevalence
  calibration, not a recent-infection target.
- `baselineRecentLTBIProportion`: inherited progression-state fraction, not a
  mutually exclusive recent-acquisition prevalence target.

The new pure calibration does not use inherited disease-progression ORs as
acquisition multipliers. The inherited semantics are documented here as risks
for future integration, but not changed.

## Active-TB observation schema

The optional active-TB observation schema supports one or more rows with:

- start date/year and end date/year;
- observed active-TB case count;
- population denominator;
- optional person-years;
- denominator type;
- whole-population versus screened-subgroup scope;
- prevalent/baseline, screen-detected or follow-up incident classification;
- passive notification, active screening, prevalence survey or combined
  ascertainment;
- pulmonary/all/unknown active-TB classification;
- source, review status, notes and uncertainty fields.

Validation requires positive denominators, non-negative counts, coherent dates,
valid enumerated meanings, source and review status. Display rates per 100,000
population and per 100,000 person-years are calculated without discarding the
underlying numerator, denominator or observation period.

Prevalent active TB at or near baseline, active TB detected during screening
and incident active TB during later follow-up remain distinct. The schema does
not infer one from another and does not fit progression hazards in Milestone 1.

## Identifiability

Recent and remote-only TBI targets identify the two infection-hazard scales
only under the specified historical hazard shapes and supplied age
distribution.

Active-TB observations do not directly identify TBI prevalence. TBI targets do
not by themselves identify disease-progression hazards. Baseline/prevalent
active TB cannot automatically be treated as future incident progression.
Risk-factor progression effects do not identify acquisition effects.
Uncertainty in age distribution and target estimates should eventually be
propagated.

Future progression calibration would require combinations such as:

- explicit recent and remote-only TBI prevalence targets;
- active-TB observations classified by timing and ascertainment;
- fixed or estimated progression-hazard structure by infection state;
- fixed test/screening detection assumptions for screen-detected active TB;
- explicit handling of prevalent baseline active TB versus incident follow-up;
- assumptions about risk-factor effects on progression, separate from
  acquisition.

## Architecture audit

Existing recent/remote and early/late implementations:

- `engine/apy/ltbi_state.py` defines `continuous_markov_recent_remote`, a
  recent-to-remote latent progression compartment with a five-year mean
  residence time.
- `engine/apy/infection_history.py` defines an experimental historical
  infection-pressure pathway. It derives a recent fraction among prevalent
  LTBI under scenario trajectory assumptions. It is not exposed in the
  reference-only Streamlit workflow.
- `engine/apy/calibration.py` preserves MATLAB-v9-compatible calibration and
  calibrates progression hazards. `matlab_v9_implicit_early_late` is a
  compatibility semantic, not measured recent infection.
- `engine/apy/simulation.py` draws baseline recent/remote progression state
  and active-TB times. The MATLAB-v9 compatibility branch starts all infected
  people in the early state and transitions those not active by the screening
  window to remote.
- `engine/apy/expected_value.py` mirrors the deterministic pathway and uses
  recent/remote progression-state probabilities for event ledgers.

Retired pathway:

- The old experimental Streamlit controls are sanitized out by
  `app/state.py::sanitize_reference_only_state`.
- Standard pages reject unsupported result metadata and do not expose the
  experimental trajectory controls.
- This milestone does not reactivate that pathway.

Compatibility and caches:

- Runner calibration cache: `engine/apy/runner.py::_CALIBRATION_CACHE`.
- Reference calibration artifact cache:
  `engine/apy/calibration_policy.py::_REFERENCE_ARTIFACT_CACHE`.
- Experimental infection-history caches:
  `engine/apy/infection_history.py::_CALIBRATION_CACHE` and
  `_AGE_RISK_GRID_CACHE`.
- Frozen reference loader cache: `engine/apy/frozen_reference.py::_load_payload`.
- Streamlit backend resources are cached in `app/state.py`.

Scenario/configuration schemas:

- `engine/apy/config.py`, `engine/apy/scenario.py`, `app/parameter_workspace.py`
  and validation modules preserve the current v9 scenario structure.
- `baselineRecentLTBIProportion` remains an inherited progression-state field.
  It is not repurposed for the new mutually exclusive recent-acquisition
  target in Milestone 1.

Population generation:

- `engine/apy/cohort.py::draw_base_population` draws age, sex, BCG and risk
  factors, then uses the inherited infection probability and disease
  multiplier functions.
- New pure calibration functions accept only explicit ages and weights.

Event ledger and economics:

- `engine/apy/event_ledger.py`, `engine/apy/event_ledger_economics.py`,
  `engine/apy/frozen_reference.py`, `engine/apy/economics.py` and DALY-related
  code are unchanged.
- Frozen reference loaders and recalculation logic remain tied to the SA Health
  compatibility identifiers.

Dynamic model:

- `engine/dynamic/*` and integration comparison contracts are outside
  Milestone 1 and unchanged.

## Planned integration

Milestone 2 should add explicit configuration fields for recent and
remote-only targets, then wire the calibrated hazards into population
generation and runner state assignment without changing the frozen SA Health
workflow. It should include a documented infection-time draw, event-ledger
metadata updates, deterministic/stochastic reconciliation and migration tests.

## Limitations and unanswered decisions

- Progression hazards are not calibrated here.
- The runner is not yet able to distinguish first infection, most recent
  infection and recent reinfection reset mechanisms.
- No risk-factor acquisition effects are included beyond age/time alive.
- Exact limiting targets may imply infinite hazards.
- Active-TB observations are validated for future use but not fitted.
- Parameter uncertainty is not propagated.
- The new pathway is not decision-ready and is not for denying care.
