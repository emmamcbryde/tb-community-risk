# Explicit recent/remote TBI v2 branch status

Branch: `feature/explicit-recent-remote-tbi-v2`

Frozen source release:

- branch: `release/sa-health-apy-he-v1.0.0`
- tag: `sa-health-apy-he-v1.0.0`
- commit: `03cc16e4a52a55e10019dab16bd4f599571c4b90`

Milestone 1 purpose: create a reviewable mathematical foundation for
explicit, mutually exclusive recent and remote-only TB infection calibration
without changing the SA Health Streamlit release, runner, event ledger,
economics, DALYs, MATLAB code, dynamic model, frozen artifacts, deployment,
release branch or release tag.

## Implemented in Milestone 1

New pure module:

- `engine/apy/explicit_recent_remote_tbi.py`

New tests:

- `tests/test_explicit_recent_remote_tbi.py`

New documentation:

- `docs/explicit_recent_remote_tbi_milestone1_spec.md`
- `docs/explicit_recent_remote_tbi_v2_status.md`

The pure module has no Streamlit, economics, DALY, runner, event-ledger,
MATLAB or dynamic-model dependency. It implements:

- age-specific recent and remote exposure durations;
- constant-window cumulative hazards;
- mutually exclusive recent, remote-only and uninfected probabilities;
- population-weighted expected prevalences;
- deterministic feasibility assessment;
- deterministic bisection calibration for `lambda_recent` and
  `lambda_remote`;
- structured calibration result serialization;
- infection-time specification helpers;
- active-TB observation validation for future calibration targets.

Selected identifiers:

- `explicit_recent_remote_tbi_foundation_v1`
- `explicit_recent_remote_tbi_history_v1`
- `constant_hazard_age_window_recent_remote_v1`
- `recent_remote_tbi_history_contract_v1`
- `explicit_recent_remote_tbi_calibration_v1`
- `active_tb_observation_targets_v1`

## Audit findings to preserve

The inherited code already has recent/remote-like language, but it is not the
new mutually exclusive recent-acquisition and remote-only prevalence
calibration:

- `engine/apy/ltbi_state.py` contains the older
  `continuous_markov_recent_remote` progression-state model.
- `engine/apy/infection_history.py` contains an experimental historical
  infection-pressure diagnostic that derives recent fraction among prevalent
  LTBI under assumed trajectory scenarios.
- `engine/apy/calibration.py`, `engine/apy/simulation.py` and
  `engine/apy/expected_value.py` preserve MATLAB-v9-compatible and
  progression-state behavior.
- `app/state.py` sanitizes unsupported experimental state from the standard SA
  Health workflow.
- `engine/apy/frozen_reference.py` remains the loader for frozen SA Health
  stochastic reference outputs and economics.

Risk-factor semantics remain separated:

- `infOR` is used by inherited infection-prevalence compatibility calibration.
- `diseaseOR` / `disOR` are applied as progression hazard multipliers despite
  odds-ratio naming.
- The new Milestone 1 calibration uses age/time alive only and does not use
  disease-progression ORs as acquisition multipliers.

## Not done in Milestone 1

Do not assume the new module is wired into the model. It is intentionally not
connected to:

- Streamlit UI;
- `run_expected_value`;
- `run_replicates`;
- population generation;
- event ledger generation;
- health economics;
- DALYs;
- MATLAB;
- dynamic transmission model;
- frozen reference loaders.

No release branch or tag should be moved. No deployment should be updated.

## Future Milestone 2 recommendation

Milestone 2 should add configuration fields for the two mutually exclusive
population targets:

- recently infected within five years;
- remote infection only.

It should then integrate the calibrated hazards into population generation and
state assignment with explicit infection-time draws, update metadata and cache
keys, and add deterministic and stochastic reconciliation tests. It should keep
the SA Health compatibility workflow isolated unless a reviewed migration plan
explicitly changes that release behavior.

The model remains for planning and sequencing, not for denying care.
