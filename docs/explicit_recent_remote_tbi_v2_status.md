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

Milestone 2A purpose: add isolated configuration, deterministic assignment,
stochastic state assignment and infection-time timing helpers without wiring
the pathway into any runner, UI, event ledger, economics, DALY, MATLAB,
frozen-reference or dynamic-transmission code.

Milestone 2B purpose: add isolated natural-history progression specification
and pure expected-value mathematics for the explicit baseline states without
wiring the pathway into the runner, UI, event ledger, intervention logic,
economics, DALYs, MATLAB, frozen-reference loading or dynamic-transmission
code.

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
- `explicit_recent_remote_tbi_config_v1`
- `explicit_recent_remote_tbi_assignment_v1`
- `explicit_recent_remote_tbi_progression_v1`
- `constant_recent_window_hazard_v1`
- `constant_remote_window_hazard_v1`

## Implemented in Milestone 2A

The new work remains in:

- `engine/apy/explicit_recent_remote_tbi.py`
- `tests/test_explicit_recent_remote_tbi.py`
- `docs/explicit_recent_remote_tbi_milestone1_spec.md`
- `docs/explicit_recent_remote_tbi_v2_status.md`

Implemented 2A pieces:

- isolated `ExplicitRecentRemoteConfig` contract with enabled status, recent
  and remote-only total-population targets, window settings, hazard-shape
  identifiers, source, reference year, review status, notes and contract IDs;
- deterministic canonical JSON serialization and SHA-256 hashing for that
  contract;
- explicit rejection of incompatible retired/frozen identifiers in the new
  config contract;
- calibration/window plumbing so feasibility, calibration and population
  prevalence helpers use supplied recent-window and remote-cap values;
- deterministic assignment summaries for expected recent, remote-only,
  uninfected, prior-remote-exposure and prior-remote-plus-recent counts and
  proportions;
- age-specific deterministic outputs for recent prevalence in the total age
  group and recent fraction among people with TBI;
- stable inverse-CDF sampling of time since the most recent infection event
  conditional on at least one event in the eligible window;
- stochastic assignment requiring an explicit seed or RNG, using no global RNG
  state;
- effective-state precedence: any recent infection yields `recent`, remote
  event(s) without recent infection yields `remote_only`, neither yields
  `uninfected`;
- prior remote exposure retained for recent reinfection without double-counting
  remote-only prevalence;
- `remainingEarlyRiskYears = max(0, recentWindowYears -
  timeSinceMostRecentInfection)`, with zero remaining early-risk years for
  remote-only and uninfected states;
- documented baseline active-TB sequencing proposal that keeps prevalent
  active TB, screen-detected active TB and incident active TB distinct.

## Implemented in Milestone 2B

New pure module:

- `engine/apy/explicit_recent_remote_progression.py`

Implemented 2B pieces:

- resolved recent-window consistency so the recent/remote classification
  boundary and remaining early-risk duration use the same `recentWindowYears`;
- recorded `age85PlusMax` and age-support provenance in the isolated
  configuration contract and deterministic hash;
- added age-support calibration sensitivity diagnostics showing the inherited
  `85+` expansion cap can affect fitted remote hazard;
- added pure prospective progression functions for cumulative hazard,
  piecewise hazard, survival probability, cumulative incidence and exact
  quantile inversion;
- added expected progression-event aggregation across weighted strata with
  explicit ascertainment probability;
- classified active-TB observation rows as baseline/prevalence,
  screen-detected, prospective incident, retrospective notification/incidence,
  or mixed/insufficiently defined;
- added a pure expected-count diagnostic only for genuinely prospective
  incident observations beginning at model baseline;
- documented that one aggregate active-TB count generally cannot identify both
  early and late progression hazards;
- listed future calibration policies without selecting a hidden default;
- documented historical notification limitations and the inherited absence of
  competing mortality in no-transmission progression;
- added worked synthetic diagnostic tables for reviewable progression
  probabilities.

The pure progression equations are:

```text
recent:      A(t) = m [lambda_E min(t, r) + lambda_L max(0, t-r)]
remote_only: A(t) = m lambda_L t
uninfected:  P(T <= t) = 0 under the current no-new-infection pathway
P(T <= t) = 1 - exp[-A(t)]
```

where `r = max(0, recentWindowYears - timeSinceMostRecentInfection)`.

## Age representation audit

The inherited Python APY path represents ages for calibration and assignment
as exact integer ages with probabilities:

- source age tables may be broad bands;
- `expand_age_distribution_table` spreads each band uniformly across integer
  ages in that band;
- `85+` expands through `age85PlusMax`, default 89;
- stochastic cohort generation samples from `exactAgeValues` and
  `exactAgeProb`;
- deterministic expected-value strata iterate the same exact-age support.

Milestone 2A assignment functions accept exact ages and weights and do not
use band midpoints silently. If a future integration receives only broad bands,
it must call or document an explicit within-band expansion before calibration
and assignment.

Milestone 2B makes the open-ended cap explicit in the new config provenance.
The default `age85PlusMax = 89` is an inherited implementation choice, not a
reviewed maximum age. A synthetic all-85-plus sensitivity changed fitted remote
hazard from `0.00666571` at support 85-89 to `0.00643236` at support 85-95
for the same 40% remote-only target.

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
- The inherited default active-TB calibration target is `10/770` and
  `earlyLateRatio` sets `lambdaLate = lambdaEarly / earlyLateRatio`; those
  are not automatically reused by the explicit pathway.
- `app/state.py` sanitizes unsupported experimental state from the standard SA
  Health workflow.
- `engine/apy/frozen_reference.py` remains the loader for frozen SA Health
  stochastic reference outputs and economics.

Risk-factor semantics remain separated:

- `infOR` is used by inherited infection-prevalence compatibility calibration.
- `diseaseOR` / `disOR` are applied as progression hazard multipliers despite
  odds-ratio naming.
- default normalized disease multipliers are jointly multiplied and can reach
  `2916` if all flags are present; no cap was found;
- The new Milestone 1 calibration uses age/time alive only and does not use
  disease-progression ORs as acquisition multipliers.
- The Milestone 2B pure progression functions accept a generic non-negative
  multiplier but do not endorse odds ratios as hazard multipliers.

## Active-TB and identifiability status

Milestone 2B keeps baseline active TB conceptually separate from TBI. Future
integration must assign baseline/prevalent active TB before latent-state
assignment; those people must not enter the ordinary TBI preventive-treatment
cascade or be counted simultaneously in TBI prevalence totals.

Prospective incident active-TB observations can be used by the new pure
expected-count diagnostic only when the observation window begins at model
baseline and the baseline population composition is defined. Retrospective
notification counts remain validation data unless a separate retrospective
population reconstruction is implemented.

No production progression-calibration policy is chosen. Supported future
policy options are external hazards, one supplied hazard plus one fitted
hazard, externally fixed early-to-late ratio with common scale, multiple
targets in a joint likelihood, or external validation only.

Competing mortality is not present in the inherited no-transmission
progression path. This likely overstates 20-year prospective progression in
older groups; future integration should allow an explicit survival function.

## Not done in Milestones 1, 2A and 2B

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

This remains true after Milestone 2A. The new assignments are not consumed by:

- `engine/apy/runner.py`;
- `engine/apy/simulation.py`;
- `engine/apy/expected_value.py`;
- event-ledger generation;
- health economics;
- DALYs;
- Streamlit pages;
- frozen-reference loading;
- MATLAB;
- `engine/dynamic/*`.

This also remains true after Milestone 2B. The new progression module is not
imported or consumed by runner, UI, event-ledger, economics, DALY,
frozen-reference, MATLAB or dynamic-model code.

No release branch or tag should be moved. No deployment should be updated.

## Future Milestone 2C recommendation

Milestone 2C should choose a reviewed progression-calibration policy or
explicit external progression hazards before connecting the new baseline
states to the APY runner. The next integration milestone should still preserve
the frozen SA Health compatibility workflow, update cache keys and metadata
only for the new pathway, and keep event-ledger, economics and DALY integration
separate unless explicitly scoped.

The model remains for planning and sequencing, not for denying care.
