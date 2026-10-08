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

Milestone 2C purpose: add explicit progression-calibration policy contracts,
prospective target eligibility, one-parameter fixed-ratio calibration,
likelihood helpers, identifiability enforcement and risk-factor safety
diagnostics without wiring the pathway into runner, UI, population generation,
interventions, event ledgers, economics, DALYs, MATLAB, frozen-reference
loading or dynamic-transmission code.

Milestone 2D purpose: correct calibration-audit issues and perform a focused
evidence review before runner integration. This separates denominator
concepts, corrects competing-mortality mathematics, makes ascertainment
identifiability explicit, strengthens practical identifiability checks,
recommends `none` as the initial production risk-factor policy, and documents
candidate natural-history parameterisations without wiring the pathway into
runner, UI, population generation, interventions, event ledgers, economics,
DALYs, MATLAB, frozen-reference loading or dynamic-transmission code.

Milestone 3A purpose: implement the approved central time-since-infection
progression curve and survivor-conditioned future-risk mathematics as isolated
pure code. This retains five-year recent infection classification, uses no
risk-factor progression multipliers in the central policy, keeps active-TB
observations validation-only by default, supports reset and no-reset
reinfection policies, and leaves competing mortality as an explicit external
input. It still does not wire the pathway into runner, UI, population
generation, interventions, event ledgers, economics, DALYs, MATLAB,
frozen-reference loading or dynamic-transmission code.

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
- `explicit_recent_remote_tbi_progression_calibration_policy_v1`
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

## Implemented in Milestone 2C

New work remains isolated in:

- `engine/apy/explicit_recent_remote_progression.py`
- `tests/test_explicit_recent_remote_tbi.py`
- `docs/explicit_recent_remote_tbi_milestone1_spec.md`
- `docs/explicit_recent_remote_tbi_v2_status.md`

Implemented 2C pieces:

- versioned progression-calibration policy contract with deterministic JSON
  serialization and hashing;
- Policy A `external_progression_hazards_v1`: supplied early/remote hazards,
  no fitting;
- Policy B `fixed_early_remote_ratio_fit_scale_v1`: `lambda_L=k` and
  `lambda_E=Rk`, with one fitted scale only;
- Policy C `validation_only_v1`: no fitting for retrospective, prevalence,
  screen-detected or unsuitable observations;
- Policy D `joint_early_remote_hazards_v1`: specified but unavailable unless
  sensitivity-rank identifiability criteria pass;
- prospective observation eligibility diagnostics with human-readable reasons;
- expected incident-case calculations retaining numerator, denominator,
  horizon and explicit ascertainment;
- deterministic bounded bisection for Policy B, with infeasible targets
  rejected rather than clipped;
- binomial and Poisson log-likelihood helpers requiring an explicit
  observation-model identifier;
- identifiability checks showing one aggregate target cannot identify both
  early and remote hazards, and duplicate-composition targets do not create
  false identifiability;
- competing-mortality interface accepting no-mortality mode, external
  constant death hazards or validated external survival curves;
- baseline active-TB strata excluded from prospective incident calibration;
- risk-factor multiplier policies: `none`,
  `reviewed_hazard_multipliers`, and
  `legacy_or_as_hazard_diagnostic_only`;
- inherited OR-as-hazard multiplication exposed only as a diagnostic, with
  warnings/blocking review status when configured thresholds are exceeded;
- concentration diagnostics reporting expected-case shares from the highest
  1%, 5% and 10% of eligible population weight;
- synthetic worked examples for eligible, ineligible, zero, impossible,
  incomplete-ascertainment, survival and high-multiplier cases.

Worked diagnostic highlights:

- external hazards `lambda_E=0.02`, `lambda_L=0.002`: expected cases
  `2.379733`, no fitting;
- fixed ratio `R=10`: fitted `lambda_L=0.00167858`,
  `lambda_E=0.01678576`, achieved cases `2.000000`;
- moderate reviewed multiplier example: fitted `lambda_L=0.00091672`;
- legacy OR-as-hazard diagnostic: fitted `lambda_L=0.00033747`, with the
  highest 1% by weight contributing `0.832978` of expected cases;
- impossible high target: rejected as `infeasible_above_achievable_range`,
  maximum achievable cases `300`;
- zero-case target: fitted early and remote hazards are both zero.

## Implemented in Milestone 2D

New work remains isolated in:

- `engine/apy/explicit_recent_remote_tbi.py`
- `engine/apy/explicit_recent_remote_progression.py`
- `tests/test_explicit_recent_remote_tbi.py`
- `docs/explicit_recent_remote_tbi_milestone1_spec.md`
- `docs/explicit_recent_remote_tbi_v2_status.md`

Implemented 2D pieces:

- active-TB observations now preserve `sourcePopulationDenominator` alongside
  the existing `populationDenominator`;
- optional numerator-composition fields are validated:
  `numeratorIncludesBaselineActiveTB`,
  `numeratorIncludesPrevalentCases`, and `baselineActiveTBCount`;
- prospective expected-case outputs now separately report
  `sourcePopulationDenominator`, `baselineActiveTBCount`,
  `prospectiveAtRiskPopulation`, `personTimeAtRisk`, and
  `tbiEligiblePopulation`;
- baseline active-TB strata remain excluded from TBI/prospective progression
  calculations, but the observed source denominator is not changed;
- fitting is blocked when baseline active TB is present and numerator
  composition is ambiguous or explicitly includes baseline/prevalent disease;
- competing mortality now uses the cause-specific cumulative-incidence
  integral `integral S_TB(t) S_D(t) h_TB(t) dt`;
- constant external death hazards are integrated analytically across the
  early/late TB-hazard boundary;
- external survival curves are validated to start at one, be finite, remain
  within `[0,1]`, be non-increasing and cover the horizon, then integrated
  numerically;
- scalar horizon survival is rejected for production-capable expected-case
  calculations unless explicitly labelled as a non-production approximation;
- Policy B now requires externally fixed ascertainment `q` in `(0,1]` with
  source and review status; `q=1` must explicitly state complete
  ascertainment;
- a diagnostic records that `q` and progression scale `k` are confounded in
  `E[C] = q E[C_true(k)]` if `q` is not fixed externally;
- Policy D practical-identifiability diagnostics now include sensitivity rank,
  singular values, condition number, composition contrast and notes that
  profile-likelihood/boundary review remains required before production;
- near-collinear target compositions fail practical-identifiability review
  even when floating-point rank is two;
- the recommended initial production risk-factor application policy is
  `none`;
- `legacy_or_as_hazard_diagnostic_only` remains available for audit but is
  blocked from production/default runner policy selection;
- the evidence table and candidate natural-history parameterisations were
  added to the scientific specification.

Follow-up after internet restoration on 2026-10-08:

- the focused evidence-review section was rechecked against live sources and
  augmented with extracted notes on WHO TPT guidance, time-since-infection
  decline, diabetes, smoking, alcohol, CKD/dialysis, chronic airway disease,
  child close-exposure evidence and mortality table sources;
- this follow-up changed documentation only and did not alter pure functions,
  tests, runner, UI, ledger, economics, DALY, MATLAB, frozen-reference or
  dynamic-model code.

Decision-support follow-up after Milestone 2D:

- `docs/explicit_recent_remote_tbi_parameter_decision_dossier.md` converts the
  evidence review into a numerical parameter decision dossier for natural
  history, calibration policy, risk-factor policy, reinfection, age effects,
  competing mortality and active-TB observation use;
- the dossier recommends retaining the five-year recent-infection
  classification with two-year and continuous-decline sensitivity analyses,
  using externally supplied hazards plus validation-only active-TB comparisons
  for initial production integration, and using `none` as the central
  risk-factor multiplier policy;
- the dossier also records candidate conservative, central and
  higher-progression parameterisations and the decisions requiring user
  approval before runner integration;
- this follow-up is documentation-only and does not approve any production
  defaults, alter code or connect the pathway to the runner.

Runner integration remains blocked pending explicit user decisions on:

- natural-history shape and whether a two-phase approximation is acceptable;
- recent-window classification and sensitivity analyses;
- progression-calibration policy and target eligibility;
- ascertainment assumptions;
- risk-factor multiplier policy;
- reinfection reset handling;
- competing-mortality source;
- active-TB observation use.

## Implemented in Milestone 3A

New work remains isolated in:

- `engine/apy/explicit_recent_remote_progression.py`
- `tests/test_explicit_recent_remote_tbi.py`
- `docs/explicit_recent_remote_tbi_milestone1_spec.md`
- `docs/explicit_recent_remote_tbi_v2_status.md`

Implemented 3A pieces:

- versioned progression-curve contract
  `explicit_recent_remote_tbi_progression_curve_v1`;
- central curve identifier
  `central_piecewise_time_since_infection_progression_v1`;
- conservative and higher-progression sensitivity identifiers
  `conservative_piecewise_time_since_infection_progression_v1` and
  `higher_progression_piecewise_time_since_infection_progression_v1`;
- deterministic canonical JSON serialization and SHA-256 hashing for curve
  contracts;
- cumulative-risk to cumulative-hazard conversion with
  `H(t) = -log[1-F(t)]`;
- strict validation of increasing time anchors, non-decreasing risks and
  hazards, non-negative segment hazards and explicit post-final-anchor hazard;
- central cumulative-risk anchors:
  `0.038` at 1 year, `0.050` at 2 years, `0.066` at 5 years,
  `0.072` at 10 years and `0.079` at 25 years;
- transformed central cumulative-hazard anchors:
  `0.0387408283`, `0.0512932944`, `0.0682788408`,
  `0.0747235462` and `0.0822952427`;
- derived central segment hazards:
  `0.0387408283`, `0.0125524661`, `0.0056618488`,
  `0.0012889411` and `0.0005047798` per year;
- central post-final-anchor hazard `0.0005047798` per year, explicitly
  extrapolated from the 10-25 year cumulative-hazard slope;
- survivor-conditioned future risk:
  `P(T <= t | T > s) = 1 - exp{-[H(s+t)-H(s)]}`;
- segment exposure diagnostics and instantaneous segment-hazard lookup;
- reset policy `recent_reinfection_resets_progression_clock_v1`;
- no-reset sensitivity policy
  `recent_reinfection_does_not_reset_progression_clock_v1`;
- pure curve expected-event aggregation with
  `riskFactorProgressionPolicy = none`;
- regression coverage proving inherited risk-factor flags do not alter central
  curve results;
- validation-only comparison of eligible prospective active-TB observations
  against a fixed curve, preserving source denominator, at-risk population,
  observation horizon, ascertainment and optional likelihood;
- curve competing-risk integration with a supplied external survival curve or
  external constant death hazard, without horizon-level survival
  multiplication;
- worked central diagnostic table showing that a person infected 4.9 years
  before baseline has 0.686% five-year survivor-conditioned future risk, not a
  fresh 6.6% five-year risk.

Milestone 3A did not approve risk-factor multipliers, did not fit the curve to
active-TB observations, did not bundle mortality data, and did not touch the
inherited MATLAB-v9 compatibility path.

Milestone 2D retained the risk diagnostic that inherited OR-labelled factors
can multiply to `2916`, and retained the synthetic example in which the
highest 1% of population weight contributes about `0.832978` of expected
cases. This remains a warning that shrinking the fitted baseline hazard does
not correct misspecified relative-risk concentration.

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

## Not done in Milestones 1, 2A, 2B, 2C and 2D

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

This also remains true after Milestone 2C. The policy/fitting helpers are
pure diagnostics and are not imported or consumed by runner, UI, population
generation, intervention, event-ledger, economics, DALY, frozen-reference,
MATLAB or dynamic-model code.

This also remains true after Milestone 2D. Denominator, ascertainment,
competing-mortality, identifiability, evidence-review and risk-policy audit
corrections are still pure explicit-pathway helpers and documentation only.

This also remains true after Milestone 3A. Time-since-infection progression
curves, conditional future-risk calculations, reinfection-clock policies,
curve validation outputs and curve mortality integration remain isolated pure
helpers and are not imported or consumed by the runner, UI, population
generation, interventions, event ledger, economics, DALY, frozen-reference,
MATLAB or dynamic-model code.

No release branch or tag should be moved. No deployment should be updated.

## Future runner-integration milestone recommendation

The next milestone should begin runner integration behind explicit
new-pathway metadata and cache keys, using the 3A central curve and
validation-only active-TB policy exactly as approved. It should keep the
inherited MATLAB-v9 compatibility pathway available and unchanged, and it
should not proceed to event-ledger, economics or DALY integration unless those
surfaces are explicitly scoped.

The model remains for planning and sequencing, not for denying care.
