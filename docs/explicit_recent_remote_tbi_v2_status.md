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

Milestone 3B purpose: add a guarded, isolated agent-based runner integration
for the explicit recent/remote pathway. The pathway is activated only by
`analysisPathway = explicit_recent_remote_tbi_v2` plus an enabled
`explicitRecentRemoteTBI` runner contract. It calibrates recent and
remote-only infection hazards from mutually exclusive total-population
targets, assigns baseline states in the stochastic cohort, and samples
untreated active-TB times from the approved survivor-conditioned
time-since-infection curve. It does not add Streamlit controls, report-facing
outputs, event-ledger rows, event-ledger economics, DALYs, MATLAB changes,
dynamic-model changes or active-TB progression calibration.

Milestone 3C purpose: harden the shared runner dispatch and add isolated
non-economic epidemiological intervention outputs for the explicit pathway.
This milestone records the inherited screening/TPT semantics in a versioned
policy contract, uses paired baseline untreated natural histories for
comparator/intervention comparisons, strengthens deterministic-stochastic
reconciliation diagnostics, and keeps event ledgers, economics, DALYs,
Streamlit controls, reports, MATLAB and the dynamic model out of scope.

Milestone 3D purpose: validate explicit intervention outputs and add the
first isolated Streamlit workflow for explicit recent/remote TBI
epidemiology. The new workflow is a separate development entry point
(`explicit_tbi_ui/app.py`, isolated in Milestone 3F) and is not imported by the frozen Streamlit
application. It adds no costable event ledger, economics, DALYs, QALYs,
deployment, merge or tag.

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
  `recent_reinfection_does_not_reset_remote_clock_v1`;
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

## Milestone 3B runner integration status

Milestone 3B intentionally changes the boundary for `run_replicates` only.
The new integration files are:

- `engine/apy/explicit_recent_remote_runner.py`;
- guarded imports in `engine/apy/calibration_policy.py`;
- guarded imports and metadata/cache-key additions in `engine/apy/runner.py`;
- guarded stochastic state and active-time assignment in
  `engine/apy/simulation.py`;
- an explicit rejection guard in `engine/apy/frozen_reference.py`.

Activation rules:

- top-level `analysisPathway` must be
  `explicit_recent_remote_tbi_v2`;
- nested `explicitRecentRemoteTBI.enabled` must be true;
- nested `pathwayIdentifier` must also be
  `explicit_recent_remote_tbi_v2`;
- inconsistent or unknown pathway provenance fails before calibration;
- the presence of the new fields without positive pathway selection does not
  silently activate the pathway;
- the frozen SA Health reference loader rejects explicit-pathway configs.

Runner behaviour:

- baseline demographic and risk-factor population generation still uses the
  inherited `draw_base_population` machinery for age, sex, BCG and risk
  factors;
- explicit-pathway calibration then replaces inherited infection prevalence
  assignment with mutually exclusive `recent`, `remote_only` and `uninfected`
  states from calibrated recent and remote hazards;
- baseline/prevalent active TB is sequenced first when configured, reported
  separately, and excluded from TBI states and the latent-treatment cascade;
- prior remote exposure is retained for recent reinfection audit fields;
- central progression uses
  `central_piecewise_time_since_infection_progression_v1`;
- untreated active-TB times are sampled by inverse cumulative hazard from
  `H(s+t)-H(s)`, so future risk is survivor-conditioned on being disease-free
  at baseline;
- isolated deterministic review helpers integrate over exact age and
  conditional infection time with Gauss-Legendre quadrature rather than using
  mean time since infection;
- isolated stochastic natural-history helpers require an explicit seed or RNG
  and reconcile to deterministic expectations within predeclared Monte Carlo
  tolerances in the focused tests;
- central risk-factor progression policy is `none`, and inherited disease ORs
  are not applied to explicit active-TB progression;
- active-TB observations remain validation-only and do not fit the curve;
- competing mortality remains `not_modelled` unless a later milestone supplies
  a reviewed mortality interface to the runner.
- explicit-pathway `run_replicates` returns a structured
  unavailable event-ledger marker, updated in Milestone 3C to
  `explicit_recent_remote_tbi_v2_event_ledger_out_of_scope_milestone_3c`,
  rather than a costable event-ledger bundle.

Still not done after Milestone 3B:

- no Streamlit controls or report-facing outputs expose the pathway;
- `run_expected_value` and report-facing expected-value outputs are not
  integrated;
- no event-ledger rows, event-ledger economic interpretation or DALY outputs
  are added for explicit-pathway outputs;
- no health economics or DALY code consumes explicit assignments;
- no MATLAB code is changed;
- no dynamic-transmission model code is changed;
- no frozen SA Health artifacts are regenerated;
- no release branch or tag is moved.

No release branch or tag should be moved. No deployment should be updated.

## Superseded 3B runner-integration recommendation

Milestone 3C supersedes the original 3B next-step recommendation by adding the
isolated epidemiological intervention-output contract and the shared-runner
hardening described below. Streamlit controls, report-facing outputs,
event-ledger rows, economics and DALYs still require separate scoped
milestones.

The model remains for planning and sequencing, not for denying care.

## Implemented in Milestone 3C

Milestone 3C retains the explicit pathway identifier
`explicit_recent_remote_tbi_v2` and central natural-history identifier
`central_piecewise_time_since_infection_progression_v1`. The shared runner
dispatch remains fail-closed:

- configurations without `analysisPathway = explicit_recent_remote_tbi_v2`
  remain on the inherited path;
- `explicitRecentRemoteTBI.enabled = true` without the explicit pathway
  identifier raises an error;
- unknown pathway identifiers, including `continuous_markov_recent_remote`,
  are rejected;
- the frozen SA Health reference eligibility guard continues to reject the
  explicit pathway;
- explicit and compatibility results carry distinct metadata and cache-key
  inputs.

The explicit runner contract now also records:

- `tptEffectPolicyIdentifier =
  effective_tpt_prevents_future_progression_v1`;
- `commonRandomNumberContractIdentifier =
  paired_baseline_untreated_progression_shared_v1`;
- deterministic hashes for both contracts.

Inherited intervention semantics audit:

- screening selection uses the inherited strategy and coverage;
- screening time is sampled uniformly over the configured screening window;
- test sensitivity applies to infected screened people, including those
  already active at screening;
- test specificity applies to uninfected screened people;
- baseline active TB is excluded before ordinary latent-TB TPT eligibility;
- active TB present before or at screening enters work-up/exclusion and does
  not receive ordinary latent-TB TPT;
- treatment starts at screening time in the inherited stochastic engine;
- treatment stop/effective time is `screenTime + regimenMonths/12 *
  doseFractionTaken`;
- a successful full or partial TPT efficacy draw prevents only active TB whose
  untreated event time is after the treatment stop/effective time;
- incomplete or ineffective treatment retains the untreated
  survivor-conditioned progression clock;
- uninfected false positives can contribute to screening, test-positive,
  treatment-start and completion counts, but receive no epidemiological TB
  prevention benefit;
- historical recent/remote infection state is retained after effective TPT.

The selected explicit TPT-effect policy is
`effective_tpt_prevents_future_progression_v1`. It does not sterilise or
rewrite historical infection prevalence, does not grant retrospective
protection from baseline, and does not apply any inherited OR-labelled
progression multipliers. This is a runner-integration working policy, not a
new scientific proof of TPT mechanism.

The common-random-number contract shares baseline demographics, baseline
active-TB assignment, recent/remote/uninfected assignment, infection clocks,
untreated active-TB progression potential and mortality histories when
supplied. Intervention-specific draws include screening, testing, treatment
initiation, discontinuation/completion and TPT efficacy.

New non-economic explicit runner output:

- top-level `explicitEpidemiology` is returned only for the explicit pathway;
- comparator and intervention summaries include population, baseline active
  TB, recent TBI, remote-only TBI, uninfected, prior-remote-plus-recent,
  screened, test-positive, true-positive, false-positive, TPT starts,
  completions, effective TPT, work-up active TB, timing categories for active
  TB and direct active TB cases averted;
- outputs are stratified by age group and effective state where available;
- the case-averted definition is "direct active TB cases averted among the
  modelled population";
- no costs, DALYs, QALYs or event-ledger rows are generated.

Event-ledger protection:

- explicit-pathway `eventLedger.available` remains false;
- the unavailable ledger marker is
  `explicit_recent_remote_tbi_v2_event_ledger_out_of_scope_milestone_3c`;
- its validation status is false, so event-ledger economics rejects it as not
  authoritative.

Deterministic intervention helper:

- `deterministic_explicit_recent_remote_intervention_expectations` integrates
  over exact age, infection state, conditional infection time and screening
  time with Gauss-Legendre quadrature;
- it records
  `gauss_legendre_infection_and_screening_time_quadrature_v1`;
- it does not replace nonlinear timing distributions with a single mean;
- it applies effective TPT only to future progression after the treatment
  stop/effective time;
- mortality remains `not_modelled` unless an external survival curve or death
  hazard is supplied to the isolated helper.

Reconciliation hardening:

- the deterministic-stochastic reconciliation helper now supports multiple
  prespecified seeds and reports Monte Carlo standard error, standardised
  difference, the acceptance rule and pass/fail for each metric;
- reconciliation rows include baseline active TB, recent, remote-only,
  uninfected, prior-remote-plus-recent and incident TB horizons;
- stochastic natural-history helpers now expose age-specific state and
  age-specific incident-TB diagnostic rows.

Still out of scope after Milestone 3C:

- Streamlit controls and report-facing outputs;
- event-ledger generation for the explicit pathway;
- health economics, DALYs or QALYs for the explicit pathway;
- MATLAB changes;
- dynamic-transmission model changes;
- frozen artifact regeneration or release-pointer changes.

## Implemented in Milestone 3D

Milestone 3D adds a standalone development workflow:

- entry point: originally `explicit_tbi_app.py`; moved to `explicit_tbi_ui/app.py` in Milestone 3F;
- support module: `app/explicit_tbi_workflow.py`;
- tests: `tests/test_explicit_tbi_streamlit_app.py`.

The frozen application entry point remains `streamlit_app.py`. Its navigation
still lists only the existing frozen workflow pages under `pages/`; it does
not import, link to or navigate to the explicit application.

The new workflow pages are:

1. `Set up`;
2. `Run analysis`;
3. `Results`;
4. `Evidence and technical information`.

The workflow is labelled as a development/research version and uses only
`analysisPathway = explicit_recent_remote_tbi_v2`. It rejects incompatible
result objects and does not consume MATLAB-v9 compatibility, frozen-reference
or retired `continuous_markov_recent_remote` results.

Set up controls:

- total modelled population;
- analysis horizon;
- baseline active TB as a state separate from TBI;
- editable recent and remote-only mutually exclusive total-population TBI
  targets;
- screening coverage and screening period;
- test, regimen, TPT initiation, completion and efficacy;
- deterministic preview versus stochastic analysis, with repetitions and seed
  shown only for stochastic mode;
- optional active-TB observation fields for validation only.

Validation behavior:

- recent, remote-only and baseline active-TB proportions must be within
  `[0,1]`;
- their sum must not exceed the total population;
- infeasible recent/remote calibration is blocked before analysis;
- observed active-TB denominators, counts, durations and ascertainment
  probabilities are validated;
- active-TB observations do not recalibrate infection hazards or progression
  anchors.

State and stale-result behavior:

- material epidemiological edits invalidate completed results;
- active-TB observation edits update validation displays only and do not change
  fitted recent/remote hazards, baseline states, comparator TB, intervention
  TB or cases averted;
- `Restore working defaults` is the single reset action and clears stale
  results.

Milestone 3D also closes the intervention-output reconciliation gap:

- `stochastic_explicit_recent_remote_intervention_summary` provides an
  isolated stochastic intervention summary for reconciliation, including
  optional external competing mortality;
- `reconcile_deterministic_stochastic_explicit_intervention` compares
  deterministic expectations with stochastic means using prespecified
  Monte Carlo standard errors or empirical simulation uncertainty;
- reconciliation rows cover eligible population, screened, true positives,
  false positives, active TB detected during work-up, TPT starts,
  completions, effective TPT, comparator TB, intervention TB, direct cases
  averted, timing categories and recent/remote cases averted.

Defects found while closing 3D:

- the first isolated stochastic intervention helper omitted the conventional
  `recentLTBIAtBaseline` and `remoteLTBIAtBaseline` audit fields before
  building summaries; this was corrected;
- the first helper incorrectly limited `effectiveTPT` to people who would
  develop TB by the analysis horizon. It now matches the inherited semantics:
  effective TPT is counted among infected people whose untreated progression
  time is after the treatment effective time, while cases averted remain
  horizon-limited.

Still out of scope after Milestone 3D:

- exposing the pathway through the frozen Streamlit application;
- event-ledger generation or costing for the explicit pathway;
- economics, DALYs, QALYs or ICERs;
- MATLAB changes;
- dynamic-transmission model changes;
- frozen artifact regeneration or release-pointer changes;
- deployment, merge or tag.

## Milestone 3E scientific validation and usability review

Milestone 3E resolves a reporting inconsistency from the Milestone 3D
summary. The reconciliation table with population `1200`, screened about
`720` and true-positive latent tests about `358` was not generated from the
ordinary working defaults. It was a high-prevalence intervention stress-test
configuration:

- population `1200`;
- recent TBI target `20%`;
- remote-only TBI target `30%`;
- screening coverage `60%`;
- IGRA sensitivity `100%`;
- IGRA specificity `98%`;
- treatment initiation `100%`;
- completion `90%`;
- adverse-event stop `0%`;
- full-course efficacy `90%`.

The ordinary working defaults remain population `1500`, recent TBI `2%`,
remote-only TBI `8%`, baseline active TB `0%`, screening coverage `30%`,
IGRA, 3HP, 20-year horizon, reset reinfection policy, no competing mortality
and the current cascade defaults. Under those defaults the deterministic
expected screened count is `450`, because screening coverage is programme
coverage, not annual coverage.

Authoritative coverage meaning:

```text
The proportion of the eligible population screened at least once over the
complete screening programme window.
```

There is no accidental repeat screening policy in the explicit workflow. A
one-year and two-year screening programme with `30%` coverage both screen
`30%` of the eligible population in expectation unless an explicitly modelled
competing event prevents screening before the sampled screening time.

Milestone 3E hardens accounting validation:

- stochastic summaries now report `notScreened`, infected and uninfected
  screened counts, latent false negatives, true negatives, TPT eligibility
  after testing/clinical assessment, and true-positive versus false-positive
  TPT starts;
- deterministic intervention expectations expose the same accounting
  proportions;
- deterministic and stochastic accounting validators check baseline-state
  totals, screened plus not-screened totals, test-result identities, TPT
  bounds, cases-averted identities, timing-component reconciliation and
  recent-plus-remote cases-averted reconciliation;
- the Streamlit workflow raises an error rather than displaying results if
  those identities fail.

The workflow configuration and provenance downloads now record:

- `screeningCoveragePolicyIdentifier =
  programme_coverage_at_least_once_v1`;
- the coverage definition above;
- pathway, runner contract, infection targets, fitted hazards, progression
  curve, reinfection policy, intervention policy, population, horizon,
  analysis method, repetitions and seed through the existing configuration
  JSON.

Performance review found that the pure deterministic helper remains
configurable, but the Streamlit "Quick deterministic preview" now uses
`16` infection-time and `12` screening-time Gauss-Legendre points. Against
the higher `32 x 24` diagnostic on the working defaults, the largest
difference across screened, true-positive latent, comparator TB, intervention
TB and direct cases averted was about `0.00025` people in a population of
`1500`, while runtime fell from about `65` seconds to about `17` seconds on
the review machine. Publication-quality deterministic checks can still call
the pure helper with tighter quadrature.

Milestone 3E does not add economics, DALYs, QALYs, ICERs, a costable event
ledger, MATLAB changes, dynamic-model changes, deployment, merge or tag.

## Milestone 3F isolated UI and presentation

Starting audit: repository root `tb-community-risk-tbi-calibration`, branch
`feature/explicit-recent-remote-tbi-v2`, HEAD
`4142c23a5f37e0bf9480c84f0858548002df161d`, upstream
`origin/feature/explicit-recent-remote-tbi-v2`, clean working tree and no
ahead/behind commits. Python was `3.12.4`; installed Streamlit was `1.32.0`,
while `requirements.txt` pins `1.39.0`. Streamlit 1.39 was installed only into
the ignored `abm/output/milestone3f/streamlit139` validation directory. The
user's installed environment and repository-wide Streamlit configuration were
not changed.

### Reproduced navigation leak and root cause

The original command was actually started:

```powershell
python -m streamlit run explicit_tbi_app.py
```

The live server's rendered messages combined the explicit opening content and
four-option workflow radio with a registry of 12 automatically discovered
pages. This included Scenario, Economics, Decision Analysis, Dynamic Model,
Compare, Dynamic ABM Compare and Technical Settings. Evidence was saved in
`abm/output/milestone3f/navigation-leak.json`.

The root entry point was adjacent to the legacy `pages/` directory. Streamlit
1.32's actual `source_util.get_pages()` scans
`Path(main_script_path).parent / "pages"`. A custom sidebar radio does not
replace that registry. The explicit entry did not call `st.navigation`, and
1.32 does not provide that API. The separate frozen entry point's call to
`st.navigation` is not executed when starting the explicit entry point.

This is ordinary legacy multipage discovery caused by the entry point's
location. It is not duplicate execution from a collision between two
`st.navigation` calls. Streamlit's [pages directory documentation](https://docs.streamlit.io/develop/concepts/multipage-apps/pages-directory)
describes automatic discovery beside the entry point. Its
[`st.navigation` documentation](https://docs.streamlit.io/develop/api-reference/navigation/st.navigation)
states that calling that API switches away from directory discovery. The
structural fix does not rely on that newer API or on hiding the sidebar.

### Final entry point and navigation

The former root file was moved to `explicit_tbi_ui/app.py`; the root file no
longer exists. The single supported explicit-app command is:

```powershell
python -m streamlit run explicit_tbi_ui\app.py
```

The isolated directory has no `pages/` subdirectory. Its server registry
contains only the isolated entry script, with exactly four internal workflow
views rendered in the sidebar:

1. Set up
2. Run analysis
3. Results
4. Evidence and technical information

The legacy `pages/` directory was neither renamed nor moved. A live direct
request for the Economics route produces `page_not_found`; no legacy page
script is registered by the isolated app. There is no separate Scenario view.

### Session boundary and deliberate non-costability

`app/explicit_tbi_session.py` owns versioned `explicit_tbi_ui_v1:` state and
widget keys. It rejects unmarked/incompatible namespaces, malformed result
provenance, non-finite input state, retired explicit-app keys and foreign
session state. A contaminated session opens a clean Set up workflow with no
old results or economics. Ordinary reruns and navigation retain valid explicit
inputs and results; material input edits mark a stored run stale.

`pages/4_Economics.py` adds an early guard only for an explicit session or an
explicit result put directly into the legacy results slot. It clears stale
legacy result/economics slots and stops before loading defaults, frozen
economics, backend calculations or economics displays. Clean frozen sessions
follow the existing code path. Tests exercise both namespace contamination and
direct legacy-slot contamination.

Every explicit workflow result carries the existing unavailable, invalid and
non-costable event-ledger marker. No costable ledger, economics, DALYs, QALYs
or ICER calculation was added.

### Presentation and precision

The title is now `TB screening and preventive treatment model`. A small caption
states development status and the planning/sequencing boundary. The first view
is Set up with `Input known TBI prevalence` and the requested recent and
remote-only labels. Both inputs visibly use `% of total population`.

`app/explicit_tbi_formatting.py` centralises display-only formatting:

- prevalence, coverage and cascade probabilities: one decimal percentage;
- baseline active TB: one decimal, allowing a meaningful second decimal;
- counts: thousands separators; whole people and individual simulated counts
  remain integers;
- expected values, stochastic means and simulation bounds: at most two decimal
  places, dropping unnecessary zero decimals;
- fractional years remain deliberately supported;
- fitted hazards, residuals, long hashes and diagnostic precision appear only
  under Evidence and technical information or technical downloads.

| Previous presentation | Ordinary presentation after 3F |
| --- | --- |
| Recent input `2.000000`, model caption `0.02000000` | `2.0%` with population unit |
| Remote-only input `8.000000`, model caption `0.08000000` | `8.0%` with population unit |
| Baseline active TB `0.000000` | `0.0%` |
| Expected screened `450.000000` | `450` |
| Comparator TB `2.218736` | `2.22` expected cases |
| Direct cases averted `0.277932` | Approximately `0.28` expected direct active TB cases averted |

Percentage callbacks write to model proportions only when the user edits a
widget. Display rounding never feeds back into configuration, health
signatures, hashes, calibration or calculations. Tests retain custom
full-precision proportions across reruns and navigation. All ordinary tables
are formatted copies; JSON/CSV downloads retain numeric full precision. The
full technical JSON includes configuration, results, calibration, marker and
timings. Validation CSV exports the stored run's validation payload.

Deterministic summaries use Expected value and N/A simulation intervals.
Stochastic summaries use Mean, low/high 95% simulation intervals and the
number of simulated populations from the stored run. Individual replicate
tables are explicitly identified. Stored population/horizon values describe
stale results rather than incorrectly using newly edited inputs.

The unused future currency formatter rounds to whole AUD with thousands
separators, consistent negatives and ties away from zero: `AUD 95,052`,
`-AUD 92,370`, `AUD 0`. The same presentation contract can later format an
authorised arithmetic ICER per DALY. No economic output is connected to it.

### Runtime investigation

The live audit `scripts/validate_explicit_tbi_ui.py` starts a fresh localhost
server, captures real rendered messages, dispatches two deterministic clicks
with unchanged defaults, captures stage timings and verifies the dispatch
counter advances exactly once for each click. Result-only navigation and
reruns do not dispatch analysis. AppTest additionally spies on the real
deterministic quadrature helper to verify a single call from a real click.

The first installed-version audit (Streamlit 1.32.0) measured:

| Measurement (seconds) | Cold click | Warm unchanged click |
| --- | ---: | ---: |
| Click through complete server rerun | 36.09 | 41.30 |
| Calibration | 0.03 | 0.06 |
| Deterministic intervention quadrature | 32.41 | 37.12 |
| Result construction, including baseline/age summaries | 3.52 | 3.96 |
| Dispatch wrapper | 35.96 | 41.16 |
| Rerun/serialization/transport above dispatch | 0.13 | 0.15 |
| Results-only rerun | 0.73 | 0.29 |

Cold initial entry import/render took 4.23 seconds. The quadrature settings
remain exactly 16 infection-time by 12 screening-time points; the engine and
existing accuracy contract are unchanged. A previously repeated baseline
natural-history calculation is now reused for age summaries, without changing
its values. Calibration is resolved once inside the analysis dispatch rather
than once in validation and again for the same analysis.

Quadrature dominates these measurements. Directory discovery registers pages;
it does not execute every page or repeat the deterministic calculation.
Duplicate analysis, navigation dispatches, large calibration costs and a
53-second cold import are excluded by the measured counters/timings. A warm
explicit rerun recomputes on a fresh click; it is not a result-cache hit. The
earlier approximately 17-20-second benchmark and the user's approximately
53-second observation have no matched CPU/load trace, so their exact difference
cannot honestly be assigned to hardware alone. The measured variation is in
the numerical computation, not a second analysis or leaked economics.

The final audit, after regression-test processes had finished, used the pinned
Streamlit 1.39.0 in a fresh server:

| Measurement (seconds) | Cold click | Warm unchanged click |
| --- | ---: | ---: |
| Click through complete server rerun | 22.17 | 26.50 |
| Process CPU time inside analysis | 20.81 | 24.86 |
| Calibration | 0.02 | 0.02 |
| Deterministic intervention quadrature | 20.30 | 24.95 |
| Result construction | 1.77 | 1.45 |
| Dispatch wrapper | 22.09 | 26.43 |
| Rerun/serialization/transport above dispatch | 0.08 | 0.07 |
| Results-only rerun | 0.42 | 0.18 |

Initial cold entry import/render took 2.51 seconds. The counter was 1 after the
cold click and 2 after the warm click. All four result tables were rendered;
the Economics route was rejected. The earlier 1.32 audit took place during
active verification work, so these are not controlled measurements of a
Streamlit-version effect. The difference is consistent with numerical workload
and execution-environment variation; it does not prove a specific hardware
cause. Full evidence is in `isolated-1.39.0-audit.json` in the milestone output
directory.

### Rendered validation and frozen boundary

AppTest covers four-only navigation, absence of Scenario/economics/dynamic/
comparison routes, the opening labels/units, baseline precision, all ordinary
result tables, stochastic summaries, invalid totals, single dispatch, reruns,
state contamination and direct legacy economics rejection. Machine-readable
exports and future currency formatting have separate precision tests.

No interactive Chrome or in-app browser surface was available in this session.
The real server audit therefore complements AppTest with rendered protobuf
elements and Arrow tables, not HTTP 200 or source inspection. Streamlit 1.39
AppTest returns an empty tree for the existing frozen `st.navigation` router;
the frozen-router test uses the live server instead of accepting that empty
tree as verification.

The frozen live server renders its original nine labels: Set up, Run Analysis,
Results, Health Economics, Explore Decisions, Evidence & Assumptions,
Technical settings, Dynamic transmission model and Model comparison. Its
Results view renders populated frozen tables. The frozen-reference loader
suite checks the existing 2,000 paired replicates and pinned headline outputs.

Final checks against the finished source:

| Check | Outcome | Test runtime |
| --- | --- | ---: |
| Installed Streamlit 1.32 isolation/formatting AppTest suite | 19 passed; 1 frozen-router version skip | 48.91 seconds |
| Pinned Streamlit 1.39 complete explicit workflow, isolation and frozen-reference suite | 41 passed | 466.19 seconds |
| Fresh 1.39 live isolated server and two real deterministic clicks | Passed; one dispatch per click | 22.17 / 26.50 second clicks |
| Live frozen 1.39 navigation and populated Results | Passed; original nine entries | Initial render 3.95 seconds; Results 9.35 seconds |
| Python compile and Git whitespace checks | Passed | Lightweight |

The complete supported-version command was
`python scripts/validate_explicit_tbi_ui.py --streamlit-target abm/output/milestone3f/streamlit139 --apptest`.
The installed-version command was
`python -m pytest tests/test_explicit_tbi_ui_isolation.py -q`.
Warnings were existing dependency deprecations. The version skip is exercised
successfully by the supported-version live frozen-router test. Population input
also explicitly renders `1,500`, accepts comma-grouped whole counts and rejects
fractional people; invalid drafts remain blocked across navigation.

`streamlit_app.py`, repository Streamlit configuration, MATLAB, the dynamic
model, all engine formulas, scientific defaults, progression curve, TPT policy,
frozen report artifacts, reference data and MATLAB reference artifacts remain
unchanged from the milestone starting commit. Local and remote release branch
and peeled tag remain at `03cc16e4a52a55e10019dab16bd4f599571c4b90`;
the annotated tag object remains `e99d36725dfc3ff40edc0ed68127e0b6cd927da5`.
No deployment, merge or tag operation was performed.

Changed files: the relocated `explicit_tbi_ui/app.py`, its README,
`app/explicit_tbi_workflow.py`, `app/explicit_tbi_session.py`,
`app/explicit_tbi_formatting.py`, the explicit-only guard in
`pages/4_Economics.py`, `scripts/validate_explicit_tbi_ui.py`,
`tests/test_explicit_tbi_streamlit_app.py`,
`tests/test_explicit_tbi_ui_isolation.py`, root README and the three existing
explicit-pathway specification/status/evidence documents.

## Next milestone recommendation

The next milestone should review the explicit development workflow and decide
whether to keep it as a separate research app, broaden rendered validation,
or start a carefully guarded report-facing integration plan. Event-ledger
costing, health economics and DALYs should remain later milestones after the
epidemiological workflow is reviewed.
