# Explicit recent/remote TBI calibration - v2 foundation specification

This branch starts a new model-development line derived from the frozen SA
Health Streamlit release. It does not modify the submitted SA Health release
and does not expose a new pathway in the standard Streamlit workflow.

The Milestone 1 scope is limited to audit, mathematical specification, pure
calibration functions and tests. The epidemiological runner, event ledger,
health economics, DALYs, MATLAB, dynamic transmission model, frozen artifacts,
deployment, release branch and release tag are intentionally unchanged.

Milestone 2A extends the isolated foundation with a configuration contract,
deterministic assignment summaries, stochastic state and infection-time
assignment helpers, and a documented baseline active-TB sequencing proposal.
It still does not connect the new pathway to the runner, Streamlit UI,
active-TB progression, active-TB calibration, interventions, event ledger,
health economics, DALYs, MATLAB or the dynamic-transmission model.

Milestone 2B adds pure natural-history progression mathematics for the new
baseline states and expected-value diagnostics for prospective active-TB
targets. It still does not connect the pathway to the runner, Streamlit UI,
event ledger, intervention logic, economics, DALYs, MATLAB,
frozen-reference loading or the dynamic-transmission model.

Milestone 2C adds explicit progression-calibration policy contracts,
eligibility checks, one-parameter fixed-ratio calibration, likelihood helpers,
identifiability enforcement and risk-factor safety diagnostics. It remains
pure code only and still does not connect the pathway to the runner, Streamlit
UI, population generation, interventions, event ledger, economics, DALYs,
MATLAB, frozen-reference loading or the dynamic-transmission model.

Milestone 2D corrects calibration-audit issues before any runner integration:
denominator concepts are separated, competing mortality is treated as a
cause-specific competing-risk integral rather than horizon survival
multiplication, ascertainment is required as an externally fixed assumption
for Policy B, practical identifiability checks are strengthened, the initial
production risk-factor policy is `none`, and candidate natural-history
parameterisations are documented from a focused evidence review. It remains
pure code only and still does not connect the pathway to the runner,
Streamlit UI, population generation, interventions, event ledger, economics,
DALYs, MATLAB, frozen-reference loading or the dynamic-transmission model.

## Selected identifiers

- Analysis basis: `explicit_recent_remote_tbi_foundation_v1`
- Natural-history semantics: `explicit_recent_remote_tbi_history_v1`
- Calibration method: `constant_hazard_age_window_recent_remote_v1`
- Infection-history contract: `recent_remote_tbi_history_contract_v1`
- Calibration contract: `explicit_recent_remote_tbi_calibration_v1`
- Active-TB observation schema: `active_tb_observation_targets_v1`
- Configuration contract: `explicit_recent_remote_tbi_config_v1`
- Assignment contract: `explicit_recent_remote_tbi_assignment_v1`
- Progression contract: `explicit_recent_remote_tbi_progression_v1`
- Progression-calibration policy contract:
  `explicit_recent_remote_tbi_progression_calibration_policy_v1`
- Recent hazard shape: `constant_recent_window_hazard_v1`
- Remote hazard shape: `constant_remote_window_hazard_v1`

These are new identifiers. They do not reuse or reactivate the retired
`continuous_markov_recent_remote` pathway and do not replace the frozen
`sa_health_matlab_v9_compatibility_reference` / `matlab_v9_implicit_early_late`
metadata.

## Definitions

At baseline, each person has one effective infection state:

1. `uninfected`: no infection in either exposure window.
2. `remote_only`: at least one infection acquired five or more years before
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

## Age representation in this repository

The inherited source tables can contain broad age bands. The Python APY data
loader does not use band midpoints for cohort generation or expected-value
calculation. Instead:

- `engine/apy/age_distribution.py::expand_age_distribution_table` parses each
  age band and spreads the band weight uniformly across exact integer ages in
  the band.
- The open-ended `85+` band is expanded through the configured
  `age85PlusMax`, which defaults to age 89.
- `engine/apy/cohort.py::draw_base_population` samples from
  `pars["exactAgeValues"]` using `pars["exactAgeProb"]`.
- `engine/apy/expected_value.py::_build_strata` iterates the same exact-age
  support and probabilities.
- The explicit recent/remote config now records `age85PlusMax` and
  age-support provenance, and those fields participate in deterministic
  configuration hashing.

Milestone 2A assignment functions therefore accept exact ages and probabilities
and use the same representation for calibration and assignment. They do not
silently substitute age-band midpoints. If a future caller supplies only broad
bands, that caller must first apply an explicit within-band expansion rule
with labelled provenance.

The default `age85PlusMax = 89` is an inherited implementation choice, not an
evidence-based biological maximum. A synthetic sensitivity using all population
weight in the open-ended age group showed that changing the support from
85-89 to 85-95 changed the fitted remote hazard for a 40% remote-only target
from `0.00666571` to `0.00643236` per year while preserving the target. This
demonstrates that the cap can matter and should remain explicit.

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

Milestone 2A made the recent-window duration and remote-history cap explicit
arguments to feasibility assessment, calibration and population-prevalence
helpers. The defaults remain five years and 100 years.

## Configuration contract

The isolated configuration contract is
`explicit_recent_remote_tbi_config_v1`. It contains:

- enabled/disabled status;
- recent TBI target as a total-population proportion;
- remote-only TBI target as a total-population proportion;
- recent-window duration, default five years;
- remote-history cap, default 100 years;
- recent and remote hazard-shape identifiers;
- target source;
- target reference year;
- review status;
- notes;
- calibration, analysis-basis, natural-history, method and infection-history
  contract identifiers.

Configuration JSON is serialized canonically with sorted keys and compact
separators, then hashed with SHA-256 for deterministic provenance. Enabled
configuration requires a target source and review status. The contract rejects
retired or incompatible identifiers rather than silently interpreting them as
the new pathway.

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

Milestone 1 specified timing; Milestone 2A implements isolated timing helpers
but still does not wire them into the runner.

For a constant Poisson hazard over an eligible window of length `L`, the time
since the most recent event conditional on at least one event has density:

```text
f(t | at least one event) = lambda * exp(-lambda * t) / (1 - exp(-lambda * L))
0 <= t <= L
```

The inverse-CDF used for sampling is:

```text
t = -log1p(-q * (1 - exp(-lambda * L))) / lambda
```

The implementation uses `expm1` and `log1p` for numerical stability. Small
hazards approach uniform timing over the eligible window. Large hazards
concentrate the most recent event close to baseline.

Conditional on the effective state:

- `recent`: draw the most recent infection event within
  `[0, min(recentWindowYears, age)]` years before baseline from the
  constant-hazard distribution truncated to the available recent window.
- `remote_only`: draw the most recent remote-window infection event within
  `[recentWindowYears, min(remoteHistoryCapYears, age)]` years before baseline
  from the constant-hazard distribution truncated to the remote window,
  conditional on no recent event.
- both remote and recent: classify effective state as `recent`, retain
  `priorRemoteExposure=true`, and treat recent reinfection as resetting the
  higher-progression-risk state in future integration unless a later scientific
  decision chooses a different mechanism.

The effective time since infection is the most recent event governing the
effective state. For recent-with-prior-remote exposure, the older remote event
is retained only as an auditable history flag and is not used as the effective
progression clock.

The isolated assignment output also calculates remaining early-risk duration
using the same authoritative recent window that defines recent versus remote
infection:

```text
remainingEarlyRiskYears =
    max(0, recentWindowYears - timeSinceMostRecentInfection)
```

Recent states therefore receive between zero and five remaining early-risk
years under the default five-year window, or the corresponding amount under a
reviewed non-default window. Remote-only states receive zero. Uninfected states
have no infection time and zero remaining early-risk years in the current
serializable output. This avoids giving every recently infected person a fresh
five-year high-risk period at baseline, and it prevents the recent/remote
classification boundary from diverging from the early-risk duration.

Audit finding: the inherited runner currently uses a baseline recent flag plus
a Markov recent-to-remote progression compartment. It does not distinguish
first infection, most recent infection and reinfection reset mechanisms. That
mechanism remains unchanged in this milestone.

## Deterministic and stochastic assignment

Milestone 2A adds pure deterministic helpers that return expected counts and
proportions for:

- `recent`;
- `remote_only`;
- `uninfected`;
- prior remote exposure;
- prior remote plus recent reinfection;
- age-specific state distributions;
- age-specific recent prevalence in the total age group;
- age-specific recent fraction among people with TBI;
- expected remaining early-risk duration among recently infected people.

The deterministic path uses no random draws. When supplied calibrated hazards,
its population totals reconcile to the requested recent and remote-only
targets within numerical tolerance.

The stochastic helpers require an explicit seed or random generator. They use
local `numpy.random.Generator` state and do not touch NumPy global RNG state.
They draw mutually exclusive effective states, draw infection times only for
people whose effective state requires an infection clock, retain prior remote
exposure for recent reinfection, and avoid draws for structurally impossible
events.

## Pure prospective progression

Milestone 2B adds `engine/apy/explicit_recent_remote_progression.py`, a pure
module with no Streamlit, runner, event-ledger, economics, DALY, MATLAB or
dynamic-model dependency.

For a recently infected person with time since most recent infection `s`,
recent-window duration `W`, remaining early-risk duration
`r = max(0, W - s)`, early progression hazard `lambda_E`, remote progression
hazard `lambda_L` and multiplier `m_i`, the prospective cumulative hazard is:

```text
A_i(t) = m_i [lambda_E min(t, r) + lambda_L max(0, t-r)]
P_i(T <= t) = 1 - exp[-A_i(t)]
```

For remote-only infection:

```text
A_i(t) = m_i lambda_L t
```

For uninfected people under the current no-new-infection pathway:

```text
P_i(T <= t) = 0
```

The pure module supports cumulative hazard, piecewise instantaneous hazard,
survival probability, cumulative incidence by horizon, expected events across
weighted strata and exact quantile inversion for future stochastic sampling.
It does not generate stochastic active-TB times.

The cumulative hazard is continuous at the end of remaining early-risk time,
although the instantaneous hazard may change from `lambda_E` to `lambda_L`.

## Progression-calibration policies

Milestone 2C defines four explicit progression-calibration policies:

- Policy A, `external_progression_hazards_v1`: externally supplied early and
  remote hazards, units, source, reference population, review status and
  notes. No fitting occurs.
- Policy B, `fixed_early_remote_ratio_fit_scale_v1`: externally supplied
  ratio `R` with provenance, then `lambda_L = k` and `lambda_E = Rk`; only
  the common scale `k` is fitted.
- Policy C, `validation_only_v1`: observations are compared with supplied
  hazards but no parameter is fitted. Retrospective notifications, prevalence
  observations, screen-detected disease and mixed/insufficient observations
  route here.
- Policy D, `joint_early_remote_hazards_v1`: specified but unavailable unless
  explicit identifiability criteria pass. It is not a fallback optimiser.

These policies deliberately do not reuse the inherited `earlyLateRatio`,
MATLAB-v9 compatibility assumptions or the `10/770` active-TB target.

## Progression target eligibility

A target is eligible for prospective progression calibration only when it has:

- a defined baseline population and compatible non-baseline-active-TB strata;
- a positive denominator;
- a prospective incident observation window beginning at model baseline;
- no baseline prevalent cases mixed into the numerator;
- explicit ascertainment probability `q`;
- explicit observation duration.

Baseline prevalence, cross-sectional screening prevalence, retrospective
notifications, mixed prevalent/incident counts, undefined denominators,
unclear ascertainment and missing periods are ineligible for fitting and are
routed to validation-only. Zero-case prospective incident targets are eligible
and fit to zero progression hazards under Policy B.

## Expected incident cases and fitting

For eligible prospective targets:

```text
E[C] = sum_i w_i [1 - exp{-A_i(T)}] q_i
```

where `w_i` is a count or population weight, `T` is the observation horizon,
and `q_i` is explicit ascertainment. Uninfected people contribute zero under
the current no-new-infection pathway. Baseline/prevalent active-TB strata are
excluded from prospective TBI progression denominators and reported as
excluded weight.

Policy B uses deterministic bounded bisection for `k`. It returns requested
cases, denominator, horizon, ascertainment, supplied ratio, fitted remote
hazard, derived early hazard, achieved cases, residual, convergence status,
feasibility status, warnings, provenance and contract version. Targets above
the achievable range are rejected rather than clipped.

## Likelihood options

Milestone 2C implements likelihood helpers for future uncertainty analysis:

- `binomial_person_event_v1` when the denominator is persons with at most one
  relevant event;
- `poisson_count_v1` for count/person-time rare-event settings where
  appropriate.

The observation-model identifier must be supplied explicitly. The code does
not choose between binomial and Poisson from numerical values alone.
Recurrent disease, migration, changing denominators or changing ascertainment
can violate both simple observation models.

## Identifiability enforcement

The pure identifiability checks enforce:

- one aggregate target cannot identify both `lambda_E` and `lambda_L`;
- fixed ratio plus scale is structurally one-dimensional;
- externally supplied hazards require no fitting;
- two or more targets do not automatically imply identifiability;
- joint fitting requires target sensitivity vectors with rank two;
- duplicate recent/remote composition has rank one and is insufficient.

Policy D is returned as unavailable unless the rank check passes. No fragile
two-parameter optimiser is implemented.

## Risk-factor safety diagnostics

Milestone 2C adds explicit risk-factor application policies:

1. `none`: all progression multipliers are one.
2. `reviewed_hazard_multipliers`: only effects explicitly reviewed as hazard
   multipliers are included.
3. `legacy_or_as_hazard_diagnostic_only`: reproduces inherited OR
   multiplication for diagnostic comparison, labelled scientifically
   provisional and not a production default.

Diagnostics report individual factor effects, effect-measure labels, combined
multiplier, log contributions, maximum possible multiplier, weighted
multiplier distribution, review-threshold warnings/blocking status,
progression probabilities and strata dominating expected cases. Multipliers
are not silently capped. Review thresholds are safety rules, not biological
evidence.

The inherited default OR-labelled effects can still combine to `2916` under
legacy diagnostic multiplication.

## Corrected competing-mortality interface

Milestone 2D removes scalar horizon survival from production-capable
expected-case calculations. A scalar `S_D(T)` cannot reconstruct the
competing-risk cumulative incidence curve and must not be used as:

```text
P(TB by T) * S_D(T)
```

For a cause-specific TB progression hazard `h_TB(t)` and external competing
death survival `S_D(t)`, prospective TB cumulative incidence is:

```text
F_TB(T) = integral_0^T S_TB(t) S_D(t) h_TB(t) dt
```

where:

```text
S_TB(t) = exp[-integral_0^t h_TB(u) du]
```

The pure functions now support:

- no competing mortality, recorded as `competingMortality = not_modelled`;
- an externally supplied constant death hazard, integrated analytically across
  the early/late TB-hazard boundary;
- an externally supplied survival curve, validated to start at one, remain
  finite, lie in `[0,1]`, be non-increasing and cover the horizon, then
  integrated numerically by trapezoid quadrature.

Scalar horizon survival can only be used behind an explicit
`allow_scalar_survival_approximation=True` flag and is labelled
`scalar_horizon_survival_nonproduction_approximation`. It is not used by
Milestone 2D tests, examples or fitting policies.

For a constant TB hazard `h_TB` and constant death hazard `h_D`, the analytic
test identity is:

```text
F_TB(T) = h_TB / (h_TB + h_D) * [1 - exp{-(h_TB + h_D)T}]
```

For piecewise early/late TB hazards the same expression is evaluated by
segment, carrying forward TB survival and death survival at each segment
start. Including competing mortality cannot increase expected TB cases.

## Inherited progression audit

Current inherited progression behavior remains unchanged:

- `engine/apy/timing.py` contains scientifically reusable pure piecewise
  early/late survival helpers for a fixed early-period duration.
- `engine/apy/ltbi_state.py` defines the retired
  `continuous_markov_recent_remote` latent-state model, including a default
  recent-to-remote transition rate of `1/5` per year.
- `engine/apy/calibration.py` calibrates infection prevalence and active-TB
  progression to inherited targets. The default active-TB target is `10/770`
  at the active-TB calibration horizon. `earlyLateRatio` must be at least one;
  `lambdaLate = lambdaEarly / earlyLateRatio`.
- `engine/apy/calibration.py::MATLAB_V9_IMPLICIT_EARLY_LATE` preserves MATLAB
  v9 compatibility by treating all infected people as beginning in the early
  state for the calibration window. This is compatibility behavior, not the
  new explicit recent/remote pathway.
- `engine/apy/simulation.py::_draw_ltbi_state_history` stochastically draws
  inherited recent/remote latent state, recent-to-remote transition time and
  untreated active-TB time. In the MATLAB-v9 branch, people not active during
  the screening window transition to remote at the screening-window boundary.
- `engine/apy/expected_value.py` uses inherited recent/remote probabilities
  and `mixed_baseline_survival`/`mixed_baseline_event_between` for deterministic
  event calculations.

Reusable for the new pathway: pure hazard/survival algebra where the state and
remaining early-risk time are explicit. Unsafe to reuse silently:
`continuous_markov_recent_remote`, inherited `baselineRecentLTBIProportion`,
MATLAB-v9 compatibility semantics and automatic reuse of the `10/770`
active-TB target.

## Progression multiplier crosswalk

The inherited engine currently applies disease-risk effects by multiplying
progression hazards. The parameter names and default normalized values are:

| parameter | label | recorded measure | current engine use |
| --- | --- | --- | --- |
| `disOR["MJ"]` | Other current risk-factor prevalence | OR disease | hazard multiplier |
| `disOR["contact"]` | Contact history | OR disease | hazard multiplier |
| `disOR["renal"]` | Renal impairment | OR disease | hazard multiplier |
| `disOR["diabetes"]` | Diabetes | OR disease | hazard multiplier |
| `disOR["smoking"]` | Smoking | OR disease | hazard multiplier |
| `disOR["cld"]` | Chronic lung disease | OR disease | hazard multiplier |
| `disOR["alcohol"]` | Alcohol/drug exposure | OR disease | hazard multiplier |

The defaults loaded through the normal Python APY configuration are:
`MJ=3.0`, `contact=5.0`, `renal=3.6`, `diabetes=3.0`, `smoking=2.0`,
`cld=3.0`, `alcohol=3.0`. They are multiplied jointly, giving a possible
combined multiplier of `2916` if all flags are present. No cap was found. The
input fields and source CSV label these as odds ratios; their use as hazard
multipliers is an inherited modelling choice and is not scientifically
endorsed by Milestone 2B.

The new pure progression functions accept a generic non-negative multiplier
`m_i`; they do not infer, cap or reinterpret risk-factor ORs.

## Risk-factor separation

Milestones 1, 2A, 2B and 2C do not use non-age risk factors in the new
acquisition calibration. Age affects exposure only through time alive.

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
not infer one from another and does not fit progression hazards in Milestones
1, 2A or 2B.

Milestone 2B classifies each observation row for future use as one of:

- `baseline_prevalence_target`;
- `screen_detected_disease_target`;
- `prospective_incident_disease_target`;
- `retrospective_notification_incidence_target`;
- `mixed_or_insufficiently_defined`.

Only a genuinely prospective incident target can be passed to the pure
expected-count diagnostic function. That requires the baseline population to
be defined, the observation window to begin at model baseline, the
recent/remote/uninfected composition to be supplied externally, ascertainment
to be complete or explicitly parameterised, and competing mortality to be
explicitly addressed. The expected count calculation is:

```text
E[C] = sum_i w_i P_i(T <= T_obs) q_i
```

where `q_i` is an ascertainment probability if used. The function retains the
observed numerator, denominator, person-years if supplied and observation
period; it does not reduce the row to an annual rate and discard those
components.

## Baseline active-TB sequencing proposal

Future integration must assign baseline active TB separately from TBI state.
A person with prevalent active TB at baseline should not simultaneously enter
the TBI preventive-treatment cascade. The proposed sequencing is:

1. Assign or import baseline active-TB status first, using a prevalence or
   observation model that is explicit about ascertainment.
2. Exclude baseline active-TB cases from the baseline TBI preventive-treatment
   cascade while retaining them for active-TB burden reporting.
3. Among people without baseline active TB, assign `recent`, `remote_only` or
   `uninfected` using the explicit recent/remote TBI pathway.
4. Treat active TB detected during screening as a screening outcome, distinct
   from prevalent active TB assigned before screening and distinct from later
   incident active TB.
5. Treat incident active TB after baseline as a future progression outcome
   conditional on infection state, progression hazards and interventions.

The optional active-TB observation rows remain data for future calibration and
review. Milestone 2C uses eligible rows only in pure diagnostic expected-count
and fixed-ratio calculations; it does not fit or alter runner progression
hazards.

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

One aggregate active-TB count generally cannot identify both `lambda_E` and
`lambda_L`. Milestone 2B therefore exposes identifiability diagnostics but
does not choose a hidden progression-calibration policy.

Potential future progression-calibration policies are:

1. Both hazards externally supplied.
2. One hazard externally supplied and the other fitted.
3. Early-to-late hazard ratio fixed externally and a common scale fitted.
4. Multiple sufficiently informative targets used in a joint likelihood.
5. Active-TB data used only for external validation.

Future progression calibration would require combinations such as:

- explicit recent and remote-only TBI prevalence targets;
- active-TB observations classified by timing and ascertainment;
- fixed or estimated progression-hazard structure by infection state;
- fixed test/screening detection assumptions for screen-detected active TB;
- explicit handling of prevalent baseline active TB versus incident follow-up;
- assumptions about risk-factor effects on progression, separate from
  acquisition.

The inherited `10/770` active-TB target is not automatically reused for the
new explicit recent/remote pathway.

## Historical observations

A retrospective notification count from before model baseline is not
automatically equivalent to prospective active TB arising from the present
baseline cohort. Reasons include population turnover, migration, depletion
through prior disease, prior preventive or active-TB treatment, mortality,
changing ascertainment, changing infection pressure and changing denominators.
Such rows are classified as validation data unless a separate retrospective
population reconstruction is implemented.

## Competing mortality

The inherited no-transmission APY natural-history progression path does not
apply competing all-cause mortality before active-TB progression. Mortality and
TB fatality appear later in burden/economic calculations, not as a competing
risk in the untreated progression-time generation or deterministic
expected-value progression functions.

Over a long horizon such as 20 years, omitting competing mortality will tend
to overstate prospective progression events among older age groups, because it
leaves people at risk for active TB after they might otherwise have died. The
importance depends on age structure, comorbidity, follow-up duration and local
mortality. Milestone 2B designs the pure expected-count interface so a future
survival function can be incorporated, but it does not invent mortality data.

## Worked progression diagnostics

The following synthetic diagnostics use illustrative hazards
`lambda_E = 0.04` per year and `lambda_L = 0.004` per year. They are not
reviewed defaults.

| state/scenario | multiplier | P(1y) | P(2y) | P(5y) | P(10y) | P(20y) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| recent, infected 0.5y ago | 1 | 0.039211 | 0.076884 | 0.166399 | 0.182905 | 0.214944 |
| recent, infected 0.5y ago | 4 | 0.147856 | 0.273851 | 0.517126 | 0.554251 | 0.620158 |
| recent, infected 2.5y ago | 1 | 0.039211 | 0.076884 | 0.104166 | 0.121905 | 0.156335 |
| recent, infected 2.5y ago | 4 | 0.147856 | 0.273851 | 0.355964 | 0.405479 | 0.493383 |
| recent, infected 4.9y ago | 1 | 0.007571 | 0.011533 | 0.023324 | 0.042663 | 0.080201 |
| recent, infected 4.9y ago | 4 | 0.029943 | 0.045340 | 0.090081 | 0.160039 | 0.284233 |
| remote-only | 1 | 0.003992 | 0.007968 | 0.019801 | 0.039211 | 0.076884 |
| remote-only | 4 | 0.015873 | 0.031493 | 0.076884 | 0.147856 | 0.273851 |
| uninfected | 1 | 0 | 0 | 0 | 0 | 0 |
| uninfected | 4 | 0 | 0 | 0 | 0 | 0 |

The same pure module includes a synthetic cost-invariance diagnostic:
changing a cost or intervention parameter does not change natural-history
progression probabilities because the progression functions do not read
economic or intervention inputs.

## Worked progression-calibration diagnostics

Milestone 2C synthetic examples use a 2026 one-year prospective incident
target with two observed cases in a denominator of 1,000 unless otherwise
specified. These are not reviewed defaults.

| example | result |
| --- | --- |
| external hazards `lambda_E=0.02`, `lambda_L=0.002` | expected cases `2.379733`; no fitting |
| fixed ratio `R=10` | fitted `k=lambda_L=0.00167858`, `lambda_E=0.01678576`, achieved cases `2.000000` |
| no multipliers | fitted `k=0.00167858`, maximum individual probability `0.016646` |
| moderate reviewed multiplier example | fitted `k=0.00091672`, maximum individual probability `0.018167` |
| legacy OR-as-hazard diagnostic | fitted `k=0.00033747`, top 1% by weight contributes `0.832978` of expected cases |
| retrospective notification | routed to validation-only as `retrospective_notification_incidence_target` |
| baseline prevalence | routed to validation-only as `baseline_prevalence_target` |
| impossible high target | rejected as `infeasible_above_achievable_range`; maximum achievable cases `300` |
| zero-case target | fitted `lambda_E=lambda_L=0` |
| incomplete ascertainment `q=0.5` | fitted `k=0.00338139`; ascertainment retained explicitly |
| external constant death hazard `0.01` | fitted `k=lambda_L=0.00168702`, `lambda_E=0.01687023`; mortality mode recorded as external constant mortality hazard |

The multiplier examples show that calibration can shrink the fitted baseline
hazard when high multipliers are present. The concentration diagnostics are
therefore part of the calibration output and should not be suppressed.

## Milestone 2D calibration audit corrections

Denominator semantics are now explicit:

- `sourcePopulationDenominator`: denominator recorded in the observed data.
  It is preserved and is not silently rewritten.
- `baselineActiveTBCount`: people already in the baseline/prevalent
  active-TB state. They are excluded from TBI states and prospective
  progression calculations.
- `prospectiveAtRiskPopulation`: non-baseline-active-TB population eligible
  to contribute a first prospective incident active-TB event.
- `personTimeAtRisk`: supplied person-years when present, otherwise a
  labelled calculation from prospective at-risk population and horizon.
- `tbiEligiblePopulation`: non-baseline-active-TB population potentially
  entering future TBI screening/TPT sequencing. This is separate from the
  observed source denominator.

If baseline active TB is present and the observed numerator does not state
whether baseline/prevalent cases are included, fitting is blocked. If the
numerator explicitly includes baseline/prevalent active TB, prospective
incident calibration is also blocked. Example: a source row with 10 cases
among 1,500 people and two baseline active-TB cases still reports
`sourcePopulationDenominator = 1500`; the at-risk population is reported
separately as 1,498, and fitting proceeds only if the numerator explicitly
excludes baseline/prevalent disease.

Ascertainment is now a fitted-policy requirement, not a hidden default. Policy
B requires `q` to be fixed externally, in `(0,1]`, sourced and assigned review
status. If `q=1`, the source must state that complete ascertainment is an
explicit assumption. With `E[C]=q E[C_true(k)]`, one aggregate count cannot
identify both `q` and the progression scale `k`; lower `q` can be offset by a
higher fitted scale.

Policy D now requires structural and practical identifiability. Diagnostics
include sensitivity/Jacobian rank, singular values, condition number,
recent/remote composition contrast, and flags stating that profile-likelihood
width and parameter-boundary effects remain required before production joint
estimation. The condition-number and composition-contrast thresholds are
numerical review rules, not biological evidence. Near-collinear targets do
not pass simply because floating-point rank is two.

The initial recommended production risk-factor policy is:

```text
none
```

`legacy_or_as_hazard_diagnostic_only` remains available for audit comparisons
but is blocked from being selected as a production/default runner policy. The
diagnostic still shows inherited factors can combine to `2916`, and the
synthetic high-multiplier example still has the highest 1% by population
weight contributing about `0.832978` of expected cases. Calibrating a lower
baseline hazard does not fix misspecified relative-risk concentration.

## Focused evidence review

The following evidence table is for scientific review. It does not create
production defaults and does not convert odds ratios, relative risks or
cumulative risks into hazards without a future documented transformation.

| topic | candidate evidence | population/setting/design | measure and outcome | applicability and limitations | status |
| --- | --- | --- | --- | --- | --- |
| Progression by time since infection | Behr et al. 2018, *BMJ*; Menzies et al. / Campbell et al. time-since-infection syntheses; classic tuberculin-converter cohorts including Sutherland, Ferebee and Comstock | Cohort/synthesis of tuberculosis after infection or conversion, mostly non-contemporary settings | Absolute risk and hazard/risk concentrated soon after infection | Supports high early risk after recent infection, but estimates vary by age, setting, comorbidity and treatment era | reviewed candidate |
| Higher-risk period duration | WHO TB preventive treatment guidance and cohort syntheses commonly emphasise recent infection/contact, especially within two years, while this model uses a five-year recent window | Programmatic guidance plus cohort evidence | Eligibility/risk period rather than direct constant-hazard estimate | Five years is a modelling window; two-year and continuous-decline alternatives should be sensitivity analyses | sensitivity only |
| Early versus remote hazards/ratio | Natural-history model syntheses including Vynnycky/Fine-style lifetime risk work and recent time-since-infection reviews | Model synthesis from historical cohorts | Early:remote relationship inferred from time-varying risk | A fixed two-phase ratio is an approximation; evidence may favor a declining continuous hazard | unresolved |
| Reinfection | Molecular epidemiology and high-incidence cohort studies such as Verver et al.; systematic reviews of recurrent TB distinguish relapse/reinfection | High-burden settings, often recurrence after disease rather than asymptomatic baseline TBI | Reinfection contribution and recurrence risk | Does not directly prove that recent reinfection resets latent progression risk in this model; reset assumption remains review-required | unresolved |
| Age effects | WHO guidance for child contacts; individual-participant meta-analyses of child household contacts; natural-history syntheses | Children and household contacts, mixed settings | Disease risk after infection/contact by age | Stronger support for young-child risk; adult age effects need careful mortality/comorbidity separation | reviewed candidate for age sensitivity |
| Diabetes | Jeon and Murray 2008 systematic review/meta-analysis; later updates | Observational studies, global | Relative risk/OR for active TB among people with diabetes | Suitable for risk review, but not automatically a hazard multiplier in this model | reviewed candidate, transformation unresolved |
| Renal disease | CKD/dialysis TB risk systematic reviews and national guidance identifying dialysis/renal failure as high risk | CKD and dialysis populations | Relative risks/incidence ratios for active TB | High-risk group but effect size depends on dialysis, transplant, setting and screening | reviewed candidate/sensitivity |
| Smoking | Bates et al. 2007 and later meta-analyses | Observational studies | Relative risk/OR for infection, disease and mortality outcomes | Smoking may affect infection and progression; do not conflate acquisition and progression | sensitivity only until pathway-specific effect reviewed |
| Harmful alcohol/drug exposure | Lonnroth et al. 2008 systematic review/meta-analysis and WHO risk-factor discussions | Observational studies | RR/OR for active TB and heavy alcohol use | Confounding and exposure definition vary; drug exposure not identical to alcohol | sensitivity only |
| Close contact | Household-contact systematic reviews/IPD meta-analyses and WHO preventive-treatment guidance | Household/close contacts | Incident TB and infection risk after exposure | Contact is an acquisition/exposure marker and may also proxy recent infection; should not be reused as progression HR without review | unsuitable as generic progression multiplier |
| Chronic lung disease | COPD/chronic airway disease meta-analyses | Observational cohorts/case-control | RR/OR/HR for active TB among COPD/chronic lung disease | Potential progression and detection bias; local disease definitions needed | sensitivity only |
| Competing all-cause mortality | Australian Bureau of Statistics life tables, Australian Government Actuary tables, WHO life tables | National life tables | Age/sex all-cause death rates or survival | Use as external mortality input only; do not invent mortality in this milestone | reviewed candidate source |

Candidate sources reviewed or queued for production review include:

- WHO consolidated TB preventive-treatment guidance, 2020 and second edition
  2024:
  https://www.who.int/publications/i/item/9789240096196
- Behr, Edelstein and Ramakrishnan 2018, *BMJ*, "Revisiting the timetable of
  tuberculosis": https://pubmed.ncbi.nlm.nih.gov/30139910/
- Menzies et al. 2018, progression assumptions in TB transmission models:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6070419/
- Time-since-infection natural-history synthesis for the United States:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC7707158/
- Jeon and Murray 2008, diabetes and active TB systematic review:
  https://doi.org/10.1371/journal.pmed.0050152
- Bates et al. 2007, tobacco smoke and TB systematic review/meta-analysis:
  https://pubmed.ncbi.nlm.nih.gov/17325294/
- Lonnroth et al. 2008, alcohol and TB systematic review:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC2533327/
- Martinez et al. 2020, child close-exposure individual-participant
  meta-analysis: https://pmc.ncbi.nlm.nih.gov/articles/PMC7289654/
- CKD without kidney failure TB-risk systematic review/meta-analysis:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10573716/
- CKD/dialysis TB-incidence systematic review/meta-analysis:
  https://pubmed.ncbi.nlm.nih.gov/35609860/
- Chronic airway disease and active TB systematic review:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC9070518/
- ABS Life expectancy/life tables 2021-2023:
  https://www.abs.gov.au/statistics/people/population/life-expectancy/2021-2023
- WHO Global Health Observatory life tables:
  https://www.who.int/data/gho/data/themes/mortality-and-global-health-estimates/ghe-life-expectancy-and-healthy-life-expectancy

Live source verification was repeated on 2026-10-08 after network access was
restored. Additional extracted review notes:

- WHO 2024 TPT guidance states that roughly one fourth of the world population
  is estimated to have been infected with TB bacilli, and about 5-10% of
  infected people develop TB disease in their lifetime; it frames TPT around
  groups at highest risk of progression rather than a single universal hazard.
- The United States time-since-infection synthesis estimated that, for a newly
  infected adult without other progression risk factors, progression rates
  decline from about 38 to 0.38 per 1,000 person-years between the first and
  25th year since infection, with 25-year cumulative risk about 7.9%. This
  supports declining risk since infection and cautions against treating a
  five-year window as biologically flat.
- The late-reactivation systematic review found declining TB rates over time,
  reaching approximately 200 cases per 100,000 person-years or less by the
  fifth year in eligible untreated cohorts, with limited evidence beyond ten
  years.
- Jeon and Murray reported cohort-study RR 3.11 (95% CI 2.27-4.26) for
  diabetes and active TB; this is a relative-risk estimate for disease, not an
  acquisition parameter and not automatically a hazard multiplier.
- Bates et al. reported smoking associations with TB infection, pulmonary
  disease and mortality; summarized TB disease RRs were about 2.3-2.7 in the
  JAMA Internal Medicine summary, reinforcing that smoking may mix acquisition
  and progression pathways.
- Lonnroth et al. and later alcohol meta-analyses support increased TB risk
  for heavy alcohol use or alcohol use disorder, but exposure definitions vary.
- CKD/dialysis review evidence identifies high TB incidence in CKD populations,
  with pooled TB incidence about 3,718 per 100,000 and higher estimates for
  hemodialysis and peritoneal dialysis groups; applicability depends on CKD
  stage and care setting.
- Chronic airway disease evidence shows COPD-associated incident TB hazard
  ratios ranging from 1.44 to 3.14 across high-income cohort studies, with
  substantial heterogeneity and limited high-burden-country evidence.
- The child close-exposure individual-participant meta-analysis included
  137,647 exposed children from 46 cohorts and supports age/contact-specific
  review; it does not justify reusing a generic adult contact OR as a
  progression hazard multiplier.
- ABS 2021-2023 life tables are a suitable Australian all-cause mortality
  source for future competing-risk inputs; WHO GHO life tables provide a
  consistent cross-country alternative. No mortality table is adopted as a
  default in Milestone 2D.

## Candidate natural-history parameterisations

These are proposals for review, not implemented production defaults.

| candidate | recent-window definition | early/remote progression | age modification | reinfection treatment | evidence quality and limitations |
| --- | --- | --- | --- | --- | --- |
| Conservative | Five-year recent window retained for compatibility with TBI target definitions | Lower early hazard, low remote hazard, weak early:remote contrast | none initially; age used only through infection history until reviewed | recent reinfection reset retained as assumption | least aggressive; may understate short-term disease among very recent infections |
| Central | Five-year recent window, with diagnostics also reporting two-year sensitivity | Fixed early:remote ratio plus fitted scale or externally supplied hazards from reviewed synthesis | consider child/older-age sensitivity only after mortality and age-specific evidence review | recent reinfection reset flagged for scientific sign-off | compatible with current two-phase functions but evidence may favor continuous decline |
| Higher-progression sensitivity | Two-year high-risk core inside five-year recent state, or higher early hazard over remaining recent window | higher early hazard and/or higher early:remote ratio; remote hazard externally bounded | optional higher risk in young children and clinically reviewed comorbidity strata | reset assumption retained but stress-tested | useful for sensitivity; not a default without stronger local evidence |

If the evidence review concludes that risk declines continuously with time
since infection, the current two-phase constant-hazard approximation should
be treated as a simplification and compared against a continuous-time hazard
shape before production integration.

## Active-TB observation examples after 2D

| example | fitting permitted? | policy/use | denominator retained | at-risk population | required information |
| --- | --- | --- | --- | --- | --- |
| Baseline prevalence found by screening | no | validation only or baseline active-TB assignment model | source denominator retained | not prospective TBI progression | prevalence ascertainment and baseline active-TB sequencing |
| Incident cases prospectively observed after baseline | yes, if eligible | external hazards or Policy B/other reviewed policy | source denominator retained | non-baseline-active-TB population | explicit q, numerator excludes baseline disease, composition and horizon |
| Retrospective notifications over previous years | no | validation only unless retrospective reconstruction exists | source denominator retained | not present baseline cohort | migration, mortality, turnover, ascertainment and historical infection pressure |
| Mixed prevalent and incident numerator | no | validation only or split after source review | source denominator retained | ambiguous | numerator decomposition |
| Changing population denominator with changed detections | rerun required | calibration/validation rerun with both old and new source rows visible | each source denominator retained | recalculated from modelled baseline composition | source numerator, denominator, period and q for each row |
| Incomplete ascertainment | yes only if q fixed | same eligible policy with q retained | source denominator retained | non-baseline-active-TB population | q, source and review status |

When stakeholders change the population denominator and observed/expected TB
detections change, calibration must rerun. The original source numerator,
denominator and observation period remain visible; the model does not overwrite
the observed denominator with a derived at-risk count.

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
- New pure calibration and assignment functions accept explicit ages and
  weights. They are not consumed by inherited population generation yet.

Event ledger and economics:

- `engine/apy/event_ledger.py`, `engine/apy/event_ledger_economics.py`,
  `engine/apy/frozen_reference.py`, `engine/apy/economics.py` and DALY-related
  code are unchanged.
- Frozen reference loaders and recalculation logic remain tied to the SA Health
  compatibility identifiers.

Dynamic model:

- `engine/dynamic/*` and integration comparison contracts are outside
  Milestones 1 and 2A and unchanged.

## Planned integration

Milestones 1, 2A, 2B and 2C have added explicit calibration, isolated
assignment helpers, pure prospective progression mathematics and policy-level
progression-calibration diagnostics. Future integration must still decide how
to connect these pieces to the APY runner while preserving the frozen SA
Health compatibility workflow. That work should include cache-key updates,
metadata propagation, baseline active-TB sequencing and migration tests before
any event-ledger, economics or DALY integration.

## Limitations and unanswered decisions

- Progression hazards are calibrated only in the isolated diagnostic Policy B
  fixed-ratio scale helper; no production policy is wired into the runner.
- The runner is not yet able to distinguish first infection, most recent
  infection and recent reinfection reset mechanisms.
- Recent reinfection resetting the higher-risk clock is a modelling assumption
  requiring scientific review.
- No risk-factor acquisition effects are included beyond age/time alive.
- Exact limiting targets may imply infinite hazards.
- Active-TB observations are validated, classified and optionally used in pure
  diagnostic expected-count/fixed-ratio calculations, but not fitted in the
  main model.
- Competing mortality is supported only through an external death hazard or a
  validated external survival curve; no mortality data are invented.
- Parameter uncertainty is not propagated.
- The new pathway is not decision-ready and is not for denying care.
