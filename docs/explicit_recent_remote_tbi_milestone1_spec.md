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

Milestones 1, 2A and 2B do not use non-age risk factors in the new acquisition
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
review. Milestone 2B does not use them to fit progression hazards.

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

Milestones 1, 2A and 2B have added explicit calibration, isolated assignment
helpers and pure prospective progression mathematics. Future integration must
still decide how to connect these pieces to the APY runner while preserving
the frozen SA Health compatibility workflow. That work should include reviewed
progression-calibration policy, cache-key updates, metadata propagation,
baseline active-TB sequencing and migration tests before any event-ledger,
economics or DALY integration.

## Limitations and unanswered decisions

- Progression hazards are not calibrated here.
- The runner is not yet able to distinguish first infection, most recent
  infection and recent reinfection reset mechanisms.
- Recent reinfection resetting the higher-risk clock is a modelling assumption
  requiring scientific review.
- No risk-factor acquisition effects are included beyond age/time alive.
- Exact limiting targets may imply infinite hazards.
- Active-TB observations are validated and classified for future use but not
  fitted.
- Competing mortality is not yet incorporated into prospective progression.
- Parameter uncertainty is not propagated.
- The new pathway is not decision-ready and is not for denying care.
