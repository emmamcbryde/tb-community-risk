# Background exposure to TB infection: catalytic extension specification

Status: **designed, not implemented in the model. Requires scientific review.** This document
specifies an optional extension. It changes no calculations. The code added so far is a
configuration schema (`engine/profiles/background_exposure.py`, section 9) and the isolated
first-infection mathematics (`engine/profiles/background_exposure_hazard.py`, section 3a). Only
tests import them; nothing in the engine, event ledger, economics or interface reads them.

Background infection pressure is an externally supplied exposure hazard. It is not inferred automatically from reported active-TB incidence.

Related: `no_transmission_model_scope.md` (model boundary and applicability) and
`recent_remote_infection_history.md` (the existing experimental historical-hazard pathway).

## 1. Purpose and core distinction

The current model assigns infection status at baseline and simulates no infection afterwards. In
settings with non-negligible local exposure, that can overstate the benefit of preventive
treatment (treated people may be reinfected) and understate future TB risk.

The extension adds an **exogenous force of infection**, a catalytic-model input:

* an externally supplied annual infection hazard λ(t, a) may generate new infection in
  uninfected people and reinfection in previously infected people;
* λ does **not** depend on the number of infectious people generated within the model;
* screening, preventive treatment or any other intervention must **not** feed back into λ;
* settings where that feedback matters (interventions that change transmission) belong in the
  Starsim workstream.

**Terminology.** In ordinary interface text, use "background exposure to TB infection". Use
"catalytic model" and "exogenous force of infection" only in technical documentation. Do not
call it "transmission" or "importation" (section 6).

## 2. Future-exposure modes

| Mode | Schema value | Definition | Status |
| --- | --- | --- | --- |
| No ongoing exposure | `none` | λ(t, a) = 0 exactly for all t after baseline. Reproduces the current model and its results. | Current behaviour; the default. |
| Constant background exposure | `constant` | λ(t, a) = λ₀, or an age-specific λ(a), constant in calendar time. | Designed. |
| Changing background exposure | `time_series` | λ(t) from a validated, user-supplied annual series. | Designed. |

**Units.** λ is a hazard (rate) in **infections per person-year**, for susceptible people. It is not
a probability, and not TB disease incidence per 100,000. The schema rejects other units, and
negative or non-finite values. A hazard is a rate, so values above 1 per person-year are
mathematically valid and accepted. They raise a **non-blocking** plausibility warning, because they
usually indicate a unit error. The warning threshold (1 per person-year, an annual infection
probability of about 63%) is provisional and not evidence-based. It is not a validity limit.

**Time-series rules** (adopted; implemented in the mathematics, section 3a):
* Each row gives a constant hazard over the half-open calendar year [y, y + 1). Rows are ordered
  by year, and duplicate years are rejected.
* A missing year stays *missing*. It is never read as zero.
* The series must cover every year that a requested interval touches. A gap, a year before the
  first row or a year after the last row is a validation error. There is **no automatic
  interpolation or extrapolation**: values are never carried forward, backfilled, interpolated
  across missing years, or replaced with zero. Any explicit fill rule would be a separate,
  recorded scientific decision; none is provided.
* Anchoring engine time 0 to calendar time is a wiring decision, made when the mathematics is
  connected.
* The series is **never inferred automatically from WHO TB disease incidence**. Incidence reflects
  past infection, progression and detection; it is not a hazard of infection
  (`INCIDENCE_TO_INFECTION_POLICY` in `engine/model_scope.py`). The schema rejects the
  incidence quantity and WHO-snapshot provenance.

## 3. Mathematical definition

Hazard λ(t, a) is in infections per person-year at calendar time t and age a. For a person who
is susceptible throughout the interval [t, t + Δt):

```text
P(infection during [t, t + Δt)) = 1 - exp( - ∫_t^{t+Δt} λ(s, a + s - t) ds )
```

For a constant hazard this reduces to `1 - exp(-λ Δt)`. The annual probability `p = 1 - exp(-λ)`
must never be entered in place of λ, or vice versa. They differ materially at high exposure
(λ = 0.05 gives p = 0.0488; λ = 0.5 gives p = 0.393).

For a piecewise-constant annual series, the cumulative hazard over [t₁, t₂) is

```text
H(t₁, t₂) = Σ_y λ_y × | [y, y + 1) ∩ [t₁, t₂) |
```

and the time of first infection is sampled by inversion: draw U ~ Uniform(0, 1), then find the
smallest τ with H(t₁, t₁ + τ) = -log(U).

**Baseline prevalence from a historical hazard.** For a person of age a at baseline time t₀,
ignoring differential mortality, self-clearance and reinfection:

```text
P(infected by age a at time t₀) = 1 - exp( - ∫_0^a λ_hist(t₀ - u, a - u) du )
```

Here u is time before baseline (the lookback), and a - u is the person's age at that time.
Because of these simplifications the formula gives "ever infected", not "currently infected
and at risk of progression". The existing experimental pathway (`engine/apy/infection_history.py`,
documented in `recent_remote_infection_history.md`) implements a parametric version:
`λ(a, τ) = s·exp(-g τ)·(a - τ + 0.5)^γ`.

**Identifiability.** Cross-sectional LTBI prevalence by age, which is the typical data, does not
uniquely identify the historical trajectory λ_hist(t, a). Period effects (calendar time) and
age effects trade off along each birth cohort's diagonal. Many different trajectories fit the
same prevalence curve, and they imply different recent-versus-remote mixes and hence different
progression risk. A trajectory is therefore a scenario assumption unless it is informed by other
data (repeated cross-sections, tuberculin surveys, cohort data, or an explicit prior). The
existing pathway already fixes the slope g as a scenario for this reason.

## 3a. Adopted mathematical conventions (implemented and validated, not connected)

`engine/profiles/background_exposure_hazard.py` provides pure, deterministic functions. They draw
no random numbers.

| Function | Meaning |
| --- | --- |
| `HazardSchedule.none()`, `.constant(λ)`, `.annual_series({year: λ})`; `schedule_from_config(config)` | Executable hazard λ(t). |
| `hazard_at(schedule, t)` | λ(t). |
| `cumulative_hazard(schedule, t0, t1)` | H(t0, t1) = ∫ λ(s) ds. |
| `infection_probability(schedule, t0, t1)`, `probability_from_cumulative_hazard(H)` | P = 1 - exp(-H), computed as `-expm1(-H)`. |
| `exponential_threshold(U)`, `first_infection_time(schedule, t0, t1, uniform=… or threshold=…)` | Inversion for one person. Returns a typed result: `InfectionOccurs(time)` or `NoInfectionInInterval(cumulative_hazard, residual_threshold)`. |
| `first_infection_times(schedule, t0, t1, thresholds)`, `infection_probabilities(...)` | Vectorised forms with the same semantics. |
| `requires_random_draw(schedule)` | False for mode `none`. |
| `plausibility_warnings(schedule)` | The non-blocking warnings described in section 2. |

Conventions:

* **Hazard, not probability.** λ is in infections per susceptible person-year. For constant λ over
  Δt, P = 1 - exp(-λΔt). λ is never used as a probability. Hazards must be finite and
  non-negative, and there is no upper limit.
* **Time.** Calendar time is in decimal years. Intervals are half-open [t0, t1) with t0 ≤ t1. A
  reversed interval or a non-finite time is an error. A zero-duration interval gives H = 0 and
  P = 0, needs no coverage, and gives no infection.
* **Piecewise-constant series.** Year y applies on [y, y + 1). An interval needs the years from
  floor(t0) to ceil(t1) - 1, so an interval ending exactly at y + 1 does not need year y + 1.
  Coverage gaps are errors (section 2).
* **Inversion.** The caller supplies U in the **open** interval (0, 1), so 0, 1, values outside
  [0, 1] and non-finite values are rejected. Alternatively the caller supplies E = -log U directly
  (finite, E ≥ 0). Infection occurs in [t0, t1) **if and only if E < H(t0, t1)**, at
  τ = inf{s : H(t0, s) > E}: the first time at which accumulated hazard exceeds E, which skips
  zero-hazard years. If E ≥ H, the result is "no infection in this interval" and carries the
  residual threshold E - H. Passing that residual to the next interval gives the same infection
  time as a single call over the whole interval. With this convention
  P(infection) = P(E < H) = 1 - exp(-H) exactly.
* **Numerics.** H is summed with `math.fsum`. P uses `expm1`, which is accurate for very small H
  and exactly 1.0 for very large H. An overflowing H is an error.
* **Scope: first infection only**, for people susceptible at t0. Reinfection (of remotely infected
  people, after preventive treatment or after active-TB treatment), partial immunity, repeated
  screening and progression after a new infection are **not implemented and not decided**
  (section 4).
* **Baseline LTBI prevalence remains a separate input.** Baseline infection status is set by the
  existing model. The prospective hazard applies only after time zero and is not used to generate
  or recalibrate baseline prevalence. Historical catalytic calibration is a later milestone.
* **Age-specific hazards: specified, not yet executable.** `background_exposure_v1` accepts
  contiguous integer age bands [ageLower, ageUpper) in constant mode. `schedule_from_config` raises
  `AgeSpecificHazardNotExecutableError` for such a configuration rather than ignoring the bands.
  Executing them needs decisions on the age representation at t0 (exact or completed years),
  whether birthdays split exposure at the exact crossing time, and how age-band boundaries
  interact with calendar-year boundaries. It also depends on decision 13 (acquisition versus
  progression modifiers).
* **No endogenous feedback.** λ is a fixed input. It does not depend on infectious people in the
  model or on any intervention.
* **Mode `none`.** λ = 0, H = 0, P = 0 and no infection. `requires_random_draw` is False, and the
  inversion functions *reject* a variate for mode `none`, so a caller cannot consume a draw. The
  `none` configuration and its hash are unchanged. A constant hazard of 0 is a value, distinct from
  `none`: it requires a draw under this API. Whether that draw comes from a separate exposure
  stream (decision 17) is a wiring decision.

**Validation results** (`tests/test_background_exposure_hazard.py`):

* **Analytical:** exact or tight-tolerance agreement with hand-calculated values:
  * H, P and inverted times for constant hazards;
  * piecewise series within one year, across two and several years, with fractional ends and at
    exact year boundaries, including additivity and chained residual thresholds;
  * a zero-hazard year skipped by inversion.
* **Numerical:**
  * relative error of P below 1e-15 at H = 1e-18 to 1e-8 (the naive formula loses about six
    significant digits at H = 1e-12);
  * P = 1.0 exactly at H = 5 × 10⁴;
  * 200-year series summed exactly;
  * monotone in hazard and duration, and bounded in [0, 1];
  * overflow rejected.
* **Validation:** errors for:
  * negative, NaN and ±∞ hazards;
  * duplicate, non-integer or missing years, and intervals outside coverage;
  * reversed intervals and non-finite times;
  * U ∈ {0, 1}, values outside [0, 1], and non-finite or negative thresholds.
* **Statistical:** fixed seed, n = 200,000 caller-supplied uniforms per scenario, and a tolerance
  of 4 standard errors fixed in advance. Each scenario compares the proportion infected and the
  first-infection-time CDF at 20%, 50% and 80% of the interval.

  | Scenario | Analytical P | Empirical | z | Max abs(z), CDF points |
  | --- | --- | --- | --- | --- |
  | constant 0.01 over 10 years | 0.09516 | 0.09455 | -0.93 | 1.39 |
  | constant 0.1 over [2020.25, 2023.75) | 0.29531 | 0.29603 | +0.70 | 0.87 |
  | constant 2.0 over half a year | 0.63212 | 0.63104 | -1.00 | 0.75 |
  | series 0.02, 0, 0.15, 0.05 over [2020.4, 2023.6) | 0.17469 | 0.17367 | -1.21 | 1.06 |
  | series 0.02, 0.05, 0.10, 0.01 over [2020, 2024) | 0.16473 | 0.16516 | +0.52 | 0.46 |

* **Compatibility and isolation:**
  * the `none` configuration and hash are identical to commit `9003146`;
  * a reference caller consumes zero draws for `none`, and the generator state is unchanged;
  * the module contains no random-number calls;
  * importing the runner, simulation, expected-value, event-ledger, economics, engine-mapping and
    general-application modules does not load it;
  * no engine, app, UI or page file imports it;
  * it imports neither Starsim nor Streamlit.

## 4. Scientific decisions required

Proposed positions are given for review. None is adopted until review, except two that are
adopted for the isolated mathematics only (section 3a):
* decisions 1 and 2: the prospective hazard is separate, and baseline prevalence stays a direct
  input;
* decision 12: age-specific execution is deferred.

| # | Decision | Proposed position |
| --- | --- | --- |
| 1 | **Historical versus prospective hazard.** | Keep them separate. λ_hist (before t₀) only establishes baseline infection status. λ_fut (after t₀) drives new infection during follow-up. Supplying λ_fut must not change baseline prevalence, and vice versa. |
| 2 | **Baseline prevalence: direct input or generated?** | It stays a direct, calibrated input by default, as now. Generating it from λ_hist is an alternative mode, never both at once, so exposure is not counted twice. |
| 3 | **Recent versus remote after incident infection.** | An incident infection enters the recent (early, higher-progression) state at its infection time, then follows the same recent-to-remote transition as baseline recent infection. |
| 4 | **Early higher-risk state: duration and interpretation.** | Keep the two concepts distinct. The current natural history uses an early state with a 5-year mean residence time. The 2-year "recent" window of the infection-history pathway is a timing label only (`recent_remote_infection_history.md`). Choose one definition for incident infection and document it. |
| 5 | **Reinfection while remotely infected.** | Candidate policies: (a) none; (b) reinfection returns the person to the recent state, with partial protection (a relative hazard of reinfection disease below 1; published estimates vary, for example Andrews et al., Clin Infect Dis 2012, reported about 79% lower risk among people already infected, to be verified in review); (c) full susceptibility. Implement (a) as a reference, make (b) the proposed primary, and use (c) as a bound. |
| 6 | **Reinfection after completed preventive treatment.** | Treated people become susceptible to λ_fut, with the same partial-protection policy as in decision 5. Preventive treatment does **not** protect against future infection. |
| 7 | **Reinfection after treatment for active TB.** | Out of scope for the first implementation unless disease after treatment is modelled. If added: susceptible to λ_fut, with a separate relapse-versus-reinfection distinction. |
| 8 | **Partial immunity after infection or disease.** | Represented only through the relative reinfection hazard in decisions 5 to 7. It is not implied by BCG status unless separately specified. |
| 9 | **Mechanism of preventive treatment.** | The current engine represents protection as cure of infection: a completed course protects with probability equal to full efficacy, a partial course with partial efficacy (`protected_full`, `protected_partial` in `engine/apy/simulation.py`). Under exposure, "cleared, then susceptible to reinfection" differs from "reduced progression hazard", and the choice changes results. Keep "clears infection" as the proposed primary, with "reduces progression" as a sensitivity. |
| 10 | **Infection before, during or after screening.** | An infection before the screening time is latent at screening and can test positive (with the recent/remote test-accuracy rules). An infection during a preventive-treatment course is proposed to be cleared by a completed course (sensitivity: not cleared). An infection after screening is not detected by a one-time screen. |
| 11 | **Repeated screening.** | Detecting later incident infections needs repeat screening rounds. That is a separate intervention design, and not part of this extension. |
| 12 | **Age-specific hazard.** | The schema specifies contiguous age bands for the constant mode, starting at 0 with the last band open-ended. They are **not yet executable**: the calculator rejects them (section 3a). Proposed for execution: age is the person's current age, which advances during follow-up, with exposure split at the exact time a band boundary is crossed. Age-by-year hazards are not specified (the schema rejects them in time-series mode). |
| 13 | **Risk factors: acquisition versus progression.** | Current risk-factor effects act on progression. Acquisition modifiers would need separate, acquisition-specific evidence. **An odds ratio for disease progression is not an infection-acquisition hazard ratio** and must never be reused as one. |
| 14 | **Competing risks, death and ageing.** | The current individual-based simulation has no background mortality. With exposure, the time at risk depends on survival, so implementation must decide whether to add death as a competing risk (proposed, from the life table used by the economics) or to state explicitly that exposure is applied to survivors-to-horizon. Ageing must update age-specific λ. |
| 15 | **Event-ledger representation.** | Add annual events: `incident_infections`, `reinfections`, and active TB split by infection origin (baseline infection versus incident infection). Existing event names and meanings stay unchanged. |
| 16 | **DALY and economic consequences.** | No new cost items. Additional active TB from incident infection flows through the existing DALY and cost logic. Outcome labels must still say "direct, among the modelled population". |
| 17 | **Stochastic pairing.** | Common random numbers: one exposure uniform stream per person, used identically in the comparator and intervention arms. Arm differences then arise only from intervention effects (for example a person cured by preventive treatment becoming susceptible). Each reinfection draw uses its own stream index so the arms stay aligned. |
| 18 | **Deterministic expected-value approximation.** | Extend the expected-value path with the closed-form infection probabilities in section 3 over annual steps. Report its agreement with a large simulation as a validation target, not as an identity. |
| 19 | **Provenance for every exposure input.** | Every λ value records source, citation, provenance (bundled, user-defined or local upload; never a WHO incidence snapshot), review status, units and applicability notes. The configuration hash enters the result-currency hash. |

## 5. Relationship to disease-transition parameters

The exposure hazard adds new entries into the infected states. It does not alter progression
hazards, test accuracy, cascade, efficacy or costs. Any parameter that changes because exposure
is enabled (for example re-calibrated progression, if the calibration target was set assuming no
ongoing infection) must be flagged and reviewed. It must not be adjusted silently.

## 6. External infection versus imported disease versus transmission

These four mechanisms are distinct and must not be merged:

| Mechanism | What it is | How it would be represented |
| --- | --- | --- |
| **Background (exogenous) force of infection** | Infection acquired by modelled people from exposure outside the model's own dynamics | This extension: λ(t, a) acting on the susceptible and reinfectable modelled population |
| **Entry of a person already infected** | A new person joining the population with latent infection (for example a migrant) | A separate demographic entry mechanism, with its own infection-status distribution. It is not λ. |
| **Entry of a person with active TB** | A new person joining with prevalent disease | A separate demographic entry mechanism. It must never be represented as force of infection. |
| **Endogenous transmission** | Infection caused by infectious people generated within the modelled population | Not represented. It requires a transmission model (Starsim workstream). |

"Importation" must not be used as a synonym for background exposure. Population entry is not
part of this extension and would need its own specification.

## 7. Compatibility and model identity

* **Exact reproducibility.** Mode `none` must leave every current output unchanged, bit for bit
  for the stochastic engine with the same seed. No new random draws may be consumed when the
  mode is `none` (the draws would shift the random stream). The configuration hash of existing
  analyses must not change.
* **Description.** Once exposure is available, the model is described as:
  > A screening and preventive-treatment model with optional exogenous infection pressure and no
  > endogenous transmission feedback.

  `engine/model_scope.py` currently holds the description without the optional exposure
  ("... with no endogenous transmission feedback"). It is updated only when the extension is
  implemented.
* **Short name.** "No-transmission model" remains accurate in the sense that matters for
  applicability: the model has no transmission *dynamics* and no intervention feedback. With
  exposure enabled, though, people in the model do acquire infection, so a reader could take
  "no transmission" to mean "no infection after baseline". Proposal: keep
  "no-transmission model" as the workstream name. In outputs, use "no transmission feedback"
  or "no endogenous transmission", and state the exposure mode (for example "background
  exposure: none"). This is a decision for review.

## 8. Applicability

Background exposure widens the range of settings where the model may be useful. Examples are
populations with modest ongoing exposure, or long horizons where reinfection after preventive
treatment matters. It does **not** make the model suitable where intervention-induced changes in
transmission are material: the exposure is fixed and cannot respond to the intervention. The
applicability framework in `no_transmission_model_scope.md` stays advisory. No incidence value
(10, 40, 100 per 100,000 or any other) is a validated engine-selection threshold.

## 9. Configuration schema (scaffolding only)

`engine/profiles/background_exposure.py`, schema `background_exposure_v1`:

| Field | Meaning |
| --- | --- |
| `schemaVersion` | `background_exposure_v1` |
| `mode` | `none`, `constant` or `time_series` |
| `quantity` | Fixed: `exogenous_infection_hazard`. `estimated_tb_disease_incidence` is rejected. |
| `unit` | Fixed: `infections per person-year` |
| `constantHazard` | One hazard value (constant mode, all ages) |
| `ageSpecificHazards` | Contiguous bands `{ageLower, ageUpper, hazard}` from 0, last open-ended (constant mode) |
| `timeSeries` | Rows `{year, hazard}`, sorted, unique years; a hazard may be missing but not every hazard (time-series mode) |
| `source`, `citation`, `reviewStatus`, `applicabilityNotes` | Provenance for the configuration as a whole |

Each hazard value reuses the profile `ProfileValue` contract (value, state, provenance, review
status, unit, source, notes), so zero and missing are distinct and user-defined values keep their
provenance. The schema rejects negative and non-finite hazards, duplicate years and WHO-snapshot
provenance. Hazards above the provisional threshold produce non-blocking warnings
(`plausibility_warnings()`). Serialisation uses sorted keys, and the SHA-256 hash is
deterministic. Profile payloads without a `backgroundExposure` block (all existing and legacy
profiles) map to mode `none`. The population-profile contract and its hash are unchanged.
Adding the block to the profile contract, as a future `population_profile_v3`, is a later
decision.

## 10. Analytical validation plan (before implementation is used)

| Test | Pass criterion |
| --- | --- |
| Constant-hazard infection probability | Simulated proportion infected by time t matches `1 - exp(-λt)` within Monte Carlo error (for example 4 standard errors) across λ ∈ {0.001, 0.01, 0.1} and several horizons |
| Piecewise-constant series | Cumulative hazard and first-infection times match the closed form of section 3, including partial years and a change of hazard mid-follow-up |
| Zero-hazard identity | Mode `none`, and constant λ = 0, reproduce every current output exactly with the same seed; no extra random draws |
| Age-specific cumulative prevalence | Generated baseline prevalence by age matches the integral in section 3 for constant and age-banded hazards |
| Incident infection timing | The distribution of infection times among the uninfected matches the exponential (or piecewise-exponential) law |
| Recent-to-remote transition | After incident infection, residence in the early state matches the specified duration distribution |
| Reinfection policies | Each candidate policy (none, partial protection, full susceptibility) gives the analytically expected reinfection rate among the remotely infected and among people after preventive treatment |
| Conservation and exclusivity | Each person is in exactly one disease state at each time; counts are conserved annually |
| Deterministic versus simulation | The expected-value path agrees with a large simulation (for example 10⁶ person-histories) within a pre-agreed tolerance |
| Common random numbers | With the intervention switched off, arms are identical person by person; the variance of paired differences is below that of unpaired differences |
| No feedback | λ used in the intervention arm is identical to λ in the comparator arm for every person and year, whatever the intervention outcome (asserted in code and tested) |

**Cross-model comparison (future).** Compare Starsim with transmission feedback disabled, the
same exogenous λ, matched natural-history parameters, a matched starting population (age,
infection status, recent/remote mix) and agreed outcome definitions (cumulative active TB, direct
cases averted, incident infections, reinfections). No equivalence is claimed until this
comparison has been run and meets pre-agreed tolerances.

## 11. Open decisions (summary)

1. Everything in section 4, especially decisions 5, 9, 13 and 14.
2. Whether any explicit, recorded fill rule for missing time-series years should ever be offered.
   At present there is none, and gaps are errors.
3. The provisional plausibility-warning threshold on λ (non-blocking).
4. How age-specific hazards execute: age representation, birthdays and band crossings
   (decision 12).
5. Short name and output wording once exposure is enabled (section 7).
6. The source of any default λ. None is proposed: the default remains mode `none`.
