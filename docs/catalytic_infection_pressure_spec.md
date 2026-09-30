# Background exposure to TB infection: catalytic extension specification

Status: **designed, not implemented. Requires scientific review.** This document specifies an
optional extension. It changes no calculations. The only code added is a non-calculating
configuration schema (`engine/profiles/background_exposure.py`, section 9), which nothing in the
engine or the ordinary interface reads.

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
a probability, and not TB disease incidence per 100,000. The schema rejects other units and
rejects values above 1 per person-year as probable unit errors. (This bound is provisional: a
hazard of 1 corresponds to an annual infection probability of about 63%.)

**Time-series rules**, to be confirmed in review. The schema already enforces the first two.
* Each row gives a constant hazard over the calendar year [y, y + 1). Rows are ordered by year,
  and duplicate years are rejected.
* A missing year stays *missing*. It is never read as zero. Proposed handling: implementation
  refuses to run when a missing year falls inside the analysis window, unless the user chooses
  one of two explicit rules, and the choice is recorded in provenance:
  * *carry forward*: use the last observed value;
  * *linear interpolation of the hazard* between observed years.
* Years before the first row or after the last row, within the analysis horizon, need an
  explicit extrapolation rule: hold the last value constant (proposed default) or refuse to run.
  The rule is recorded.
* The series must cover the full follow-up horizon after the rule is applied. The analysis start
  year anchors engine time 0 to calendar time.
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

## 4. Scientific decisions required

Proposed positions are given for review. None is adopted until review.

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
| 12 | **Age-specific hazard.** | Supported for the constant mode through contiguous age bands starting at 0, with the last band open-ended. Age-by-year hazards are not specified (the schema rejects them in time-series mode). Age is the person's current age, which advances during follow-up. |
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
provenance. The schema rejects negative, non-finite and implausibly large hazards, duplicate
years, and WHO-snapshot provenance. Serialisation uses sorted keys, and the SHA-256 hash is
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
2. The time-series missing-year, interpolation and extrapolation rules (section 2).
3. The plausibility bound on λ.
4. Short name and output wording once exposure is enabled (section 7).
5. The source of any default λ. None is proposed: the default remains mode `none`.
