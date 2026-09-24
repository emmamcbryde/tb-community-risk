# Dynamic TB transmission model: readiness specification for country-profile integration

Status: Designed, not implemented. Audit of current code; no formulas changed.

Scope: this document audits the Python dynamic model on branch
`feature/general-community-tb-model` and specifies a later integration
milestone. Citations are `file:line` against the current working tree. Where
the code does not specify something it is marked "not represented in current
code". Note: branch `codex-supervisor-test` contains unmerged commits
(`0364d32` "two-epoch beta calibration", `3ef127f` "secular decline
calibration") and tests `tests/test_dynamic_model.py` and
`tests/test_dynamic_ui_calibration.py` that are not on this branch (only stale
`.pyc` files remain in `tests/__pycache__/`). They are out of scope here.

## 1. Files audited

- Engine: `engine/dynamic/dynamic_model.py` (`simulate_dynamic`, line 4);
  `engine/dynamic/exec_dynamic.py` (`run_dynamic_model`, line 3, a pass-through).
- Calibration and UI: `ui/dynamic_ui.py` (`render_dynamic_ui`, line 478), used by
  `pages/5_Dynamic_Model.py` and `ui/app.py`. `ui/update.py` is an older
  near-duplicate (same constants, its own `render_dynamic_ui` at line 408) and is
  not imported by any page.
- Age helpers: `ui/dynamic_age_distribution.py`.
- Infection backcast: `engine/infection_backcast.py`.
- Output contract and comparison: `engine/integration/dynamic_output_contract_v9.py`,
  `engine/integration/compare_dynamic_abm_v9.py`,
  `engine/integration/risk_factor_crosswalk_v9.csv`,
  `docs/apy_dynamic_metric_alignment_v9.md`.
- Legacy static model (not used by the dynamic model): `engine/model.py`,
  `engine/intervention.py`, `engine/params.py`.
- Profile and incidence contracts: `engine/profiles/population_profile.py`,
  `engine/profiles/engine_mapping.py`, `engine/who_incidence/trend.py`.
- Tests: `tests/test_dynamic_output_contract_v9.py`,
  `tests/test_dynamic_age_distribution.py`,
  `tests/test_dynamic_abm_compare_page_helpers.py`,
  `tests/test_compare_dynamic_abm_v9.py`. None of these call `simulate_dynamic`,
  `calibrate_beta_and_ltbi_scale`, `refine_beta_random_walk` or
  `engine/infection_backcast.py`. The transmission engine and its calibration
  therefore have no unit tests on this branch.

## 2. Current model structure

### 2.1 Compartments

Five aggregate (not age-structured) compartments, counts of people
(`dynamic_model.py:196-200`):

| State | Meaning |
|---|---|
| `S` | Uninfected (also receives people cleared by LTBI treatment) |
| `L_fast` | Recent latent infection |
| `L_slow` | Remote latent infection |
| `I` | Active TB disease (infectious; undetected and pre-treatment) |
| `R` | Detected/treated; absorbing, non-infectious |

### 2.2 Flows (`dynamic_model.py:253-271`)

| Flow | Rate expression | Parameters (code values) |
|---|---|---|
| S -> L_fast (infection) | `lambda_t * S` | `lambda_t = beta_now * I / N` (line 253) |
| L_fast -> L_slow (stabilisation) | `omega * L_fast` | `omega = 1.0 / 5.0` per year (line 94) |
| L_fast -> I (fast progression) | `sigma_fast_eff * L_fast` | `sigma_fast = 0.01` per year (line 92) x `risk_multiplier` |
| L_slow -> I (reactivation) | `sigma_slow_eff * L_slow` | `sigma_slow = 0.001` per year (line 93) x `risk_multiplier` |
| L_fast -> S, L_slow -> S (LTBI treatment) | `tau_now * L` | `tau_at_time` (line 210) |
| I -> R (detection plus treatment) | `gamma_now * I` | `gamma_at_time` (line 217), from `delta_pre`, `delta_post` |

Not represented in current code: reinfection of `L_fast`, `L_slow` or `R`
(only `S` is infected); relapse from `R`; TB-specific mortality;
background mortality; self-cure; loss to follow-up; treatment failure; drug
resistance; age- or sex-specific progression; progression dependence on time
since infection other than the two-state fast/slow split.

### 2.3 Demography

Closed population. There are no births, deaths or ageing. `N` is recomputed
each step as the compartment sum (line 231) and is conserved by the ODEs except
where the `max(.,0)` clipping (lines 281-285) truncates an overshoot (see 2.9).
Age enters only through the initial seeding (`age_counts`, `ltbi_ever`,
`ltbi_recent`, lines 174-179), after which age information is discarded.
`age_counts` is a dict single-year age -> count, built in
`ui/dynamic_ui.py:644-648` from an OWID CSV (`load_population_data`, line 125,
default `data/population_age_latest.csv`), a manual 5-year table expanded
uniformly within bands (`expand_five_year_age_distribution`,
`ui/dynamic_age_distribution.py:56`), an upload, or a uniform 0-100 default.
Population size is the UI input `Population size` (default `10000`,
`ui/dynamic_ui.py:492`).

### 2.4 Transmission

- Force of infection: `lambda_t = beta_now * I[i-1] / N` (line 253),
  frequency-dependent, homogeneous mixing, no contact matrix, no age mixing, no
  differential infectiousness (e.g. smear status), no external (imported)
  infection.
- `beta`: scalar `params.get("beta", 0.0)` (line 54), or piecewise-constant per
  model year from `params["beta_series"]` indexed by `floor(t)` and clamped to
  the last element (lines 243-251). Units: effective infections per prevalent
  case per year.
- Latent people are fully protected from reinfection (no partial-immunity
  parameter). Not represented in current code.

### 2.5 Detection

A single exit `I -> R` at rate `gamma`, with `delta_pre = params.get("delta_pre",
12.0 / max(pre_det_months, 0.1))` and `pre_det_months` default `12.0`
(lines 78-80). The UI fixes `BASELINE_DIAG_MONTHS = 12.0`
(`ui/dynamic_ui.py:41`), i.e. `delta_pre = 1.0` per year. Under intervention,
`gamma` ramps linearly from `delta_pre` to `delta_post` over `rollout_years`
(lines 217-223); `delta_post = 12 / post_det_months` with
`post_det_months = max(12 * (1 - diag_reduction_pct/100), 0.1)`
(`ui/dynamic_ui.py:1007-1011`, slider default `50`). Because `I` has no other
exit, every incident case is eventually detected: the implied long-run case
detection ratio is 1 and detection only changes the delay. Notifications are
not output (the `gamma * I` flow is not accumulated).

### 2.6 Treatment of disease and infectiousness during treatment

Treatment is instantaneous on entry to `R`; `R` is non-infectious and absorbing.
Infectiousness during treatment, treatment duration, outcomes and post-TB
states: not represented in current code.

### 2.7 Fast/slow latency

Two-state latency (2.2). Implied probability of progression via the fast route
before stabilisation is `sigma_fast / (sigma_fast + omega)` = about 4.8%
(times `risk_multiplier`). All constants are hard-coded in `simulate_dynamic`
and cannot be overridden via `params`.

### 2.8 Intervention representation

- LTBI test-and-treat (lines 127-141, 210-215): `coverage_total = ltbi_coverage
  * eff`, where `eff` comes from a local table `regimen_efficacy = {"None": 0.0,
  "1HP": 0.90, "3HP": 0.90, "4R": 0.80, "6H": 0.70, "9H": 0.75}` (line 127).
  `annual_fraction_treated = coverage_total / rollout_years` is applied as a
  continuous rate `tau` on both latent compartments while
  `time <= rollout_years`, moving people to `S` (reinfectable). Because it is a
  rate on a depleting stock, the fraction cleared over rollout is
  `1 - exp(-coverage_total)`, not `coverage_total` (e.g. 0.36 not 0.45 for
  coverage 0.5 and 3HP).
- `testing_method` is passed in `params` (`ui/dynamic_ui.py:1043`) but is never
  read by `simulate_dynamic`: no test sensitivity, specificity, false
  positives, treatment completion, adverse events or screening of `S` or `I`.
- The efficacy table differs from `engine/intervention.py:10` `REGIMENS`
  (e.g. `3HP` `0.92`, `4R` `0.90`, `6H` `0.65`, `9H` `0.85`, with separate
  `completion`), which is used only by the legacy static model.
- Diagnosis improvement: the `gamma` ramp in 2.5.

### 2.9 Risk factors

`risk_multiplier = prod_k ((1 - p_k) + p_k * RR_k)` (lines 117-119), with
hard-coded `RR = {"smoker": 1.5, "alcohol": 2.0, "diabetes": 3.0, "renal": 2.5,
"HIV_treated": 4.0, "HIV_untreated": 10.0}` (line 99) and `p_k` from the six
required `*_pct` params (lines 29-39, 65-70). The multiplier applies to both
progression rates; it assumes independent, multiplicatively combined factors
and a population-averaged (unstratified) hazard. It does not affect
infection, detection or infectiousness. With UI defaults (30, 15, 10, 5, 3, 3
per cent; `ui/dynamic_ui.py:499-504`) it equals about 2.36. Effect-measure type
is not represented. `risk_factor_crosswalk_v9.csv` maps four individual-engine
factors to dynamic keys with `comparability` `partial`; `close_contact`,
`marijuana`, `chronic_lung_disease` are individual-engine only and the two HIV
inputs are dynamic only. The crosswalk is documentation; no code reads it.

### 2.10 Time step and solver

Forward Euler, `dt = 0.1` years (line 192), `n_steps = round(years / dt)`.
Annual incidence accumulates the left-endpoint progression flow
`new_cases_rate * dt` into `annual_incidence[year_idx]` (lines 258-262).
Observed defects (verified by a scratch run, no repository changes):

1. Explicit Euler is unstable when `gamma * dt > 1`, i.e. `delta_post > 10`
   (`post_det_months < 1.2`, `diag_reduction_pct > 90`). At
   `diag_reduction_pct = 100`, `delta_post = 120`; `I` is clipped to 0 while `R`
   gains `12 * I`, and `N` grew from 10100 to 10390 over 20 years.
2. `year_idx` is clamped to `len(beta_series) - 1` (line 250) before incidence
   is accumulated, so if `beta_series` is shorter than `years` all later
   incidence is added to the last beta year and later years read 0. Current
   UI calls always pass a series of length `years`, so this is latent.
3. `years` is effectively integer: `annual_incidence` has length `int(years)`.

### 2.11 Initial conditions and backcasting

Default seeding (lines 168-187): `S0 = sum pop*(1-ltbi_ever)`,
`L_fast0 = sum pop*ltbi_recent`, `L_slow0 = sum pop*(ltbi_ever-ltbi_recent)`,
`I0 = N * initial_incidence_per_100k/1e5 * pre_det_months/12` taken from `S0`,
`R0 = 0`. Stitch mode (lines 148-166): `params["initial_state"]` supplies all
five states; `S` absorbs any mismatch with `sum(age_counts)`.

LTBI by age comes from `engine/infection_backcast.py`:

- `calc_ari_from_incidence` (line 5): Styblo-style
  `ARI_t = inc_t * f_ssplus / K * adjustment`, `f_ssplus=0.6`, `K=5000.0`,
  clamped to [0,1].
- `infection_prob_by_age_split` (line 45), `window_recent=5`:
  `ever[a] = 1 - prod_{k<a}(1-ARI_{-k})`; `remote[a]` = infected more than 5
  years ago and not in the last 5; `recent[a] = ever - remote`. No mortality,
  no clearance, no age dependence of ARI; ARI before the supplied history
  repeats the oldest value; infants (`a <= 0`) are uninfected. Recent includes
  reinfection of previously infected people.
- `compute_ltbi_from_inc_hist` (`ui/dynamic_ui.py:247`) applies
  `ARI_FLOOR = 1e-6` and shifts the history by `shift_years`, but applies it to
  the current age structure (the age structure at the start of the fit window
  is not reconstructed).
- `build_incidence_history` (`ui/dynamic_ui.py:162`) produces the history from
  a pattern (`Constant`, `Falling 3%/year`, `Rising 3%/year`) or an uploaded
  `year, incidence` CSV with geometric-trend extrapolation, 3-year centred
  smoothing for synthetic patterns, and `INCIDENCE_FLOOR = 0.1`.

Consequence: the pre-window infection history is set by assuming ARI is
proportional to disease incidence (constant `adjustment`). This is the
assumption that `engine/who_incidence/trend.py:117`
`INCIDENCE_TO_INFECTION_POLICY` says must not be made without an explicit
linkage model. The WHO snapshot and population profile are not wired to the
dynamic model (`engine/profiles/engine_mapping.py:14-15`).

### 2.12 Outputs

`simulate_dynamic` returns `time`, `annual_incidence_time`,
`annual_incidence` (counts per year interval), `annual_prevalence_I`, the five
state trajectories at `dt` resolution and `final_state` (lines 305-316). The UI
builds `df_future` (`ui/dynamic_ui.py:1068-1089`) with per-arm incidence per
100k, annual and cumulative counts and cases averted, adds placeholder bounds
`DUMMY_CI_PCT = 20.0` (line 62; the code gives +/-20 per cent although the
comment and tooltip say +/-10 per cent), and packages it with
`build_dynamic_results_bundle_v9` (`dynamic_output_contract_v9.py:90`,
`CONTRACT_VERSION = "dynamic_output_contract_v9"`,
`MODEL_VERSION = "dynamic_python_v1"`). `compare_dynamic_abm_v9`
(`compare_dynamic_abm_v9.py:221`) aligns six metrics (`ALIGNED_METRICS`,
line 9) and lists structurally non-comparable ones (line 33). Not output:
notifications, deaths, prevalence of infection by age, FOI, TPT courses,
screening counts.

## 3. Current calibration

Calibration exists, in the UI layer, not in `engine/` (`ui/dynamic_ui.py`).

- Target: annual disease incidence per 100k, point values only, over the last
  `CALIB_YEARS_FIT = 20` years of the constructed history (year-end alignment,
  lines 296-301). Source is a synthetic pattern or a user CSV; WHO bounds and
  the WHO snapshot are not used. No notifications, prevalence or LTBI targets.
- Stage 1, `calibrate_beta_and_ltbi_scale` (line 272): grid of
  `ARI_ADJ_GRID_POINTS = 21` values over `ARI_ADJ_BOUNDS = (0.05, 5.0)`; for
  each, `minimize_scalar(..., method="bounded")` over `BETA_BOUNDS = (0.01,
  50.0)` minimising RMSE (31-point grid fallback without SciPy).
- Stage 2, `refine_beta_random_walk` (line 356): one log-beta per year
  (20 parameters), L-BFGS-B (`maxiter` 120), objective
  `mean(((pred-obs)/obs_scale)^2) + BETA_RW_WEIGHT * mean((diff(x)/sigma_rw)^2)`
  with `BETA_RW_PCT = 10`, `sigma_rw = log(1.1)`, `BETA_RW_WEIGHT = 0.005`,
  `obs_scale = max(mean(obs), 1)`. `ari_adjustment` is held at the stage-1
  value.
- Fixed during calibration: `sigma_fast`, `sigma_slow`, `omega`,
  `risk_multiplier`, `delta_pre = 1.0`, no intervention.
- Projection (lines 996-1052): starts from the calibrated `final_state` with
  constant `beta_forward = beta_series_hat[-1]`.
- Uncertainty propagation: absent. Single point fit; the projection band is
  the fixed +/-20 per cent placeholder. No posterior, no profile likelihood, no
  fit diagnostics other than RMSE.

## 4. Identifiability of the current structure

| Data | Separately identifiable | Confounded / not identifiable |
|---|---|---|
| Incidence alone (current) | Overall level of progression flow; with fixed natural history, roughly the pair (`ari_adjustment`, mean `beta`), via early-window reactivation versus later transmission-driven contribution | `beta` vs `gamma` (only `beta / gamma` matters near equilibrium, and `gamma` is fixed); `beta(t)` year-to-year (20 parameters for 20 points, determined by the smoothness penalty, and the progression kernel lags FOI changes by years); `ari_adjustment` vs `risk_multiplier` vs `sigma_*` (all scale incidence) |
| Incidence plus notifications | Under current structure notifications equal lagged incidence (CDR = 1), so they add information only on `gamma` timing | CDR itself (needs a non-detected exit such as death or self-cure, or a reporting fraction) |
| Plus age-specific LTBI prevalence | `ari_adjustment` (cumulative FOI) separately from progression scale | Timing of FOI within a birth cohort's lifetime |
| Plus prevalence survey | Duration `1/gamma` hence `beta` level | - |
| Plus child incidence or genotypic clustering | Recent-FOI share (recent vs remote contribution) | - |

In low-incidence settings the reactivation term from `L_slow0` dominates for
decades, so incidence carries little information on current `beta`; `beta`
estimates near `BETA_BOUNDS` limits should be treated as unidentified.

## 5. Inputs, risk factors and economics integration

### 5.1 Country-profile mapping (`engine/profiles/population_profile.py`)

| Profile field | Dynamic model input | Gap |
|---|---|---|
| `population_size` (`ProfileValue`) | `age_counts` total | Direct |
| `age_distribution` (`AgeBand`: `label`, `lower_age`, `upper_age`, `proportion`; line 241) | `age_counts` dict by single year | Needs band expansion rule (as `expand_five_year_age_distribution`) and open upper band handling |
| `ltbi_prevalence` (scalar `ProfileValue`) | `ltbi_ever`, `ltbi_recent` by age | Profile has one scalar; engine needs age-specific ever and recent. Requires a documented infection-history model (e.g. `docs/recent_remote_infection_history.md`) fitted to the scalar |
| `incidence` (`IncidenceData`, series of `IncidencePoint` `year`, `estimate`, `lower`, `upper`; lines 268, 291; measure `estimated_tb_disease_incidence`) | calibration target | Bounds unused today; must feed the likelihood (section 6) |
| `trend` (`TrendSpec`, line 364) | none | Descriptive only; must not set FOI trend |
| `risk_factors` (`RiskFactor`, line 396: `prevalence`, `effect_estimate`, `effect_measure` in `RR`/`HR`/`OR`, `affected_transition` default `"progression_to_disease"`, `prevalence_bounds`, `effect_bounds`) | six fixed `*_pct` keys with hard-coded `RR` | Needs a generic list of (prevalence, multiplier) pairs; measure type must be recorded and any conversion reviewed (OR not converted today, `engine_mapping.py` `EFFECT_INTERPRETATION`) |
| not in profile | notifications, CDR prior, TB mortality, background mortality, births | Must be added to the profile or supplied as reviewed priors |

### 5.2 Risk-factor integration requirements

Replace the hard-coded `RR` dict with profile risk factors; keep the product
form only as an explicit, documented assumption (independence); report which
factors apply to which transition; apply the same `engine_key` semantics as the
individual-based engine so that the crosswalk can be generated from code rather
than maintained as a CSV. Stratified compartments (risk group x latency) are
required if TPT targeting by risk group is to be modelled.

### 5.3 Economics integration

The individual-based engine produces an event ledger
(`engine/apy/event_ledger.py:17`, `EVENT_LEDGER_CONTRACT_VERSION =
"ltbi_screening_event_ledger_v3"`, `EVENT_DEFINITIONS` line 26) that is costed
by `engine/apy/event_ledger_economics.py` (`COMPONENT_SPECS` line 56,
`DALY_COMPONENTS` line 83). To reuse that economics layer the dynamic model
would need to produce, per arm and per model year (same `YEAR_BIN_CONVENTION`):
`population`, `screened`, `test_positive_total`, `false_positive`,
`tpt_started_total`, `tpt_started_false_positive`, `tpt_completed_total`,
`tpt_adr_stop_total`, `active_tb_cases`, `active_tb_cases_prevented`, plus
TB deaths and person-years with TB for YLL/YLD. Currently only
`active_tb_cases` (as `annual_incidence`) exists; screening of `S`, test
accuracy, completion and ADRs are not represented, and TB mortality is not
represented, so YLL cannot be computed. Transmission-averted cases
(secondary benefit) would be a new ledger event and must be labelled
separately from direct prevention, since the ledger scope statement
(`DIRECT_EFFECTS_SCOPE_STATEMENT`) covers direct effects only.

## 6. Using WHO incidence without equating incidence and infection trends

Causal chain, with the quantity observed at each step:

1. FOI `lambda(t)` (unobserved) acts on susceptibles -> new infections.
2. Infections enter recent latency; they stabilise to remote latency
   (`omega`). Stock of remote infection reflects decades of past FOI.
3. Progression: incidence
   `inc(t) = m * [sigma_f * Lf(t) + sigma_s * Ls(t)]`, where `Lf` is a
   convolution of recent FOI and `Ls` of long-past FOI; `m` is the risk
   multiplier. Disease incidence is therefore a lagged, smoothed,
   risk-weighted transform of the FOI history plus a slowly decaying
   reactivation reservoir.
4. Disease incidence (WHO estimate, with asymmetric bounds) is what the
   profile holds.
5. Detection: notifications `= CDR(t) * inc(t)` delayed by `1/gamma`;
   detection also feeds back to FOI through `I` duration.

A falling incidence can arise from falling FOI, a shrinking remote reservoir
(cohort replacement), falling risk-factor prevalence, or (for notifications)
changes in detection. The trends separate only with additional information:

- Cumulative FOI: age-specific LTBI prevalence (IGRA/TST surveys) or ARI surveys.
- Recent versus remote share: incidence in children under 5 or 15, genotypic
  clustering, or incidence by time since migration.
- Duration and CDR: prevalence surveys, notification-to-incidence ratio,
  patient and health-system delay studies.
- Risk-factor contribution: time series of risk-factor prevalence (profile
  `RiskFactor.prevalence` is a single value today).
- Natural history: `sigma_f`, `sigma_s`, `omega` from cohort literature, with
  priors rather than fixed constants.

Required constraint form: the FOI trajectory must be a free (smooth, prior-
regularised) function, never set from the incidence slope; WHO incidence enters
only through an observation model on `inc(t)`. A likelihood consistent with
`engine/who_incidence/trend.py` (lines 26-35) treats bounds as 2.5 and 97.5
per cent points of a split-normal on the log scale; because WHO estimates are
model-derived and serially correlated, the likelihood should include a
between-year correlation or common-shift term, and the bounds must not be
treated as a posterior on the model's own quantity.

## 7. Candidate calibration formulations

### (a) Joint transmission and detection fitted to incidence plus notifications

- Data: WHO incidence with bounds; national notifications by year (not in the
  profile today); ideally age-specific notifications.
- Identifiable: `beta(t)` level and `CDR(t)` / `gamma(t)` from the ratio of
  notifications to incidence, given fixed natural history.
- Assumptions: notifications are a known fraction of detections (no
  over-reporting, no reporting of prevalent or relapse cases outside the model);
  WHO incidence independent of notifications. The latter is false for many
  countries, where WHO incidence is itself derived from notifications and an
  assumed CDR, which makes the fit circular.
- Advantages: uses routinely available data; constrains detection, which
  drives intervention effects on transmission.
- Failure modes: requires a structural change (non-detected exit, e.g. TB
  death and self-cure) because current CDR is 1; double counting of the same
  information; trade-off between `beta` and `gamma` near equilibrium;
  notification artefacts (COVID-era 2020-2022) read as epidemiology.
- Cost: moderate (two time-varying functions; optimisation over about 2 x 20
  parameters, more for uncertainty).
- Streamlit suitability: medium; needs a notification data source and
  per-country review of WHO method.

### (b) Transmission fitted to incidence, detection externally specified

- Data: WHO incidence with bounds; reviewed prior for CDR or diagnostic delay;
  optional LTBI prevalence.
- Identifiable: smooth `beta(t)` (few parameters, e.g. level and slope or a
  low-order spline) and a scale for pre-window infection history, conditional
  on detection and natural-history priors.
- Assumptions: detection prior is correct; WHO bounds summarise uncertainty in
  incidence; natural-history constants transferable.
- Advantages: no circularity with notifications; close to the current code
  (replace RMSE by WHO-bound likelihood, replace Styblo-driven history with a
  parametric FOI history); fast; explicit about what is assumed.
- Failure modes: projections of detection interventions are only as good as the
  prior; incidence alone weakly identifies current `beta` in low-incidence
  settings (section 4), so posteriors may be prior-dominated; this must be
  shown, not hidden.
- Cost: low (a few parameters; sampling-importance-resampling or short MCMC
  feasible in seconds to minutes).
- Streamlit suitability: high.

### (c) Latent state-space observation model

- Data: WHO incidence with bounds, notifications, optional prevalence and LTBI
  surveys, each with its own observation model.
- Structure: latent process for FOI (e.g. log random walk) and disease
  generation; separate observation equations for WHO estimate, notifications
  (with CDR process) and surveys; parameters estimated by particle MCMC or
  Kalman-type approximations on a linearised model.
- Identifiable: whatever the data jointly support, with posterior
  quantification of non-identified directions.
- Assumptions: observation-error models, process-noise scale, priors.
- Advantages: principled separation of disease generation from observation;
  handles missing years and disruptions; honest uncertainty.
- Failure modes: sensitivity to process-noise and prior choices; convergence
  problems; hard to explain; slow.
- Cost: high (minutes to hours per country).
- Streamlit suitability: low for interactive use; feasible as offline
  precomputed country fits loaded as snapshots.

## 8. Recommended staged approach

Stage 0 (prerequisites, no formula changes to disease equations without review):
add engine unit tests for `simulate_dynamic` (mass conservation, zero-beta,
steady state); fix or guard the Euler instability (reject or sub-step when
`gamma*dt > 1`) and the `beta_series` length clamp; move calibration out of
`ui/dynamic_ui.py` into `engine/dynamic/`. Validation: tests pass;
conservation to 1e-9 relative over 50 years for all UI slider extremes.

Stage 1 (current, descriptive): WHO incidence and trend are shown and not used
by the dynamic model. Validation: existing tests; UI labels state that dynamic
projections use user-supplied incidence and placeholder bounds.

Stage 2 (formulation b): parametric FOI history (for example log-linear or
low-order spline in calendar time, independent of the incidence slope),
explicit priors on `delta_pre` (or CDR), `sigma_*`, `omega`, and risk
multipliers from the profile; split-normal log-scale likelihood from WHO bounds
with a correlation option; posterior by sampling-importance-resampling or MCMC;
projections from posterior draws replace `DUMMY_CI_PCT`. Validation:
- Synthetic recovery: simulate from known FOI path and detection, generate
  WHO-like bounds, refit; true parameters inside 95 per cent intervals in about
  95 per cent of replicates; report where posteriors equal priors.
- Posterior predictive checks: simulated incidence inside WHO bounds for the
  nominal share of years; no systematic residual trend.
- Prior sensitivity: projections under alternative detection priors reported.
- Trend-independence check: calibrated FOI trend is reported separately from
  the incidence APC and is not constrained to equal it.
- No-transmission limit: with `beta = 0`, no intervention and matched LTBI
  recent/remote seeding and risk multiplier, cumulative dynamic cases match the
  individual-based engine `cumulative_baseline_active_tb_cases` within Monte
  Carlo error, once natural-history parameters and mortality are aligned
  (currently they differ; the comparison must first document the mapping).
- Intervention equivalence: in the same limit, TPT cases averted match the
  individual-based engine for the same coverage, efficacy and completion.

Stage 3 (add notifications, formulation a elements): add a non-detected exit
(TB death, self-cure) and notification output; include notifications with a
separate observation model; flag countries where WHO incidence is derived from
notifications. Validation: synthetic recovery of CDR; posterior predictive
checks on notifications and incidence jointly; COVID-era sensitivity
(excluding 2020-2022).

Stage 4 (formulation c, offline): state-space fits run offline and shipped as
versioned, checksummed snapshots (same pattern as the WHO snapshot). Validation:
agreement with stage 2/3 posteriors where data are informative; coverage tests
on synthetic data with process noise; run-time and convergence diagnostics
(R-hat, effective sample size) recorded in the snapshot manifest.

Throughout: outputs are for planning and sequencing, not for denying care; any
display of calibrated FOI or `beta` must show its uncertainty and prior
dependence.

## 9. Open scientific decisions requiring review

1. Natural-history parameters: keep `sigma_fast = 0.01`, `sigma_slow = 0.001`,
   `omega = 0.2` or adopt literature-based priors consistent with the
   individual-based engine; whether progression should depend on age.
2. Reinfection and partial protection of latent and treated people.
3. Whether to add demography (births, ageing, background and TB mortality);
   required for horizons over about 10 years and for YLL.
4. Pre-window infection history: replacement for the Styblo rule
   (`K=5000.0`, `f_ssplus=0.6`) and how to fit it to scalar LTBI prevalence.
5. Detection prior: source of CDR or delay per country; whether `gamma`
   varies over calendar time.
6. WHO likelihood: treatment of asymmetric bounds, serial correlation, and
   countries whose estimates are notification-derived.
7. Handling of 2020-2022 disruption in both incidence and notifications.
8. Risk factors: independence assumption, OR-to-RR/HR conversion, whether
   factors act on transmission or detection as well as progression, and
   stratification for targeted TPT.
9. TPT representation: rate versus pulse, test accuracy and false positives,
   completion, and reconciliation of the two regimen efficacy tables
   (`dynamic_model.py:127` versus `engine/intervention.py:10`).
10. Mixing: homogeneous versus age-structured or setting-specific contacts,
    and imported infection for migrant-dominated epidemics.
11. How dynamic projections relate to the individual-based engine in the
    application: comparison only, or a combined estimate of direct plus
    transmission effects feeding economics.
12. Acceptance criteria for synthetic recovery and posterior predictive checks
    before any calibrated output is shown to users.
