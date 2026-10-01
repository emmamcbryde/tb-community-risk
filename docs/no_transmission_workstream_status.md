# No-transmission workstream status

This document is replaced at the end of each milestone. It does not keep a history; `git log` does.

| Field | Value |
| --- | --- |
| Date | 2026-10-01 |
| Branch | `feature/generic-no-transmission-model` |
| Milestone | Pure mathematics for exogenous background infection pressure (isolated, not connected) |
| Commit | The commit that adds this file (see `git log -1 -- docs/no_transmission_workstream_status.md`) |
| Starting point | `9003146` (previous milestone implementation `ab9ceca`) |
| Merged, tagged or deployed | No |

## Objective

Implement and validate the pure mathematics for optional exogenous background TB infection
pressure. Keep it entirely separate from the epidemiological engine, event ledger, economics and
ordinary Streamlit interface.

## Scientific decisions applied

These are adopted for the isolated mathematics only. Details are in section 3a of
`catalytic_infection_pressure_spec.md`.

1. **Baseline LTBI prevalence remains a separate input.** The prospective hazard applies after time
   zero and does not generate or recalibrate baseline prevalence.
2. **Hazard, not probability.** λ is in infections per susceptible person-year, and
   P = 1 - exp(-H), computed as `-expm1(-H)`.
3. **Modes:** `none`, `constant` and piecewise-constant annual `time_series`. There are no
   qualitative presets.
4. **Time series:**
   * year y applies on [y, y + 1);
   * an interval needs the years from floor(t0) to ceil(t1) - 1;
   * gaps and uncovered years are errors. Nothing is carried forward, backfilled, interpolated or
     set to zero.
5. **Hazard range:** finite and non-negative, with no mathematical maximum. A hazard above
   1 per person-year produces a non-blocking plausibility warning. That threshold is provisional
   and not evidence-based.
   * **This changed the schema:** `background_exposure_v1` previously *rejected* such values. That
     rejection contradicted this decision, so it is now a warning (`plausibility_warnings()`).
   * The change is more permissive only. Every previously valid payload, and its hash, is
     unchanged.
6. **First infection only**, for people susceptible at t0. Nothing else is implemented.
7. **Mode `none`:** H = 0, P = 0, no event and no random draw.

**Age-specific hazards are specified but not yet executable.** `schedule_from_config` raises
`AgeSpecificHazardNotExecutableError`; the bands are never silently ignored.

## Implementation completed

* `engine/profiles/background_exposure_hazard.py` (new; pure; no random-number calls):
  * `HazardSchedule` (`none`, `constant`, `annual_series`) and `schedule_from_config`;
  * `hazard_at` and `cumulative_hazard`;
  * `infection_probability` and `probability_from_cumulative_hazard`;
  * `exponential_threshold` and `first_infection_time`, with a typed result: `InfectionOccurs` or
    `NoInfectionInInterval`, which carries the residual threshold for chaining intervals;
  * the vectorised `first_infection_times` and `infection_probabilities`;
  * `requires_random_draw` and `plausibility_warnings`.
* **Inversion convention.**
  * The caller supplies U ∈ (0, 1), open, or E = -log U (finite, E ≥ 0).
  * Infection occurs in [t0, t1) if and only if E < H(t0, t1), at τ = inf{s : H(t0, s) > E}.
* `engine/profiles/background_exposure.py`:
  * the hard ceiling above 1 is replaced by `plausibility_warnings()`;
  * the docstring is updated.
* Docs:
  * `catalytic_infection_pressure_spec.md` has a new section 3a (API, conventions and validation
    results), the adopted time-series rules, the hazard-versus-probability and warning text, and
    the age-specific status;
  * it states: "Background infection pressure is an externally supplied exposure hazard. It is not
    inferred automatically from reported active-TB incidence."

Nothing is connected to cohort generation, infection or disease outcomes, the event ledger,
economics, Streamlit, MATLAB, frozen artifacts or Starsim.

## Validation performed

* **New file `tests/test_background_exposure_hazard.py`** (29 tests). Its coverage:
  * mode `none`, including a zero-draw proof;
  * constant hazards and piecewise-constant series;
  * numerical behaviour and input validation;
  * statistical checks;
  * isolation.
* **Statistical checks:**
  * fixed seed, n = 200,000 caller-supplied variates per scenario, and a tolerance of 4 standard
    errors fixed in advance;
  * five scenarios, each checking the proportion infected and the infection-time CDF at three
    points;
  * worst |z| = 1.39.
* **Mode `none` draws nothing.**
  * `requires_random_draw` is False.
  * The inversion functions reject a variate for `none`.
  * A reference caller with a counting generator makes 0 calls, and the bit-generator state is
    unchanged.
  * With numpy and Python RNG entry points patched to raise, every function runs, and the global
    RNG states are unchanged.
  * The source contains no RNG usage.
* **Not wired in.**
  * In a subprocess, importing the runner, simulation, expected-value, event-ledger and economics
    modules, the frozen reference, the engine mapping and the general-application modules does not
    load the new module.
  * A source scan finds no import of it in `engine/`, `app/`, `ui/`, `pages/`, `general_pages/`,
    `adapters/` or the root scripts.
* **Mode `none` unchanged:** its configuration dict and hash (`39cac6a7…`) are identical to those
  computed with the schema at `9003146`.
* **Focused and guardrail tests**, bounded by `timeout -k 30 900`: 66 passed, 561 subtests, 19 s.
  The files are:
  * the hazard and schema tests;
  * `test_no_transmission_scope_docs.py` and `test_general_scope_terminology.py`;
  * `test_frozen_release_integrity.py` and `test_frozen_reference_loader.py`;
  * `test_apy_source_formatting.py`.
* **Broader regression**, bounded by `timeout -k 30 1500`: `test_general_*.py`,
  `test_results_page_presentation.py` and `test_legacy_static_ui.py`. Result: 115 passed,
  90 subtests, 102 s.
* **Frozen outputs:** the frozen-release integrity checks pass, covering the reference-artifact
  SHA-256 values, the release tag and branch, the frozen headline, and the release-file diff.
* **Not run:**
  * MATLAB and the full suite;
  * `test_sa_health_reference_package.py::...test_rendered_health_economics_widgets_recalculate_without_changing_health`.
    That test is the documented **pre-existing** 90 s AppTest timeout: it reproduced on a clean
    `7920889` worktree (1 failed in 379.87 s). It was deliberately not rerun, and remains
    classified as pre-existing.

## Unresolved issues

* **Scientific decisions** (spec section 4):
  * reinfection (of remotely infected people, after preventive treatment, after active-TB
    treatment);
  * partial immunity;
  * preventive-treatment semantics (clears infection or reduces progression);
  * progression after an incident infection (recent state, duration);
  * mortality and ageing, including competing risks;
  * repeated screening;
  * historical catalytic calibration of baseline prevalence;
  * deriving scenarios from external epidemiological data (never automatically from WHO incidence).
* **Execution of age-specific hazards:**
  * age representation at t0;
  * birthday and band-crossing splits;
  * interaction with calendar years;
  * acquisition versus progression modifiers.
* Whether any explicit, recorded fill rule for missing years should ever be offered.
* The provisional plausibility-warning threshold.
* Constant λ = 0 versus `none` when wired: draw-stream design (decision 17).
* Anchoring engine time 0 to calendar time.
* The pre-existing SA Health Economics AppTest timeout, which is not fixed here.
* The applicability framework and citations still need scientific review. The future of
  `engine/dynamic/` is also undecided.

## Recommended next scientific decision

**Reinfection and the preventive-treatment mechanism**, decided together (decisions 5, 6 and 9):

* whether completed preventive treatment clears infection and returns the person to the
  susceptible pool;
* the relative hazard of reinfection disease for previously infected people.

Under prospective exposure these choices determine who is "susceptible at t0" for the
first-infection mathematics, and they change the direction and size of the benefit of preventive
treatment. They should be settled before any wiring.
