# General community TB application - Milestone 2 status

Branch `feature/general-community-tb-model`. The frozen release (`sa-health-apy-he-v1.0.0`,
commit `03cc16e`) and `streamlit_app.py` are unchanged. Nothing has been merged, tagged or
deployed.

Run the general application with `streamlit run general_app.py`.

## Status legend

| Label | Meaning |
| --- | --- |
| Implemented and validated | Built, tested, and behaviour verified in the rendered interface |
| Implemented but provisional | Built and tested; scientific choices await review |
| Designed but not implemented | Specified in code interfaces or documents only |
| Blocked by external evidence | Needs data, licence confirmation or scientific input |
| Retained for backward compatibility | Kept unchanged so the frozen or earlier behaviour still works |

## Components

| Component | Status | Documentation |
| --- | --- | --- |
| Environment check and timeout root cause (duplicated calibration memoised) | Implemented and validated | `environment_and_reproducibility.md` |
| Complete WHO 2025-round incidence snapshot (217 areas, 2000-2024) and importer | Implemented and validated | `who_snapshot_schema.md`, `who_data_update_guide.md` |
| WHO licence and attribution | Implemented; one question open | `who_licence_and_attribution.md` |
| Population profile v2 (provenance, national population, bounds, v1 migration) | Implemented and validated | `population_profile_schema.md` |
| Log-linear and penalised smooth trends, propagated WHO-bound uncertainty, COVID handling, diagnostics | Implemented but provisional | `trend_estimation_methods.md` |
| State-space trend | Designed but not implemented | `trend_estimation_methods.md` |
| Country preview, explicit apply and conflict choices | Implemented and validated | this document |
| Local or subnational incidence upload | Implemented and validated | `local_incidence_upload_guide.md` |
| Risk-factor table (bounds, evidence fields, CSV) | Implemented and validated | `risk_factor_schema.md` |
| Effect-measure warnings and crosswalk | Implemented and validated | `risk_factor_schema.md` |
| Effect-measure conversion policy | Designed but not implemented | `risk_factor_schema.md` |
| Engine semantics (hazard multipliers, calibration) | Retained for backward compatibility | `risk_factor_schema.md` |
| Country risk-factor prevalence | Blocked by external evidence | `risk_factor_schema.md` |
| Incidence driving infection pressure or the dynamic model | Designed but not implemented | `dynamic_model_readiness_spec.md` |
| Hash-based result currency, cost-only economics, paired ICER plane | Implemented and validated | this document |
| Lazy provenance package | Implemented and validated | this document |
| Performance measurements | Implemented and validated | `performance_benchmark.md` |

## Workflow decisions (conservative defaults)

* **Selection does not apply.** Choosing a country only previews it. Applying it
  needs "Apply selected country incidence data", after a table of the fields that
  will change. User-supplied incidence or a subnational location needs an explicit
  keep-or-replace choice.
* **What country data change.** Only location identity, the incidence series and
  bounds, national population context, provenance and trend settings change. LTBI
  prevalence, age, test accuracy, the cascade, costs, DALYs and risk factors never
  change.
* **Trend defaults.** The default trend is log-linear over the last 10 years, with
  all WHO estimates included. Uncertainty uses 1,000 split-normal draws with seed
  20250101 and independent years.
* **Descriptive only.** Incidence data are descriptive: the configuration link
  records `incidenceUsedByEngine: false`, and results say they are not a
  country-specific estimate.
* **Result currency.** Results are current while the epidemiological
  configuration hash is unchanged. Trend and incidence edits do not invalidate
  them. Cost edits invalidate only the economics, and recalculating economics
  reuses the completed epidemiological run.
* **Accidental long runs.** A stochastic run estimated at more than 3 minutes
  needs a confirmation tick. The run button is disabled when results are already
  current, and identical inputs reuse cached results (last 3 per session).

## Rendered verification (Chrome, Streamlit 1.54, 2026-09-24/25)

| Check | Result |
| --- | --- |
| Fresh session, demonstration profile | Pass |
| Australia selected but not applied | Pass: summary, chart with legend and COVID shading, trend; profile unchanged |
| Australia applied | Pass: field-change table; population stays at the demonstration default; national population shown as context |
| South Africa applied over Australia | Pass: no conflict prompt (both WHO) |
| Country with unusual data | Pass: PRK shows "WHO publishes no incidence estimates", with no apply action |
| Local incidence CSV upload | Pass: preview, apply, "Local estimate" labels, recorded as user-supplied |
| 5-, 10- and 15-year log-linear trends; smooth; COVID years included and excluded | Pass (ZAF 10 years: -9.3%/yr, -13.4% to -4.3%; 5 years: -8.8%, -19.2% to +2.6%; smooth: -8.6% with short-series and boundary warnings) |
| Excluding 2020-2022 from a 5-year period | Pass: "Trend not estimated: only 2 year(s)" |
| User override after country application | Pass: "✎ User-defined" |
| Conflict with local data when applying a country | Pass: choice required; apply disabled until chosen |
| Restore demonstration defaults | Pass: country, local data and overrides cleared |
| No-risk-factor profile | Pass |
| Deterministic run through Results and Health economics | Pass: N/A intervals; not-country-specific notice |
| Stochastic configuration without launching | Pass: estimate shown; confirmation required; Run disabled |
| Duplicate-run prevention | Fixed during review (Run is disabled immediately after completion) |
| Provenance package preparation | Pass (lazy; download button appears only after preparation) |
| Frozen release from a clean tag extract | Pass: original wording unchanged |
| Narrow (phone) width | Not verified: the browser window could not be resized in this environment |

Defects found and fixed during the review:
* trend metrics truncated;
* "WHO" lower-cased;
* repeated odds-ratio notes;
* COVID shading unlabelled;
* the country preview lost on page navigation;
* the snapshot ID over-long in source labels;
* local-data labels mentioning WHO;
* excluded years not visually distinguished;
* a possible duplicate run straight after completion.

**Sessions.** Refreshing the browser or typing a page URL starts a new Streamlit
session, and unsaved inputs are lost. Navigate with the sidebar, and use
*Download population profile (JSON)* to keep a profile. Loading a saved profile
in the interface is recommended for Milestone 3.

## Scientific decisions still required

1. Linkage from estimated disease incidence to infection pressure: which quantity is
   calibrated, with what detection model (see `dynamic_model_readiness_spec.md`).
2. The trend method for any modelling use: log-linear or smooth, the default
   period, and COVID-era handling.
3. Whether WHO bounds should be treated as independent between years, and whether
   they are a suitable likelihood.
4. Effect-measure conversion: whether to convert ORs and RRs, and with what
   baseline risk; how to treat mutually unadjusted estimates.
5. Country-specific LTBI prevalence, age distribution, test accuracy, cascade and
   costs, and their evidence-review process.
6. A licensed source of country risk-factor prevalence.
7. Confirmation from WHO on redistributing the transformed snapshot.

## Recommended Milestone 3

* **Stage 0 of the dynamic readiness plan:** tests for the dynamic engine, and
  fixes for the Euler step-size instability and `beta_series` indexing found in
  the audit, with numerical equivalence documented.
* **Stage 2 calibration prototype:** formulation (b), transmission fitted to WHO
  incidence with detection specified externally, run on synthetic-recovery tests
  before any country use.
* **Evidence register** for country LTBI prevalence and risk-factor prevalence,
  with licence tracking.
* **Profile save and load** in the interface, with snapshot-ID verification on
  reload.
