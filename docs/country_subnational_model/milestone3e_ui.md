# Milestone 3E isolated scientific-review interface

The software review interface is implemented; scientific central parameters and
production incidence inference are **not approved**. This consumes retained 2A–3D
contracts rather than changing their semantics. The 3D package remains internal
evidence synthesis/model checking, not external validation.

## Architecture and launch

Entry: `country_model_review_ui/app.py`. Pure adapters: `state.py`, `model.py`,
`forward.py`, `gates.py`. Display-only formatting/export boundary:
`presentation.py`. Seven explicit views reside in `views/`:

1. Set up population
2. Baseline TBI
3. Risk factors
4. Natural history and scientific gates
5. Run diagnostic
6. Results and validation
7. Evidence and unresolved decisions

There is no adjacent `pages/` directory. Streamlit discovers only this entry
point; the seven workflow views are an explicit radio allow-list. No root page,
deployed entry point, frozen SA Health, APY, explicit-TBI, generic no-transmission,
Starsim or legacy economics/results module is imported.

Exact local command, from repository root:

```powershell
python -m streamlit run country_model_review_ui/app.py --server.address 127.0.0.1 --server.enableCORS true --server.enableXsrfProtection true --browser.gatherUsageStats false
```

Default port 8501; `--server.port 8765` selects the review's development port.
No deployment, merge or tag is involved. Global configuration and shared Python
packages are unchanged; the pinned 1.39 test copy lives in ignored workspace
output, separate from the installed 1.32 environment.

## Population and source semantics

WHO 2025 country/territory catalogue and incidence, ABS June 2024 ERP and UN WPP 2024
medium demography use the existing offline frozen loaders. The UI retains
snapshot IDs, input/scientific hashes, estimate/projection flags, source vintage,
missingness and published-both-sex reconciliation. Sex populations are retained
separately; combined display is derived only. Source rounding discrepancies are
not forced to zero.

2024 is the only supported baseline, not a disguised latest-data lookup. WPP's
complete snapshot validation is explicitly loaded and cached; ordinary edits
do not repeatedly deserialize it. PRK/missing incidence is visible and blocked,
with no positive replacement or inherited calibration.

Custom/subnational populations require their own geography namespace and counts,
or explicit selection of national age/sex proportions as a diagnostic proxy.
National incidence remains context only. Explicit source scaling preserves raw
source counts independently and records any 80-digit arithmetic residual.

Subgroups require unique stable IDs, positive reconciling denominators and an
explicit disjoint partition. They may have separate age/sex CSVs or an explicitly
selected parent-proportion policy. Their age/sex sums must reconcile to the parent;
roundoff residuals are retained, never silently removed. Add/duplicate/remove/edit
controls act on review session data only. Duplicates receive new IDs and a missing
denominator requiring a decision, avoiding apparent double counting.

Current residence, birth, previous residence, exposure geography, arrival date/year
and residence duration are distinct metadata. Missing values remain explicit.
These do not reconstruct migration or automatically choose historical exposure.
The selected pressure proxy applies only under the explicit proxy selection;
different exposure histories need reviewed source choices/direct profiles, not
an assumption that residence equals acquisition geography.

## User-TBI route

The default functional fork is observed/user-specified TBI. An overall input
alone produces a direct total profile with unavailable recency. Explicitly
selecting catalytic age/recency derivation fits the retained one common scale
per subgroup to the covered population. Smoothed WHO incidence is an infection-
pressure **shape proxy**, not an estimated force of infection. Earliest normalized
shape is held constant before the observed series; reconstruction covers the
complete finite-age lifetime. No post 2024 WHO extrapolation is used.

Direct age or age/sex CSVs use exact single-year or explicit band allocation
policies. Overlap fails and gaps remain unresolved. Supplied values and notation
remain in the source CSV/receipts; model comparisons are separate. Direct recent,
remote-only and never-infected states must reconcile exactly. Prior remote plus
recent is a subset of recent and never another mutually exclusive state.

Total-only overrides never imply recency. The optional
`model_conditional_recency_allocation_v1` explicitly allocates the supplied total
using model-conditional composition and labels it hybrid/model-derived. Direct
states retain their values. A separate opt-in permits conditional model infection-
time distributions for forward diagnosis; a prevalence or mean does not invent
an infection-time distribution.

Assay positivity is distinct from latent infection. Overall positivity requires
fixed assay/product, threshold, Se, Sp and provenance. The observation equation is
`positive=Se*TBI+(1-Sp)*(1-TBI)`. Performance is not fitted withalpha. Uniform
performance across included ages is an explicit diagnostic assumption; pediatric,
BCG and immunosuppression transport is not approved. Direct age/state CSVs require
latent TBI, not mixed/untransformed positivity. Infeasible values fail without
clipping or fabricated counts.

Subgroup outputs combine expected counts, not unweighted prevalence means.
Missing or incidence-selected subgroup profiles retain their population and
block complete mixture results. A missing subgroup value never inherits a
global prevalence silently. Covered native mixtures and unresolved parent rows
remain separately represented in the UI envelope.

## Risk and association

The retained eleven-factor catalogue preserves definitions and legacy effects.
Inputs distinguish measured zero, missing, excluded, not applicable, not reviewed
and definition mismatch. Direct age/sex and explicit band/uniform policies are
supported. Child exposure applicability is separate from adult effect evidence.
Current smoking is not heavy smoking; BMI is not percent ideal body weight.
Definition mismatches remain unavailable, not silently mapped.

Central combination is `highest_applicable_reviewed_rr_v1`. There is no additive,
multiplicative or extra multimorbidity central effect. Legacy effects require an
explicit sensitivity selection and remain unreviewed as conditional hazards.
RR 74 displays the existing 94.31% one-year warning at selection, preview and affected
results. No arbitrary cap is introduced. An RR 1 reference must be explicitly
selected and cannot coexist with unacknowledged entered risk records.

Within age/sex/subgroup cells, `L=max(0,p+q-1)`, `U=min(p,q)`, `J0=p*q` and the
normalized controlc interpolates fromJ 0 toL orU. Achieved phi is separately
reported and unavailable for degenerate margins. The control is not Pearson
correlation. Factor marginals, union, TBI, recency and population remain unchanged.
CSV overrides use retained precedence; ambiguity and infeasible inputs fail.

## Forward diagnostics and units

The fast default is deterministic **instantaneous first-episode case flow**.
Retained 2C conditional last-event CDF masses are integrated against 3B piecewise
hazards. Uniform single-year age integration remains explicit; there are no age
midpoints or mean-clock substitutions. Native compressed 3A state/RR counts supply
the highest multiplier. This adapter is independently checked against 3B GL 8
conditional supports. It does not recalibrate TBI to incidence.

Adult surrogate, remote-tail extension, baseline-active-TB exclusion and mortality
must be explicitly selected for sensitivity. Children and unknown 100+ ages remain
unresolved. No extrapolation retains positive unresolved remote exposure mass;
it does not set its hazard to zero. **Any** positive unresolved population makes
the full-population rate unavailable, including very small positive counts.

Optional one-year prospective output calls unchanged 3B cohort/forward functions,
using 1-year maximum infection-time quadrature spacing and.25-year time steps.
A 50000-support-node estimate blocks excessive expansion before construction.
First disease and competing death deplete the risk set. No future infection is
assumed in this UI's baseline-cohort mode. Imported disease, prior TB, recurrence,
historical depletion, migration and endogenous transmission remain unrepresented.

ABS mortality is a validation sensitivity only and cannot apply to another
country. Custom Australian use requires explicit proxy selection. Future survival
does not change instantaneous baseline hazard; it affects prospective outcomes.
Progression/tail comparisons and baseline mortality invariance are labelled
scenarios, not confidence intervals. Long-horizon country projections and automatic
prospective multi-scenario expansion are not exposed.

WHO annual new/recurrent country incidence is displayed separately as **Observed
context**, not a matched first-episode target. Different populations, time measures
and denominators remain visible. Disagreement never modifies supplied TBI. The
fictional Viet Nam 65.33% counterfactual is not used by this interface.

## Diagnostic inversion bundle

Public mode stays `unavailable_pending_progression_and_risk_forward_operator`.
The experimental route is disabled by default. A separate uploaded JSON bundle
may describe an explicitly reviewed diagnostic configuration; no example or
legacy calibration is supplied automatically. Its maximum size is 2 MB.

Required top-level fields:

- `review_binding_hash`: technical hash of the current material inputs, excluding
  the bundle itself to avoid a circular reference;
- `component_reviews`: all 15 items shown in the UI, each with exactly resolved/
  reviewed status plus nonempty source ID, reviewer and rationale;
- `mapping_inputs`: complete JSON constructor inputs for retained 3C
  `IncidenceMapping`, including immutable history, operator and reviewed within-cell
  age support; maximum 300 strata;
- `observation`: retained `country_tb.active_tb_incidence_observation.v1`, matching
  population, episode, acquisition, period, denominator and instantaneous measure;
- `unresolved_threshold_approval`, `scientific_review`, `concentration_policy`:
  retained 3C contracts, not checkbox approvals;
- `assay_decision`: reviewed decision and source;
- `numerical_protocol`: reviewed source, policy
  `retained_3C_bounded_grid_all_roots_and_refinement_v1` and an explicit positive
  `absolute_refinement_tolerance` no looser than native 3C.

The UI checks source population/age/sex/subgroup counts, risk marginals, definitions
and requested association against selected inputs. It refuses invented birthdays,
scope boundaries cutting source cells, missing prevalences, stale bindings and
incompatible targets. An upload cannot relabel the four-point surrogate as an
approved country default or authorize production inference.

First, **Check scientific and numerical gates** runs native scientific pre-fit
checks and bounded attainable-range/continuity/derivative/all-root/refinement
checks. This adopts no TBI profile. Only a current positive clearance enables
**Run guarded incidence diagnostic**. The backend rechecks it independently;
session flags/direct navigation cannot bypass it. Nonunique, infeasible, flat or
failed-refinement mappings remain blocked. Output status is diagnostic-only,
not production approval. A fully specified synthetic uploaded fixture verifies
known-scale recovery without adopting empirical country assumptions.

## Upload schemas, precision and provenance

Examples: [population](../../country_model_review_ui/examples/population.csv),
[age/sex TBI](../../country_model_review_ui/examples/tbi_age_sex.csv),
[direct states](../../country_model_review_ui/examples/tbi_states.csv),
[risk](../../country_model_review_ui/examples/risk.csv),
[association](../../country_model_review_ui/examples/association.csv).
UTF-8 CSV, exact documented headers, at most 5000 rows/2 MB. Probabilities use 0–1;
percentage-labelled UI inputs use 0–100. Open upper age is blank. Numerics are
preserved as source strings; percentage conversion does not truncate to the
default 28-digit Decimal context.

| Upload | Required columns |
| --- | --- |
|Population|age_lower, age_upper_exclusive, sex, population, source, review_status, geography_id, baseline_year |
|Subgroup population|subgroup_id plus all population columns |
|TBI|record_id, subgroup_id, age_lower, age_upper_exclusive, sex, any_tbi, recent, remote_only, never_infected, prior_remote_plus_recent, allocation_policy, source, review_status, geography_id, baseline_year |
|Risk|record_id, subgroup_id, factor_id, age_lower, age_upper_exclusive, sex, state, prevalence, allocation_policy, supplied_exposure_definition, source, review_status, reason |
|Association|control_id, scope, subgroup_id, sex, age_lower, age_upper_exclusive, control, source, review_status |

Every run retains geography, population definition, snapshot IDs, baseline source,
original TBI/risk inputs, association, natural/pediatric/tail/reset/mortality/assay
decisions, migration/importation/recurrence metadata, numerical configuration,
native scientific hashes, application version and contract versions.9 newUI
envelope families are in [review_ui.schema.json](../../schemas/country_tb/v1/review_ui.schema.json).
Retained contracts are unchanged.

Material configuration hashes invalidate prior results and prepared downloads.
Display-only sex aggregation does not change scientific inputs. Stale results
are retained for audit but hidden as current outputs. Exports are prepared only
on request; JSON contains complete source/native contracts, gates, unresolved
population and full precision. CSV includes profile, risk/association and run
provenance fields. Rounding is display-only. No Excel workbook is generated.

## Reviewer questions

Review source/population compatibility, within-band allocation, subgroup partition
and migration metadata; assay interpretation; age/recency conditioning; factor
definitions and evidence; extreme multipliers; pediatric/tail/mortality gaps;
episode/period compatibility; unresolved-population/refinement thresholds; and
which empirical datasets should be acquired next. UI usability/testing does not
resolve these scientific decisions or substitute for external review.
