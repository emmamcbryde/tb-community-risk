# Milestone 3C diagnostics and validation

## Git and recovery state

The requested worktree started clean on
`feature/country-subnational-tb-model-v1` at
`213cacda333c5dbc334e6982891aa357bc81ab46`, tracking the same origin branch,
ahead/behind0/0. Read-only remote inspection confirmed agreement before edits.
No surviving unfinished country work was overwritten. All protected worktrees,
refs and frozen bytes were fingerprinted before edits and compared again.

The first independent3B reproduction took55.836seconds:

```
python -m engine.country_tb.preview_progression --output abm/output/country_tb_m3c_reproduction_20261030 --mc-samples 4000
```

Its immutable export and every requested number matched; see the
[scientific review](milestone3c_scientific_review.md). No3B source correction
was needed. Additional numerical corrections were confined to the new3C
adapter: nonnegative Frechet-table interpolation near degenerate marginals,
stable tiny positive risk union, structured unresolved-kernel refusal, and
refinement of future-pressure time steps as well as infection-time support.

## Focused tests

```
python -m pytest tests/country_tb/test_incidence_inversion.py tests/country_tb/test_progression_forward.py tests/country_tb/test_risk_allocation.py -q --durations=10
```

**202 passed in118.86seconds** (119.857seconds including process startup).
After the new3C adapter's final stability/gate refinements, the unchanged132
retained3A/3B tests remain accepted; final3C coverage was rerun separately:

```
python -m pytest tests/country_tb/test_incidence_inversion.py -q --durations=5
```

**100 passed in8.03seconds** (8.968seconds including startup). All model/schema/
test bytes were fingerprinted before/after this run and remained unchanged.
This includes scalar compression, near-zero/near-one probabilities, rare-group
dominance, all-root rejection, false gate resealing,100+ retention, threshold
approval, mortality snapshot requirements and both numerical refinement grids.

```
python -m pytest tests/country_tb/test_baseline_sources.py -q --durations=5 -k "no_fallback or overall_known_scale or age_specific_one_scale or direct_profile_original or direct_recent_remote or subgroup_metadata or incomplete_subgroup or open_tail_direct_value or explicit_conditional_recency"
```

**12 passed,62 deselected in6.89seconds** (7.926seconds including startup).
Overall/age/assay anchors, exact direct profiles, recency allocation, subgroup
weighting and incomplete denominators remain functional. No public incidence
source bypass exists. No full unbounded country suite was run.

## Diagnostic results

```
python -m engine.country_tb.preview_incidence_inversion --output abm/output/country_tb_m3c_diagnostics_20261030
```

Outputs are immutable under `abm/output/`; nothing is exposed in production UI.
Final export completed in315.964seconds. Scientific hash:
`8eb35e535f9ecce7b88ef3fa002a4f586f9283cc7cfa46826bc0a011e21e816c`.
Offline validation passed all six contract families and22 scientific hash
occurrences(21 unique) in0.175seconds. Final tested/exported model/schema bytes
remain unchanged. The [machine validation receipt](validation_receipt_3c.json)
records exact commands, startup-inclusive runtimes, source checksums and export
checksum. Canonical artifacts are UTF-8: readers on Windows must specify UTF-8
rather than the machine's cp1252 text default.

Known fictional scale0.017 recovers as0.01699999999837018. The actual
non-monotonic operator example finds0.0446287103 and0.3218875825 for target800;
it adopts neither. Flat/saturated maps, impossible targets and unresolved
kernels are rejected. Boundary zero is explicit rather than clipped.

The real WHO Viet Nam target remains **validation only** and the country gate
is blocked. Its annual new/recurrent observation is not relabelled as the
instantaneous first-episode rate predicted by the scan. The source, uncertainty,
period and denominator remain attached. On the declared illustrative search
0<=alpha<=0.2, forward rates range0–1278.7903/100,000PY and cross182 between
0.018 and0.020. No empirical fitted country scale is adopted.

The **separately fictional closed population with Viet Nam-shaped historical
pressure** has a numerical root near0.01819928768, total TBI65.3314%, recent
8.4174%, remote-only56.9140%. These are mathematical counterfactuals, not
Viet Nam estimates or an approved WHO calibration. Full age/sex results and
derivative/condition diagnostics are in the executable export.

| Fictional scenario | Conditional alpha | Conditional total TBI |
| --- | ---: | ---: |
| Fixed illustrative remote hazard | 0.01819929 | 65.3314% |
| Constant final-segment tail | 0.01782389 | 64.6287% |
| Declining positive plateau | 0.01852166 | 65.9219% |
| Declining zero plateau | 0.01869405 | 66.2328% |
| No remote estimate | Unavailable | Unavailable |
| Half progression hazard | 0.04066451 | 88.4173% |
| Double progression hazard | 0.00847920 | 40.3335% |
| Union association−1 | 0.02182274 | 71.3379% |
| Union association+1 | 0.01592776 | 60.8140% |
| No risk factors in fictional scope | 0.02260947 | 72.4759% |
| Hypothetical target112 | 0.01061803 | 47.2810% |
| Hypothetical target264 | 0.02787403 | 78.8701% |
| Hypothetical20% imported/recurrent adjustment | 0.01417307 | 56.8486% |
| Bounded backward trend sensitivity | 0.01825414 | 64.3801% |
| Lower WHO shape, same synthetic target | 0.01882562 | 52.2095% |
| Upper WHO shape, same synthetic target | 0.01770493 | 72.6064% |

Bounds and scenario spreads are not TBI confidence intervals. Adult/child
clinical alternatives are not approved by fitting them. Instantaneous mortality
and future-pressure sensitivity is exactly zero; prospective1-year diagnostics
separately show depletion/mortality/future effects and refuse unadopted roots
when integration refinement fails. Failed prospective fits are retained with
their numerical differences, never forced through a relaxed tolerance.

All three prospective demonstrations remain `blocked_postfit_gate`: known-scale,
mortality0.08 and future-pressure0.05 scenarios have refinement differences
0.012016,0.035076,0.006083/100,000PY respectively, above their0.001323 allowance.
Their roots are unadopted; these are resolution failures, not accepted
prospective incidence calibrations. Finer prospective numerical work is an
explicit remaining dependency. The constant-hazard analytical prospective
kernel and existing3B deterministic/Monte Carlo checks pass independently.

The representative fictional multimorbidity example has baseline rate328.0513
and case fractions49.3900%,78.6097%,89.1604% in the top1%,5%,10% of expected
risk. Omitting transplant,HIV,silicosis,haemodialysis or head/neck cancer changes
the rate respectively to318.1682,295.0464,308.8180,295.8633,318.3463. Each RR's
population, case contribution and conditional refit are exported. A separate
maximal-association fixture puts all TBI/cases inside an unreviewed0.1% transplant
group: top1% case share100%, and the scientific gate refuses it. No arbitrary cap.

Frozen demographic denominators (female+male, no forced reconciliation to
published both-sex rounding):

| Source population | Complete population |0–14 unresolved progression population |100+ unresolved age population |
| --- | ---: | ---: | ---: |
| Australia ABS2024 |27,204,809 |4,794,286(17.6229%) |6,588(0.024216%) |
| Viet Nam WPP2024 |100,987,682 |23,453,385(23.2240%) |7,139(0.007069%) |

These are demographic coverage audits, not incidence estimates. Adult scientific
approval is also absent; the remaining population is not implicitly approved.
Australian national inversion additionally lacks migration/importation/prior-TB
structure. User TBI and explicit subgroup forward validation remain preferred.

## Integrity and limits

The [protection manifest](protected_fingerprints_3c.json) records unchanged
protected HEAD/status/diff/index/dirty-file bytes, all protected refs,128 frozen
files and55 retained2A–3B files. Explicit-TBI remains six modified/two untracked.
Frozen SA Health branch/peeled tag remain at
`03cc16e4a52a55e10019dab16bd4f599571c4b90`; annotated tag object remains
`e99d36725dfc3ff40edc0ed68127e0b6cd927da5`.

Compilation, local documentation links, six-family offline JSON-schema exports,
canonical hashes and `git diff --check` are final checks. Changed paths are
restricted to new country-engine files, new country tests/schema and country
review documentation. No production UI, APY, explicit-TBI, Starsim, MATLAB,
economics, DALY/QALY, intervention or frozen-release path is modified.

Actual unresolved-threshold approval is pending; candidate0%,0.1%,1% ceilings do
not constitute approval. Country science gates remain blocked independently of
that choice. Recommend external evidence/observation validation and specific
configuration approval, not production UI integration.
