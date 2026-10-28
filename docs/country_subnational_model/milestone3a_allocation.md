# Milestone 3A: baseline population risk allocation

This pure, offline layer allocates expected baseline population counts. It is for
planning and sequencing, not denying care. It does not calculate active TB,
interventions, costs or health outcomes. No existing application imports it.

## Retained baseline semantics

The 2C/2D implementation and contract files are retained byte for byte. Infection
events remain a repeated-event non-homogeneous Poisson process. The most recent
event determines recency; a recent reinfection overrides remote-only status.
`prior_remote_plus_recent` is an auditable subset of `recent`, never an added
mutually exclusive state. This allocator neither creates events nor recalibrates
the one-scale baseline. Most-recent-event summaries remain in the embedded source.

Source selection remains explicit. User anchors, direct profiles and
population-weighted subgroup mixtures retain their existing behavior. Incidence
inference still returns
`unavailable_pending_progression_and_risk_forward_operator`. PRK remains blocked;
the 2C Australia/Viet Nam scale-0.005/year examples remain diagnostics.

## Candidate catalogue and definition compatibility

The following values transcribe [Menzies et al., 2008, Table 2, page 500](https://www.tstin3d.com/thinking3D_2008.pdf).
See the retained [evidence audit](tstin3d_evidence.md). They are candidate
progression relative risks, not approved contemporary causal hazard ratios.
Each record retains a definition ID, source/reference chain, vintage, outcome,
age/sex applicability, review status, compatibility notes and scientific hash.

| Stable ID | Catalogue exposure definition | Candidate RR |
| --- | --- | ---: |
| cxr_granuloma | Granulomatous abnormality on chest radiography | 2 |
| heavy_smoking | Cigarette smoking more than one pack/day | 3 |
| underweight | Body weight less than 90% of ideal body weight | 3 |
| diabetes | Diabetes mellitus; record type, glycaemic control and duration when available | 3.6 |
| tnf_inhibitor | TNF inhibitor treatment; record agent, dose and duration; legacy estimate originated with infliximab | 4 |
| cxr_fibronodular | Fibronodular chest radiographic disease consistent with past TB | 6 |
| head_neck_carcinoma | Carcinoma of the head or neck; record site and treatment status | 16 |
| haemodialysis | Renal failure requiring maintenance haemodialysis | 25.3 |
| silicosis | Diagnosed pulmonary silicosis | 30 |
| hiv | HIV infection; record stage, CD4 and antiretroviral treatment when known | 50 |
| transplant_immunosuppression | Transplantation requiring immunosuppressive treatment | 74 |

Audited transcription and contemporary effect approval are separate. The default
`reviewed_effects_only_v1` leaves an unreviewed positive-exposure RR unresolved.
An explicit `diagnostic_provisional_adult_effects_v1` policy permits these numbers
only in an explicitly diagnostic allocation. Reviewed age/sex-specific effects
can be supplied with evidence, vintage, review and provenance.

Exact catalogue definitions are required. Current smoking, BMI below 18.5, all
CKD, all cancer and all immunosuppression are not automatic substitutes. A mismatch
is unavailable. An explicit reviewed definition mapping retains both definitions,
rationale and provenance and emits a warning; it does not approve the effect.
HIV treatment era, severity and prior TB treatment require later review.

Recency and close contact are not ordinary factor IDs. Fixed 0.1% progression,
lifetime-to-age-80 calculations, diagnostic TST interpretation, BCG/NTM
adjustments, Styblo conversion and raw additive TSTin3D scoring are excluded.

## Prevalence, exact age and sex

The caller declares the analysis factor scope. Inputs distinguish `value`,
`missing`, `not_applicable`, `excluded`, `not_reviewed`, and
`unavailable_due_to_definition_mismatch`. A measured zero is `value: 0`.
Missing or unreviewed values are not zero-filled. Explicit exclusion and
non-applicability need a reason. The union concerns the declared scope after
explicit exclusions, not every possible progression factor.

Each record requires an allocation policy:

| Policy | Meaning |
| --- | --- |
| `direct_age_sex_values_v1` | Exact single-year values; no inferred open-tail age |
| `piecewise_constant_within_supplied_band_v1` | Explicit constant prevalence in a supplied age/sex band |
| `uniform_within_eligible_strata_v1` | Explicit overall 15+ prevalence, constant across eligible strata |
| `reviewed_age_sex_pattern_scaled_to_total_v1` | Reviewed relative weights scaled to exactly the supplied overall 15+ marginal |

Both-sex records apply the same exposure prevalence to female and male cells.
Separate records preserve sex-specific exposure patterns. Overlapping records,
out-of-range probabilities and boundaries cutting source demographic cells fail.
A reviewed uploaded table uses the same records; uploading does not imply an
allocation rule. Pattern scaling uses population weights, requires complete
eligible-cell coverage, and fails if a resulting prevalence exceeds one. It
never caps a value or converts a percentage automatically.

Single-year cells remain exact. Clinical groups are 0–4, 5–14 and 15+; demographic
display bands remain 0–4, 5–14, 15–44, 45–64 and 65+. Cells crossing a clinical
boundary remain unresolved. Display totals report crossing cells separately.
Sex changes exposure allocation only; no diagnostic or treatment rule is added.

Below 15, absent adult-factor inputs are
`not_applicable_pending_pediatric_review`. They are not numerical zero. Even if
child exposure is directly supplied, an adult RR is not applied automatically.
Reviewed pediatric effects can use the same explicit effect contract.

The 100+ denominator is always retained. Unresolved baseline TBI remains
unresolved. A reviewed direct open-cell TBI value, or a reviewed within-cell age
distribution upstream, can resolve it. An explicit band prevalence can cover
100+ without inventing an exact age when the adult effect is constant. No
arbitrary terminal age or midpoint is used.

## Overlap and highest RR

The central overlap policy is
`conditionally_independent_risk_factors_within_age_sex_stratum_v1`, before TBI
coupling. It gives `q = 1 - product(1-p_j)`, expected factor count `sum(p_j)`,
and an exact factor-count distribution. Multimorbidity is probability of two or
more factors. A separate `reviewed_joint_risk_factor_distribution_v1` accepts a
complete reviewed pattern table, checks its marginals exactly and invents no
pairwise correlations.

The central progression policy is `highest_applicable_reviewed_rr_v1`:
`max(1, RR_j for present factors)`. There is no extra multimorbidity effect.
The combination API returns all present factors, every tied highest contributor,
and effect provenance. Compressed output reports possible responsible IDs and
nonexclusive highest-contributor probabilities, counting all ties.

Multiplicative RR and additive-excess RR are explicit calculation-only
sensitivities. They cannot enter population allocation. Cap behavior is explicit:
no cap, or an explicitly reviewed sensitivity cap. The central maximum is not
capped. These are relative-risk calculations, not progression rates.

Independent allocation uses dynamic programming on highest RR and factor count,
at most 144 states for eleven factors, without enumerating membership patterns.
Reviewed joint data can explicitly supply at most 2^11 = 2,048 patterns.

## Residual TBI association

Association acts within each age/sex/subgroup cell. With TBI prevalence `p`,
union prevalence `q`, `L=max(0,p+q-1)`, `U=min(p,q)`, `J0=p*q`:

```
J = J0 + c*(J0-L)  when -1 <= c < 0
J = J0             when c = 0
J = J0 + c*(U-J0)  when  0 < c <= 1
```

The required input `c` is a normalized feasible residual association control,
not Pearson correlation. Achieved phi is reported separately as
`(J-p*q)/sqrt(p*(1-p)*q*(1-q))`; it is unavailable for degenerate marginals.
Precedence is exact-age/sex, exact-age/both, age-band/sex, age-band/both,
sex-only, then explicit global control. Ambiguous controls fail and selected
control provenance is retained.

For every risk-positive pattern, `P(TBI|R+) = J/q`; for risk-negative people,
`P(TBI|R-) = (p-J)/(1-q)`. Empty groups are handled without division by zero.
The conditional pattern distribution within R+ stays the same in both TBI
groups. Each factor marginal, union, TBI marginal and population are preserved
exactly. Factor-specific joint probabilities change mechanically with the union
control; no preferential diabetes/smoking association is inferred. Independence
before coupling does not imply independence conditional on TBI afterwards.

Within TBI, recent, remote-only and prior-remote-plus-recent proportions are
unchanged. The subset is never added again. Source total-only TBI with recency
unavailable stays unavailable. A source decimal recency closure residual, if
present in an inherited hybrid profile, is retained and reported; source values
are not repaired and complete stratification is unavailable.

## Outputs, mixtures and limitations

Outputs retain the complete baseline source, profiles, policies, review evidence,
hashes, exact expected counts, factor marginals, count categories, effective-RR
distribution and state-specific factor present/absent counts. Unknown components
have explicit coverage denominators. A known full TBI marginal may be reported
when risk allocation is incomplete, but `complete_result` is null. Full-population
TBI is null when baseline TBI coverage is incomplete.

Subgroups are allocated separately. Mixtures combine counts afterwards, never
averaged RRs, prevalence inputs or controls. Aggregate phi separates residual
within-cell covariance from structural age/sex/subgroup covariance. Incomplete
subgroups stay visible. No migration history, selective mortality or historical
exposure location is inferred.

Exact rational probabilities/counts retain numerator/denominator pairs; 80-digit
decimal displays are convenience representations. Canonical serialization and
SHA-256 hashes exclude source receipt ordering, retain scientific assumptions,
and are independent of row order and caller decimal context. Source records
preserve all supplied decimal precision.

Stochastic generation is deferred. The embedded future generator contract
requires an explicit RNG/seed, declared factor-pattern sampling, union-conditional
TBI, source-conditional recency, auditable subset membership, all highest-RR ties,
and predeclared Monte Carlo tolerances. Infection times must come from the original
reviewed distribution, not band midpoints. Parameter, source-data, scenario and
stochastic uncertainty remain distinct; no confidence interval is invented.

See [contracts](milestone3a_contracts.md), [diagnostics](milestone3a_examples.md),
[validation](milestone3a_validation.md) and [decisions](decision_register_3a.md).
