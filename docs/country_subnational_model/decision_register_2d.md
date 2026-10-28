# Milestone 2D decision register

The model supports planning and sequencing, not denying care. This register
supersedes earlier source-choice plans while preserving the 2C catalytic
mathematics and its explicit scientific limitations.

## Approved and implemented

| Decision | Implementation |
| --- | --- |
| Two explicit primary baseline-TBI modes | Exactly one versioned source mode; unknown/mixed modes fail |
| User-specified prevalence supported | Overall and age/sex calibration anchors plus source-preserving direct profiles |
| Subgroup-specific inputs | Each disjoint subgroup can have its own source, reviewed history and common within-subgroup scale |
| No automatic fallback | Blocked, unreviewed or incompatible sources retain unresolved population |
| No active-incidence inference before approved forward operator | Structured unavailable result and null TBI estimate, regardless of supplied incidence/history |
| Total-only direct profiles do not automatically imply recency | Unavailable by default; explicit hybrid policy required for conditional allocation |
| Direct recent/remote input retained | All three exclusive probabilities must reconcile; no recalculation from a comparison model |
| Population-weighted subgroup combination | Counts and positive reconciled denominators determine weights |
| Current residence differs from possible historical exposure | Distinct nullable metadata; reviewed proxy-application policy required for cross-geography history |
| Non-overlap must be supported | Reviewed disjoint partition plus complete membership checks when person IDs exist |
| Open tails stay visible | No exact age/terminal-age invention; direct values or reviewed age weights can resolve 100+ |
| Interface specification only | Exact two source labels and subgroup controls documented; incidence remains non-operational |

The 2C repeated Poisson event process is retained: previous infection does not
prevent reinfection, the most recent event determines recency, recent
reinfection overrides remote-only, and prior remote exposure remains an
included auditable subset. Total TBI calibration remains one-scale within
each subgroup. Progression is unconnected. The 0.005/year Australia/Viet Nam
outputs remain diagnostic; PRK remains blocked and its missing incidence is
not replaced. Frozen 2A/2C artifacts are not reinterpreted as country estimates.

## Deferred

- Risk-factor population allocation and TBI-risk-factor association.
- Approved recent/remote and age-specific progression forward operator.
- Active-incidence fitting and ascertainment/notification transformation.
- Migration, imported-disease and infection-selective mortality modelling.
- Endogenous transmission.
- Intervention, screening/treatment and economic simulation, DALYs/QALYs.
- Production UI, navigation and integration with any existing model runner.

## Remaining scientific/data review

Review population and survey representativeness, assay/threshold performance,
within-age and terminal-age distributions, risk-factor definitions/effects,
progression natural history, imported disease and historical exposure.
Review evidence that named population categories actually form a disjoint
partition; aggregate labels cannot prove it. Differing subgroup age grids need
reviewed demographic rebinning weights before combination. Complex survey
covariance and coherent mixture uncertainty propagation remain unimplemented.
The inherited 2C bounded multi-observation optimizer's basin coverage and
compatibility thresholds still require scientific review.

## Recommended next milestone: risk-factor population allocation

Use user-entered risk-factor prevalence, reviewed age- and sex-specific
patterns, and explicit TBI-risk-factor association to allocate risk-factor
states to the agreed subgroup population. When more than one reviewed risk
factor applies, use the **highest applicable reviewed RR**. Do not add RRs,
multiply RRs or multiply inherited ORs. Keep association and the risk summary
rule explicit and independently auditable. Do not introduce progression,
incidence inference or intervention effects as a side effect of allocation.

This is a recommendation only. Milestone 2D does not start that work.
