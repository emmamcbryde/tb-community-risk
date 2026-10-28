# Milestone 2B/2C decision register

Updated 2026-10-28 for the explicit Milestone 2C scientific correction.
See [implementation](milestone2c_calibration.md) and
[recovery/validation](milestone2c_validation.md). The original 2B milestone
sequence below remains historical context.

| Status | Decision |
| --- | --- |
| Approved, retained | 2024 baseline; WHO report vintage 2025 with final observation 2024; separate immutable future reports; no automatic migration |
| Approved, retained | Separate sex strata, no initial sex difference in diagnosis/treatment effects, source hierarchy and medium future demography |
| Approved, retained | WPP 2024 remains projection; ABS June 2024 estimate; no double-counted both-sex rows, manufactured sex split or automatic successor-state mapping |
| Approved, retained | Fine internal age resolution, display 0â€“4/5â€“14/15â€“44/45â€“64/65+; missing PRK incidence remains missing |
| Approved, retained | Incidence alone cannot identify infection scale/TBI; no Styblo; no implicit equality of incidence and pressure trends |
| Implemented 2C central, review required before inference | Repeated non-homogeneous Poisson infection events; same pressure after previous infection; no protection or clearance; most recent event determines recency; recent reinfection overrides remote-only; prior remote exposure remains auditable |
| Provisional central, review needed | One calibrated pressure scale with fixed history/assay/exposure first; no free progression scale simultaneously |
| Implemented 2C central, provisional science | Hold earliest normalized pressure shape constant before first observation, covering the full lifetime; bounded backward log-linear and reviewed supplied histories are sensitivity policies |
| Implemented 2C diagnostic, provisional science | Five-year centred local-linear OLS of positive log incidence, truncated three-point edges, normalized at 2024; raw/bounds/residual/edge/COVID diagnostics retained; no force-of-infection claim |
| Required sensitivity | Constant vs bounded early-trend backcast vs historical/anchor-informed scenarios; recent window 2 vs 5 years; adult/child age exposure; assay and survival assumptions |
| Required sensitivity | WHO coherent bounds, trend window/COVID policy, mapping/kernel, mortality selection, migration and data vintage |
| Implemented 2C clock, provisional science | Most recent reinfection resets progression clock; no connection to active-TB progression. Total TBI cannot independently identify protection, clearance, clock effects or infection-time detail |
| Scientific review | Interpretation of retained/viable infection versus persistent assay positivity; evidence for departures from no protection/no clearance; assay-specific performance and thresholds |
| Scientific review | Progression by age at infection, current age and time since infection; separate 0â€“4, 5â€“14, 15+ evidence with finer ages where supported; HIV and disease/treatment exits |
| Scientific review | Historical disease observation/detection/relapse kernel and denominator consistency; interpretation of WHO uncertainty/covariance |
| New data | Reviewed TBI surveys with assay thresholds, age/sex, geography/year, sampling/weights and uncertainty; historical infection information |
| New data | Sex/age mortality tables, migration/origin/residence exposure histories where material, open-tail age weights, historical denominators |
| New data | Usable PRK incidence or independently supplied pressure/anchor; no invention from neighbouring countries |
| Deferred dynamic | Endogenous infectiousness, contact mixing, transmission, partial immunity/reinfection feedback and indirect effects |
| Deferred later work | Medical-risk allocation, interventions, economics, DALYs/QALYs, production UI and all engine integration |

Milestone 2C should validate catalytic calibration before allocating risk factors:

1. Review the target infection state, exact baseline timing, first/last-event
   semantics and diagnostic exclusions. Choose one explicit synthetic history
   and one-scale parameter bounds.
2. Implement reviewed observation ingestion and validation, retaining counts,
   assay/threshold, design and original-source hashes. Start with synthetic
   perfect-assay data and independently checkable roots.
3. Validate one-scale likelihood/objective fitting, infeasible targets,
   uncertainty coverage and exact-age anchors; report rank/condition and
   boundary diagnostics. Test age-cell integration and open-tail limitations.
4. Add reviewed assay operators and one real survey pilot with sex strata and
   explicit backcast scenarios. If no defensible anchor exists, remain in
   scenario mode. PRK remains blocked for incidence-informed inference.
5. Specify/review the minimum natural-history disease operator before claiming
   WHO incidence calibration; then compare the parsimonious pressure-shape
   candidates, using independent/historical evidence and held-out checks.
6. Freeze a reviewed age-profile artifact and infection-time semantics usable
   by both engines. Do not yet allocate risk factors or run either engine.

Acceptance requires reproducible synthetic recovery, honest uncertainty and
incompatibility reporting, source-preserving overrides, sex/count
reconciliation and unchanged frozen/protected components. A successful
optimizer alone is insufficient.


## Milestone 2C completion decisions

| Status | Decision |
| --- | --- |
| Implemented | One shared non-negative alpha; supplied scale, one overall bounded root, actual-count binomial likelihood or explicitly selected prevalence/uncertainty WLS |
| Implemented | Deterministic most-recent-event CDF/inverse, moments/quantiles, recent and remote-only summaries; explicit RNG validation only |
| Implemented | Reviewed direct latent-TBI or fixed binary Se/Sp observation operators; no joint fitting of assay performance |
| Implemented | Uniform attained age within single-year cells; reviewed alternatives for finite broader cells; no display midpoint |
| Implemented | Unreviewed open tails remain unestimated with full population reconciliation and exclusion fraction; reviewed discrete age weights supported |
| Implemented | Separate source-preserving total-only override has null recency/time composition and a new scientific hash |
| Implemented | Separate WHO-bound/backcast, assay and synthetic-anchor uncertainty scenarios; no TBI confidence-interval label |
| Diagnostic only | Australia and Viet Nam use invented alpha .005/year and t0=2024.5; no empirical anchor adopted; PRK remains blocked |
| Required before country inference | Review infection target, historical exposure/migration, within-year/open-tail ages, assay threshold/performance, survey weights/design, dependencies and empirical anchors |
| Remaining method review | Multi-observation objective basin coverage, covariance/design-based objectives, compatibility thresholds and independent held-out observations |
| Deferred | Risk allocation/association, progression, interventions, screening/treatment, economics, DALYs/QALYs, endogenous transmission and all production/runner integrations |

Candidate surveys are listed in [examples](milestone2c_examples.md) without
adopting their estimates. The central migration assumption is an explicit
limitation, not a claim that migration is absent. Source-rounding differences
remain visible. Protected worktrees/refs and frozen 2A inputs remain unchanged.
