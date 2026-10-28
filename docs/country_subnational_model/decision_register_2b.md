# Milestone 2B decision register and narrow Milestone 2C

| Status | Decision |
| --- | --- |
| Approved, retained | 2024 baseline; WHO report vintage 2025 with final observation 2024; separate immutable future reports; no automatic migration |
| Approved, retained | Separate sex strata, no initial sex difference in diagnosis/treatment effects, source hierarchy and medium future demography |
| Approved, retained | WPP 2024 remains projection; ABS June 2024 estimate; no double-counted both-sex rows, manufactured sex split or automatic successor-state mapping |
| Approved, retained | Fine internal age resolution, display 0–4/5–14/15–44/45–64/65+; missing PRK incidence remains missing |
| Approved, retained | Incidence alone cannot identify infection scale/TBI; no Styblo; no implicit equality of incidence and pressure trends |
| Provisional central, review needed | First-infection absorbing diagnostic for software validation only; no country inference from its exclusions |
| Provisional central, review needed | One calibrated pressure scale with fixed history/assay/exposure first; no free progression scale simultaneously |
| Provisional central, review needed | Explicit constant pre-2000 pressure as a transparent reference scenario when no historical evidence; never silently selected or described as historically true |
| Provisional central, review needed | Low-dimensional log-linear and COVID-aware descriptors compared under reviewed natural history before an empirical pilot; no production method selected |
| Required sensitivity | Constant vs bounded early-trend backcast vs historical/anchor-informed scenarios; recent window 2 vs 5 years; adult/child age exposure; assay and survival assumptions |
| Required sensitivity | WHO coherent bounds, trend window/COVID policy, mapping/kernel, mortality selection, migration and data vintage |
| Scientific review | Target viable/retained/ever infection vs assay positivity; clearance vs reversion; first vs last-event recency; reinfection susceptibility and clock reset |
| Scientific review | Progression by age at infection, current age and time since infection; separate 0–4, 5–14, 15+ evidence with finer ages where supported; HIV and disease/treatment exits |
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
