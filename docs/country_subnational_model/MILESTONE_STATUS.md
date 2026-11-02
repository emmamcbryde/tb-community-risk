# Country and subnational TB milestone status

Updated 2026-11-02: **Milestone 3E implemented: isolated scientific-review UI ready for recheck.**

Launch the separate entry from repository root:

```powershell
python -m streamlit run country_model_review_ui/app.py --server.address 127.0.0.1 --server.enableCORS true --server.enableXsrfProtection true --browser.gatherUsageStats false
```

The seven-view interface is isolated under `country_model_review_ui/`, with no
adjacent `pages/` directory and no legacy/deployed navigation integration.
See [architecture and contracts](milestone3e_ui.md),
[validation and performance](milestone3e_validation.md),
[decisions](decision_register_3e.md) and the [application guide](../../country_model_review_ui/README.md).

The operational route is observed/user-specified 2024 TBI: overall total-only or
explicit catalytic age/recency derivation, direct age/sex and recent/remote CSVs,
subgroup-specific sources and population-weighted mixtures. Original values,
assay positivity versus latent infection, and model comparisons remain distinct.
Total-only input has unavailable recency unless explicit hybrid allocation is
selected. Missing subgroup inputs never inherit an overall value silently.

Risk inputs retain zero, missing, excluded, not-applicable and unreviewed states,
exact definitions, age/sex structure and provenance. Central combination is the
highest applicable reviewed RR. Legacy effects are sensitivity-only and adult
effects do not automatically apply to children. Normalized residual TBI/any-risk
association preserves marginals and reports achieved phi separately.

Deterministic instantaneous sensitivity results and optional bounded one-year
first-episode projections are available from supplied baseline populations.
Children, unknown 100+ ages, missing clocks/effects and unsupported remote tails
remain unresolved. No extrapolation never means zero progression. Positive
unresolved coverage prevents complete-population rates. Competing deaths and
first disease deplete the prospective risk set. No long-horizon country claim
is made. Material edits mark results and prepared exports stale.

Production incidence inference remains
`unavailable_pending_progression_and_risk_forward_operator`.
An explicitly enabled, separately uploaded, source-bound diagnostic configuration
must pass every scientific gate, then attainable-range/unique-root/refinement
checks. UI flags and direct navigation cannot bypass these. No real country
configuration or central parameter set is approved. Synthetic known-scale
validation does not approve empirical WHO inversion.

Milestone 3D was internal evidence synthesis/model checking, not external
validation. The four-point adult curve remains sensitivity-only; pediatric
curves, central remote tail, conditional risk effects, mortality adoption and
migration/episode reconciliation need scientific review. Its observed comparison
RMSEs remain 8.25 points (Viet Nam 2020), 2.71 (2021), 1.18 (US-born NHANES) and 22.14
(cross-campaign transport). Contact samples do not validate national inference.
The fictional Viet Nam 65.33% result remains fictional and is not used by the UI.

All retained 2A through 3D implementations/contracts, frozen artifacts and protected
worktrees/refs remain unchanged. No production entry point, APY, explicit-TBI,
Starsim/Summer, MATLAB, intervention or economics work is integrated. No merge,
tag or deployment is part of 3E. Full precision is retained in lazy JSON/CSV
review exports; no Excel workbook is added.

Next: reviewer usability/scientific-input recheck, followed by targeted matched
data acquisition and the pending component decisions. Production inference,
production UI, interventions, economics, DALYs/QALYs and endogenous transmission
remain deferred. This model is for planning and sequencing, not denying care.
