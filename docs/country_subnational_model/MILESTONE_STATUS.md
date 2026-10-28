# Country and subnational TB milestone status

Updated 2026-10-28 for Milestone 2D. Historical milestone statuses remain in Git.

Two explicit baseline-TBI modes, reviewed/user-specified TBI pathways,
subgroup composition, population-weighted mixtures and a pure future-interface
specification are implemented in the isolated country-model layer.

Active-TB-incidence inference deliberately returns
`unavailable_pending_progression_and_risk_forward_operator`. It produces no
TBI estimate or temporary progression-model calculation. Observed/user inputs
support one-scale overall/age anchors, direct age/sex profiles, direct exclusive
recent/remote states, and total-only overrides with recency unavailable unless
explicit conditional allocation is selected.

Subgroup residence and exposure metadata are distinct. Mixtures require
reviewed disjoint membership, aligned age/sex grids and reconciled populations.
Unresolved groups and 100+ cells remain in every denominator; full-population
TBI is unavailable when coverage is incomplete. Risk-factor fields remain null
placeholders. All demonstrations use invented values, not country estimates.

Milestone 2C's repeated Poisson infection events, most-recent-event recency,
remote-history audit and one-scale calibration remain unchanged. Its Australia/
Viet Nam 0.005/year results remain diagnostic; PRK remains blocked. No
progression, intervention, migration, transmission or economics model is added.

Read [behavior](milestone2d_sources_subgroups.md), [contracts](milestone2d_contracts.md),
[validation/results](milestone2d_validation.md) and [decisions](decision_register_2d.md).
Protected worktrees/refs, frozen data, production UI and all runners remain
unchanged. Recommended next work is reviewed risk-factor population allocation;
it is not started automatically. This model is for planning and sequencing,
not denying care.
