# Country and subnational TB milestone status

Updated 2026-10-28 for Milestone 3A. Historical milestone statuses remain in Git.

The isolated country-model layer now allocates user-specified progression-factor
prevalences within age/sex/subgroup strata, retains explicit definition and review
status, handles overlap and multimorbidity, selects the highest applicable RR,
and applies an explicit feasible residual association between any TBI and the
declared non-recency risk-factor union. All marginals and TBI recency are retained.

Active-TB-incidence inference deliberately returns
`unavailable_pending_progression_and_risk_forward_operator`. It produces no
TBI estimate or temporary progression-model calculation. Observed/user inputs
support one-scale overall/age anchors, direct age/sex profiles, direct exclusive
recent/remote states, and total-only overrides with recency unavailable unless
explicit conditional allocation is selected.

Subgroups are allocated separately and combined by population counts. Unresolved
groups and 100+ cells remain in denominators; full-population TBI is unavailable
when its coverage is incomplete. Adult candidate effects are not applied to
children automatically. Provisional legacy effects require explicit diagnostic
use; the default effect policy is reviewed-only. All demonstrations use invented
values, not country estimates. The 2D source contracts retain their null risk
placeholder; the new allocation artifact references the unchanged baseline.

Milestone 2C's repeated Poisson infection events, most-recent-event recency,
remote-history audit and one-scale calibration remain unchanged. Its Australia/
Viet Nam 0.005/year results remain diagnostic; PRK remains blocked. No
progression, intervention, migration, transmission or economics model is added.

Read [allocation](milestone3a_allocation.md), [contracts](milestone3a_contracts.md),
[diagnostics](milestone3a_examples.md), [validation](milestone3a_validation.md)
and [decisions](decision_register_3a.md).
Protected worktrees/refs, frozen data, production UI and all runners remain
unchanged. Recommended next work is the age-, time-since-infection- and
risk-specific active-TB progression forward operator; it has not been started.
Stochastic individual generation remains deferred. This model is for planning and sequencing,
not denying care.
