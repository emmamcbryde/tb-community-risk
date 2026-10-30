# Country and subnational TB milestone status

Updated 2026-10-30 for Milestone 3C. Historical milestone statuses remain in Git.

Independent3B diagnostics reproduce. New offline diagnostic inversion checks
scientific gates, attainable range, monotonicity/all roots, sensitivity and
numerical refinement. Production availability remains false. Adult surrogate
is sensitivity-only; pediatric curves and remote tail remain external-required;
legacy RR hazard use and mortality/observation/migration assumptions are not
approved for country inversion. Only explicitly fictional closed-population
examples fit a scale. Actual unresolved-threshold approval remains pending.
See [3C scientific review](milestone3c_scientific_review.md),
[inversion contracts](milestone3c_inversion.md) and [decisions](decision_register_3c.md).

The isolated layer now projects first active-TB episodes after a supplied baseline,
using attained age, time since most recent infection and the highest applicable
RR through cumulative hazard. It supports competing mortality and explicit
exogenous future infection/reinfection, with deterministic expected values and a
continuous event-time Monte Carlo validator. All retained 2C–3A source contracts
and implementations are unchanged; source marginals are preserved at initialization.

Active-TB-incidence inference deliberately returns
`unavailable_pending_progression_and_risk_forward_operator`. Forward code does
not approve inverse calibration. The3C review withholds country approval;
external evidence and observation reconciliation remain required. Observed/user inputs
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
intervention, migration, transmission or economics model is added. The new forward
operator is a separate explicit consumer, not an integration into those old paths.

Pediatric curves remain external-required. The five-knot adult candidate is not
adopted unqualified; a four-knot source-audited surrogate is diagnostic only.
Central tail policy is no extrapolation. No reviewed national mortality snapshot
has been adopted. Missing components block complete mathematical projections;
no output is labelled a complete country estimate.

Read [evidence](milestone3b_evidence.md), [mathematics](milestone3b_forward.md),
[contracts](milestone3b_contracts.md), [diagnostics](milestone3b_examples.md),
[validation](milestone3b_validation.md) and [decisions](decision_register_3b.md).
Protected worktrees/refs, frozen data, production UI and all runners remain
unchanged. Recommended next work is independent incidence reconciliation,
identifiability analysis and approval or rejection of incidence-based TBI inversion;
it has not been started. A full individual factor-membership population generator
remains deferred. This model is for planning and sequencing,
not denying care.
