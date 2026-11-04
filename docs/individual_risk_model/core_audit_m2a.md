# Independent Milestone 1 core audit

Reviewed actual implementation at commit8b533f81a6a1435bed11f10db52c6d8d916224a1,
not its completion report. Read contracts, observation likelihoods, missing-data
enumeration, hazard segments/competing death, factor policies, ranking, exports
and canonical hashing. New failing regression cases were written and run before
introducing reference parameters: five failed and one offline/product-identity
case passed. After focused corrections, all six plus the retained69 M1 cases
passed (75 tests in4.55s).

## Findings and corrections

1. `cxr=other` matched neither known radiographic factor and could receive RR1.
   Unclassified findings now require a reviewed effect or explicit exclusion.
2. Evidence checks inspected only top-level component statuses. Nested synthetic
   curves/effects could survive relabelling their containers reviewed. Recursive
   checks preserve nested evidence status and reject historical-only adoption in
   hazard configurations. Adoption metadata is not confused with evidence status.
3. Conditional-progression rankings tested overlap using absolute-risk intervals.
   People with different infection probabilities could have disjoint absolute
   ranges but identical conditional ranges. Separate conditional missing-scenario
   fields now drive overlap for that ranking metric; unavailable bounds stay visible.
4. Tie sorting stringified identifiers, so integer1 and text"1" could depend on
   input order. Canonical typed IDs now provide stable ordering. Source-row indices
   remain audit locators rather than invented person IDs.
5. Predictions included parameter/prior hashes but not the whole parameter record.
   The complete source definition is now retained, permitting reconstruction. Package,
   formulation, source version and explicit policy overrides add reproducibility.

## Arithmetic confirmed

Infection likelihoods update latent probability rather than identifying positivity
with infection. Conditional progression uses segment cumulative-hazard increments,
not RR times probability. Expected absolute risk integrates `p_s*r_s`; infection-
weighted conditional risk preserves the product identity under joint priors.
Highest/additive-excess/product policies, factor ties, unknown/absent/excluded
states, explicit dependence, first-event depletion and competing death were
reviewed and retained. No shared numerical curves were tuned for historical fit.

No prediction path imports a networking client or makes a web request. Tests deny
socket connections for both formulations. Source browsing and PDF inspection are
separate development activities; ordinary offline predictions do not fetch evidence.

The current model remains a baseline-infection model: no future acquisition,
reinfection, endogenous transmission, migration/changing exposure, recurrence or
future changes in factors. Missing history and unsupported pediatric effects are
not silently approximated. These are scientific limitations, not audit successes
that establish clinical validity.
