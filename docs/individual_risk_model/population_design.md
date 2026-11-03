# Future synthetic-population design: contract only

`SyntheticPopulationContract` accepts a population ID, source-bearing marginal
inputs, subgroup definitions, an explicit dependence policy, optional correlation
constraints, a weighted-record policy and a future sampling seed. Generation is
not implemented. Its status is always `design_only_not_implemented`.

Future marginal inputs may cover exact age/sex, historical incidence, reviewed
TBI anchors, diabetes, smoking, HIV, renal disease, underweight and additional
factors. Every marginal needs geography/subgroup, year, definition, denominator,
uncertainty, review status and provenance. Missingness cannot be encoded as zero.
National figures must not be silently assigned to migrant or clinical subgroups.

Marginal prevalences do not identify the joint distribution. An independence
assumption must be a visible selected policy, not a generator default. Reviewed
joint tables and explicit correlation constraints should take precedence.
Inconsistent or infeasible constraints must fail without clipping.

Proposed sequence for a later milestone:

1. Construct source-preserving weighted seed records with stable IDs and subgroup
   membership. Weights reconcile population, age and sex denominators.
2. Fit approved marginals using iterative proportional fitting/raking with
   documented structural zeros, convergence rules and feasibility checks.
3. Add reviewed joint/correlation constraints and assess whether multiple joint
   distributions remain possible. Report sensitivity rather than claiming a
   uniquely identified population.
4. Only if individual sampling is needed, use an explicit RNG/seed, retain weights
   and the full generation contract, and validate Monte Carlo marginals under
   prespecified tolerances. Never use a global RNG or hide stochastic uncertainty
   within parameter uncertainty.
5. Predict and rank with the same person-record core, retaining every marginal's
   source and every generation assumption in scientific hashes.

WHO disease incidence does not directly identify TBI prevalence. Any future
incidence-to-infection input requires separately reviewed infection-pressure
history, progression, observation, mortality and migration models. This milestone
adds neither that conversion nor endogenous transmission. No Starsim or other
dynamic-model workstream is modified.
