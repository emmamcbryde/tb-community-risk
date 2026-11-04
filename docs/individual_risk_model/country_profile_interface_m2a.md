# Future frozen country-profile interface: specification only

Proposed contract ID: `individual_tb.frozen_country_profile.v1`. No repository,
download, profile loader or population generation is implemented in2A.

A future immutable profile must include profile ID/version, geography code/level,
source years, observation/reference year, estimate/projection status, population
definition, subgroup membership and exposure/residence distinctions, licence,
source citation/URL/access date, raw-file checksum and scientific hash.

Population rows require exact age or explicit lower/upper-exclusive bounds,
open-age indicator, sex, population count, source unit/rounding and coverage.
Female and male records remain separate; published both-sex totals reconcile
rather than replace them. Open ages require reviewed internal age weights or stay
unresolved; display midpoints are not calculation ages.

Factor rows require stable factor ID, exact definition, age range, sex, subgroup,
prevalence, numerator/denominator when available, uncertainty bounds and type,
source year, sampling frame, review status, missingness and per-input provenance.
HIV requires treatment/context definitions; diabetes type/ascertainment, smoking
current/former/intensity and >1-pack/day distinction, and underweight metric and
threshold remain explicit. Renal disease is not automatically haemodialysis.
Unknown, measured zero, not applicable, excluded, invalid and unavailable values
are distinct. A missing cell does not inherit a national or opposite-sex value.

Candidate allocation policy:
`conditional_independence_within_exact_age_sex_and_subgroup_v1`. This is a proposed
initial assumption requiring explicit selection, not an identified joint
population. Marginal prevalences do not determine correlations or higher-order
overlap. Reviewed joint data and constraints should supersede it.

Under the highest-RR rule, positive clustering can make effects redundant;
independence may spread elevated risk over more people and yield higher aggregate
risk. In that sense it is plausibly conservative for some burden assessments.
It is **not guaranteed conservative** in every population or for individual ranks,
association with infection, age structure, mortality or intervention comparisons.
The direction must be assessed under the full reviewed joint structure.

Future inputs may include source-preserving historical incidence and TBI anchors,
but WHO incidence must not be translated directly to TBI without a separately
reviewed infection-pressure/natural-history/observation model. Weighted records,
raking/IPF, explicit correlations and seeded stochastic sampling remain future
work described in the [population design](population_design.md).
