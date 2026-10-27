# Country and subnational TB model architecture review

Milestone 2A implements the frozen incidence and demographic data layer. The model is for planning and sequencing, not for denying care. The original Milestone 1 architecture remains below as the scientific boundary; no projection engine or intervention is connected.

Current review: [data contracts](data_layer_contracts.md), [duplicate forensics](duplicate_forensics.md), [source policy](demographic_source_policy.md), [annual updates](annual_data_updates.md), [worked examples](worked_examples_2a.md), and [current milestone status](MILESTONE_STATUS.md).

Review the [architecture](architecture.md), then the [infection identifiability note](incidence_to_infection.md), [association mathematics](association.md), and [versioned contracts](schemas.md). The [implementation plan](implementation_plan.md) lists decisions required before executable modelling.

| Deliverable | Location |
| --- | --- |
| Comprehensive architecture | [architecture.md](architecture.md) |
| Candidate progression evidence and exclusions | [tstin3d_evidence.md](tstin3d_evidence.md) |
| WHO repository and installed snapshot audit | [who_2025_audit.md](who_2025_audit.md) |
| Demographic loader audit | [demography_reuse_audit.md](demography_reuse_audit.md) |
| Existing code reuse assessment | [reuse_matrix.md](reuse_matrix.md) |
| Incidence and infection identifiability | [incidence_to_infection.md](incidence_to_infection.md) |
| Feasible association scale | [association.md](association.md) |
| Seven proposed schemas and identifiers | [schemas.md](schemas.md), [machine schema](../../schemas/country_tb/v1/contracts.schema.json) |
| Staged implementation | [implementation_plan.md](implementation_plan.md) |
| Replaceable milestone status | [MILESTONE_STATUS.md](MILESTONE_STATUS.md) |
| Starting Git state and protection checks | [git_audit.md](git_audit.md) |

The isolated mathematics and validation prototype under engine/country_tb imports only the Python standard library. It is not imported by existing apps, runners or engines. The schema example is invented software test data, not a country evidence default. No existing model code was copied or cherry-picked.
