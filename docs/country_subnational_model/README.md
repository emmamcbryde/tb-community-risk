# Country and subnational TB model architecture review

Milestone 3A implements isolated risk-factor prevalence allocation, multimorbidity and explicit within-age/sex TBI association. Active-TB-incidence inference remains deliberately blocked pending an approved progression and risk forward operator. The model is for planning and sequencing, not denying care. All demonstrations remain synthetic/diagnostic; no production application or runner is connected.

Current 3A review: [allocation and evidence](milestone3a_allocation.md), [contracts](milestone3a_contracts.md), [invented diagnostics](milestone3a_examples.md), [validation](milestone3a_validation.md), and [decisions and next milestone](decision_register_3a.md).

Retained 2D review: [source and subgroup behavior](milestone2d_sources_subgroups.md), [executable contracts](milestone2d_contracts.md), [validation and synthetic results](milestone2d_validation.md), and [approved/deferred decisions](decision_register_2d.md).

Retained 2C review: [scientific implementation](milestone2c_calibration.md), [executable contracts](milestone2c_contracts.md), [diagnostic results](milestone2c_examples.md), [recovery and validation](milestone2c_validation.md), and [updated decisions](decision_register_2b.md).

Retained 2B review: [source recheck](frozen_recheck_2b.md), [scientific specification](calibration_specification.md), [mapping candidates](mapping_candidates.md), [backcast sensitivities](backcast_synthetic_examples.md), [proposed contracts](calibration_contracts.md), [decision register and Milestone 2C](decision_register_2b.md), and [validation](milestone2b_validation.md).

Retained 2A layer: [data contracts](data_layer_contracts.md), [duplicate forensics](duplicate_forensics.md), [source policy](demographic_source_policy.md), [annual updates](annual_data_updates.md), [worked examples](worked_examples_2a.md), and [current milestone status](MILESTONE_STATUS.md).

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
