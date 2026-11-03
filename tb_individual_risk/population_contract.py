"""Design contract only. No synthetic population generator in Milestone 1."""
from dataclasses import dataclass, field, asdict
from .contracts import ValidationError
from .provenance import seal

@dataclass(frozen=True)
class SyntheticPopulationContract:
    population_id: str
    marginals: dict
    provenance: dict
    joint_policy: str = "external_required"
    subgroup_definitions: tuple = ()
    correlation_constraints: tuple = ()
    weighted_record_policy: str = "external_required"
    sampling_seed: object = None
    schema_id: str = "individual_tb.synthetic_population_specification.v1"
    status: str = "design_only_not_implemented"

    def __post_init__(self):
        if self.schema_id != "individual_tb.synthetic_population_specification.v1" or self.status != "design_only_not_implemented":
            raise ValidationError("population", "Unsupported contract; generation is not implemented.")
        if not self.population_id or self.joint_policy not in {"external_required", "explicit_independence", "reviewed_joint_distribution"}:
            raise ValidationError("population", "Stable ID and explicit joint policy required.")
        if set(self.marginals) != set(self.provenance) or any(not s for s in self.provenance.values()):
            raise ValidationError("provenance", "Each marginal requires its own source.")

    def record(self):
        return seal(asdict(self))
