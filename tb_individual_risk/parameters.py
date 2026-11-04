"""Explicit parameter objects; no clinical defaults or embedded imputation."""
from copy import deepcopy
from dataclasses import dataclass, field, asdict
from .contracts import ValidationError, finite, FIELDS
from .provenance import seal

EXTERNAL = "external_parameter_required"
STATUSES = {EXTERNAL, "reviewed", "synthetic_not_clinical", "sensitivity_only"}
COMPONENTS = ("infection_probability", "tst_interpretation", "igra_interpretation",
              "baseline_progression_hazard", "recent_infection_progression",
              "remote_infection_progression", "age_effects", "radiographic_effects",
              "medical_risk_effects", "competing_mortality", "prediction_horizons")
POLICIES = {"highest_applicable_rr", "additive_excess_rr", "multiplicative_rr"}

@dataclass(frozen=True)
class Parameter:
    name: str
    status: str = EXTERNAL
    source: str = EXTERNAL
    config: dict = field(default_factory=dict)
    schema_id: str = "individual_tb.parameter.v1"

    def __post_init__(self):
        if self.status not in STATUSES or not self.source:
            raise ValidationError(self.name, "Parameter requires source and recognized evidence status.")
        if self.schema_id != "individual_tb.parameter.v1":
            raise ValidationError("schema_id", "Unsupported parameter contract.")
        if not isinstance(self.config, dict):
            raise ValidationError(self.name, "Parameter config must be a mapping.")

@dataclass(frozen=True)
class PopulationPrior:
    source: str = EXTERNAL
    status: str = EXTERNAL
    distributions: dict = field(default_factory=dict)
    bounds: dict = field(default_factory=dict)
    dependence_policy: str = "external_required"
    joint_states: tuple = ()
    joint_bound_states: tuple = ()
    schema_id: str = "individual_tb.population_prior.v1"

    def __post_init__(self):
        if self.schema_id != "individual_tb.population_prior.v1" or self.status not in STATUSES or not self.source:
            raise ValidationError("priors", "Invalid prior version/status/source.")
        if self.dependence_policy not in {"external_required", "explicit_independence", "explicit_joint"}:
            raise ValidationError("priors", "Unknown dependence policy.")
        if self.distributions and self.dependence_policy != "explicit_independence":
            raise ValidationError("priors", "Marginal distributions require explicit independence; otherwise supply a joint table.")
        if (self.joint_states or self.joint_bound_states) and self.dependence_policy != "explicit_joint":
            raise ValidationError("priors", "Joint supports require explicit_joint policy.")
        if self.status == EXTERNAL and (self.distributions or self.bounds or self.joint_states or self.joint_bound_states):
            raise ValidationError("priors", "External-required priors cannot supply numerical assumptions.")
        for name, rows in self.distributions.items():
            if not rows or any(set(r) != {"value", "probability"} for r in rows):
                raise ValidationError(name, "Distribution needs value/probability rows.")
            weights = [finite(r["probability"], name, maximum=1) for r in rows]
            if abs(sum(weights) - 1) > 1e-12:
                raise ValidationError(name, "Prior probabilities must sum to one; no clipping.")
        for name, values in self.bounds.items():
            if not values:
                raise ValidationError(name, "Bounds need prespecified nonempty states.")
        for rows in (self.joint_states, self.joint_bound_states):
            if rows:
                names = set(rows[0]["values"])
                if not names or any(set(r["values"]) != names for r in rows):
                    raise ValidationError("priors", "Joint rows must cover the same explicit fields.")
        if self.joint_states:
            weights = [finite(r["probability"], "joint prior", maximum=1) for r in self.joint_states]
            if abs(sum(weights)-1) > 1e-12:
                raise ValidationError("priors", "Joint probabilities must sum to one.")

@dataclass(frozen=True)
class ModelParameters:
    components: dict = field(default_factory=lambda: {n: Parameter(n) for n in COMPONENTS})
    population_prior: PopulationPrior = field(default_factory=PopulationPrior)
    combination_policy: str = "highest_applicable_rr"
    classification: str = "external_required"
    max_scenarios: int = 8192
    extreme_multiplier_warning: float = 20.0
    schema_id: str = "individual_tb.model_parameters.v1"

    def __post_init__(self):
        if self.schema_id != "individual_tb.model_parameters.v1" or set(self.components) != set(COMPONENTS):
            raise ValidationError("parameters", "Exactly the versioned component families are required.")
        if self.combination_policy not in POLICIES or self.classification not in {"external_required", "reviewed", "synthetic_not_clinical", "sensitivity_only"}:
            raise ValidationError("parameters", "Unknown policy/classification; no fallback.")
        if type(self.max_scenarios) is not int or not 1 <= self.max_scenarios <= 1000000:
            raise ValidationError("max_scenarios", "Explicit finite enumeration budget required.")
        finite(self.extreme_multiplier_warning, "warning threshold", minimum=1)
        for name, component in self.components.items():
            if not isinstance(component, Parameter) or component.name != name:
                raise ValidationError(name, "Component identity mismatch.")
        statuses = [p.status for p in self.components.values()] + [self.population_prior.status]
        def evidence_states(value):
            if isinstance(value,dict):
                for key,item in value.items():
                    if isinstance(item,str) and (key=="status" or key=="adoption" and item=="historical_only"):yield item
                    yield from evidence_states(item)
            elif isinstance(value,(list,tuple)):
                for item in value:yield from evidence_states(item)
        nested=list(evidence_states({k:asdict(v) for k,v in self.components.items()}))
        if "historical_only" in nested:
            raise ValidationError("parameters","Historical-only evidence cannot enter the hazard formulation; use the separate reference package.")
        statuses += nested
        if "synthetic_not_clinical" in statuses and self.classification != "synthetic_not_clinical":
            raise ValidationError("classification", "Synthetic parameters cannot be relabelled reviewed or clinical.")
        if self.classification == "reviewed" and any(s != "reviewed" for s in statuses):
            raise ValidationError("classification", "Reviewed set cannot contain unresolved/sensitivity components.")

    def record(self):
        return seal(asdict(self))

    @classmethod
    def from_dict(cls, raw):
        data = deepcopy(raw)
        if "scientific_hash" in data and seal(data)["scientific_hash"] != data["scientific_hash"]:
            raise ValidationError("scientific_hash", "Parameter scientific hash mismatch.")
        data.pop("scientific_hash", None)
        data["components"] = {k: Parameter(**v) for k,v in data["components"].items()}
        data["population_prior"] = PopulationPrior(**data["population_prior"])
        return cls(**data)

    def component(self, name):
        return self.components[name]
