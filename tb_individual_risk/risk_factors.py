"""Explicit multimorbidity policies. No default numerical clinical effects."""
import math
from .contracts import finite, ValidationError
from .parameters import POLICIES, EXTERNAL, STATUSES
from .infection_probability import Unavailable,require_component


def combine_rr(effects, policy="highest_applicable_rr"):
    if policy not in POLICIES:
        raise ValidationError("combination_policy", "Unknown combination policy; no fallback.")
    values = {k: finite(v, k, minimum=1) for k,v in effects.items()}
    if policy == "highest_applicable_rr":
        rr = max([1.0] + list(values.values()))
    elif policy == "additive_excess_rr":
        try:rr = 1 + math.fsum(v-1 for v in values.values())
        except OverflowError as exc:raise Unavailable("Additive multiplier overflows numeric range; no cap applied.") from exc
    else:
        rr = math.prod(values.values())
    if not math.isfinite(rr):
        raise Unavailable("Effective multiplier overflows numeric range; no cap or policy switch applied.")
    highest = max([1.0] + list(values.values()))
    return {"effective_rr": rr, "policy": policy, "factors_present": sorted(values),
            "highest_contributors": sorted(k for k,v in values.items() if v == highest),
            "individual_rr": values, "cap": None}


def factor_definitions(parameters):
    definitions = {}
    for name in ("medical_risk_effects", "radiographic_effects"):
        component = parameters.component(name)
        for factor_id, spec in component.config.get("factors", {}).items():
            if factor_id in definitions:
                raise ValidationError("factors", "Duplicate factor ID across effect components.")
            definitions[factor_id] = (component, spec)
    return definitions


def risk_factors(values, parameters):
    for name in ("medical_risk_effects","radiographic_effects"):
        cfg=require_component(parameters.component(name))
        if "factors" not in cfg:raise Unavailable(name+": explicit included/excluded factor scope required")
        if cfg.get("cap") is not None:raise ValidationError("cap","Multiplier caps are not implemented; no cap silently applied.")
    effects = {}
    evidence = {}
    age = values.get("current_age")
    for factor_id, (component, spec) in factor_definitions(parameters).items():
        field = spec["field"]
        value = values.get(field)
        if value is None:
            raise Unavailable(field + ": missing factor state requires prior or bound")
        present = value == spec.get("present_value", True)
        if not present:
            continue
        for context,expected in spec.get("required_context",{}).items():
            if values.get(context) is None or values[context]!=expected:
                raise Unavailable(factor_id+": treatment/exposure context missing or incompatible with effect definition")
        if spec.get("status",EXTERNAL) not in STATUSES or spec.get("status",EXTERNAL)==EXTERNAL or spec.get("rr") is None:
            raise Unavailable(factor_id+": reviewed numerical effect external_parameter_required")
        if "age_min" not in spec or "age_max" not in spec:
            raise Unavailable(factor_id+": explicit age applicability external_parameter_required")
        if age is None:
            raise Unavailable("current_age required for factor applicability")
        if not float(spec.get("age_min", 0)) <= float(age) < float(spec.get("age_max", math.inf)):
            raise Unavailable(factor_id + ": not_applicable_pending_age_specific_review; adult RR not applied to children")
        if component.status == EXTERNAL or spec.get("status", EXTERNAL) == EXTERNAL or spec.get("rr") is None:
            raise Unavailable(factor_id + ": external_parameter_required")
        if not spec.get("source") or not spec.get("definition"):
            raise Unavailable(factor_id + ": effect provenance/definition required")
        effects[factor_id] = spec["rr"]
        evidence[factor_id] = {k: spec.get(k) for k in ("definition", "status", "source", "effect_measure")}
    result = combine_rr(effects, parameters.combination_policy)
    result["evidence"] = evidence
    result["warnings"] = []
    if result["effective_rr"] >= parameters.extreme_multiplier_warning:
        result["warnings"].append("Extreme multiplier: scientific plausibility requires review; no cap imposed.")
    if sum(float(v) >= parameters.extreme_multiplier_warning for v in effects.values()) > 1:
        result["warnings"].append("Multiple extreme effects present; policy sensitivity may dominate the risk.")
    if effects:
        result["warnings"].append("Relative risks used as hazard multipliers are an explicit approximation, not verified hazard ratios.")
    return result
