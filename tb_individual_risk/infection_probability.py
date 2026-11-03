"""Latent infection probability is distinct from assay positivity."""
import math
from .contracts import finite, ValidationError
from .parameters import EXTERNAL

class Unavailable(ValueError):
    pass


def require_component(component):
    if component.status == EXTERNAL:
        raise Unavailable(component.name + ": external_parameter_required")
    return component.config


def required_fields(parameters, values):
    cfg = parameters.component("infection_probability").config
    mode = cfg.get("mode")
    if mode == "supplied_probability":
        return {"tbi_probability"}
    if mode != "bayes_assays":
        return set()
    required = {"pretest_tbi_probability"}
    for assay in cfg.get("assays", []):
        result = values.get(assay + "_result")
        # Recorded TST induration can resolve a missing binary interpretation.
        if not (assay == "tst" and result is None and values.get("tst_mm") is not None):
            required.add(assay + "_result")
        acfg = parameters.component(assay + "_interpretation").config
        if result != "not_done":
            if acfg.get("by_bcg"):
                required.add("bcg")
            if acfg.get("product"):
                required.add("igra_product")
    return required


def infection_probability(values, parameters):
    cfg = require_component(parameters.component("infection_probability"))
    mode = cfg.get("mode")
    if mode == "supplied_probability":
        if values.get("tbi_probability") is None:
            raise Unavailable("tbi_probability: explicit evidence probability or prior required")
        p = finite(values["tbi_probability"], "tbi_probability", maximum=1)
        return {"probability": p, "method": mode, "updates": [],
                "assumptions": ["Supplied latent-infection probability is already conditional on the stated evidence; tests are not counted again."]}
    if mode != "bayes_assays":
        raise Unavailable("No implemented reviewed infection-probability model selected.")
    if values.get("pretest_tbi_probability") is None:
        raise Unavailable("pretest_tbi_probability: external probability/prior required")
    p = finite(values["pretest_tbi_probability"], "pretest_tbi_probability", maximum=1)
    assays = cfg.get("assays", [])
    if len(assays) != len(set(assays)) or any(a not in {"tst", "igra"} for a in assays):
        raise ValidationError("assays", "Use distinct supported assays.")
    if len(assays) > 1 and cfg.get("joint_policy") != "conditionally_independent_given_infection":
        raise Unavailable("Joint assays require an explicit reviewed joint-likelihood policy; no automatic independence.")
    # Log likelihoods avoid underflow and preserve deterministic p=0/1.
    li, lu = 0.0, 0.0
    trace = []
    for assay in assays:
        result = values.get(assay + "_result")
        if result == "not_done":
            trace.append({"assay": assay, "result": result, "update": "none_explicit_not_done"})
            continue
        acfg = require_component(parameters.component(assay + "_interpretation"))
        if assay == "tst" and values.get("tst_mm") is not None:
            if "positive_threshold_mm" not in acfg:
                raise Unavailable("TST threshold external_parameter_required")
            threshold = finite(acfg["positive_threshold_mm"], "TST threshold")
            interpreted = "positive" if float(values["tst_mm"]) >= threshold else "negative"
            if result is not None and result != interpreted:
                raise ValidationError("tst_result", "Recorded interpretation contradicts the selected TST threshold.")
            result = interpreted
        if result not in {"positive", "negative"}:
            raise Unavailable(assay + ": unknown/indeterminate/borderline requires an explicit observation model")
        if acfg.get("product") and values.get("igra_product") != acfg["product"]:
            raise Unavailable("IGRA product missing or incompatible; products are not interchangeable.")
        accuracy = acfg
        if acfg.get("by_bcg"):
            if values.get("bcg") is None:
                raise Unavailable("BCG status required by this TST observation model")
            accuracy = acfg["by_bcg"].get(str(values["bcg"]).lower(), {})
        if "sensitivity" not in accuracy or "specificity" not in accuracy:
            raise Unavailable(assay + " sensitivity/specificity external_parameter_required")
        se = finite(accuracy["sensitivity"], "sensitivity", maximum=1)
        sp = finite(accuracy["specificity"], "specificity", maximum=1)
        a, b = (se, 1-sp) if result == "positive" else (1-se, sp)
        li += math.log(a) if a else -math.inf
        lu += math.log(b) if b else -math.inf
        trace.append({"assay": assay, "result": result, "sensitivity": se, "specificity": sp,
                      "likelihood_given_infection": a, "likelihood_given_no_infection": b})
    a = (math.log(p) if p else -math.inf) + li
    b = (math.log1p(-p) if p < 1 else -math.inf) + lu
    scale = max(a, b)
    if scale == -math.inf:
        raise Unavailable("Evidence has zero probability under the selected observation model.")
    posterior = math.exp(a-scale)/(math.exp(a-scale)+math.exp(b-scale))
    return {"probability": posterior, "method": mode, "updates": trace,
            "assumptions": ["Fixed externally supplied assay accuracy; not estimated with infection probability."] +
            (["Assays conditionally independent given latent infection: explicit approximation."] if len(assays)>1 else [])}
