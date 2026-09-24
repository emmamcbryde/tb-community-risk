"""Effect-measure semantics of the current analysis engine, and warnings.

No effect measure is converted in this milestone. ``ConversionPolicy.NONE`` is the
only active policy; the other policies describe designs awaiting scientific review
and raise if selected.

Current engine semantics (engine/apy/cohort.py, engine/apy/calibration.py):

* Progression to disease: each person's early and late progression hazards are
  ``lambda * prod(effect_k ** flag_k)`` over the disease-effect factors (cannabis,
  contact, renal, diabetes, smoking, chronic lung disease, alcohol/drugs). The
  declared measure type is ignored: RR, HR and OR values are all used as hazard
  multipliers. There is no cap. The baseline hazard ``lambda`` is recalibrated so that
  the population-average risk of active TB within the calibration horizon matches
  the target, so the effects redistribute risk between people rather than add
  to it.
* Infection: separate infection effects (cannabis, contact, renal, from the bundled
  input table) enter a complementary log-log model,
  ``P(infected) = 1 - exp(-H(age) * prod(OR_k ** flag_k))``. They therefore also act
  as multipliers on the cumulative infection hazard. They are not editable in the
  general interface.
* Joint effects multiply (no interaction terms), so correlated factors (for example
  smoking and chronic lung disease) may be double counted if each estimate was
  unadjusted for the others.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any

from engine.profiles.population_profile import EffectMeasure, PopulationProfile, ValueState


EXTREME_COMBINED_MULTIPLIER = 20.0


class ConversionPolicy(str, Enum):
    NONE = "none"
    OR_TO_RR_BASELINE_RISK = "or_to_rr_baseline_risk"
    RR_TO_HR_CONSTANT_HAZARD = "rr_to_hr_constant_hazard"


ACTIVE_CONVERSION_POLICY = ConversionPolicy.NONE
POLICY_STATUS = {
    ConversionPolicy.NONE: "Active: effect estimates are used exactly as entered.",
    ConversionPolicy.OR_TO_RR_BASELINE_RISK: "Designed, inactive: RR = OR / (1 - p0 + p0 * OR) needs a baseline risk p0; awaiting review.",
    ConversionPolicy.RR_TO_HR_CONSTANT_HAZARD: "Designed, inactive: HR from RR under a constant-hazard assumption over a stated horizon; awaiting review.",
}

ENGINE_INTERPRETATION = {
    EffectMeasure.HR: "Used directly as a progression-hazard multiplier (consistent with its definition).",
    EffectMeasure.RR: "Used as a progression-hazard multiplier without conversion; approximates a hazard ratio only when risks are small.",
    EffectMeasure.OR: "Used as a progression-hazard multiplier without conversion; overstates the effect when outcomes are common.",
    None: "No measure type recorded; the value is used as an unspecified hazard multiplier.",
}
SCIENTIFIC_SUPPORT = {
    EffectMeasure.HR: "Supported",
    EffectMeasure.RR: "Approximate",
    EffectMeasure.OR: "Approximate; review",
    None: "Unclear",
}


@dataclass(frozen=True)
class EffectWarning:
    code: str
    message: str
    risk_factor: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {"code": self.code, "message": self.message, "riskFactor": self.risk_factor}


def convert_effect(value: float, measure: EffectMeasure | None, policy: ConversionPolicy = ACTIVE_CONVERSION_POLICY, **_: Any) -> float:
    """Return the value the engine will use. Only the no-conversion policy is active."""
    if policy is ConversionPolicy.NONE:
        return value
    raise NotImplementedError(f"Conversion policy {policy.value!r} is designed but inactive pending scientific review.")


def enabled_engine_factors(profile: PopulationProfile) -> list:
    return [
        factor
        for factor in profile.risk_factors
        if factor.enabled and factor.engine_key is not None and factor.prevalence.state is ValueState.VALUE
    ]


def effect_warnings(profile: PopulationProfile) -> list[EffectWarning]:
    """Warnings about how effect estimates are applied; they never change results."""
    warnings: list[EffectWarning] = []
    factors = enabled_engine_factors(profile)
    for factor in factors:
        if factor.effect_measure is EffectMeasure.OR:
            warnings.append(
                EffectWarning(
                    "or_as_hazard_multiplier",
                    f"{factor.label}: an odds ratio is applied as a progression-hazard multiplier without conversion.",
                    factor.label,
                )
            )
        elif factor.effect_measure is None and factor.effect_estimate.state is ValueState.VALUE:
            warnings.append(
                EffectWarning("missing_measure_type", f"{factor.label}: effect-measure type is missing.", factor.label)
            )
    if len(factors) > 1:
        warnings.append(
            EffectWarning(
                "multiplied_factors",
                f"{len(factors)} risk factors are combined by multiplying their effects; if the estimates are not mutually "
                "adjusted, correlated factors may be double counted.",
            )
        )
    combined = 1.0
    for factor in factors:
        if factor.effect_estimate.value is not None:
            combined *= factor.effect_estimate.value
    if factors and combined > EXTREME_COMBINED_MULTIPLIER:
        warnings.append(
            EffectWarning(
                "extreme_combined_multiplier",
                f"A person with all enabled risk factors would have a {combined:,.0f}-fold progression hazard.",
            )
        )
    return warnings


def crosswalk_rows(profile: PopulationProfile) -> list[dict[str, str]]:
    rows = []
    for factor in profile.risk_factors:
        measure = factor.effect_measure
        rows.append(
            {
                "Risk factor": factor.label,
                "Reported measure": measure.value if measure else "Not recorded",
                "Engine interpretation": ENGINE_INTERPRETATION[measure] if factor.engine_key else "Not used by the current engine.",
                "Affected quantity": "Early and late progression hazards (infection to disease)" if factor.engine_key else "None",
                "Scientific support": SCIENTIFIC_SUPPORT[measure] if factor.engine_key else "Not applicable",
                "Joint effects": "Multiplied with other enabled factors; no interaction terms",
                "Cap applied": "No",
                "Source estimate adjusted": "Not recorded" if not factor.notes else factor.notes,
            }
        )
    return rows
