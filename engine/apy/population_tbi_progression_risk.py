"""Aggregate scenario adjustment, separate from TBI prevalence and the base curve.

This is not a risk-factor OR product, fitted target, or independently estimated
relative risk. The runner's general default is standard; the isolated UI chooses
the provisional risk-enriched working scenario explicitly.
"""
from __future__ import annotations

from math import isfinite
from numbers import Real
from decimal import Decimal
from typing import Any, Mapping

CONTRACT_VERSION = "population_tbi_progression_risk_adjustment_v1"
ADJUSTMENT_FIELD = "populationTBIProgressionRiskAdjustment"
MODES = ("standard", "risk_enriched", "custom")


def positive_multiplier(value: Any) -> float:
    if isinstance(value, bool) or not isinstance(value, (Real, Decimal)):
        raise ValueError("Progression-risk multiplier must be an explicitly supplied positive finite number.")
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("Progression-risk multiplier must be a positive finite number.") from exc
    if not isfinite(number) or number <= 0.0:
        raise ValueError("Progression-risk multiplier must be positive and finite; values are not clipped.")
    return number


def build_progression_risk_adjustment(
    mode: str = "standard", *, multiplier: Any = None,
    source: str | None = None, review_status: str | None = None, notes: str = "",
) -> dict[str, Any]:
    if mode not in MODES:
        raise ValueError(f"Progression-risk mode must be one of {MODES}.")
    resolved = positive_multiplier(multiplier) if mode == "custom" else {"standard": 1.0, "risk_enriched": 8.0}[mode]
    if mode != "custom" and multiplier is not None and positive_multiplier(multiplier) != resolved:
        raise ValueError(f"{mode} progression risk must resolve to {resolved:g}.")
    payload = {
        "contractVersion": CONTRACT_VERSION,
        "version": 1,
        "mode": mode,
        "multiplier": resolved,
        "source": source if source is not None else ("Central curve reference scenario" if mode == "standard" else "Working scenario assumption"),
        "reviewStatus": review_status if review_status is not None else ("reference_scenario" if mode == "standard" else "provisional_requires_scientific_review"),
        "notes": notes,
    }
    return validate_progression_risk_adjustment(payload)


def validate_progression_risk_adjustment(payload: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(payload, Mapping):
        raise ValueError("Progression-risk adjustment must be a versioned mapping.")
    if payload.get("contractVersion") != CONTRACT_VERSION or payload.get("version") != 1:
        raise ValueError("Unsupported population TBI progression-risk adjustment contract/version.")
    mode = payload.get("mode")
    if mode not in MODES:
        raise ValueError("Unsupported progression-risk mode.")
    multiplier = positive_multiplier(payload.get("multiplier"))
    if mode != "custom" and multiplier != {"standard": 1.0, "risk_enriched": 8.0}[mode]:
        raise ValueError("Selected progression-risk mode and resolved multiplier disagree.")
    for field in ("source", "reviewStatus"):
        if not isinstance(payload.get(field), str) or not payload[field].strip():
            raise ValueError(f"Progression-risk {field} must be a non-empty string.")
    if not isinstance(payload.get("notes"), str):
        raise ValueError("Progression-risk notes must be a string.")
    return {**dict(payload), "multiplier": multiplier}


def adjustment_for_runner(runner_config: Mapping[str, Any]) -> dict[str, Any]:
    # Pre-3G explicit configurations remain the unchanged standard scenario.
    if ADJUSTMENT_FIELD not in runner_config:
        return build_progression_risk_adjustment()
    return validate_progression_risk_adjustment(runner_config[ADJUSTMENT_FIELD])
