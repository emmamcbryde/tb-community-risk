from __future__ import annotations

import math
from copy import deepcopy
from typing import Any

import numpy as np


SUPPORTED_PARTIAL_MODES = {"none", "linear", "threshold80", "knots"}


def default_regimen_library(short_mode: str | None = "threshold80") -> dict[str, dict[str, Any]]:
    return {
        "r3HP": _make_regimen(
            "3HP",
            3,
            12,
            0.80,
            0.05,
            0.85,
            short_mode,
            [0, 10 / 12, 11 / 12, 1.00],
            [0, 0.00, 0.85, 0.85],
            11 / 12,
        ),
        "r4R": _make_regimen(
            "4R", 4, 120, 0.78, 0.02, 0.85, short_mode, [0, 0.80, 1.00], [0, 0.00, 0.85], 0.80
        ),
        "r3HR": _make_regimen(
            "3HR", 3, 90, 0.80, 0.03, 0.85, short_mode, [0, 0.80, 1.00], [0, 0.00, 0.85], 0.80
        ),
        "r6H": _make_regimen(
            "6H", 6, 180, 0.67, 0.03, 0.65, "knots", [0, 0.50, 1.00], [0, 0.30, 0.65], math.nan
        ),
        "r9H": _make_regimen(
            "9H",
            9,
            270,
            0.60,
            0.04,
            0.85,
            "knots",
            [0, 90 / 270, 180 / 270, 1.00],
            [0, 0.30, 0.65, 0.85],
            math.nan,
        ),
    }


def regimen_field_name(label: str) -> str:
    normalised = str(label).strip().upper()
    mapping = {"3HP": "r3HP", "4R": "r4R", "3HR": "r3HR", "6H": "r6H", "9H": "r9H"}
    if normalised not in mapping:
        raise ValueError(f"Unknown regimen label: {label}")
    return mapping[normalised]


def get_regimen_from_library(label: str, library: dict[str, dict[str, Any]]) -> dict[str, Any]:
    field = regimen_field_name(label)
    if field not in library:
        raise ValueError(f"Regimen {label} not found in regimenLibrary.")
    return deepcopy(library[field])


def merge_regimen_library(
    library: dict[str, dict[str, Any]], custom_library: dict[str, dict[str, Any]] | None
) -> dict[str, dict[str, Any]]:
    merged = deepcopy(library)
    if custom_library:
        merged.update(deepcopy(custom_library))
    return merged


def apply_regimen_overrides(reg: dict[str, Any], config: dict[str, Any]) -> dict[str, Any]:
    out = deepcopy(reg)
    if config.get("regimenPComplete") is not None:
        out["pComplete"] = config["regimenPComplete"]
    if config.get("regimenADRstop") is not None:
        out["pADRstop"] = config["regimenADRstop"]
    if config.get("regimenEffFull") is not None:
        out["effFull"] = config["regimenEffFull"]
    if config.get("pCompleteTPT") is not None:
        out["pComplete"] = config["pCompleteTPT"]
    if config.get("pSterilise") is not None:
        out["effFull"] = config["pSterilise"]
    if config.get("tptDuration") is not None:
        out["months"] = 12 * float(config["tptDuration"])
        out["targetDoses"] = round(30 * out["months"])
    return out


def validate_regimen(reg: dict[str, Any]) -> dict[str, Any]:
    errors = []
    required = [
        "label",
        "months",
        "targetDoses",
        "pComplete",
        "pADRstop",
        "effFull",
        "partialMode",
        "partialDoseKnots",
        "partialEffKnots",
        "partialThreshold",
    ]
    for field in required:
        if field not in reg:
            errors.append(f'Regimen definition is missing field "{field}".')
    if errors:
        return {"isValid": False, "errors": errors}

    if reg["months"] <= 0 or reg["targetDoses"] <= 0:
        errors.append("Regimen months and targetDoses must be > 0.")
    if not 0 <= reg["pComplete"] <= 1:
        errors.append("Regimen pComplete must be in [0,1].")
    if not 0 <= reg["pADRstop"] <= 1:
        errors.append("Regimen pADRstop must be in [0,1].")
    if not 0 <= reg["effFull"] <= 1:
        errors.append("Regimen effFull must be in [0,1].")
    if reg["pComplete"] > (1 - reg["pADRstop"] + 1e-12):
        errors.append(
            "Regimen pComplete exceeds the maximum possible completion after "
            "accounting for treatment-limiting ADR."
        )
    if str(reg["partialMode"]).lower() not in SUPPORTED_PARTIAL_MODES:
        errors.append("Regimen partialMode is not supported.")
    if not reg["partialDoseKnots"] or not reg["partialEffKnots"]:
        errors.append("Regimen partialDoseKnots and partialEffKnots cannot be empty.")
    elif len(reg["partialDoseKnots"]) != len(reg["partialEffKnots"]):
        errors.append("Regimen partialDoseKnots and partialEffKnots must have the same length.")
    return {"isValid": not errors, "errors": errors}


def regimen_partial_efficacy(
    reg: dict[str, Any], dose_fraction: float | list[float]
) -> float | list[float]:
    is_scalar = np.isscalar(dose_fraction)
    frac = np.clip(np.asarray(dose_fraction, dtype=float), 0, 1)
    mode = str(reg["partialMode"]).strip().lower()
    if mode == "none":
        eff = np.zeros_like(frac, dtype=float)
    elif mode == "linear":
        eff = float(reg["effFull"]) * frac
    elif mode == "threshold80":
        threshold = reg.get("partialThreshold", math.nan)
        if math.isnan(float(threshold)):
            threshold = 0.80
        eff = float(reg["effFull"]) * (frac >= float(threshold))
    elif mode == "knots":
        eff = np.interp(frac, reg["partialDoseKnots"], reg["partialEffKnots"])
        eff = np.clip(eff, 0, float(reg["effFull"]))
    else:
        raise ValueError(f"Unknown regimen partialMode: {reg['partialMode']}")
    return float(eff) if is_scalar else eff.tolist()


def _make_regimen(
    label: str,
    months: float,
    target_doses: float,
    p_complete: float,
    p_adr_stop: float,
    eff_full: float,
    partial_mode: str | None,
    partial_dose_knots: list[float],
    partial_eff_knots: list[float],
    partial_threshold: float,
) -> dict[str, Any]:
    return {
        "label": label,
        "months": months,
        "targetDoses": target_doses,
        "pComplete": p_complete,
        "pADRstop": p_adr_stop,
        "effFull": eff_full,
        "partialMode": partial_mode,
        "partialDoseKnots": partial_dose_knots,
        "partialEffKnots": partial_eff_knots,
        "partialThreshold": partial_threshold,
    }
