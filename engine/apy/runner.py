from __future__ import annotations

import json
from copy import deepcopy
from typing import Any

import numpy as np
import pandas as pd

from engine.apy.calibration import calibrate_from_config
from engine.apy.config import normalise_config
from engine.apy.regimen import (
    apply_regimen_overrides,
    default_regimen_library,
    get_regimen_from_library,
    regimen_partial_efficacy,
    validate_regimen,
)
from engine.apy.simulation import simulate_one_cohort
from engine.apy.summary import summarise_numeric_rows
from engine.apy.validation import validate_config


MODEL_VERSION = "python_apy_v9_port"
BACKEND = "python"
_CALIBRATION_CACHE: dict[str, dict[str, Any]] = {}


def run_replicates(
    config: dict[str, Any] | None = None,
    n_reps: int | None = None,
    seed: int | None = None,
    n: int | None = None,
    keep_example_cohort: bool = True,
) -> dict[str, Any]:
    cfg = validate_config(normalise_config(config))
    if n is not None:
        cfg["N"] = int(n)
    if n_reps is not None:
        cfg["nReps"] = int(n_reps)
    if seed is not None:
        cfg["seed"] = int(seed)
    if int(cfg["N"]) <= 0:
        raise ValueError("N must be > 0.")
    if int(cfg["nReps"]) <= 0:
        raise ValueError("nReps must be > 0.")

    calibration_with_parameters = _cached_calibration_from_config(cfg)
    pars = calibration_with_parameters["parameters"]
    calibration = {
        key: value
        for key, value in calibration_with_parameters.items()
        if key != "parameters"
    }
    reg = _build_regimen(cfg)
    parent_rng = np.random.default_rng(int(cfg["seed"]))
    replicate_seeds = parent_rng.integers(0, np.iinfo(np.uint32).max, size=int(cfg["nReps"]))

    raw_rows = []
    example_cohort = None
    for idx, replicate_seed in enumerate(replicate_seeds, start=1):
        out = simulate_one_cohort(
            int(cfg["N"]),
            pars,
            reg,
            calibration,
            cfg,
            np.random.default_rng(int(replicate_seed)),
            return_cohort=keep_example_cohort and idx == 1,
        )
        row = {"rep": idx, "seed": int(replicate_seed)}
        row.update(out["raw"])
        raw_rows.append(row)
        if idx == 1:
            example_cohort = out["cohort"]

    raw = pd.DataFrame(raw_rows)
    summary = summarise_numeric_rows(raw)
    return {
        "raw": raw,
        "summary": summary,
        "exampleCohort": example_cohort,
        "parameters": pars,
        "calibration": calibration,
        "strategy": build_strategy_metadata(cfg, reg),
        "interfaceConfig": cfg,
        "modelVersion": MODEL_VERSION,
        "backend": BACKEND,
    }


def run_scenario(config: dict[str, Any] | None = None) -> dict[str, Any]:
    return run_replicates(config=config)


def run_scenario_with_do_nothing(
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    from engine.apy.natural_history import run_do_nothing_summary
    from engine.apy.results_bundle import build_results_bundle

    results = run_scenario(config)
    do_nothing = run_do_nothing_summary(results)
    bundle = build_results_bundle(results, do_nothing=do_nothing)
    return {
        "results": results,
        "doNothing": do_nothing,
        "bundle": bundle,
    }


def build_strategy_metadata(config: dict[str, Any], reg: dict[str, Any]) -> dict[str, Any]:
    return {
        "testType": str(config["testType"]).upper(),
        "screeningStrategy": str(config["screeningStrategy"]).lower(),
        "regimen": reg["label"],
        "screenCoverage": float(config["screenCoverage"]),
        "screenWindow": float(config["screenWindow"]),
        "followHorizon": float(config["followHorizon"]),
        "pStartTPT": _default_if_empty(config.get("pStartTPT"), 0.85),
        "regimenPComplete": reg["pComplete"],
        "regimenADRstop": reg["pADRstop"],
        "regimenEffFull": reg["effFull"],
        "N": int(config["N"]),
        "nReps": int(config["nReps"]),
        "seed": int(config["seed"]),
        "regimenMonths": reg["months"],
        "targetDoses": reg["targetDoses"],
        "partialMode": reg["partialMode"],
        "partialDoseFractionADR": _default_if_empty(
            config.get("partialDoseFractionADR"), 0.30
        ),
        "partialDoseFractionOther": _default_if_empty(
            config.get("partialDoseFractionOther"), 0.60
        ),
        "partialEfficacyAt50pct": regimen_partial_efficacy(reg, 0.50),
        "partialEfficacyAt80pct": regimen_partial_efficacy(reg, 0.80),
        "partialEfficacyAt100pct": regimen_partial_efficacy(reg, 1.00),
    }


def _build_regimen(config: dict[str, Any]) -> dict[str, Any]:
    library = default_regimen_library(config.get("partialShortCourseMode") or "threshold80")
    reg = get_regimen_from_library(config["regimen"], library)
    reg = apply_regimen_overrides(reg, config)
    validation = validate_regimen(reg)
    if not validation["isValid"]:
        raise ValueError("; ".join(validation["errors"]))
    return reg


def _default_if_empty(value: Any, default: float) -> float:
    if value is None or value == []:
        return float(default)
    return float(value)


def _cached_calibration_from_config(config: dict[str, Any]) -> dict[str, Any]:
    key = json.dumps(_calibration_key_payload(config), sort_keys=True, default=str)
    if key not in _CALIBRATION_CACHE:
        _CALIBRATION_CACHE[key] = calibrate_from_config(config)
    return deepcopy(_CALIBRATION_CACHE[key])


def _calibration_key_payload(config: dict[str, Any]) -> dict[str, Any]:
    fields = [
        "csvFile",
        "ltbiPrevalence",
        "activeTBPrevalence",
        "targetAgeOR",
        "ageDistributionFile",
        "ageDistributionTable",
        "age85PlusMax",
        "ageDistributionSheet",
        "riskPrev",
        "diseaseOR",
        "screenWindow",
        "earlyLateRatio",
    ]
    return {field: config.get(field) for field in fields}
