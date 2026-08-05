from __future__ import annotations

import json
from copy import deepcopy
from typing import Any

import numpy as np
import pandas as pd

from engine.apy.calibration import calibrate_from_config
from engine.apy.config import normalise_config
from engine.apy.event_ledger import (
    annual_records_from_event_times,
    make_bundle as make_event_ledger_bundle,
    metadata_from_config as event_metadata_from_config,
    zero_comparator_wide,
)
from engine.apy.scenario import DEFAULT_COMPARATOR, DEFAULT_INTERVENTION
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
    event_total_rows = []
    annual_frames = []
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
        _append_agent_based_ledger_rows(
            event_total_rows,
            annual_frames,
            cfg,
            idx,
            int(replicate_seed),
            out["eventLedgerData"],
        )
        if idx == 1:
            example_cohort = out["cohort"]

    raw = pd.DataFrame(raw_rows)
    annual = pd.concat(annual_frames, ignore_index=True) if annual_frames else pd.DataFrame()
    event_ledger = make_event_ledger_bundle(
        metadata=event_metadata_from_config(
            cfg,
            model_type="agent_based",
            backend=BACKEND,
            model_version=MODEL_VERSION,
        ),
        replicate_totals_wide=pd.DataFrame(event_total_rows),
        annual_events=annual,
    )
    summary = summarise_numeric_rows(raw)
    return {
        "raw": raw,
        "summary": summary,
        "eventLedger": event_ledger,
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


def _append_agent_based_ledger_rows(
    total_rows: list[dict[str, Any]],
    annual_frames: list[pd.DataFrame],
    config: dict[str, Any],
    replicate_id: int,
    replicate_seed: int,
    ledger_data: dict[str, Any],
) -> None:
    base = _ledger_base(config, replicate_id, replicate_seed)
    population = int(config["N"])
    follow_horizon = float(config["followHorizon"])
    comparator_active_times = ledger_data["comparatorActiveTimes"]
    comparator = zero_comparator_wide(
        base,
        population=population,
        active_tb_cases=len(comparator_active_times),
    )
    intervention = {**base, "arm": "intervention"}
    intervention.update({key: float(value) for key, value in ledger_data["interventionTotals"].items()})
    total_rows.extend([comparator, intervention])

    comparator_events = {
        "population": np.zeros(population),
        "eligible_population": np.zeros(population),
        "active_tb_cases": comparator_active_times,
    }
    intervention_events = {
        "population": np.zeros(population),
        "eligible_population": np.zeros(population),
        **ledger_data["interventionEventTimes"],
    }
    annual_frames.append(
        annual_records_from_event_times(
            base={**base, "arm": "comparator"},
            event_times=comparator_events,
            follow_horizon=follow_horizon,
        )
    )
    annual_frames.append(
        annual_records_from_event_times(
            base={**base, "arm": "intervention"},
            event_times=intervention_events,
            follow_horizon=follow_horizon,
        )
    )


def _ledger_base(config: dict[str, Any], replicate_id: int, replicate_seed: int | None) -> dict[str, Any]:
    scenario = config.get("scenario") if isinstance(config.get("scenario"), dict) else {}
    return {
        "contractVersion": "ltbi_screening_event_ledger_v1",
        "scenarioId": config.get("scenarioLabel") or scenario.get("scenarioName"),
        "scenarioVersion": scenario.get("scenarioVersion", config.get("configVersion")),
        "populationPresetId": config.get("populationPresetId"),
        "modelType": "agent_based" if replicate_seed is not None else "expected_value",
        "backend": BACKEND,
        "modelVersion": MODEL_VERSION,
        "arm": "",
        "comparator": DEFAULT_COMPARATOR,
        "intervention": DEFAULT_INTERVENTION,
        "replicateId": int(replicate_id),
        "pairedReplicateId": int(replicate_id),
        "replicateSeed": replicate_seed,
        "valueType": "simulated_count" if replicate_seed is not None else "expected",
        "screeningWindow": float(config["screenWindow"]),
        "followUpHorizon": float(config["followHorizon"]),
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
