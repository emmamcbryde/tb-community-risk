from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from typing import Any, Callable

import pandas as pd

from adapters.serialization import to_json_like
from engine.apy.config import normalise_config
from engine.apy.evidence import assess_apy_reference_readiness
from engine.apy.event_ledger import EVENT_LEDGER_CONTRACT_VERSION
from engine.apy.event_ledger_economics import (
    HEALTH_ECONOMICS_CONTRACT_VERSION,
    PRIMARY_DISCOUNT_RATE,
    run_event_ledger_health_economics,
)
from engine.apy.expected_value import run_expected_value
from engine.apy.runner import MODEL_VERSION, EXPECTED_VALUE_MODEL_VERSION, run_replicates
from engine.apy.validation import validate_config


DECISION_ANALYSIS_CONTRACT_VERSION = "apy_clinician_decision_analysis_v1"


_AllowedAdapter = Callable[[dict[str, Any], Any], None]


def run_scenario_comparison(
    base_config: dict[str, Any],
    economics_config: dict[str, Any] | None,
    scenario_definitions: list[dict[str, Any]],
    *,
    model_type: str,
    n_reps: int | None = None,
    master_seed: int | None = None,
) -> dict[str, Any]:
    """Run readiness-aware APY strategy comparisons over validated adapters."""

    if model_type not in {"expected_value", "agent_based"}:
        raise ValueError("model_type must be 'expected_value' or 'agent_based'.")
    if not scenario_definitions:
        raise ValueError("At least one scenario definition is required.")

    base = normalise_config(base_config)
    scenario_results = []
    for index, definition in enumerate(scenario_definitions):
        cfg, changed, inherited = build_scenario_config(base, definition)
        if n_reps is not None:
            cfg["nReps"] = int(n_reps)
        if master_seed is not None:
            cfg["seed"] = int(master_seed)
        cfg["scenarioLabel"] = definition.get("label") or definition["scenarioId"]
        cfg.setdefault("scenario", {})
        if isinstance(cfg["scenario"], dict):
            cfg["scenario"]["scenarioName"] = definition["scenarioId"]
        model_result = _run_model(cfg, model_type)
        econ_result = _run_economics_if_possible(model_result, economics_config, cfg)
        readiness = assess_apy_reference_readiness(cfg, economics_config or {})
        scenario_results.append(
            {
                "scenarioId": definition["scenarioId"],
                "label": definition.get("label", definition["scenarioId"]),
                "description": definition.get("description", ""),
                "baseScenarioId": definition.get("baseScenarioId", base.get("scenarioLabel", "base")),
                "changedFields": changed,
                "unchangedInheritedFields": inherited,
                "modelType": model_type,
                "valueType": "expected" if model_type == "expected_value" else "simulated_count",
                "evidenceReadiness": readiness,
                "provisional": _is_provisional(model_result, econ_result, readiness),
                "eventLedgerVersion": EVENT_LEDGER_CONTRACT_VERSION,
                "economicResultsVersion": (
                    econ_result.get("contractVersion") if isinstance(econ_result, dict) else None
                ),
                "configurationHash": config_hash(cfg),
                "modelVersion": EXPECTED_VALUE_MODEL_VERSION if model_type == "expected_value" else MODEL_VERSION,
                "seed": cfg.get("seed") if model_type == "agent_based" else None,
                "nReps": cfg.get("nReps") if model_type == "agent_based" else None,
                "eventLedger": model_result["eventLedger"],
                "economics": econ_result,
                "metrics": _scenario_metrics(model_result["eventLedger"], econ_result),
                "cohortFingerprint": _cohort_fingerprint(model_result),
            }
        )

    comparisons = _paired_strategy_comparisons(scenario_results)
    validation = _validate_scenario_comparison(scenario_results, comparisons)
    return {
        "contractVersion": DECISION_ANALYSIS_CONTRACT_VERSION,
        "analysisType": "scenario_comparison",
        "modelType": model_type,
        "metadata": {
            "scenarioCount": len(scenario_results),
            "eventLedgerContractVersion": EVENT_LEDGER_CONTRACT_VERSION,
            "economicContractVersion": HEALTH_ECONOMICS_CONTRACT_VERSION,
        },
        "scenarios": scenario_results,
        "scenarioSummaries": [_summary_record(row) for row in scenario_results],
        "pairedComparisons": comparisons,
        "validation": validation,
    }


def build_scenario_config(
    base_config: dict[str, Any],
    definition: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    if "scenarioId" not in definition:
        raise ValueError("Scenario definition requires scenarioId.")
    cfg = normalise_config(deepcopy(base_config))
    changes = deepcopy(definition.get("changes") or {})
    changed: dict[str, Any] = {}
    for field, value in changes.items():
        adapter = _ALLOWED_SCENARIO_ADAPTERS.get(field)
        if adapter is None:
            raise ValueError(f"Unsupported or unvalidated scenario change: {field}")
        adapter(cfg, value)
        changed[field] = value
    cfg = validate_config(normalise_config(cfg))
    inherited = {
        field: cfg.get(field)
        for field in [
            "populationPresetId",
            "N",
            "testType",
            "regimen",
            "screeningStrategy",
            "screenCoverage",
            "screeningWindowYears",
            "followUpHorizonYears",
            "pStartTPT",
        ]
        if field not in changed
    }
    return cfg, changed, inherited


def config_hash(config: dict[str, Any]) -> str:
    payload = json.dumps(to_json_like(config), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _run_model(config: dict[str, Any], model_type: str) -> dict[str, Any]:
    if model_type == "expected_value":
        return run_expected_value(config)
    return run_replicates(config, keep_example_cohort=True)


def _run_economics_if_possible(
    model_result: dict[str, Any],
    economics_config: dict[str, Any] | None,
    scenario_config: dict[str, Any],
) -> dict[str, Any] | None:
    if economics_config is None:
        return None
    econ = deepcopy(economics_config)
    econ.setdefault("scenarioConfig", {})
    econ["scenarioConfig"].update(
        {
            "testType": scenario_config.get("testType"),
            "regimen": scenario_config.get("regimen"),
        }
    )
    try:
        return run_event_ledger_health_economics(model_result, econ)
    except Exception as exc:
        return {
            "available": False,
            "error": str(exc),
            "contractVersion": HEALTH_ECONOMICS_CONTRACT_VERSION,
        }


def _scenario_metrics(ledger: dict[str, Any], economics: dict[str, Any] | None) -> dict[str, Any]:
    metrics = _event_metrics(ledger)
    if economics and economics.get("available") is True and isinstance(economics.get("replicateResults"), pd.DataFrame):
        reps = economics["replicateResults"]
        primary = reps[reps["discountRate"] == PRIMARY_DISCOUNT_RATE]
        complete = primary[primary["economicPairComplete"].astype(bool)] if not primary.empty else primary
        if not complete.empty:
            metrics.update(
                {
                    "comparatorCost": _mean(complete, "comparatorCost"),
                    "interventionCost": _mean(complete, "interventionCost"),
                    "incrementalCost": _mean(complete, "incrementalCost"),
                    "comparatorDALYs": _mean(complete, "comparatorDALYs"),
                    "interventionDALYs": _mean(complete, "interventionDALYs"),
                    "dalysAverted": _mean(complete, "dalysAverted"),
                    "icerClassification": _first_non_empty(complete, "classification"),
                    "interpretableICER": _primary_icer(economics),
                    "nmb": _mean(complete, "netMonetaryBenefit"),
                    "probabilityPositiveNMB": _probability_positive_nmb(economics),
                }
            )
        else:
            metrics.update(
                {
                    "comparatorCost": None,
                    "interventionCost": None,
                    "incrementalCost": None,
                    "comparatorDALYs": None,
                    "interventionDALYs": None,
                    "dalysAverted": None,
                    "icerClassification": "incomplete / not calculated",
                    "interpretableICER": None,
                    "nmb": None,
                    "probabilityPositiveNMB": None,
                }
            )
    return metrics


def _event_metrics(ledger: dict[str, Any]) -> dict[str, Any]:
    totals = ledger.get("replicateTotals")
    frame = totals if isinstance(totals, pd.DataFrame) else pd.DataFrame(totals or [])
    out: dict[str, Any] = {}
    for event in [
        "screened",
        "test_positive_total",
        "true_positive_recent",
        "true_positive_remote",
        "false_positive",
        "tpt_started_total",
        "tpt_completed_total",
        "infection_effectively_treated_total",
        "active_tb_cases_prevented",
    ]:
        out[event] = _mean_event(frame, "intervention", event)
    out["comparator_active_tb"] = _mean_event(frame, "comparator", "active_tb_cases")
    out["intervention_active_tb"] = _mean_event(frame, "intervention", "active_tb_cases")
    return out


def _paired_strategy_comparisons(scenarios: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if len(scenarios) < 2:
        return []
    base = scenarios[0]
    rows = []
    for other in scenarios[1:]:
        rows.append(
            {
                "baseScenarioId": base["scenarioId"],
                "comparisonScenarioId": other["scenarioId"],
                "modelType": other["modelType"],
                "pairedCohortFingerprintMatch": (
                    base.get("cohortFingerprint") == other.get("cohortFingerprint")
                    if other["modelType"] == "agent_based"
                    else None
                ),
                "deltaActiveTBCasesPrevented": _subtract(
                    other["metrics"].get("active_tb_cases_prevented"),
                    base["metrics"].get("active_tb_cases_prevented"),
                ),
                "deltaIncrementalCost": _subtract(
                    other["metrics"].get("incrementalCost"),
                    base["metrics"].get("incrementalCost"),
                ),
                "deltaDALYsAverted": _subtract(
                    other["metrics"].get("dalysAverted"),
                    base["metrics"].get("dalysAverted"),
                ),
                "deltaNMB": _subtract(other["metrics"].get("nmb"), base["metrics"].get("nmb")),
            }
        )
    return rows


def _validate_scenario_comparison(
    scenarios: list[dict[str, Any]],
    comparisons: list[dict[str, Any]],
) -> dict[str, Any]:
    errors = []
    warnings = []
    for scenario in scenarios:
        ledger_validation = scenario["eventLedger"].get("validation") or {}
        if ledger_validation.get("isValid") is not True:
            errors.append({"field": scenario["scenarioId"], "message": "Event ledger is not valid."})
        economics = scenario.get("economics")
        if economics and economics.get("available") is True:
            econ_validation = economics.get("validation") or {}
            if econ_validation.get("structurallyValid") is not True:
                errors.append({"field": scenario["scenarioId"], "message": "Economic result is not structurally valid."})
            if econ_validation.get("economicallyComplete") is not True:
                warnings.append({"field": scenario["scenarioId"], "message": "Economic result is incomplete."})
    for row in comparisons:
        if row.get("modelType") == "agent_based" and row.get("pairedCohortFingerprintMatch") is False:
            errors.append(
                {
                    "field": "pairedCohortFingerprint",
                    "message": "Stochastic strategy comparison does not preserve the same baseline cohort fingerprint.",
                    "comparisonScenarioId": row.get("comparisonScenarioId"),
                }
            )
    return {"isValid": not errors, "errors": errors, "warnings": warnings}


def _summary_record(scenario: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "screened",
        "test_positive_total",
        "true_positive_recent",
        "true_positive_remote",
        "false_positive",
        "tpt_started_total",
        "tpt_completed_total",
        "infection_effectively_treated_total",
        "active_tb_cases_prevented",
        "comparator_active_tb",
        "intervention_active_tb",
        "incrementalCost",
        "dalysAverted",
        "icerClassification",
        "interpretableICER",
        "nmb",
        "probabilityPositiveNMB",
    ]
    return {
        "scenarioId": scenario["scenarioId"],
        "label": scenario["label"],
        "modelType": scenario["modelType"],
        "valueType": scenario["valueType"],
        "provisional": scenario["provisional"],
        "configurationHash": scenario["configurationHash"],
        **{key: scenario["metrics"].get(key) for key in keys},
    }


def _cohort_fingerprint(result: dict[str, Any]) -> str | None:
    cohort = result.get("exampleCohort")
    if not isinstance(cohort, pd.DataFrame) or cohort.empty:
        return None
    columns = [
        "id",
        "ageYears",
        "BCG",
        "MJ",
        "contact",
        "renal",
        "diabetes",
        "smoking",
        "chronicLungDisease",
        "alcoholDrugs",
        "infected",
        "recentLTBIAtBaseline",
        "remoteLTBIAtBaseline",
        "tActiveUntreated",
    ]
    present = [col for col in columns if col in cohort.columns]
    payload = cohort[present].to_json(orient="split", double_precision=12)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _is_provisional(
    model_result: dict[str, Any],
    econ_result: dict[str, Any] | None,
    readiness: dict[str, Any],
) -> bool:
    ledger_meta = model_result.get("eventLedger", {}).get("metadata", {})
    return bool(
        ledger_meta.get("ltbiStateProvisional")
        or (econ_result or {}).get("metadata", {}).get("isProvisional")
        or not readiness.get("overallClinicianReady")
    )


def _mean_event(frame: pd.DataFrame, arm: str, event: str) -> float | None:
    if frame.empty:
        return None
    rows = frame[(frame["arm"] == arm) & (frame["eventName"] == event)]
    if rows.empty:
        return None
    return float(pd.to_numeric(rows["value"], errors="coerce").mean())


def _mean(frame: pd.DataFrame, column: str) -> float | None:
    if column not in frame or frame.empty:
        return None
    values = pd.to_numeric(frame[column], errors="coerce").dropna()
    return None if values.empty else float(values.mean())


def _first_non_empty(frame: pd.DataFrame, column: str) -> Any:
    if column not in frame or frame.empty:
        return None
    for value in frame[column]:
        if value not in (None, ""):
            return value
    return None


def _primary_icer(economics: dict[str, Any]) -> float | None:
    summaries = economics.get("summaries")
    frame = summaries if isinstance(summaries, pd.DataFrame) else pd.DataFrame(summaries or [])
    rows = frame[
        (frame.get("discountRate") == PRIMARY_DISCOUNT_RATE)
        & (frame.get("metric") == "primaryICER_ratioOfMeans")
    ] if not frame.empty else frame
    if rows.empty:
        return None
    value = rows.iloc[0].get("mean")
    return None if pd.isna(value) else float(value) if value is not None else None


def _probability_positive_nmb(economics: dict[str, Any]) -> float | None:
    summaries = economics.get("summaries")
    frame = summaries if isinstance(summaries, pd.DataFrame) else pd.DataFrame(summaries or [])
    rows = frame[
        (frame.get("discountRate") == PRIMARY_DISCOUNT_RATE)
        & (frame.get("metric") == "probabilityPositiveNMB_fixedParameterSimulation")
    ] if not frame.empty else frame
    if rows.empty:
        return None
    value = rows.iloc[0].get("mean")
    return None if pd.isna(value) else float(value) if value is not None else None


def _subtract(a: Any, b: Any) -> float | None:
    try:
        if a is None or b is None or pd.isna(a) or pd.isna(b):
            return None
        return float(a) - float(b)
    except (TypeError, ValueError):
        return None


def _set(cfg: dict[str, Any], field: str, value: Any) -> None:
    cfg[field] = value


def _set_test(cfg: dict[str, Any], value: Any) -> None:
    test = str(value).upper()
    if test not in {"IGRA", "TST"}:
        raise ValueError("test must be IGRA or TST.")
    cfg["testType"] = test


def _set_regimen(cfg: dict[str, Any], value: Any) -> None:
    regimen = str(value).upper()
    if regimen not in {"3HP", "4R", "3HR", "6H", "9H"}:
        raise ValueError("Unsupported preventive regimen.")
    cfg["regimen"] = regimen


def _set_probability(cfg: dict[str, Any], field: str, value: Any) -> None:
    value = float(value)
    if not 0 <= value <= 1:
        raise ValueError(f"{field} must be in [0,1].")
    cfg[field] = value


def _set_test_characteristic(cfg: dict[str, Any], field: str, value: Any) -> None:
    _set_probability(cfg, field, value)
    cfg.setdefault("testCharacteristics", {})
    cfg["testCharacteristics"].setdefault("IGRA", {})
    cfg["testCharacteristics"].setdefault("TST", {})
    if field == "testSensitivity":
        cfg["testCharacteristics"]["IGRA"]["sensitivity"] = float(value)
    elif field == "testSpecificity":
        cfg["testCharacteristics"]["IGRA"]["specificity"] = float(value)
    elif field == "tstSensitivity":
        cfg["testCharacteristics"]["TST"]["sensitivity"] = float(value)
    elif field == "tstSpecificityNoBCG":
        cfg["testCharacteristics"]["TST"]["specificity"] = float(value)
    elif field == "tstSpecificityBCG":
        cfg["testCharacteristics"]["TST"]["specificityBCG"] = float(value)


def _set_time(cfg: dict[str, Any], field: str, value: Any) -> None:
    value = float(value)
    if value <= 0:
        raise ValueError(f"{field} must be > 0.")
    cfg[field] = value
    if field == "screeningWindowYears":
        cfg["screenWindow"] = value
    if field == "followUpHorizonYears":
        cfg["followHorizon"] = value


def _set_strategy(cfg: dict[str, Any], value: Any) -> None:
    strategy = str(value).lower()
    if strategy not in {"random", "ltbi", "cure", "prevent"}:
        raise ValueError("Unsupported screening strategy.")
    cfg["screeningStrategy"] = strategy


_ALLOWED_SCENARIO_ADAPTERS: dict[str, _AllowedAdapter] = {
    "test": _set_test,
    "testType": _set_test,
    "preventiveRegimen": _set_regimen,
    "regimen": _set_regimen,
    "screeningStrategy": _set_strategy,
    "targetingObjective": _set_strategy,
    "screeningCoverage": lambda cfg, value: _set_probability(cfg, "screenCoverage", value),
    "screenCoverage": lambda cfg, value: _set_probability(cfg, "screenCoverage", value),
    "screeningWindowYears": lambda cfg, value: _set_time(cfg, "screeningWindowYears", value),
    "followUpHorizonYears": lambda cfg, value: _set_time(cfg, "followUpHorizonYears", value),
    "treatmentInitiation": lambda cfg, value: _set_probability(cfg, "pStartTPT", value),
    "pStartTPT": lambda cfg, value: _set_probability(cfg, "pStartTPT", value),
    "testSensitivity": lambda cfg, value: _set_test_characteristic(cfg, "testSensitivity", value),
    "testSpecificity": lambda cfg, value: _set_test_characteristic(cfg, "testSpecificity", value),
    "tstSensitivity": lambda cfg, value: _set_test_characteristic(cfg, "tstSensitivity", value),
    "tstSpecificityNoBCG": lambda cfg, value: _set_test_characteristic(cfg, "tstSpecificityNoBCG", value),
    "tstSpecificityBCG": lambda cfg, value: _set_test_characteristic(cfg, "tstSpecificityBCG", value),
    "regimenPComplete": lambda cfg, value: _set_probability(cfg, "regimenPComplete", value),
    "regimenADRstop": lambda cfg, value: _set_probability(cfg, "regimenADRstop", value),
    "regimenEffFull": lambda cfg, value: _set_probability(cfg, "regimenEffFull", value),
}
