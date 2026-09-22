from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import pandas as pd

from engine.apy.costing import normalise_cost_table
from engine.apy.event_ledger import (
    EVENT_LEDGER_CONTRACT_VERSION,
    event_definitions_frame,
)
from engine.apy.event_ledger_economics import (
    COST_COMPONENTS,
    DALY_COMPONENTS,
    HEALTH_ECONOMICS_CONTRACT_VERSION,
    ID_COLS,
    INTERPRETABLE_ICER_CLASSIFICATION,
    PROGRAM_COMPONENTS,
    PRIMARY_DISCOUNT_RATE,
    _build_assumptions,
    _cost_lookup,
    _resolve_daly_inputs,
    _resolve_threshold,
    _summary_rows as _economics_summary_rows_from_result,
    classify_incremental_result,
)
from engine.apy.infection_history import configure_compatibility_reference_assumptions
from engine.apy.summary import empirical_quantile, summarise_numeric_rows


FROZEN_REFERENCE_CONTRACT_VERSION = "sa_health_frozen_reference_streamlit_v1"
SUPPORTED_ANALYSIS_BASIS = "sa_health_matlab_v9_compatibility_reference"
SUPPORTED_NATURAL_HISTORY_SEMANTICS = "matlab_v9_implicit_early_late"
FROZEN_REFERENCE_PACKAGE_ID = (
    "sa_health_apy_matlab_v9_compatible_working_reference_igra_3hp_prevent_30pct"
)
FROZEN_REFERENCE_REPS = 2000
FROZEN_REFERENCE_SEED = 1

_ROOT = Path(__file__).resolve().parents[2]
_DATA_DIR = _ROOT / "app" / "reference_data"
FROZEN_REFERENCE_MANIFEST_PATH = _DATA_DIR / "sa_health_frozen_reference_manifest.json"
FROZEN_REFERENCE_CONFIG_PATH = _DATA_DIR / "sa_health_frozen_reference_config.json"
FROZEN_REFERENCE_ECONOMICS_CONFIG_PATH = _DATA_DIR / "sa_health_frozen_reference_economics_config.json"
FROZEN_REFERENCE_TOTALS_PATH = _DATA_DIR / "sa_health_frozen_event_ledger_totals.csv.gz"
FROZEN_REFERENCE_ANNUAL_PATH = _DATA_DIR / "sa_health_frozen_event_ledger_annual.csv.gz"
FROZEN_REFERENCE_PRIMARY_REPLICATES_PATH = (
    _DATA_DIR / "sa_health_frozen_primary_economic_replicates.csv.gz"
)
FROZEN_REFERENCE_PRIMARY_ANNUAL_ECONOMICS_PATH = (
    _DATA_DIR / "sa_health_frozen_primary_economic_annual_by_arm.csv.gz"
)


_CONFIG_MATCH_KEYS = (
    "N",
    "activeTBCalibrationHorizonYears",
    "activeTBPrevalence",
    "age85PlusMax",
    "analysisMethod",
    "baselineRecentLTBIProportion",
    "diseaseOR",
    "earlyLateRatio",
    "earlyProgressionPeriodYears",
    "followHorizon",
    "followUpHorizonYears",
    "ltbiPrevalence",
    "nReps",
    "naturalHistorySemantics",
    "pStartTPT",
    "partialDoseFractionADR",
    "partialDoseFractionOther",
    "partialShortCourseMode",
    "populationPresetId",
    "recentToRemoteTransitionRatePerYear",
    "regimen",
    "regimenADRstop",
    "regimenAssumptions",
    "regimenEffFull",
    "regimenPComplete",
    "riskPrev",
    "screenCoverage",
    "screenWindow",
    "screeningStrategy",
    "screeningWindowYears",
    "seed",
    "targetAgeOR",
    "testCharacteristics",
    "testSensitivity",
    "testSpecificity",
    "testType",
    "tstSensitivity",
    "tstSpecificityBCG",
    "tstSpecificityNoBCG",
)


def is_frozen_sa_health_reference_eligible(config: dict[str, Any] | None) -> bool:
    """Return True when config exactly matches the frozen 2,000-run report scenario."""
    if not isinstance(config, dict):
        return False
    try:
        candidate = _normalised_reference_config(config)
        frozen = _normalised_reference_config(frozen_reference_config())
    except Exception:
        return False
    return _comparison_payload(candidate) == _comparison_payload(frozen)


def frozen_reference_config() -> dict[str, Any]:
    config = _read_json(FROZEN_REFERENCE_CONFIG_PATH)
    return _normalised_reference_config(config)


def frozen_reference_economics_config() -> dict[str, Any]:
    return deepcopy(_load_payload()["economicsConfig"])


def load_frozen_reference_results() -> dict[str, Any]:
    """Load the frozen SA Health reference result, event ledger and economics config.

    The returned object contains DataFrames and is intended for in-process use by
    Streamlit, not JSON export.
    """
    payload = _load_payload()
    return {
        "resultsBundle": _copy_bundle(payload["resultsBundle"]),
        "economicsConfig": deepcopy(payload["economicsConfig"]),
        "referenceEconomics": _copy_economics(payload["referenceEconomics"]),
        "manifest": deepcopy(payload["manifest"]),
    }


def validate_stochastic_replicates(
    replicate_results: Any,
    *,
    expected_reps: int,
    configuration_hash: str | None = None,
) -> dict[str, Any]:
    frame = replicate_results.copy() if isinstance(replicate_results, pd.DataFrame) else pd.DataFrame(replicate_results or [])
    if "discountProfile" in frame.columns:
        frame = frame[frame["discountProfile"].astype(str).eq("primary")].copy()
    errors: list[str] = []
    required = {
        "replicateId",
        "incrementalCost",
        "dalysAverted",
        "activeTBCasesPrevented",
        "comparatorActiveTBCases",
        "interventionActiveTBCases",
        "configurationHash",
        "economicConfigurationHash",
        "analysisBasis",
        "naturalHistorySemantics",
        "seed",
        "nReps",
        "replicateContractVersion",
        "economicPairComplete",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        errors.append(f"Missing replicate fields: {', '.join(missing)}")
    if not frame.empty and "replicateId" in frame.columns:
        ids = pd.to_numeric(frame["replicateId"], errors="coerce")
        if ids.isna().any():
            errors.append("Replicate identifiers must be present and numeric.")
        if int(ids.nunique(dropna=True)) != len(frame):
            errors.append("Replicate identifiers must be unique.")
    if len(frame) != int(expected_reps):
        errors.append(f"Expected {int(expected_reps)} primary replicate rows, found {len(frame)}.")
    for column in [
        "incrementalCost",
        "dalysAverted",
        "activeTBCasesPrevented",
        "comparatorActiveTBCases",
        "interventionActiveTBCases",
    ]:
        if column in frame.columns:
            values = pd.to_numeric(frame[column], errors="coerce")
            if values.isna().any():
                errors.append(f"{column} contains missing or non-numeric values.")
            elif not values.map(lambda item: math.isfinite(float(item))).all():
                errors.append(f"{column} contains non-finite values.")
    if "economicPairComplete" in frame.columns:
        complete = frame["economicPairComplete"].map(_boolish)
        if not bool(complete.all()):
            errors.append("All primary replicate economic pairs must be complete.")
    if configuration_hash and "configurationHash" in frame.columns:
        hashes = {str(item) for item in frame["configurationHash"].dropna().unique()}
        if hashes and hashes != {configuration_hash}:
            errors.append("Replicate configuration hashes do not match the completed analysis.")
    if "analysisBasis" in frame.columns:
        bases = {str(item) for item in frame["analysisBasis"].dropna().unique()}
        if bases and bases != {SUPPORTED_ANALYSIS_BASIS}:
            errors.append("Replicate analysis basis is unsupported or mixed.")
    if "naturalHistorySemantics" in frame.columns:
        semantics = {str(item) for item in frame["naturalHistorySemantics"].dropna().unique()}
        if semantics and semantics != {SUPPORTED_NATURAL_HISTORY_SEMANTICS}:
            errors.append("Replicate natural-history semantics are unsupported or mixed.")
    return {
        "isValid": not errors,
        "errors": errors,
        "rowCount": int(len(frame)),
        "uniqueReplicateIds": int(frame["replicateId"].nunique()) if "replicateId" in frame.columns and not frame.empty else 0,
    }


def recalculate_frozen_reference_economics(
    results_bundle: dict[str, Any] | None,
    econ_config: dict[str, Any] | None,
) -> dict[str, Any] | None:
    """Fast cost-only recalculation for the frozen SA Health reference ledger."""
    if not _is_frozen_reference_bundle(results_bundle):
        return None
    econ_config = deepcopy(econ_config or {})
    payload = _load_payload()
    if not _is_primary_default_discount(econ_config):
        return None
    if not _daly_assumptions_recalculation_equivalent(
        (econ_config or {}).get("dalyAssumptions"),
        (payload["economicsConfig"] or {}).get("dalyAssumptions"),
    ):
        return None

    ledger_metadata = deepcopy(
        (((results_bundle or {}).get("technical") or {}).get("eventLedger") or {}).get("metadata") or {}
    )
    assumptions = _build_assumptions(econ_config, ledger_metadata)
    unresolved: list[dict[str, Any]] = []
    cost_items = normalise_cost_table(econ_config.get("costItems") or [])
    costs = _cost_lookup(cost_items, econ_config, ledger_metadata, assumptions, unresolved)
    daly_inputs = _resolve_daly_inputs(assumptions["daly"], unresolved)
    threshold = _resolve_threshold(econ_config.get("threshold") or {}, assumptions, unresolved)
    annual = payload["referenceEconomics"]["annualByArm"].copy()
    if annual.empty:
        return None
    _recalculate_cost_columns(annual, costs, ledger_metadata)
    _recalculate_totals(annual)
    configuration_hash = ledger_metadata.get("configurationHash")
    replicate_results = _fast_replicate_results(
        annual,
        threshold,
        configuration_hash=str(configuration_hash or ""),
        economic_configuration_hash=_hash_json(econ_config),
    )
    summaries = _fast_summaries(replicate_results, threshold)
    validation = _fast_economic_validation(replicate_results)
    result = {
        "available": True,
        "source": "frozen_reference_vectorized_cost_recalculation",
        "contractVersion": HEALTH_ECONOMICS_CONTRACT_VERSION,
        "metadata": {
            "economicContractVersion": HEALTH_ECONOMICS_CONTRACT_VERSION,
            "eventLedgerContractVersion": ledger_metadata.get("contractVersion"),
            "scenarioId": ledger_metadata.get("scenarioId"),
            "modelType": ledger_metadata.get("modelType"),
            "valueType": "simulated_count",
            "perspective": assumptions["metadata"].get("perspective"),
            "targetCurrency": assumptions["metadata"].get("targetCurrency", "AUD"),
            "targetPriceYear": assumptions["metadata"].get("targetPriceYear", "2019"),
            "economicHorizonYears": assumptions["metadata"].get("economicHorizonYears"),
            "scopeStatement": ledger_metadata.get("scopeStatement"),
            "primaryDiscountRate": PRIMARY_DISCOUNT_RATE,
            "primaryCostDiscountRate": PRIMARY_DISCOUNT_RATE,
            "primaryHealthDiscountRate": PRIMARY_DISCOUNT_RATE,
            "isProvisional": True,
            "workingDefault": (econ_config.get("metadata") or {}).get("workingDefault", True),
            "referenceStatus": (econ_config.get("metadata") or {}).get("referenceStatus", ""),
            "conclusionPermitted": False,
        },
        "assumptions": {**assumptions, "daly": daly_inputs},
        "costItems": cost_items,
        "annualByArm": annual,
        "replicateResults": replicate_results,
        "summaries": summaries,
        "validation": validation,
        "unresolvedInputs": unresolved,
        "warnings": ["Inputs are provisional; no clinician-ready cost-effectiveness conclusion is produced."],
        "provenance": {
            "epidemiologySource": "frozen SA Health paired event ledger",
            "costSource": "vectorized recalculation from current economic costItems",
            "dalySource": "frozen report DALY assumptions; DALYs unchanged for cost-only edits",
        },
    }
    _attach_fast_legacy_compatibility_fields(result, econ_config, costs)
    result["summaryRows"] = _economics_summary_rows_from_result(result)
    result["summaryTable"] = result["summaryRows"]
    result["status"] = {
        "isComplete": bool(validation.get("economicallyComplete")),
        "missingInputs": [item.get("field") for item in unresolved],
        "notCalculated": [],
        "messages": result["warnings"],
        "validationReport": validation,
    }
    return result


@lru_cache(maxsize=1)
def _load_payload() -> dict[str, Any]:
    manifest = _read_json(FROZEN_REFERENCE_MANIFEST_PATH)
    config = frozen_reference_config()
    economics_config = _read_json(FROZEN_REFERENCE_ECONOMICS_CONFIG_PATH)
    totals = _read_ledger_csv(FROZEN_REFERENCE_TOTALS_PATH, manifest)
    annual = _read_ledger_csv(FROZEN_REFERENCE_ANNUAL_PATH, manifest)
    replicate_economics = pd.read_csv(FROZEN_REFERENCE_PRIMARY_REPLICATES_PATH)
    annual_economics = pd.read_csv(FROZEN_REFERENCE_PRIMARY_ANNUAL_ECONOMICS_PATH)
    for column in ("replicateId", "pairedReplicateId", "replicateSeed"):
        if column in replicate_economics.columns:
            replicate_economics[column] = pd.to_numeric(replicate_economics[column], errors="coerce").astype("Int64")
    for column in ("costPairComplete", "dalyPairComplete", "economicPairComplete"):
        if column in replicate_economics.columns:
            replicate_economics[column] = replicate_economics[column].map(_boolish)
    metadata = _metadata(config, manifest)
    _attach_replicate_provenance(
        replicate_economics,
        totals,
        metadata,
        economics_config,
    )
    ledger = {
        "metadata": metadata,
        "definitions": event_definitions_frame(),
        "replicateTotals": totals,
        "annualEvents": annual,
        "validation": {
            "isValid": True,
            "errors": [],
            "warnings": [],
            "source": "validated frozen SA Health reference artifact",
        },
        "summaries": _ledger_summaries(totals),
    }
    raw = _raw_from_totals(totals)
    summary = summarise_numeric_rows(raw)
    dynamic_summary = summarise_numeric_rows(_dynamic_raw_from_totals(totals))
    results_bundle = {
        "metadata": {
            "available": True,
            "modelVersion": metadata["modelVersion"],
            "backend": metadata["backend"],
            "scenarioLabel": config.get("scenarioLabel"),
            "modelType": "agent_based",
            "analysisMethod": "agent_based",
            "analysisBasis": SUPPORTED_ANALYSIS_BASIS,
            "naturalHistorySemantics": SUPPORTED_NATURAL_HISTORY_SEMANTICS,
            "nReps": FROZEN_REFERENCE_REPS,
            "seed": FROZEN_REFERENCE_SEED,
            "configurationHash": manifest.get("configurationHash"),
            "contractVersion": "apy_results_bundle_v9_python_port",
            "frozenReferenceArtifact": manifest.get("artifactContractVersion"),
        },
        "headline": {
            "available": True,
            "strategy": _strategy_metadata(config),
            "calibration": {},
            "keyMetricsRows": _summary_rows(summary, key_only=True),
            "summaryRows": _summary_rows(summary, key_only=False),
        },
        "technical": {
            "available": True,
            "interfaceConfig": config,
            "calibration": {},
            "tableMetadata": {
                "rawRows": int(len(raw)),
                "rawColumns": list(raw.columns),
                "summaryRows": int(len(summary)),
                "summaryColumns": list(summary.columns),
                "source": "frozen SA Health reference artifact",
            },
            "dynamicComparison": {
                "available": True,
                "source": "frozen technical.eventLedger",
                "population": config.get("N"),
                "followHorizon": config.get("followUpHorizonYears", config.get("followHorizon")),
                "metricRows": _summary_rows(dynamic_summary, key_only=False),
                "missingFields": [],
                "notes": "Derived from paired frozen event-ledger rows.",
            },
            "eventLedger": ledger,
        },
        "downloads": {
            "available": False,
            "source": "frozen SA Health reference artifact loaded in-app",
        },
    }
    reference_economics = {
        "available": True,
        "source": "frozen_sa_health_reference_primary_replicates",
        "contractVersion": manifest.get("healthEconomicsContractVersion"),
        "metadata": {
            "economicContractVersion": manifest.get("healthEconomicsContractVersion"),
            "eventLedgerContractVersion": EVENT_LEDGER_CONTRACT_VERSION,
            "scenarioId": metadata.get("scenarioId"),
            "modelType": "agent_based",
            "valueType": "simulated_count",
            "primaryDiscountRate": 0.03,
            "primaryCostDiscountRate": 0.03,
            "primaryHealthDiscountRate": 0.03,
            "isProvisional": True,
            "workingDefault": True,
            "referenceStatus": "frozen_sa_health_report_reference",
        },
        "replicateResults": replicate_economics,
        "annualByArm": annual_economics,
        "costItems": normalise_cost_table(economics_config.get("costItems") or []),
        "summaryRows": _economic_summary_rows(replicate_economics),
        "summaryTable": _economic_summary_rows(replicate_economics),
        "validation": validate_stochastic_replicates(
            replicate_economics,
            expected_reps=FROZEN_REFERENCE_REPS,
        ),
        "warnings": [
            "Inputs are provisional; no clinician-ready cost-effectiveness conclusion is produced."
        ],
    }
    return {
        "manifest": manifest,
        "config": config,
        "economicsConfig": economics_config,
        "resultsBundle": results_bundle,
        "referenceEconomics": reference_economics,
    }


def _normalised_reference_config(config: dict[str, Any]) -> dict[str, Any]:
    cfg = configure_compatibility_reference_assumptions(deepcopy(config))
    cfg["analysisMethod"] = str(cfg.get("analysisMethod") or "agent_based")
    if cfg["analysisMethod"] == "agent_based":
        cfg["analysisMethodLabel"] = "SA Health report analysis - 2,000 simulated communities"
    cfg["nReps"] = int(float(cfg.get("nReps") or FROZEN_REFERENCE_REPS))
    cfg["seed"] = int(float(cfg.get("seed") or FROZEN_REFERENCE_SEED))
    cfg["screeningWindowYears"] = float(cfg.get("screeningWindowYears", cfg.get("screenWindow", 2)))
    cfg["screenWindow"] = cfg["screeningWindowYears"]
    cfg["followUpHorizonYears"] = float(cfg.get("followUpHorizonYears", cfg.get("followHorizon", 20)))
    cfg["followHorizon"] = cfg["followUpHorizonYears"]
    cfg["analysisBasis"] = SUPPORTED_ANALYSIS_BASIS
    cfg["naturalHistorySemantics"] = SUPPORTED_NATURAL_HISTORY_SEMANTICS
    return cfg


def _daly_assumptions_recalculation_equivalent(
    current: dict[str, Any] | None,
    reference: dict[str, Any] | None,
) -> bool:
    """Compare only DALY inputs that change the frozen recalculation.

    The Streamlit editable workspace can round-trip source notes and citations
    through user-facing rows. Those text fields do not affect the frozen
    vectorized cost recalculation, but a byte-for-byte comparison would reject
    the fast path after a purely economic edit and force a slow generic
    recomputation. Keep the guard strict on values, inclusion flags and review
    states, while ignoring explanatory wording.
    """

    return _json_like(_daly_recalculation_payload(current or {})) == _json_like(
        _daly_recalculation_payload(reference or {})
    )


def _daly_recalculation_payload(daly: dict[str, Any]) -> dict[str, Any]:
    numeric_records = (
        "activeTBDisabilityWeight",
        "activeTBDurationYears",
        "tbCaseFatalityRisk",
        "yllPerTBDeath",
        "dalyLossPerTPTStarted",
        "dalyLossPerADRStop",
        "postTBDALYsPerActiveTBCase",
    )
    flags = (
        "method",
        "includeTPTHealthLoss",
        "includeADRHealthLoss",
        "includePostTBSequelae",
        "tptHealthLossExclusionStatus",
        "adrHealthLossExclusionStatus",
        "postTBSequelaeStatus",
    )
    payload: dict[str, Any] = {key: daly.get(key) for key in flags}
    for key in numeric_records:
        payload[key] = _assumption_recalculation_record(daly.get(key))
    return payload


def _assumption_recalculation_record(value: Any) -> Any:
    if not isinstance(value, dict):
        return value
    return {
        "value": value.get("value"),
        "status": value.get("status"),
        "provisional": value.get("provisional"),
        "unit": value.get("unit"),
    }


def _comparison_payload(config: dict[str, Any]) -> dict[str, Any]:
    return {key: _json_like(config.get(key)) for key in _CONFIG_MATCH_KEYS}


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _read_ledger_csv(path: Path, manifest: dict[str, Any]) -> pd.DataFrame:
    frame = pd.read_csv(path)
    frame["contractVersion"] = EVENT_LEDGER_CONTRACT_VERSION
    frame["scenarioId"] = "SA Health APY MATLAB-v9-compatible working-reference analysis"
    frame["scenarioVersion"] = "m1"
    frame["populationPresetId"] = "apy_demonstration"
    frame["modelType"] = "agent_based"
    frame["backend"] = "python"
    frame["modelVersion"] = "python_apy_v9_port"
    frame["comparator"] = "current practice / no additional systematic LTBI screening"
    frame["intervention"] = "targeted LTBI screening and preventive treatment"
    frame["valueType"] = "simulated_count"
    frame["screeningWindow"] = 2.0
    frame["screeningWindowYears"] = 2.0
    frame["earlyProgressionPeriodYears"] = 2.0
    frame["activeTBCalibrationHorizonYears"] = 2.0
    frame["followUpHorizon"] = 20.0
    frame["followUpHorizonYears"] = 20.0
    frame["configurationHash"] = manifest.get("configurationHash")
    frame["analysisBasis"] = SUPPORTED_ANALYSIS_BASIS
    frame["naturalHistorySemantics"] = SUPPORTED_NATURAL_HISTORY_SEMANTICS
    if "withinFollowUp" in frame.columns:
        frame["withinFollowUp"] = frame["withinFollowUp"].map(_boolish)
    for column in ("replicateId", "pairedReplicateId", "replicateSeed", "modelYear"):
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    if "value" in frame.columns:
        frame["value"] = pd.to_numeric(frame["value"], errors="coerce").fillna(0.0)
    return frame


def _attach_replicate_provenance(
    replicates: pd.DataFrame,
    totals: pd.DataFrame,
    metadata: dict[str, Any],
    economics_config: dict[str, Any],
) -> None:
    counts = _active_tb_counts_by_replicate(totals)
    if not counts.empty:
        keyed = counts.set_index("replicateId")
        replicate_ids = replicates["replicateId"]
        for column in counts.columns:
            if column == "replicateId":
                continue
            replicates[column] = replicate_ids.map(keyed[column])
    replicates["configurationHash"] = metadata.get("configurationHash")
    replicates["economicConfigurationHash"] = _hash_json(economics_config)
    replicates["analysisBasis"] = SUPPORTED_ANALYSIS_BASIS
    replicates["naturalHistorySemantics"] = SUPPORTED_NATURAL_HISTORY_SEMANTICS
    replicates["seed"] = FROZEN_REFERENCE_SEED
    replicates["nReps"] = FROZEN_REFERENCE_REPS
    replicates["replicateContractVersion"] = FROZEN_REFERENCE_CONTRACT_VERSION


def _active_tb_counts_by_replicate(totals: pd.DataFrame) -> pd.DataFrame:
    if totals.empty or not {"arm", "eventName", "replicateId", "value"}.issubset(totals.columns):
        return pd.DataFrame()
    subset = totals[totals["eventName"].astype(str).eq("active_tb_cases")].copy()
    if subset.empty:
        return pd.DataFrame()
    wide = subset.pivot_table(index="replicateId", columns="arm", values="value", aggfunc="first")
    out = pd.DataFrame(index=wide.index)
    out["comparatorActiveTBCases"] = pd.to_numeric(wide.get("comparator"), errors="coerce")
    out["interventionActiveTBCases"] = pd.to_numeric(wide.get("intervention"), errors="coerce")
    return out.reset_index()


def _metadata(config: dict[str, Any], manifest: dict[str, Any]) -> dict[str, Any]:
    return {
        "contractVersion": EVENT_LEDGER_CONTRACT_VERSION,
        "scenarioId": config.get("scenarioLabel"),
        "scenarioVersion": (config.get("scenario") or {}).get("scenarioVersion", config.get("configVersion"))
        if isinstance(config.get("scenario"), dict)
        else config.get("configVersion"),
        "populationPresetId": config.get("populationPresetId"),
        "modelType": "agent_based",
        "backend": "python",
        "analysisBasis": SUPPORTED_ANALYSIS_BASIS,
        "naturalHistorySemantics": SUPPORTED_NATURAL_HISTORY_SEMANTICS,
        "screeningWindow": config.get("screenWindow"),
        "screeningWindowYears": config.get("screeningWindowYears", config.get("screenWindow")),
        "earlyProgressionPeriodYears": config.get("earlyProgressionPeriodYears"),
        "activeTBCalibrationHorizonYears": config.get("activeTBCalibrationHorizonYears"),
        "followUpHorizon": config.get("followHorizon"),
        "followUpHorizonYears": config.get("followUpHorizonYears", config.get("followHorizon")),
        "modelVersion": "python_apy_v9_port",
        "calibrationPolicy": config.get("calibrationPolicy"),
        "scopeStatement": "Direct benefits, harms and costs only; transmission-mediated benefits are not included.",
        "comparator": "current practice / no additional systematic LTBI screening",
        "intervention": "targeted LTBI screening and preventive treatment",
        "seed": FROZEN_REFERENCE_SEED,
        "nReps": FROZEN_REFERENCE_REPS,
        "configurationHash": manifest.get("configurationHash"),
        "packageId": manifest.get("sourcePackageId") or FROZEN_REFERENCE_PACKAGE_ID,
        "frozenReferenceArtifact": manifest.get("artifactContractVersion"),
    }


def _raw_from_totals(totals: pd.DataFrame) -> pd.DataFrame:
    intervention = _wide_by_replicate(totals, "intervention")
    comparator = _wide_by_replicate(totals, "comparator")
    rows = pd.DataFrame(
        {
            "rep": intervention.index,
            "seed": intervention.get("replicateSeed"),
            "nScreened": intervention.get("screened", 0.0),
            "nTestPositiveNonActive": intervention.get("true_positive_latent", 0.0)
            + intervention.get("false_positive", 0.0),
            "nFalsePositiveTreated": intervention.get("tpt_started_false_positive", 0.0),
            "nTotalCoursesStarted": intervention.get("tpt_started_total", 0.0),
            "nTotalCoursesCompleted": intervention.get("tpt_completed_total", 0.0),
            "nADRstop": intervention.get("tpt_adr_stop_total", 0.0),
            "nCuredInfection": intervention.get("infection_effectively_treated_total", 0.0),
            "nPreventedActiveTB": intervention.get("active_tb_cases_prevented", 0.0),
            "nActiveBy20y": intervention.get("active_tb_cases", 0.0),
            "nActiveComparatorBy20y": comparator.get("active_tb_cases", 0.0),
        }
    )
    return rows.reset_index(drop=True)


def _dynamic_raw_from_totals(totals: pd.DataFrame) -> pd.DataFrame:
    intervention = _wide_by_replicate(totals, "intervention")
    comparator = _wide_by_replicate(totals, "comparator")
    comp = pd.to_numeric(comparator.get("active_tb_cases", 0.0), errors="coerce")
    inter = pd.to_numeric(intervention.get("active_tb_cases", 0.0), errors="coerce")
    prevented = pd.to_numeric(intervention.get("active_tb_cases_prevented", 0.0), errors="coerce")
    return pd.DataFrame(
        {
            "cumulative_baseline_active_tb_cases": comp,
            "cumulative_intervention_active_tb_cases": inter,
            "cumulative_cases_averted": prevented,
            "relative_reduction_cumulative_active_tb_cases": prevented / comp.replace(0, pd.NA),
        }
    )


def _wide_by_replicate(totals: pd.DataFrame, arm: str) -> pd.DataFrame:
    subset = totals[totals["arm"].astype(str).eq(arm)].copy()
    wide = subset.pivot_table(
        index="replicateId",
        columns="eventName",
        values="value",
        aggfunc="first",
    )
    seed = subset[["replicateId", "replicateSeed"]].drop_duplicates().set_index("replicateId")
    wide = wide.join(seed, how="left")
    return wide.sort_index()


def _ledger_summaries(totals: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (arm, event), group in totals.groupby(["arm", "eventName"], dropna=False):
        values = pd.to_numeric(group["value"], errors="coerce").dropna()
        rows.append(
            {
                "arm": arm,
                "eventName": event,
                "mean": float(values.mean()) if not values.empty else None,
                "median": float(values.median()) if not values.empty else None,
                "min": float(values.min()) if not values.empty else None,
                "max": float(values.max()) if not values.empty else None,
                "n": int(values.count()),
            }
        )
    return pd.DataFrame(rows)


def _summary_rows(summary: pd.DataFrame, *, key_only: bool) -> list[dict[str, Any]]:
    rows = summary.to_dict(orient="records")
    if not key_only:
        return rows
    key_metrics = {
        "nScreened",
        "nTestPositiveNonActive",
        "nTotalCoursesStarted",
        "nTotalCoursesCompleted",
        "nCuredInfection",
        "nPreventedActiveTB",
    }
    return [row for row in rows if row.get("Metric") in key_metrics]


def _economic_summary_rows(replicates: pd.DataFrame) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    metrics = [
        "comparatorCost",
        "interventionCost",
        "incrementalCost",
        "comparatorDALYs",
        "interventionDALYs",
        "dalysAverted",
        "activeTBCasesPrevented",
    ]
    rows.append(
        {
            "discountProfile": "primary",
            "discountRate": 0.03,
            "costDiscountRate": 0.03,
            "healthDiscountRate": 0.03,
            "metric": "pairedReplicateCompleteness",
            "totalPairedReplicates": int(len(replicates)),
            "completePairedReplicates": int(len(replicates)),
            "excludedPairedReplicates": 0,
            "exclusionReasons": "",
        }
    )
    mean_inc = float(pd.to_numeric(replicates["incrementalCost"], errors="coerce").mean())
    mean_daly = float(pd.to_numeric(replicates["dalysAverted"], errors="coerce").mean())
    classification = "dominant" if mean_inc < 0 and mean_daly > 0 else "not calculated"
    rows.append(
        {
            "discountProfile": "primary",
            "discountRate": 0.03,
            "costDiscountRate": 0.03,
            "healthDiscountRate": 0.03,
            "metric": "primaryICER_ratioOfMeans",
            "mean": None if classification == "dominant" else mean_inc / mean_daly,
            "classification": classification,
            "n": int(len(replicates)),
            "intervalLabel": "simulation distribution across replicates",
        }
    )
    for metric in metrics:
        values = pd.to_numeric(replicates[metric], errors="coerce").dropna()
        rows.append(
            {
                "discountProfile": "primary",
                "discountRate": 0.03,
                "costDiscountRate": 0.03,
                "healthDiscountRate": 0.03,
                "metric": metric,
                "mean": float(values.mean()),
                "classification": "",
                "n": int(values.count()),
                "intervalLabel": "simulation distribution across replicates",
                "sd": float(values.std(ddof=1)),
                "median": float(values.median()),
                "p2_5": float(empirical_quantile(values.to_numpy(), 0.025)),
                "p97_5": float(empirical_quantile(values.to_numpy(), 0.975)),
                "min": float(values.min()),
                "max": float(values.max()),
            }
        )
    return rows


def _strategy_metadata(config: dict[str, Any]) -> dict[str, Any]:
    return {
        "testType": str(config.get("testType", "")).upper(),
        "screeningStrategy": str(config.get("screeningStrategy", "")).lower(),
        "regimen": config.get("regimen"),
        "screenCoverage": config.get("screenCoverage"),
        "screeningWindowYears": config.get("screeningWindowYears"),
        "followUpHorizonYears": config.get("followUpHorizonYears"),
        "N": config.get("N"),
        "nReps": config.get("nReps"),
        "seed": config.get("seed"),
    }


def _copy_bundle(bundle: dict[str, Any]) -> dict[str, Any]:
    out = deepcopy({key: value for key, value in bundle.items() if key != "technical"})
    technical = deepcopy({key: value for key, value in bundle["technical"].items() if key != "eventLedger"})
    ledger = bundle["technical"]["eventLedger"]
    technical["eventLedger"] = {
        **deepcopy({key: value for key, value in ledger.items() if key not in {"replicateTotals", "annualEvents", "definitions", "summaries"}}),
        "replicateTotals": ledger["replicateTotals"].copy(),
        "annualEvents": ledger["annualEvents"].copy(),
        "definitions": ledger["definitions"].copy(),
        "summaries": ledger["summaries"].copy(),
    }
    out["technical"] = technical
    return out


def _copy_economics(economics: dict[str, Any]) -> dict[str, Any]:
    out = deepcopy({key: value for key, value in economics.items() if key not in {"replicateResults", "annualByArm"}})
    out["replicateResults"] = economics["replicateResults"].copy()
    if "annualByArm" in economics:
        out["annualByArm"] = economics["annualByArm"].copy()
    return out


def _is_frozen_reference_bundle(results_bundle: dict[str, Any] | None) -> bool:
    if not isinstance(results_bundle, dict):
        return False
    metadata = results_bundle.get("metadata") or {}
    technical = results_bundle.get("technical") or {}
    ledger_metadata = ((technical.get("eventLedger") or {}).get("metadata") or {}) if isinstance(technical, dict) else {}
    return (
        (metadata.get("frozenReferenceArtifact") == FROZEN_REFERENCE_CONTRACT_VERSION)
        or (ledger_metadata.get("frozenReferenceArtifact") == FROZEN_REFERENCE_CONTRACT_VERSION)
    )


def _is_primary_default_discount(econ_config: dict[str, Any]) -> bool:
    discounting = econ_config.get("discounting") or {}
    values = [
        discounting.get("primaryDisplayedRate"),
        discounting.get("selectedAnnualRate"),
        discounting.get("healthOutcomeAnnualRate"),
    ]
    for value in values:
        if value in (None, "", []):
            continue
        try:
            if abs(float(value) - PRIMARY_DISCOUNT_RATE) > 1e-12:
                return False
        except (TypeError, ValueError):
            return False
    return True


def _recalculate_cost_columns(annual: pd.DataFrame, costs: dict[str, Any], metadata: dict[str, Any]) -> None:
    q = {name: pd.to_numeric(annual.get(name, 0.0), errors="coerce").fillna(0.0) for name in [
        "screened",
        "tpt_started_total",
        "tpt_started_false_positive",
        "tpt_adr_stop_total",
        "active_tb_cases",
    ]}
    intervention = annual["arm"].astype(str).eq("intervention")
    year = pd.to_numeric(annual["modelYear"], errors="coerce").fillna(0.0)
    first_year = float(metadata.get("programRunningFirstYear") or 0)
    duration = float(metadata.get("programRunningDurationYears") or metadata.get("screeningWindowYears") or 0)
    running_qty = intervention.astype(float) * ((year >= first_year) & (year < first_year + duration)).astype(float)
    setup_qty = intervention.astype(float) * year.eq(0).astype(float)
    annual["screeningTestCost"] = _component(q["screened"], costs.get("test"))
    annual["tptRegimenCost"] = _component(q["tpt_started_total"], costs.get("regimen"))
    annual["falsePositiveIncrementalCost"] = _component(q["tpt_started_false_positive"], costs.get("false_positive"))
    annual["activeTBDiseaseCost"] = _component(q["active_tb_cases"], costs.get("active_tb"))
    annual["programSetupCost"] = _component(setup_qty, costs.get("setup"))
    annual["programRunningCost"] = _component(running_qty, costs.get("running"))
    annual["adrManagementCost"] = _component(q["tpt_adr_stop_total"], costs.get("adr"))
    annual["returnForResultsCost"] = _component(q["screened"], costs.get("return_results"))
    annual["clinicalReviewCost"] = _component(q["tpt_started_total"], costs.get("clinical_review"))
    annual["activeTBExclusionWorkupCost"] = _component(q["tpt_started_total"], costs.get("active_tb_exclusion"))
    annual["travelOutreachStaffSupportCost"] = _component(q["screened"], costs.get("travel_outreach"))


def _component(quantity: Any, unit_cost: Any) -> pd.Series:
    values = pd.Series(quantity, dtype="float64")
    if unit_cost is None:
        return values * 0.0
    return values * float(unit_cost)


def _recalculate_totals(annual: pd.DataFrame) -> None:
    for column in COST_COMPONENTS:
        if column not in annual.columns:
            annual[column] = 0.0
        annual[column] = pd.to_numeric(annual[column], errors="coerce").fillna(0.0)
    for column in DALY_COMPONENTS:
        if column not in annual.columns:
            annual[column] = 0.0
        annual[column] = pd.to_numeric(annual[column], errors="coerce").fillna(0.0)
    annual["costComplete"] = True
    annual["dalyComplete"] = True
    annual["missingCostComponents"] = ""
    annual["missingDALYComponents"] = ""
    annual["totalUndiscountedCost"] = annual[COST_COMPONENTS].sum(axis=1)
    annual["totalDiscountedCost"] = annual["totalUndiscountedCost"] * pd.to_numeric(
        annual.get("costDiscountFactor", 1.0), errors="coerce"
    ).fillna(1.0)
    annual["totalUndiscountedDALYs"] = annual[DALY_COMPONENTS].sum(axis=1)
    annual["totalDiscountedDALYs"] = annual["totalUndiscountedDALYs"] * pd.to_numeric(
        annual.get("healthDiscountFactor", 1.0), errors="coerce"
    ).fillna(1.0)


def _fast_replicate_results(
    annual: pd.DataFrame,
    threshold: float | None,
    *,
    configuration_hash: str,
    economic_configuration_hash: str,
) -> pd.DataFrame:
    """Construct paired replicate results from frozen annual economics rows.

    This mirrors the generic economics result contract, but avoids the per-row
    Python loops used by the full event-ledger path. The frozen annual artifact
    is already primary-discount-profile only and structurally complete.
    """
    included = annual[annual["includedInEconomicAnalysis"].astype(bool)].copy()
    if included.empty:
        return pd.DataFrame()
    for column in ID_COLS:
        if column not in included.columns:
            included[column] = "" if column not in {"replicateId", "pairedReplicateId", "replicateSeed"} else 0
    group_cols = ID_COLS + ["arm"]
    factor_cost = pd.to_numeric(included.get("costDiscountFactor", 1.0), errors="coerce").fillna(1.0)
    factor_health = pd.to_numeric(included.get("healthDiscountFactor", 1.0), errors="coerce").fillna(1.0)
    component_columns: list[str] = []
    for component in COST_COMPONENTS:
        values = pd.to_numeric(included.get(component, 0.0), errors="coerce").fillna(0.0)
        included[f"{component}Discounted"] = values * factor_cost
        included[f"{component}Undiscounted"] = values
        component_columns.extend([f"{component}Discounted", f"{component}Undiscounted"])
    for component in DALY_COMPONENTS:
        values = pd.to_numeric(included.get(component, 0.0), errors="coerce").fillna(0.0)
        included[f"{component}Discounted"] = values * factor_health
        included[f"{component}Undiscounted"] = values
        component_columns.extend([f"{component}Discounted", f"{component}Undiscounted"])
    programme_discounted = sum(
        pd.to_numeric(included.get(component, 0.0), errors="coerce").fillna(0.0) * factor_cost
        for component in PROGRAM_COMPONENTS
    )
    programme_undiscounted = sum(
        pd.to_numeric(included.get(component, 0.0), errors="coerce").fillna(0.0)
        for component in PROGRAM_COMPONENTS
    )
    included["totalProgrammeCostDiscounted"] = programme_discounted
    included["totalProgrammeCostUndiscounted"] = programme_undiscounted
    component_columns.extend(["totalProgrammeCostDiscounted", "totalProgrammeCostUndiscounted"])

    aggregations = {
        "totalDiscountedCost": "sum",
        "totalDiscountedDALYs": "sum",
        "totalUndiscountedCost": "sum",
        "totalUndiscountedDALYs": "sum",
        "active_tb_cases": "sum",
        "active_tb_cases_prevented": "sum",
        "infection_effectively_treated_total": "sum",
        **{column: "sum" for column in component_columns},
    }
    arms = included.groupby(group_cols, dropna=False, as_index=False).agg(aggregations)
    comp = arms[arms["arm"].astype(str).eq("comparator")].copy()
    inter = arms[arms["arm"].astype(str).eq("intervention")].copy()
    if comp.empty or inter.empty:
        return pd.DataFrame()
    pair = comp.merge(inter, on=ID_COLS, suffixes=("_comparator", "_intervention"), how="inner")
    rows = pair[ID_COLS].copy()
    rows["costPairComplete"] = True
    rows["dalyPairComplete"] = True
    rows["economicPairComplete"] = True
    rows["exclusionReasons"] = ""
    rows["comparatorCost"] = pair["totalDiscountedCost_comparator"]
    rows["interventionCost"] = pair["totalDiscountedCost_intervention"]
    rows["incrementalCost"] = rows["interventionCost"] - rows["comparatorCost"]
    rows["comparatorDALYs"] = pair["totalDiscountedDALYs_comparator"]
    rows["interventionDALYs"] = pair["totalDiscountedDALYs_intervention"]
    rows["dalysAverted"] = rows["comparatorDALYs"] - rows["interventionDALYs"]
    rows["comparatorActiveTBCases"] = pair["active_tb_cases_comparator"]
    rows["interventionActiveTBCases"] = pair["active_tb_cases_intervention"]
    rows["activeTBCasesPrevented"] = pair["active_tb_cases_prevented_intervention"]
    rows["infectionsEffectivelyTreated"] = pair["infection_effectively_treated_total_intervention"]
    rows["costPerActiveTBCasePrevented"] = _divide_series(rows["incrementalCost"], rows["activeTBCasesPrevented"])
    rows["costPerInfectionEffectivelyTreated"] = _divide_series(
        rows["incrementalCost"],
        rows["infectionsEffectivelyTreated"],
    )
    classifications = [
        classify_incremental_result(cost, dalys)
        for cost, dalys in zip(rows["incrementalCost"], rows["dalysAverted"])
    ]
    rows["classification"] = classifications
    rows["primaryICERInterpretable"] = [
        item == INTERPRETABLE_ICER_CLASSIFICATION for item in classifications
    ]
    rows["replicateICER"] = _divide_series(rows["incrementalCost"], rows["dalysAverted"])
    rows.loc[~rows["primaryICERInterpretable"].astype(bool), "replicateICER"] = None
    rows["netMonetaryBenefit"] = (
        None if threshold is None else threshold * rows["dalysAverted"] - rows["incrementalCost"]
    )
    rows["configurationHash"] = configuration_hash
    rows["economicConfigurationHash"] = economic_configuration_hash
    rows["analysisBasis"] = SUPPORTED_ANALYSIS_BASIS
    rows["naturalHistorySemantics"] = SUPPORTED_NATURAL_HISTORY_SEMANTICS
    rows["seed"] = FROZEN_REFERENCE_SEED
    rows["nReps"] = FROZEN_REFERENCE_REPS
    rows["replicateContractVersion"] = FROZEN_REFERENCE_CONTRACT_VERSION
    paired_columns: dict[str, Any] = {}
    for side, suffix in (("comparator", "comparator"), ("intervention", "intervention")):
        paired_columns[f"{side}_totalArmCostDiscounted"] = pair[f"totalDiscountedCost_{suffix}"]
        paired_columns[f"{side}_totalArmCostUndiscounted"] = pair[f"totalUndiscountedCost_{suffix}"]
        paired_columns[f"{side}_totalProgrammeCostDiscounted"] = pair[f"totalProgrammeCostDiscounted_{suffix}"]
        paired_columns[f"{side}_totalProgrammeCostUndiscounted"] = pair[f"totalProgrammeCostUndiscounted_{suffix}"]
        paired_columns[f"{side}_activeTBDiseaseCareDiscounted"] = pair[f"activeTBDiseaseCostDiscounted_{suffix}"]
        for component in COST_COMPONENTS:
            paired_columns[f"{side}_{component}Discounted"] = pair[f"{component}Discounted_{suffix}"]
            paired_columns[f"{side}_{component}Undiscounted"] = pair[f"{component}Undiscounted_{suffix}"]
        for component in DALY_COMPONENTS:
            paired_columns[f"{side}_{component}Discounted"] = pair[f"{component}Discounted_{suffix}"]
            paired_columns[f"{side}_{component}Undiscounted"] = pair[f"{component}Undiscounted_{suffix}"]
    rows = pd.concat([rows, pd.DataFrame(paired_columns, index=rows.index)], axis=1)
    return rows.astype(object).where(pd.notna(rows), None)


def _fast_summaries(replicates: pd.DataFrame, threshold: float | None) -> pd.DataFrame:
    if not isinstance(replicates, pd.DataFrame) or replicates.empty:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    summary_group_cols = ["discountProfile", "discountRate", "costDiscountRate", "healthDiscountRate"]
    metrics = [
        "comparatorCost",
        "interventionCost",
        "incrementalCost",
        "comparatorDALYs",
        "interventionDALYs",
        "dalysAverted",
        "activeTBCasesPrevented",
        "netMonetaryBenefit",
    ]
    for keys, group in replicates.groupby(summary_group_cols, dropna=False):
        profile, rate, cost_rate, health_rate = keys
        complete = group[group["economicPairComplete"].astype(bool)].copy()
        total_pairs = int(len(group))
        complete_pairs = int(len(complete))
        rows.append(
            {
                "discountProfile": profile,
                "discountRate": rate,
                "costDiscountRate": cost_rate,
                "healthDiscountRate": health_rate,
                "metric": "pairedReplicateCompleteness",
                "totalPairedReplicates": total_pairs,
                "completePairedReplicates": complete_pairs,
                "excludedPairedReplicates": total_pairs - complete_pairs,
                "exclusionReasons": "",
            }
        )
        mean_inc = float(pd.to_numeric(complete["incrementalCost"], errors="coerce").mean()) if complete_pairs else None
        mean_daly = float(pd.to_numeric(complete["dalysAverted"], errors="coerce").mean()) if complete_pairs else None
        classification = classify_incremental_result(mean_inc, mean_daly) if complete_pairs else "incomplete / not calculated"
        rows.append(
            {
                "discountProfile": profile,
                "discountRate": rate,
                "costDiscountRate": cost_rate,
                "healthDiscountRate": health_rate,
                "metric": "primaryICER_ratioOfMeans",
                "mean": _div(mean_inc, mean_daly) if classification == INTERPRETABLE_ICER_CLASSIFICATION else None,
                "classification": classification,
                "n": complete_pairs,
                "totalPairedReplicates": total_pairs,
                "completePairedReplicates": complete_pairs,
                "excludedPairedReplicates": total_pairs - complete_pairs,
                "intervalLabel": "simulation distribution across replicates",
            }
        )
        if threshold is not None and "netMonetaryBenefit" in complete:
            valid_nmb = pd.to_numeric(complete["netMonetaryBenefit"], errors="coerce").dropna()
            rows.append(
                {
                    "discountProfile": profile,
                    "discountRate": rate,
                    "costDiscountRate": cost_rate,
                    "healthDiscountRate": health_rate,
                    "metric": "probabilityPositiveNMB_fixedParameterSimulation",
                    "mean": None if valid_nmb.empty else float((valid_nmb > 0).mean()),
                    "n": int(valid_nmb.count()),
                    "numerator": int((valid_nmb > 0).sum()) if not valid_nmb.empty else 0,
                    "denominator": int(valid_nmb.count()),
                    "intervalLabel": "probability of positive NMB across finite-population simulation replicates under fixed parameter assumptions",
                }
            )
        for metric in metrics:
            if metric not in complete:
                continue
            values = pd.to_numeric(complete[metric], errors="coerce").dropna()
            if values.empty:
                continue
            rows.append(
                {
                    "discountProfile": profile,
                    "discountRate": rate,
                    "costDiscountRate": cost_rate,
                    "healthDiscountRate": health_rate,
                    "metric": metric,
                    "n": int(values.count()),
                    "mean": float(values.mean()),
                    "sd": None if len(values) == 1 else float(values.std(ddof=1)),
                    "median": float(values.median()),
                    "p2_5": float(empirical_quantile(values.to_numpy(), 0.025)),
                    "p97_5": float(empirical_quantile(values.to_numpy(), 0.975)),
                    "min": float(values.min()),
                    "max": float(values.max()),
                    "classification": "",
                    "intervalLabel": "simulation distribution across replicates",
                }
            )
    return pd.DataFrame(rows).astype(object).where(pd.notna(pd.DataFrame(rows)), None)


def _attach_fast_legacy_compatibility_fields(
    result: dict[str, Any],
    econ_config: dict[str, Any],
    costs: dict[str, Any],
) -> None:
    annual = result["annualByArm"]
    if not isinstance(annual, pd.DataFrame) or annual.empty:
        primary = pd.DataFrame()
    else:
        primary = annual.copy()
        if "discountProfile" in primary:
            primary = primary[primary["discountProfile"].astype(str).eq("primary")]
        if "includedInEconomicAnalysis" in primary:
            primary = primary[primary["includedInEconomicAnalysis"].astype(bool)]
    comp = primary[primary["arm"].astype(str).eq("comparator")] if not primary.empty else pd.DataFrame()
    inter = primary[primary["arm"].astype(str).eq("intervention")] if not primary.empty else pd.DataFrame()

    def mean_sum(frame: pd.DataFrame, column: str) -> float | None:
        if frame.empty or column not in frame:
            return None
        ids = [column for column in ("replicateId", "pairedReplicateId", "replicateSeed") if column in frame.columns]
        if not ids:
            values = pd.to_numeric(frame[column], errors="coerce").dropna()
            return None if values.empty else float(values.sum())
        grouped = pd.to_numeric(frame[column], errors="coerce").fillna(0.0).groupby(
            [frame[id_column] for id_column in ids],
            dropna=False,
        ).sum()
        return None if grouped.empty else float(grouped.mean())

    quantities = {
        "nScreened": mean_sum(inter, "screened"),
        "nTotalCoursesStarted": mean_sum(inter, "tpt_started_total"),
        "nFalsePositiveTreated": mean_sum(inter, "tpt_started_false_positive"),
        "nCuredInfection": mean_sum(inter, "infection_effectively_treated_total"),
        "nPreventedActiveTB": mean_sum(inter, "active_tb_cases_prevented"),
        "baselineActiveTBCases": mean_sum(comp, "active_tb_cases"),
        "interventionActiveTBCases": mean_sum(inter, "active_tb_cases"),
    }
    cost_values = {
        "testingCost": mean_sum(inter, "screeningTestCost"),
        "treatmentCost": mean_sum(inter, "tptRegimenCost"),
        "falsePositiveIncrementalCost": mean_sum(inter, "falsePositiveIncrementalCost"),
        "programSetupCost": mean_sum(inter, "programSetupCost"),
        "programRunningCost": mean_sum(inter, "programRunningCost"),
        "adrManagementCost": mean_sum(inter, "adrManagementCost"),
        "returnForResultsCost": mean_sum(inter, "returnForResultsCost"),
        "clinicalReviewCost": mean_sum(inter, "clinicalReviewCost"),
        "activeTBExclusionWorkupCost": mean_sum(inter, "activeTBExclusionWorkupCost"),
        "travelOutreachStaffSupportCost": mean_sum(inter, "travelOutreachStaffSupportCost"),
        "baselineTBDiseaseCost": mean_sum(comp, "activeTBDiseaseCost"),
        "interventionTBDiseaseCost": mean_sum(inter, "activeTBDiseaseCost"),
    }
    cost_values["tbDiseaseCostsAverted"] = _subtract(cost_values["baselineTBDiseaseCost"], cost_values["interventionTBDiseaseCost"])
    programme_cost_fields = [
        "testingCost",
        "treatmentCost",
        "falsePositiveIncrementalCost",
        "programSetupCost",
        "programRunningCost",
        "adrManagementCost",
        "returnForResultsCost",
        "clinicalReviewCost",
        "activeTBExclusionWorkupCost",
        "travelOutreachStaffSupportCost",
    ]
    cost_values["totalProgramCost"] = _sum_or_none([cost_values[field] for field in programme_cost_fields])
    cost_values["netCostVsBaseline"] = _subtract(cost_values["totalProgramCost"], cost_values["tbDiseaseCostsAverted"])
    result["inputs"] = deepcopy(econ_config)
    result["costNormalisation"] = deepcopy(econ_config.get("costNormalisation", {}))
    result["discounting"] = deepcopy(econ_config.get("discounting", {}))
    result["healthOutcome"] = deepcopy(econ_config.get("healthOutcome", {}))
    result["threshold"] = deepcopy(econ_config.get("threshold", {}))
    result["scopeStatement"] = result["metadata"].get(
        "scopeStatement",
        "Direct benefits, harms and costs only; transmission-mediated benefits are not included.",
    )
    result["strategy"] = {
        "testType": (econ_config.get("metadata") or {}).get("testType"),
        "regimen": (econ_config.get("metadata") or {}).get("regimen"),
    }
    result["quantities"] = quantities
    result["unitCosts"] = {
        "testPerPerson": costs.get("test"),
        "treatmentPerStarted": costs.get("regimen"),
        "falsePositiveIncrementalPerPerson": costs.get("false_positive"),
        "activeTBDiseasePerCase": costs.get("active_tb"),
        "returnForResultsPerScreened": costs.get("return_results"),
        "clinicalReviewPerTPTStarted": costs.get("clinical_review"),
        "activeTBExclusionWorkupPerTPTStarted": costs.get("active_tb_exclusion"),
        "travelOutreachStaffSupportPerScreened": costs.get("travel_outreach"),
    }
    result["costs"] = cost_values
    result["costEffectiveness"] = {
        "costPerInfectionCured": _div(cost_values["netCostVsBaseline"], quantities["nCuredInfection"]),
        "costPerTBCasePrevented": _div(cost_values["netCostVsBaseline"], quantities["nPreventedActiveTB"]),
    }
    result["_legacyCompatibilityStatus"] = {
        "missingInputs": [item["field"] for item in result.get("unresolvedInputs") or []],
        "notCalculated": [key for key, value in cost_values.items() if value is None],
        "messages": [],
        "partialCalculations": [],
    }


def _divide_series(numerator: Any, denominator: Any) -> pd.Series:
    num = pd.to_numeric(pd.Series(numerator), errors="coerce")
    den = pd.to_numeric(pd.Series(denominator), errors="coerce")
    out = pd.Series([None] * len(num), dtype="object")
    mask = num.notna() & den.notna() & den.ne(0)
    out.loc[mask] = (num.loc[mask] / den.loc[mask]).astype(float)
    return out


def _div(numerator: Any, denominator: Any) -> float | None:
    try:
        num = float(numerator)
        den = float(denominator)
    except (TypeError, ValueError):
        return None
    if den == 0:
        return None
    return num / den


def _subtract(left: Any, right: Any) -> float | None:
    if left is None or right is None:
        return None
    return float(left) - float(right)


def _sum_or_none(values: list[Any]) -> float | None:
    total = 0.0
    for value in values:
        if value is None:
            return None
        total += float(value)
    return total


def _fast_economic_validation(replicates: pd.DataFrame) -> dict[str, Any]:
    total = int(len(replicates)) if isinstance(replicates, pd.DataFrame) else 0
    complete = int(replicates["economicPairComplete"].sum()) if total and "economicPairComplete" in replicates else 0
    return {
        "isValid": bool(total and total == complete),
        "structurallyValid": bool(total and total == complete),
        "economicallyComplete": bool(total and total == complete),
        "conclusionPermitted": False,
        "totalPairedReplicates": total,
        "completePairedReplicates": complete,
        "excludedPairedReplicates": total - complete,
        "exclusionReasons": "",
        "errors": [] if total and total == complete else [{"field": "replicateResults", "message": "Incomplete frozen replicate pairs."}],
        "warnings": [],
    }


def _json_like(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_like(item) for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))}
    if isinstance(value, list):
        return [_json_like(item) for item in value]
    if isinstance(value, float):
        return round(value, 12)
    return value


def _hash_json(value: Any) -> str:
    payload = json.dumps(_json_like(value), sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _boolish(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return bool(value)
