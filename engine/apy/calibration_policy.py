from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from typing import Any

from adapters.serialization import to_json_like
from engine.apy.calibration import (
    calibrate_early_hazard,
    calibrate_from_config,
    expected_infection_prevalence_exact,
    solve_loglambda_for_gamma,
)
from engine.apy.config import normalise_config
from engine.apy.data import load_parameters_from_config
from engine.apy.ltbi_state import resolve_ltbi_state_assumptions
from engine.apy.timing import resolve_time_settings


CALIBRATION_ARTIFACT_VERSION = "apy_reference_calibration_artifact_v1"
CALIBRATION_POLICIES = {
    "full_reference_calibration",
    "infection_intercept_only",
    "progression_hazards_only",
    "none",
}
_REFERENCE_ARTIFACT_CACHE: dict[str, dict[str, Any]] = {}


def build_reference_calibration_artifact(config: dict[str, Any]) -> dict[str, Any]:
    cfg = normalise_config(config)
    key = config_hash({**cfg, "referenceCalibrationArtifact": None, "calibrationPolicy": "full_reference_calibration"})
    if key in _REFERENCE_ARTIFACT_CACHE:
        return deepcopy(_REFERENCE_ARTIFACT_CACHE[key])
    calibration = calibrate_from_config(cfg)
    ltbi_state = resolve_ltbi_state_assumptions(cfg)
    timing = resolve_time_settings(cfg)
    artifact = {
        "artifactVersion": CALIBRATION_ARTIFACT_VERSION,
        "calibrationVersion": "apy_python_calibration_v1",
        "configurationHash": config_hash(cfg),
        "infectionAgeShapeParameter": calibration["ageInfGamma"],
        "infectionIntercept": calibration["ageInfLogLambda"],
        "expectedInfPrev": calibration["expectedInfPrev"],
        "targetInfPrev": calibration["targetInfPrev"],
        "recentLTBIProportion": ltbi_state.get("baselineRecentLTBIProportion"),
        "recentLTBIProportionStatus": ltbi_state.get("baselineRecentLTBIProportionStatus") or ltbi_state.get("status"),
        "recentToRemoteTransitionModel": ltbi_state.get("transitionModel"),
        "recentToRemoteTransitionRatePerYear": ltbi_state.get("recentToRemoteTransitionRatePerYear"),
        "earlyProgressionHazard": calibration["lambdaEarly"],
        "lateProgressionHazard": calibration["lambdaLate"],
        "targetActiveAtCalibrationHorizon": calibration["targetActiveAtCalibrationHorizon"],
        "expectedActiveAtCalibrationHorizon": calibration["expectedActiveAtCalibrationHorizon"],
        "earlyProgressionPeriodYears": timing["earlyProgressionPeriodYears"],
        "activeTBCalibrationHorizonYears": timing["activeTBCalibrationHorizonYears"],
        "targetAgeOR": calibration["targetAgeOR"],
        "readinessStatus": ltbi_state.get("status"),
        "provisional": bool(ltbi_state.get("provisional")),
    }
    artifact["artifactHash"] = calibration_artifact_hash(artifact)
    _REFERENCE_ARTIFACT_CACHE[key] = deepcopy(artifact)
    return artifact


def resolve_calibration_for_config(config: dict[str, Any]) -> dict[str, Any]:
    cfg = normalise_config(config)
    policy = str(cfg.get("calibrationPolicy") or "full_reference_calibration")
    if policy not in CALIBRATION_POLICIES:
        raise ValueError(f"Unknown calibrationPolicy: {policy}")
    if policy == "full_reference_calibration":
        calibration = calibrate_from_config(cfg)
        reference = build_reference_calibration_artifact(cfg)
        return _attach_policy_metadata(
            calibration,
            policy=policy,
            reference=reference,
            retained=[],
            recalibrated=["infectionIntercept", "infectionAgeShapeParameter", "earlyProgressionHazard", "lateProgressionHazard"],
        )

    reference = cfg.get("referenceCalibrationArtifact")
    if not isinstance(reference, dict):
        reference = build_reference_calibration_artifact({**cfg, "calibrationPolicy": "full_reference_calibration"})
    pars = load_parameters_from_config(cfg)
    calibration = _calibration_from_artifact(reference, pars)
    if policy == "none":
        return _attach_policy_metadata(
            calibration,
            policy=policy,
            reference=reference,
            retained=["infectionIntercept", "infectionAgeShapeParameter", "earlyProgressionHazard", "lateProgressionHazard"],
            recalibrated=[],
        )
    if policy == "infection_intercept_only":
        target = _target_ltbi_prevalence(cfg)
        if target <= 0:
            calibration["ageInfLogLambda"] = -745.0
            calibration["expectedInfPrev"] = 0.0
            calibration["zeroInfectionPrevalence"] = True
        else:
            calibration["ageInfLogLambda"] = solve_loglambda_for_gamma(
                float(reference["infectionAgeShapeParameter"]),
                pars,
                target,
            )
            calibration["expectedInfPrev"] = expected_infection_prevalence_exact(
                calibration["ageInfLogLambda"],
                calibration["ageInfGamma"],
                pars,
            )
            calibration["zeroInfectionPrevalence"] = False
        calibration["targetInfPrev"] = target
        return _attach_policy_metadata(
            calibration,
            policy=policy,
            reference=reference,
            retained=["infectionAgeShapeParameter", "earlyProgressionHazard", "lateProgressionHazard", "recentToRemoteTransitionRatePerYear"],
            recalibrated=["infectionIntercept"],
        )
    if policy == "progression_hazards_only":
        target_active = _target_active_prevalence(cfg)
        timing = resolve_time_settings(cfg)
        ltbi = resolve_ltbi_state_assumptions(cfg)
        hazard = calibrate_early_hazard(
            pars,
            calibration["ageInfLogLambda"],
            calibration["ageInfGamma"],
            target_active,
            timing["activeTBCalibrationHorizonYears"],
            _default_if_empty(cfg.get("earlyLateRatio"), 5.0),
            timing["earlyProgressionPeriodYears"],
            float(ltbi["baselineRecentLTBIProportion"]),
            float(ltbi["recentToRemoteTransitionRatePerYear"]),
        )
        calibration["lambdaEarly"] = hazard["lambdaEarly"]
        calibration["lambdaLate"] = hazard["lambdaLate"]
        calibration["expectedActiveAtCalibrationHorizon"] = hazard["expectedActiveAtCalibrationHorizon"]
        calibration["expectedActive2y"] = hazard["expectedActiveAtCalibrationHorizon"]
        calibration["targetActiveAtCalibrationHorizon"] = target_active
        calibration["targetActive2y"] = target_active
        return _attach_policy_metadata(
            calibration,
            policy=policy,
            reference=reference,
            retained=["infectionIntercept", "infectionAgeShapeParameter"],
            recalibrated=["earlyProgressionHazard", "lateProgressionHazard"],
        )
    raise AssertionError("Unhandled calibration policy.")


def config_hash(config: dict[str, Any]) -> str:
    payload = json.dumps(to_json_like(config), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def calibration_artifact_hash(artifact: dict[str, Any]) -> str:
    payload = deepcopy(artifact)
    payload.pop("artifactHash", None)
    return config_hash(payload)


def _calibration_from_artifact(artifact: dict[str, Any], pars: dict[str, Any]) -> dict[str, Any]:
    return {
        "parameters": pars,
        "ageInfLogLambda": artifact.get("infectionIntercept"),
        "ageInfGamma": artifact.get("infectionAgeShapeParameter"),
        "expectedInfPrev": artifact.get("expectedInfPrev"),
        "expectedAgeOR": None,
        "lambdaEarly": artifact.get("earlyProgressionHazard"),
        "lambdaLate": artifact.get("lateProgressionHazard"),
        "expectedActiveAtCalibrationHorizon": artifact.get("expectedActiveAtCalibrationHorizon"),
        "expectedActive2y": artifact.get("expectedActiveAtCalibrationHorizon"),
        "targetInfPrev": artifact.get("targetInfPrev"),
        "targetAgeOR": artifact.get("targetAgeOR"),
        "targetActiveAtCalibrationHorizon": artifact.get("targetActiveAtCalibrationHorizon"),
        "targetActive2y": artifact.get("targetActiveAtCalibrationHorizon"),
        "earlyProgressionPeriodYears": artifact.get("earlyProgressionPeriodYears"),
        "activeTBCalibrationHorizonYears": artifact.get("activeTBCalibrationHorizonYears"),
        "baselineRecentLTBIProportion": artifact.get("recentLTBIProportion"),
        "recentToRemoteTransitionRatePerYear": artifact.get("recentToRemoteTransitionRatePerYear"),
        "ltbiStateAssumptionStatus": artifact.get("recentLTBIProportionStatus"),
        "zeroInfectionPrevalence": False,
    }


def _attach_policy_metadata(
    calibration: dict[str, Any],
    *,
    policy: str,
    reference: dict[str, Any],
    retained: list[str],
    recalibrated: list[str],
) -> dict[str, Any]:
    out = deepcopy(calibration)
    out["calibrationPolicy"] = policy
    out["referenceCalibrationHash"] = reference.get("artifactHash") or calibration_artifact_hash(reference)
    out["referenceCalibrationArtifactVersion"] = reference.get("artifactVersion")
    out["calibrationParametersRetained"] = retained
    out["calibrationParametersRecalibrated"] = recalibrated
    out["referenceCalibrationArtifact"] = deepcopy(reference)
    return out


def _target_ltbi_prevalence(config: dict[str, Any]) -> float:
    return _default_if_empty(config.get("ltbiPrevalence"), 47 / 624)


def _target_active_prevalence(config: dict[str, Any]) -> float:
    return _default_if_empty(config.get("activeTBPrevalence"), 10 / 770)


def _default_if_empty(value: Any, default: float) -> float:
    if value is None or value in ("", []):
        return float(default)
    return float(value)
