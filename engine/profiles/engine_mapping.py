"""Map a population profile onto the existing validated engine configuration.

Engine semantics are preserved exactly:

* Bundled demonstration values are passed as "use engine default" so that the
  engine's age-specific demonstration inputs are applied unchanged.
* A user-defined risk-factor prevalence is applied uniformly across age groups
  (the engine's existing single-value override behaviour).
* An effect estimate is applied by the engine as a multiplier on the hazard of
  progression from infection to disease, whatever its declared measure type
  (RR, HR or OR). The measure type is retained; no conversion is performed.
* A disabled or excluded risk factor, or a profile without risk factors, is run
  without stratification by that factor (prevalence 0, so no multiplier applies).
* The attached incidence series is descriptive in this milestone; it does not
  yet drive calibration, infection pressure or the dynamic model.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any

from engine.profiles.demonstration import (
    ENGINE_EFFECT_APPLICATION,
    ENGINE_KEYS,
    GENERAL_DEFAULT_STOCHASTIC_SIMULATIONS,
)
from engine.profiles.population_profile import (
    EffectMeasure,
    PopulationProfile,
    Provenance,
    ValueState,
    validate_profile,
)


GENERAL_PROFILE_LINK_VERSION = "general_profile_link_v1"
DEFAULT_ANALYSIS = {
    "analysisMethod": "agent_based",
    "nReps": GENERAL_DEFAULT_STOCHASTIC_SIMULATIONS,
    "seed": 1,
}
ANALYSIS_METHOD_LABELS = {
    "agent_based": "Stochastic analysis - simulated populations",
    "expected_value": "Deterministic expected-value preview",
}
EFFECT_INTERPRETATION = {
    EffectMeasure.HR: "Hazard ratio applied directly as a progression-hazard multiplier.",
    EffectMeasure.RR: "Risk ratio applied as a progression-hazard multiplier without conversion.",
    EffectMeasure.OR: "Odds ratio applied as a progression-hazard multiplier without conversion.",
}


class ProfileMappingError(ValueError):
    """Raised when a profile cannot yet be expressed in the engine configuration."""


def default_intervention() -> dict[str, Any]:
    config = _base_engine_config()
    return {
        "testType": config["testType"],
        "regimen": config["regimen"],
        "screenCoverage": float(config["screenCoverage"]),
        "screeningStrategy": config["screeningStrategy"],
        "screeningWindowYears": int(config["screeningWindowYears"]),
    }


def default_analysis() -> dict[str, Any]:
    return dict(DEFAULT_ANALYSIS)


def build_engine_config(
    profile: PopulationProfile,
    *,
    intervention: dict[str, Any] | None = None,
    analysis: dict[str, Any] | None = None,
    snapshot_manifest: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Return an engine configuration for ``profile``; raises on blocking issues."""
    issues = [issue for issue in validate_profile(profile) if issue["severity"] in {"error", "blocking"}]
    if issues:
        raise ProfileMappingError("; ".join(f"{issue['field']}: {issue['message']}" for issue in issues))
    if any(band.proportion.is_user_override for band in profile.age_distribution):
        raise ProfileMappingError("User-defined age distributions are not yet supported by the analysis engine.")

    config = _base_engine_config()
    intervention = {**default_intervention(), **(intervention or {})}
    analysis = {**DEFAULT_ANALYSIS, **(analysis or {})}

    population = int(profile.population_size.value)
    config["N"] = population
    config.setdefault("scenario", {})["populationSize"] = population

    if profile.ltbi_prevalence.provenance is Provenance.USER_DEFINED:
        config["ltbiPrevalence"] = profile.ltbi_prevalence.value

    risk_prev = dict(config.get("riskPrev") or {})
    disease_effects = dict(config.get("diseaseOR") or {})
    mapped_keys = set()
    for factor in profile.risk_factors:
        key = factor.engine_key
        if key is None:
            continue
        if key not in ENGINE_KEYS:
            raise ProfileMappingError(f"Unknown engine risk-factor key {key!r}.")
        mapped_keys.add(key)
        if not factor.enabled or factor.prevalence.state in {ValueState.EXCLUDED, ValueState.NOT_APPLICABLE}:
            risk_prev[key] = 0.0
            continue
        if factor.prevalence.provenance is Provenance.USER_DEFINED:
            risk_prev[key] = factor.prevalence.value
        if factor.effect_estimate.provenance is Provenance.USER_DEFINED:
            disease_effects[key] = factor.effect_estimate.value
    for key in ENGINE_KEYS:
        if key not in mapped_keys:
            risk_prev[key] = 0.0
    config["riskPrev"] = risk_prev
    config["diseaseOR"] = disease_effects
    config["scenario"].setdefault("riskFactorAssumptions", {})
    config["scenario"]["riskFactorAssumptions"]["riskPrev"] = deepcopy(risk_prev)
    config["scenario"]["riskFactorAssumptions"]["diseaseOR"] = deepcopy(disease_effects)

    config["testType"] = intervention["testType"]
    config["regimen"] = intervention["regimen"]
    config["screenCoverage"] = float(intervention["screenCoverage"])
    config["screeningStrategy"] = intervention["screeningStrategy"]
    config["screeningWindowYears"] = int(intervention["screeningWindowYears"])
    config["screenWindow"] = int(intervention["screeningWindowYears"])

    method = str(analysis["analysisMethod"])
    if method not in ANALYSIS_METHOD_LABELS:
        raise ProfileMappingError(f"Unknown analysis method {method!r}.")
    config["analysisMethod"] = method
    config["analysisMethodLabel"] = ANALYSIS_METHOD_LABELS[method]
    config["nReps"] = int(analysis["nReps"])
    config["seed"] = int(analysis["seed"])
    config["simulationMode"] = "custom"
    config["simulationModeLabel"] = f"{int(analysis['nReps']):,} simulated populations"
    config["scenarioLabel"] = f"General community TB analysis: {profile.name}"
    config["generalProfileLink"] = profile_link(profile, snapshot_manifest)
    return config


def profile_link(profile: PopulationProfile, snapshot_manifest: dict[str, Any] | None = None) -> dict[str, Any]:
    """Provenance attached to runs so saved analyses identify their exact inputs."""
    manifest = snapshot_manifest or {}
    detail = dict(profile.incidence.source_detail)
    return {
        "linkVersion": GENERAL_PROFILE_LINK_VERSION,
        "profileId": profile.profile_id,
        "profileName": profile.name,
        "profileSchemaVersion": profile.schema_version,
        "profileHash": profile.profile_hash(),
        "location": profile.location.name,
        "iso3": profile.location.iso3,
        "incidenceProvenance": profile.incidence.provenance.value,
        "incidenceSnapshotId": profile.incidence.snapshot_id,
        "incidenceSnapshotDataSha256": detail.get("dataSha256") or ((manifest.get("dataFile") or {}).get("sha256")),
        "incidenceDataHash": profile.incidence.data_hash,
        "incidenceUsedByEngine": False,
        "dataVintage": profile.data_vintage,
    }


NON_EPIDEMIOLOGICAL_KEYS = {"generalProfileLink", "scenarioLabel", "analysisMethodLabel", "simulationModeLabel"}


def epidemiological_config_hash(config: dict[str, Any]) -> str:
    """Hash of the engine inputs that determine health outcomes (labels and provenance excluded)."""
    import hashlib
    import json

    payload = {key: value for key, value in config.items() if key not in NON_EPIDEMIOLOGICAL_KEYS}
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def effect_interpretation_rows(profile: PopulationProfile) -> list[dict[str, str]]:
    """Technical statement of how each risk-factor effect is applied internally."""
    rows = []
    for factor in profile.risk_factors:
        if factor.engine_key is None:
            application = "Not used by the current analysis engine."
        elif factor.effect_measure is None:
            application = "No effect-measure type recorded."
        else:
            application = EFFECT_INTERPRETATION[factor.effect_measure]
        rows.append(
            {
                "Risk factor": factor.label,
                "Declared measure": "" if factor.effect_measure is None else factor.effect_measure.value,
                "Internal application": (
                    "Progression hazard multiplier"
                    if factor.engine_application == ENGINE_EFFECT_APPLICATION
                    else "Not applied"
                ),
                "Note": application,
            }
        )
    return rows


def _base_engine_config() -> dict[str, Any]:
    from engine.apy.working_defaults import build_unified_working_default_preset

    return deepcopy(build_unified_working_default_preset()["config"])
