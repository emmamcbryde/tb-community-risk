from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

from engine.apy.config import normalise_config
from engine.apy.data import load_parameters_from_config
from engine.apy.scenario import (
    DEFAULT_POPULATION_PRESET_ID,
    build_scenario_contract,
)


RISK_FACTOR_LABELS = {
    "female": "Female",
    "BCG": "Prior BCG vaccination",
    "contact": "Contact history",
    "smoking": "Smoking",
    "MJ": "Cannabis exposure",
    "renal": "Renal impairment",
    "diabetes": "Diabetes",
    "cld": "Chronic lung disease",
    "alcohol": "Alcohol/drug exposure",
}


def resolved_demographic_profile(config: dict[str, Any]) -> dict[str, Any]:
    """Return the APY demographic values actually consumed by model loaders."""
    cfg = normalise_config(config)
    pars = load_parameters_from_config(cfg)
    scenario = cfg.get("scenario") or {}
    source = _relative_path(pars.get("ageDistributionSource") or cfg.get("ageDistributionFile") or "")
    eligible = _eligible_population(cfg)
    return {
        "populationSize": int(cfg.get("N", 0) or 0),
        "eligiblePopulation": eligible,
        "populationPresetId": cfg.get("populationPresetId", DEFAULT_POPULATION_PRESET_ID),
        "ageDistributionSource": source or "Repository default age distribution",
        "ageEvidenceStatus": (
            "Repository APY default used by the frozen compatibility analysis; "
            "external demographic provenance has not been independently reviewed in this workflow."
        ),
        "riskFactorSource": _risk_source_label(scenario),
        "profileUse": "Used in cohort generation, LTBI probability, risk-factor assignment and targeting.",
        "editableStatus": "Editable; demographic changes require rerunning epidemiology.",
        "ageRows": _age_rows(pars),
        "broadAgeRows": _broad_age_rows(pars),
        "riskRows": _risk_rows(pars),
    }


def demographic_summary_rows(config: dict[str, Any]) -> list[dict[str, Any]]:
    profile = resolved_demographic_profile(config)
    return [
        {"Item": "Total eligible population", "Value": profile["eligiblePopulation"]},
        {"Item": "Population preset", "Value": profile["populationPresetId"]},
        {"Item": "Demographic source", "Value": profile["ageDistributionSource"]},
        {"Item": "Evidence status", "Value": profile["ageEvidenceStatus"]},
        {"Item": "How the profile is used", "Value": profile["profileUse"]},
        {"Item": "Editing status", "Value": profile["editableStatus"]},
    ]


def age_distribution_rows(config: dict[str, Any]) -> list[dict[str, Any]]:
    return list(resolved_demographic_profile(config)["ageRows"])


def broad_age_distribution_rows(config: dict[str, Any]) -> list[dict[str, Any]]:
    return list(resolved_demographic_profile(config)["broadAgeRows"])


def risk_factor_rows(config: dict[str, Any]) -> list[dict[str, Any]]:
    return list(resolved_demographic_profile(config)["riskRows"])


def demographic_profile_hash(config: dict[str, Any]) -> tuple[Any, ...]:
    profile = resolved_demographic_profile(config)
    age = tuple((row["Age group"], round(float(row["Proportion"]), 12)) for row in profile["ageRows"])
    risk = tuple((row["Risk factor"], round(float(row["Prevalence"]), 12)) for row in profile["riskRows"])
    return (
        profile["populationSize"],
        profile["eligiblePopulation"],
        profile["ageDistributionSource"],
        age,
        risk,
    )


def restore_apy_demographic_defaults(config: dict[str, Any]) -> dict[str, Any]:
    """Reset only APY population/demographic inputs from the scenario contract."""
    restored = deepcopy(config)
    scenario = build_scenario_contract(DEFAULT_POPULATION_PRESET_ID)
    risk = scenario.get("riskFactorAssumptions") or {}
    age_source = scenario.get("ageProfileSource") or {}

    restored["populationPresetId"] = scenario.get("populationPresetId", DEFAULT_POPULATION_PRESET_ID)
    restored["N"] = scenario.get("populationSize", restored.get("N"))
    restored["ageDistributionFile"] = age_source.get("path", "")
    restored["ageDistributionTable"] = []
    restored["riskPrev"] = deepcopy(risk.get("riskPrev") or {})
    if isinstance(risk.get("diseaseOR"), dict):
        restored["diseaseOR"] = deepcopy(risk["diseaseOR"])

    nested = deepcopy(restored.get("scenario") or {})
    nested["populationPresetId"] = scenario.get("populationPresetId", DEFAULT_POPULATION_PRESET_ID)
    nested["populationSize"] = scenario.get("populationSize", nested.get("populationSize"))
    nested["ageProfileSource"] = deepcopy(age_source)
    nested["riskFactorAssumptions"] = deepcopy(risk)
    nested["eligible"] = deepcopy(scenario.get("eligible") or nested.get("eligible") or {})
    nested["sourcesAndNotes"] = deepcopy(scenario.get("sourcesAndNotes") or nested.get("sourcesAndNotes") or [])
    restored["scenario"] = nested

    source_files = deepcopy(restored.get("sourceDataFiles") or {})
    if age_source.get("path"):
        source_files["ageDistributionFile"] = Path(str(age_source["path"])).name
    restored["sourceDataFiles"] = source_files
    return restored


def _age_rows(pars: dict[str, Any]) -> list[dict[str, Any]]:
    table = pars.get("ageDistributionTable") or []
    if table:
        weights = [_age_row_weight(row) for row in table]
        total = sum(weights)
        rows = []
        for row, weight in zip(table, weights):
            age_group = str(row.get("age_group") or row.get("age") or row.get("Age") or "").strip()
            if not age_group or weight <= 0:
                continue
            rows.append(
                {
                    "Age group": age_group,
                    "Proportion": weight / total if total > 0 else 0.0,
                    "Source proportion": weight,
                }
            )
        return rows

    by_age = {}
    for age, prob in zip(pars.get("exactAgeValues") or [], pars.get("exactAgeProb") or []):
        by_age[str(age)] = by_age.get(str(age), 0.0) + float(prob)
    return [{"Age group": age, "Proportion": proportion, "Source proportion": proportion} for age, proportion in by_age.items()]


def _broad_age_rows(pars: dict[str, Any]) -> list[dict[str, Any]]:
    labels = ["0-4 years", "5-14 years", "15+ years"]
    return [
        {"Age group": label, "Proportion": float(value)}
        for label, value in zip(labels, pars.get("popFrac") or [])
    ]


def _risk_rows(pars: dict[str, Any]) -> list[dict[str, Any]]:
    keys = {
        "female": "totalFemalePrev",
        "BCG": "totalBCGPrev",
        "contact": "totalContactPrev",
        "smoking": "totalCurrentSmokerPrev",
        "MJ": "totalMJPrev",
        "renal": "totalRenalPrev",
        "diabetes": "totalDiabetesPrev",
        "cld": "totalCLDPrev",
        "alcohol": "totalAlcoholPrev",
    }
    return [
        {"Risk factor": RISK_FACTOR_LABELS[key], "Prevalence": float(pars[value_key])}
        for key, value_key in keys.items()
        if value_key in pars
    ]


def _age_row_weight(row: dict[str, Any]) -> float:
    for key in ["smoothed proportion", "Smoothed proportion", "proportion", "Proportion"]:
        try:
            value = row.get(key)
            if value not in (None, ""):
                return float(value)
        except (TypeError, ValueError):
            return 0.0
    return 0.0


def _eligible_population(config: dict[str, Any]) -> int:
    n = int(config.get("N", 0) or 0)
    eligible = ((config.get("scenario") or {}).get("eligible") or {})
    if eligible.get("number") not in (None, ""):
        return int(float(eligible["number"]))
    if eligible.get("proportion") not in (None, ""):
        return int(round(n * float(eligible["proportion"])))
    return n


def _risk_source_label(scenario: dict[str, Any]) -> str:
    risk = scenario.get("riskFactorAssumptions") or {}
    source = risk.get("source") or "Repository APY risk-factor defaults"
    return str(source)


def _relative_path(value: str) -> str:
    if not value:
        return ""
    path = Path(value)
    try:
        root = Path(__file__).resolve().parents[1]
        return str(path.resolve().relative_to(root))
    except (OSError, ValueError):
        return str(value)
