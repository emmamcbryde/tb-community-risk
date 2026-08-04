from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

from engine.apy.scenario import (
    DEFAULT_POPULATION_PRESET_ID,
    config_updates_from_scenario,
    build_scenario_contract,
)


REPO_ROOT = Path(__file__).resolve().parents[2]


def build_default_config() -> dict[str, Any]:
    scenario = build_scenario_contract(DEFAULT_POPULATION_PRESET_ID)
    config = {
        "configVersion": "apy_v9_config_v1",
        "modelVersion": "v9",
        "scenarioLabel": "APY default scenario",
        "populationPresetId": DEFAULT_POPULATION_PRESET_ID,
        "useDefaults": True,
        "csvFile": "abm/default_data.csv",
        "outputDir": "abm/output",
        "sourceDataFiles": {
            "tbDataFile": "default_data.csv",
            "ageDistributionFile": "default_age_distribution.csv",
        },
        "N": 1500,
        "nReps": 2000,
        "seed": 1,
        "screenWindow": 3,
        "followHorizon": 20,
        "screenCoverage": 0.30,
        "screeningStrategy": "prevent",
        "ltbiPrevalence": None,
        "activeTBPrevalence": None,
        "targetAgeOR": 7.54,
        "ageDistributionFile": "",
        "ageDistributionTable": [],
        "age85PlusMax": 89,
        "ageDistributionSheet": 1,
        "riskPrev": {
            "MJ": None,
            "contact": None,
            "renal": None,
            "diabetes": None,
            "smoking": None,
            "cld": None,
            "alcohol": None,
            "female": None,
            "BCG": None,
        },
        "diseaseOR": {
            "MJ": None,
            "contact": None,
            "renal": 3.6,
            "diabetes": None,
            "smoking": None,
            "cld": None,
            "alcohol": None,
        },
        "testType": "IGRA",
        "testSensitivity": 0.95,
        "testSpecificity": 0.98,
        "tstSensitivity": 0.80,
        "tstSpecificityBCG": 0.55,
        "tstSpecificityNoBCG": 0.97,
        "testCharacteristics": {
            "IGRA": {
                "name": "IGRA",
                "sensitivity": 0.95,
                "specificity": 0.98,
                "units": "probability",
                "source": "Existing APY default parameter",
                "notes": "Editable; not folded into coverage or regimen effectiveness.",
                "resourceUse": {},
            },
            "TST": {
                "name": "TST",
                "sensitivity": 0.80,
                "specificity": 0.97,
                "specificityBCG": 0.55,
                "units": "probability",
                "source": "Existing APY default parameter",
                "notes": "Specificity can differ by BCG status. TST may require a return visit for reading; no extra resource-use value is invented here.",
                "resourceUse": {"returnVisitForReading": True},
            },
        },
        "regimen": "3HP",
        "regimenAssumptions": {
            "defaultRegimen": "3HP",
            "availableRegimens": ["3HP", "4R", "3HR", "6H", "9H"],
            "notes": (
                "Treatment offered, started, completed, adverse-event discontinuation, "
                "regimen efficacy and effectively treated true infection remain distinct model quantities."
            ),
        },
        "pStartTPT": None,
        "regimenPComplete": None,
        "regimenADRstop": None,
        "regimenEffFull": None,
        "partialShortCourseMode": None,
        "partialDoseFractionADR": None,
        "partialDoseFractionOther": None,
        "earlyLateRatio": None,
    }
    config.update(config_updates_from_scenario(scenario))
    return config


def normalise_config(config: dict[str, Any] | None = None) -> dict[str, Any]:
    normalised = build_default_config()
    if config:
        _deep_update(normalised, config)
    normalised["testType"] = _normalise_text(normalised.get("testType"))
    normalised["regimen"] = _normalise_text(normalised.get("regimen"))
    normalised["screeningStrategy"] = _normalise_text(
        normalised.get("screeningStrategy")
    )
    normalised["partialShortCourseMode"] = _normalise_text(
        normalised.get("partialShortCourseMode")
    )
    return normalised


def strip_empty_fields(d: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in d.items():
        if isinstance(value, dict):
            nested = strip_empty_fields(value)
            if nested:
                out[key] = nested
        elif value is None:
            continue
        elif value == "":
            continue
        elif isinstance(value, (list, tuple)) and len(value) == 0:
            continue
        else:
            out[key] = value
    return out


def resolve_repo_path(path: str | Path) -> Path:
    path_value = Path(path)
    if path_value.is_absolute():
        return path_value
    return REPO_ROOT / path_value


def _deep_update(base: dict[str, Any], overrides: dict[str, Any]) -> None:
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            _deep_update(base[key], value)
        else:
            base[key] = deepcopy(value)


def _normalise_text(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None
