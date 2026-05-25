from __future__ import annotations

from copy import deepcopy
from pathlib import Path, PureWindowsPath
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]


def build_default_config() -> dict[str, Any]:
    return {
        "configVersion": "apy_v9_config_v1",
        "modelVersion": "v9",
        "scenarioLabel": "APY default scenario",
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
        "screenWindow": 2,
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
        "testSensitivity": None,
        "testSpecificity": None,
        "tstSensitivity": None,
        "tstSpecificityBCG": None,
        "tstSpecificityNoBCG": None,
        "regimen": "3HP",
        "pStartTPT": None,
        "regimenPComplete": None,
        "regimenADRstop": None,
        "regimenEffFull": None,
        "partialShortCourseMode": None,
        "partialDoseFractionADR": None,
        "partialDoseFractionOther": None,
        "earlyLateRatio": None,
    }


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
        if not path_value.exists():
            repo_local_path = _repo_local_suffix_path(PureWindowsPath(str(path)))
            if repo_local_path is not None:
                return repo_local_path
        return path_value
    windows_path = PureWindowsPath(str(path))
    if windows_path.is_absolute():
        repo_local_path = _repo_local_suffix_path(windows_path)
        if repo_local_path is not None:
            return repo_local_path
    return REPO_ROOT / path_value


def _repo_local_suffix_path(path: PureWindowsPath) -> Path | None:
    parts = [
        part
        for part in path.parts
        if part not in {path.drive, path.root, path.anchor}
    ]
    for i in range(len(parts)):
        candidate = REPO_ROOT.joinpath(*parts[i:])
        if candidate.exists():
            return candidate
    return None


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
