from __future__ import annotations

import csv
import math
import re
from copy import deepcopy
from pathlib import Path
from typing import Any

import pandas as pd

from engine.apy.age_distribution import (
    broad_age_group_from_years,
    build_fallback_exact_age_prob,
    expand_age_distribution_table,
)
from engine.apy.config import normalise_config, resolve_repo_path


AGE_NAMES = ["0-4", "5-14", ">=15"]


def load_tb_csv(csv_file: str | Path) -> dict[str, Any]:
    rows = _read_messy_csv(csv_file)
    pars: dict[str, Any] = {
        "ageNames": AGE_NAMES.copy(),
        "popFrac": _get_row_values(rows, 0, "population_fraction", [1, 2, 3]),
        "baseInfByAge": _get_row_values(rows, 0, "infection_risk_by_age", [1, 2, 3]),
        "mjPrevByAge": _get_row_values(rows, 0, "MJ_use", [1, 2, 3]),
        "contactPrevByAge": _get_row_values(rows, 0, "contact", [1, 2, 3]),
        "renalPrevByAge": _get_row_values(rows, 0, "renal", [1, 2, 3]),
        "diabetesPrevByAge": _get_row_values(rows, 0, "diabetes", [1, 2, 3]),
        "smokingPrevByAge": _get_row_values(rows, 0, "smoking", [1, 2, 3]),
        "cldPrevByAge": _get_row_values(rows, 0, "chronic lung disease", [1, 2, 3]),
        "alcoholPrevByAge": _get_row_values(rows, 0, "alcohol/drugs", [1, 2, 3]),
        "infOR": {
            "MJ": _safe_default(_get_row_values(rows, 0, "MJ_use", 5, True), 1),
            "contact": _safe_default(_get_row_values(rows, 0, "contact", 5, True), 1),
            "renal": _safe_default(_get_row_values(rows, 0, "renal", 5, True), 1),
        },
        "disOR": {
            "MJ": _safe_default(_get_row_values(rows, 0, "MJ_use", 4, True), 1),
            "contact": _safe_default(_get_row_values(rows, 0, "contact", 4, True), 1),
            "renal": _safe_default(_get_row_values(rows, 0, "renal", 4, True), 1),
            "diabetes": _safe_default(_get_row_values(rows, 0, "diabetes", 4, True), 1),
            "smoking": _safe_default(_get_row_values(rows, 0, "smoking", 4, True), 1),
            "cld": _safe_default(_get_row_values(rows, 0, "chronic lung disease", 4, True), 1),
            "alcohol": _safe_default(_get_row_values(rows, 0, "alcohol/drugs", 4, True), 1),
        },
        "totalFemalePrev": _get_row_values(rows, 8, "female sex", 9, True),
        "totalContactPrev": _get_row_values(rows, 8, "close contact", 9, True),
        "totalCurrentSmokerPrev": _get_row_values(rows, 8, "current smoker", 9, True),
        "totalMJPrev": _get_row_values(rows, 8, "water-pipe marijuana", 9, True),
        "totalRenalPrev": _get_row_values(rows, 8, "renal disease", 9, True),
        "totalDiabetesPrev": _get_row_values(rows, 8, "diabetes", 9, True),
        "totalCLDPrev": _get_row_values(rows, 8, "chronic lung disease", 9, True),
        "totalBCGPrev": _get_row_values(rows, 8, "bcg vaccinated", 9, True),
        "comorbidityCountCats": [0, 1, 2, 3],
        "comorbidityCountProb": [
            _get_row_values(rows, 7, "proportion comorbidity", 8),
            _get_row_values(rows, 7, "proportion comorbidity", 9),
            _get_row_values(rows, 7, "proportion comorbidity", 10),
            _get_row_values(rows, 7, "proportion comorbidity", 11),
        ],
    }
    return pars


def resolve_age_distribution_file(
    age_distribution_file: str | Path | None, csv_file: str | Path
) -> Path | None:
    if age_distribution_file is not None and str(age_distribution_file).strip():
        return resolve_repo_path(age_distribution_file)

    csv_path = resolve_repo_path(csv_file)
    for name in [
        "default_age_distribution.csv",
        "default_age_distribution.xlsx",
        "age_groups_SA.xlsx",
    ]:
        candidate = csv_path.parent / name
        if candidate.is_file():
            return candidate
    return None


def load_age_distribution(
    pars: dict[str, Any],
    age_file: str | Path | None,
    age_distribution_table: Any = None,
    age85_plus_max: int = 89,
) -> dict[str, Any]:
    out = deepcopy(pars)
    if age_distribution_table is not None and _has_rows(age_distribution_table):
        table = pd.DataFrame(age_distribution_table)
        source = "provided_table"
    elif age_file is None or not Path(age_file).is_file():
        exact_ages = list(range(0, 90))
        out["exactAgeValues"] = exact_ages
        out["exactAgeProb"] = build_fallback_exact_age_prob(exact_ages, out["popFrac"])
        out["ageDistributionSource"] = ""
        out["ageDistributionTable"] = []
        return out
    else:
        path = Path(age_file)
        table = pd.read_csv(path) if path.suffix.lower() == ".csv" else pd.read_excel(path)
        source = str(path)

    exact_ages, exact_prob = expand_age_distribution_table(table, age85_plus_max)
    out["exactAgeValues"] = exact_ages
    out["exactAgeProb"] = exact_prob
    out["ageDistributionSource"] = source
    out["ageDistributionTable"] = table.to_dict(orient="records")
    return out


def refresh_total_prevalence_summaries(pars: dict[str, Any]) -> dict[str, Any]:
    out = deepcopy(pars)
    if out.get("exactAgeValues") and out.get("exactAgeProb"):
        groups = broad_age_group_from_years(out["exactAgeValues"])
        weights = out["exactAgeProb"]
        out["popFrac"] = [
            sum(w for g, w in zip(groups, weights) if g == 1),
            sum(w for g, w in zip(groups, weights) if g == 2),
            sum(w for g, w in zip(groups, weights) if g == 3),
        ]
    else:
        exact_ages = list(range(0, 90))
        groups = broad_age_group_from_years(exact_ages)
        weights = build_fallback_exact_age_prob(exact_ages, out["popFrac"])

    out["totalMJPrev"] = _weighted_age_lookup(out["mjPrevByAge"], groups, weights)
    out["totalContactPrev"] = _weighted_age_lookup(out["contactPrevByAge"], groups, weights)
    out["totalCurrentSmokerPrev"] = _weighted_age_lookup(out["smokingPrevByAge"], groups, weights)
    out["totalRenalPrev"] = _weighted_age_lookup(out["renalPrevByAge"], groups, weights)
    out["totalDiabetesPrev"] = _weighted_age_lookup(out["diabetesPrevByAge"], groups, weights)
    out["totalCLDPrev"] = _weighted_age_lookup(out["cldPrevByAge"], groups, weights)
    out["totalFemalePrev"] = _safe_default(out.get("totalFemalePrev"), 0.56)
    out["totalBCGPrev"] = _safe_default(out.get("totalBCGPrev"), 0.68)
    return out


def apply_risk_prevalence_overrides(
    pars: dict[str, Any], overrides: dict[str, Any] | None
) -> dict[str, Any]:
    out = deepcopy(pars)
    if not overrides:
        return out
    for raw_name, value in overrides.items():
        if _is_empty(value):
            continue
        name = _normalise_name(raw_name)
        if name in {"mj", "marijuana", "waterpipemarijuana", "waterpipemarijuanause"}:
            out["mjPrevByAge"] = _normalise_age_prevalence_override(value, raw_name)
        elif name in {"contact", "closecontact", "tbcontact", "closecontactwithtbcase"}:
            out["contactPrevByAge"] = _normalise_age_prevalence_override(value, raw_name)
        elif name in {"renal", "renaldisease"}:
            out["renalPrevByAge"] = _normalise_age_prevalence_override(value, raw_name)
        elif name == "diabetes":
            out["diabetesPrevByAge"] = _normalise_age_prevalence_override(value, raw_name)
        elif name in {"smoking", "currentsmoker", "eversmoker"}:
            out["smokingPrevByAge"] = _normalise_age_prevalence_override(value, raw_name)
        elif name in {"cld", "chroniclungdisease", "chroniclung"}:
            out["cldPrevByAge"] = _normalise_age_prevalence_override(value, raw_name)
        elif name in {"alcohol", "alcoholdrugs", "harmfulalcohol"}:
            out["alcoholPrevByAge"] = _normalise_age_prevalence_override(value, raw_name)
        elif name in {"female", "femalesex", "sexfemale"}:
            out["totalFemalePrev"] = _normalise_overall_prevalence_override(value, raw_name)
        elif name in {"bcg", "priorbcg", "bcgvaccinated"}:
            out["totalBCGPrev"] = _normalise_overall_prevalence_override(value, raw_name)
    return out


def apply_disease_or_overrides(
    pars: dict[str, Any], overrides: dict[str, Any] | None
) -> dict[str, Any]:
    out = deepcopy(pars)
    if not overrides:
        return out
    for raw_name, value in overrides.items():
        if _is_empty(value):
            continue
        name = _normalise_name(raw_name)
        target = {
            "mj": "MJ",
            "marijuana": "MJ",
            "waterpipemarijuana": "MJ",
            "waterpipemarijuanause": "MJ",
            "contact": "contact",
            "closecontact": "contact",
            "tbcontact": "contact",
            "closecontactwithtbcase": "contact",
            "renal": "renal",
            "renaldisease": "renal",
            "diabetes": "diabetes",
            "smoking": "smoking",
            "currentsmoker": "smoking",
            "eversmoker": "smoking",
            "cld": "cld",
            "chroniclungdisease": "cld",
            "chroniclung": "cld",
            "alcohol": "alcohol",
            "alcoholdrugs": "alcohol",
            "harmfulalcohol": "alcohol",
        }.get(name)
        if target is not None:
            out["disOR"][target] = _normalise_positive_or_override(value, raw_name)
    return out


def load_parameters_from_config(config: dict[str, Any]) -> dict[str, Any]:
    cfg = normalise_config(config)
    pars = load_tb_csv(cfg["csvFile"])
    age_file = resolve_age_distribution_file(cfg.get("ageDistributionFile"), cfg["csvFile"])
    pars = load_age_distribution(
        pars,
        age_file,
        cfg.get("ageDistributionTable"),
        int(cfg.get("age85PlusMax") or 89),
    )

    default_alcohol = [0.0, 0.0, 0.4]
    alcohol = pars.get("alcoholPrevByAge")
    if alcohol is None or all(_is_nan(v) for v in alcohol):
        pars["alcoholPrevByAge"] = default_alcohol
    else:
        pars["alcoholPrevByAge"] = [
            default_alcohol[i] if _is_nan(v) else float(v) for i, v in enumerate(alcohol)
        ]

    pars = apply_risk_prevalence_overrides(pars, cfg.get("riskPrev"))
    pars = apply_disease_or_overrides(pars, cfg.get("diseaseOR"))
    pars["totalFemalePrev"] = _safe_default(pars.get("totalFemalePrev"), 0.56)
    pars["totalBCGPrev"] = _safe_default(pars.get("totalBCGPrev"), 0.68)
    return refresh_total_prevalence_summaries(pars)


def _read_messy_csv(csv_file: str | Path) -> list[list[str]]:
    path = resolve_repo_path(csv_file)
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return [row for row in csv.reader(f)]


def _get_row_values(
    rows: list[list[str]],
    search_col: int,
    label: str,
    out_cols: int | list[int],
    allow_missing: bool = False,
) -> Any:
    row = next(
        (
            row
            for row in rows
            if search_col < len(row)
            and label.lower() in str(row[search_col]).strip().lower()
        ),
        None,
    )
    if row is None:
        if allow_missing:
            return math.nan if isinstance(out_cols, int) else [math.nan] * len(out_cols)
        raise ValueError(f'Could not find label "{label}" in column {search_col + 1}.')

    if isinstance(out_cols, int):
        return _parse_number(row[out_cols] if out_cols < len(row) else "")
    return [_parse_number(row[col] if col < len(row) else "") for col in out_cols]


def _parse_number(value: object) -> float:
    if value is None:
        return math.nan
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return float(value)
    text = str(value).strip().replace(",", "")
    if not text:
        return math.nan
    is_percent = text.endswith("%")
    if is_percent:
        text = text[:-1]
    try:
        number = float(text)
    except ValueError:
        match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text)
        if match is None:
            return math.nan
        number = float(match.group(0))
    return number / 100 if is_percent else number


def _weighted_age_lookup(values: list[float], groups: list[int], weights: list[float]) -> float:
    return sum(float(values[group - 1]) * float(weight) for group, weight in zip(groups, weights))


def _normalise_age_prevalence_override(value: Any, label: str) -> list[float]:
    values = _as_float_list(value)
    if len(values) == 1:
        values = values * 3
    if len(values) != 3 or any(_is_nan(v) or v < 0 or v > 1 for v in values):
        raise ValueError(
            f"Risk prevalence override for {label} must be a scalar or "
            "a 3-element vector [0-4, 5-14, >=15] within [0,1]."
        )
    return values


def _normalise_overall_prevalence_override(value: Any, label: str) -> float:
    values = _as_float_list(value)
    if len(values) != 1 or _is_nan(values[0]) or values[0] < 0 or values[0] > 1:
        raise ValueError(
            f"Overall prevalence override for {label} must be a scalar probability between 0 and 1."
        )
    return values[0]


def _normalise_positive_or_override(value: Any, label: str) -> float:
    values = _as_float_list(value)
    if len(values) != 1 or _is_nan(values[0]) or values[0] <= 0:
        raise ValueError(f"Disease OR override for {label} must be a positive scalar.")
    return values[0]


def _as_float_list(value: Any) -> list[float]:
    if isinstance(value, (list, tuple)):
        return [float(v) for v in value]
    return [float(value)]


def _normalise_name(value: object) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).strip().lower())


def _safe_default(value: Any, fallback: float) -> float:
    return fallback if _is_nan(value) else float(value)


def _is_nan(value: Any) -> bool:
    try:
        return math.isnan(float(value))
    except (TypeError, ValueError):
        return True


def _is_empty(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, (list, tuple)) and not value:
        return True
    if isinstance(value, (list, tuple)) and all(_is_nan(v) for v in value):
        return True
    return _is_nan(value) if not isinstance(value, (list, tuple, dict)) else False


def _has_rows(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, pd.DataFrame):
        return not value.empty
    if isinstance(value, list):
        return len(value) > 0
    return False
