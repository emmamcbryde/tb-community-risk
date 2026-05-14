from __future__ import annotations

import math
import re
from typing import Iterable

import pandas as pd


def parse_age_band(value: object, age85_plus_max: int = 89) -> tuple[int, int]:
    text = str(value).strip()
    text = text.replace("–", "-").replace("—", "-")
    text = text.replace(">=", "").replace("≥", "")
    if text.endswith("+"):
        lo = int(float(text[:-1]))
        return lo, int(age85_plus_max)
    if "-" in text:
        lo_text, hi_text = text.split("-", 1)
        return int(float(lo_text)), int(float(hi_text))
    age = int(float(text))
    return age, age


def expand_age_distribution_table(
    df: pd.DataFrame, age85_plus_max: int = 89
) -> tuple[list[int], list[float]]:
    if df.empty:
        raise ValueError("Age distribution table must not be empty.")

    lower_columns = {str(col).lower(): col for col in df.columns}
    age_col = next((col for key, col in lower_columns.items() if "age" in key), df.columns[0])
    weight_col = next(
        (col for key, col in lower_columns.items() if "smoothed" in key),
        None,
    )
    if weight_col is None:
        weight_col = next(
            (col for key, col in lower_columns.items() if "proportion" in key),
            None,
        )
    if weight_col is None:
        raise ValueError("Could not identify a proportion column in age distribution table.")

    exact_ages: list[int] = []
    exact_prob: list[float] = []
    for _, row in df.iterrows():
        weight = _parse_number(row[weight_col])
        if weight is None or math.isnan(weight):
            continue
        lo, hi = parse_age_band(row[age_col], age85_plus_max)
        ages = list(range(lo, hi + 1))
        if not ages:
            continue
        per_age = weight / len(ages)
        exact_ages.extend(ages)
        exact_prob.extend([per_age] * len(ages))

    total = sum(exact_prob)
    if total <= 0:
        raise ValueError("Age distribution probabilities must sum to a positive value.")
    exact_prob = [p / total for p in exact_prob]
    return exact_ages, exact_prob


def build_fallback_exact_age_prob(
    exact_ages: Iterable[int], pop_frac_3: Iterable[float]
) -> list[float]:
    pop_frac = list(pop_frac_3)
    exact_age_values = [int(age) for age in exact_ages]
    adult_count = sum(1 for age in exact_age_values if age >= 15)
    prob = []
    for age in exact_age_values:
        if age <= 4:
            prob.append(pop_frac[0] / 5)
        elif 5 <= age <= 14:
            prob.append(pop_frac[1] / 10)
        else:
            prob.append(pop_frac[2] / adult_count)
    total = sum(prob)
    return [p / total for p in prob]


def broad_age_group_from_years(age_years: Iterable[int | float]) -> list[int]:
    groups = []
    for age in age_years:
        age_value = float(age)
        if age_value <= 4:
            groups.append(1)
        elif age_value <= 14:
            groups.append(2)
        else:
            groups.append(3)
    return groups


def _parse_number(value: object) -> float | None:
    if value is None:
        return None
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return float(value)
    text = str(value).strip().replace(",", "")
    if not text or text.lower() == "nan":
        return None
    is_percent = text.endswith("%")
    if is_percent:
        text = text[:-1]
    try:
        number = float(text)
    except ValueError:
        match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text)
        if not match:
            return None
        number = float(match.group(0))
    return number / 100 if is_percent else number
