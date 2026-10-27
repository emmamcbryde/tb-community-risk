"""Pure validation mathematics for proposed country TB contracts.

No calibration, population generation, projection, network access or file writes.
Canonical serialization is Python JSON v1, not a claim of RFC 8785 compliance.
"""

from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
import math
from typing import Any


CONTRACT_VERSION = "1.0.0"
FACTOR_IDS = (
    "cxr_granuloma", "heavy_smoking", "underweight", "diabetes",
    "tnf_inhibitor", "cxr_fibronodular", "head_neck_carcinoma",
    "haemodialysis", "silicosis", "hiv", "transplant_immunosuppression",
)
VALUE_STATES = ("value", "missing", "excluded", "not_applicable")


def _json_value(value: Any) -> None:
    if value is None or type(value) in (str, bool, int):
        return
    if type(value) is float and math.isfinite(value):
        return
    if type(value) is list:
        for item in value:
            _json_value(item)
        return
    if type(value) is dict and all(type(key) is str for key in value):
        for item in value.values():
            _json_value(item)
        return
    raise ValueError("Canonical input must contain JSON values and string keys only.")


def canonical_bytes(payload: Any) -> bytes:
    """Hash-free analytical payload; arrays retain their specified input order."""
    _json_value(payload)
    return json.dumps(
        payload, sort_keys=True, ensure_ascii=True, allow_nan=False,
        separators=(",", ":"),
    ).encode("ascii")


def content_hash(payload: Any) -> str:
    """No hidden removal of provenance, timestamps or other fields."""
    return hashlib.sha256(canonical_bytes(payload)).hexdigest()


def _number(value: Any, name: str) -> Fraction:
    if type(value) not in (int, float) or (type(value) is float and not math.isfinite(value)):
        raise ValueError(f"{name} must be a finite JSON number, not a boolean or string.")
    return Fraction(str(value))


def _probability(value: Any, name: str) -> Fraction:
    number = _number(value, name)
    if not 0 <= number <= 1:
        raise ValueError(f"{name} must be between zero and one.")
    return number


def validate_prevalence(row: dict[str, Any]) -> float | None:
    """Validate percentages and state; return an unrounded internal proportion."""
    if row.get("factor_id") not in FACTOR_IDS:
        raise ValueError("Unknown medical factor; recency/contact are not prevalence inputs.")
    state = row.get("state")
    if state not in VALUE_STATES:
        raise ValueError("Unknown prevalence state.")
    fields = ("value_percent", "lower_percent", "upper_percent")
    if any(field not in row for field in fields):
        raise ValueError("All percentage fields must be explicit, including nulls.")
    if state != "value":
        if any(row[field] is not None for field in fields):
            raise ValueError("Non-value states must carry null percentage fields.")
        if state in ("excluded", "not_applicable") and not str(row.get("reason") or "").strip():
            raise ValueError("Exclusion or non-applicability needs a reason.")
        return None
    value = _number(row["value_percent"], "value_percent")
    if not 0 <= value <= 100:
        raise ValueError("Prevalence must be between 0% and 100%.")
    lower, upper = row["lower_percent"], row["upper_percent"]
    if (lower is None) != (upper is None):
        raise ValueError("Supply both uncertainty bounds or neither.")
    if lower is not None:
        lo, hi = _number(lower, "lower_percent"), _number(upper, "upper_percent")
        if not 0 <= lo <= value <= hi <= 100:
            raise ValueError("Bounds must enclose prevalence within 0% to 100%.")
    return float(value / 100)


def validate_age_bands(bands: list[dict[str, Any]]) -> None:
    """Require an ordered partition [0, infinity), using exclusive upper ages."""
    if not bands:
        raise ValueError("Age bands must not be empty.")
    expected = 0
    for index, band in enumerate(bands):
        lo, hi = band.get("age_lower"), band.get("age_upper_exclusive")
        if type(lo) is not int or lo != expected:
            raise ValueError("Age bands must start at zero and have no gaps or overlaps.")
        if hi is None:
            if index != len(bands) - 1:
                raise ValueError("Only the last age band may be open.")
        elif type(hi) is not int or hi <= lo:
            raise ValueError("Upper age must exceed lower age.")
        if hi is not None and index == len(bands) - 1:
            raise ValueError("The last band must preserve the open oldest-age tail.")
        expected = hi


def reconcile_age_counts(
    source: list[dict[str, Any]], target: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Aggregate whole source cells; never split a coarse/open source cell."""
    validate_age_bands(source)
    validate_age_bands(target)
    with localcontext() as context:
        context.prec = 50
        counts = [_number(row["count"], "count") for row in source]
        if any(value < 0 for value in counts) or sum(counts) <= 0:
            raise ValueError("Counts must be nonnegative with a positive total.")
        totals = [Fraction(0) for _ in target]
        for row, count in zip(source, counts):
            matches = [
                i for i, band in enumerate(target)
                if row["age_lower"] >= band["age_lower"]
                and (
                    band["age_upper_exclusive"] is None
                    or (
                        row["age_upper_exclusive"] is not None
                        and row["age_upper_exclusive"] <= band["age_upper_exclusive"]
                    )
                )
            ]
            if len(matches) != 1:
                raise ValueError("A source cell crosses a target boundary; an allocation rule is required.")
            totals[matches[0]] += count
        total = sum(counts)
        return [
            {**band, "count": float(count), "proportion": float(count / total)}
            for band, count in zip(target, totals)
        ]


def _joint(p: Fraction, q: Fraction, overlap: Fraction) -> dict[str, Any]:
    low, high = max(Fraction(0), p + q - 1), min(p, q)
    if not low <= overlap <= high:
        raise ValueError("Overlap is outside the Frechet bounds; no clipping is performed.")
    cells = {
        "tbi_and_risk": overlap, "tbi_only": p - overlap,
        "risk_only": q - overlap, "neither": 1 - p - q + overlap,
    }
    if any(cell < 0 for cell in cells.values()) or sum(cells.values()) != 1:
        raise ValueError("Invalid joint distribution.")
    variance = p * (1 - p) * q * (1 - q)
    if variance == 0:
        phi = None
    else:
        covariance = overlap - p * q
        var_d = Decimal(variance.numerator) / Decimal(variance.denominator)
        cov_d = Decimal(covariance.numerator) / Decimal(covariance.denominator)
        phi = float(cov_d / var_d.sqrt())
    return {
        "tbi_prevalence": float(p), "any_risk_prevalence": float(q),
        "minimum_overlap": float(low), "maximum_overlap": float(high),
        "overlap": float(overlap), "phi": phi,
        "phi_status": "undefined_degenerate_marginal" if phi is None else "defined",
        "cells": {key: float(value) for key, value in cells.items()},
    }


def validate_overlap(p: float, q: float, overlap: float) -> dict[str, Any]:
    with localcontext() as context:
        context.prec = 50
        return _joint(
            _probability(p, "TBI"), _probability(q, "any factor"),
            _probability(overlap, "overlap"),
        )


def feasible_association(p: float, q: float, alpha: float) -> dict[str, Any]:
    """Normalized overlap scale, never an assertion that alpha equals phi."""
    with localcontext() as context:
        context.prec = 50
        p_d, q_d = _probability(p, "TBI"), _probability(q, "any factor")
        a_d = _number(alpha, "normalized association")
        if not -1 <= a_d <= 1:
            raise ValueError("Normalized association must lie in [-1,1].")
        low, high = max(Fraction(0), p_d + q_d - 1), min(p_d, q_d)
        independent = p_d * q_d
        overlap = independent + a_d * (
            independent - low if a_d < 0 else high - independent
        )
        return {"requested_alpha": float(a_d), **_joint(p_d, q_d, overlap)}


def association_summary(strata: list[dict[str, Any]]) -> dict[str, Any]:
    """Aggregate from demographic counts, including age-induced association."""
    if not strata:
        raise ValueError("At least one stratum is required.")
    with localcontext() as context:
        context.prec = 50
        counts = [_number(row["count"], "stratum count") for row in strata]
        if any(count < 0 for count in counts) or sum(counts) <= 0:
            raise ValueError("Stratum counts must be nonnegative with positive total.")
        total = sum(counts)
        weights = [count / total for count in counts]
        results = [feasible_association(row["p"], row["q"], row["alpha"]) for row in strata]
        # Reconstruct exact rational overlap, avoiding rounded JSON results.
        p_values = [_probability(row["p"], "TBI") for row in strata]
        q_values = [_probability(row["q"], "any factor") for row in strata]
        joints = []
        for row, p, q in zip(strata, p_values, q_values):
            a = _number(row["alpha"], "normalized association")
            j0 = p * q
            joints.append(j0 + a * (j0 - max(Fraction(0), p + q - 1) if a < 0 else min(p, q) - j0))
        p_all = sum(w * p for w, p in zip(weights, p_values))
        q_all = sum(w * q for w, q in zip(weights, q_values))
        j_all = sum(w * j for w, j in zip(weights, joints))
        expected_product = sum(w * p * q for w, p, q in zip(weights, p_values, q_values))
        return {
            "strata": results, "weights": [float(w) for w in weights],
            "overall": _joint(p_all, q_all, j_all),
            "within_age_covariance": float(j_all - expected_product),
            "age_driven_covariance": float(expected_product - p_all * q_all),
        }
