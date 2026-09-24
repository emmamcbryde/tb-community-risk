"""Pure helpers translating between population profiles and editable tables."""

from __future__ import annotations

from dataclasses import replace
import math
from typing import Any

from app.general.terminology import USER_DEFINED_MARK
from engine.profiles.population_profile import (
    EffectMeasure,
    PopulationProfile,
    ProfileValue,
    ReviewStatus,
    RiskFactor,
    ValueState,
    user_value,
)


REVIEW_LABELS = {
    ReviewStatus.NOT_REVIEWED: "Not reviewed",
    ReviewStatus.REVIEWED: "Reviewed",
    ReviewStatus.NOT_REQUIRED: "Not required",
}
REVIEW_CODES = {label: code for code, label in REVIEW_LABELS.items()}
EFFECT_MEASURE_OPTIONS = [member.value for member in EffectMeasure]
RISK_FACTOR_COLUMNS = [
    "Enabled",
    "Risk factor",
    "Prevalence (%)",
    "Effect estimate",
    "Effect measure",
    "Source",
    "Review status",
    "Status",
    "Notes",
    "id",
]


def value_status_text(value: ProfileValue) -> str:
    if value.state is ValueState.MISSING:
        return "Missing"
    if value.state is ValueState.NOT_APPLICABLE:
        return "Not applicable"
    if value.state is ValueState.EXCLUDED:
        return "Excluded"
    parts = ["Zero"] if value.value == 0 else []
    parts.append("User-defined" if value.is_user_override else "Bundled")
    if value.review_status is ReviewStatus.NOT_REVIEWED:
        parts.append("not reviewed")
    return ", ".join(parts)


def risk_factor_status(factor: RiskFactor) -> str:
    if not factor.enabled:
        base = "Not included"
    elif factor.prevalence.state is ValueState.MISSING or factor.effect_estimate.state is ValueState.MISSING:
        base = "Missing value"
    elif factor.prevalence.is_zero:
        base = "Zero prevalence"
    elif factor.review_status is ReviewStatus.REVIEWED:
        base = "Reviewed"
    else:
        base = "Not reviewed"
    if factor.engine_key is None:
        base += "; not used by the current model"
    return f"{USER_DEFINED_MARK} - {base}" if factor.is_user_override else base


def risk_factor_rows(profile: PopulationProfile) -> list[dict[str, Any]]:
    rows = []
    for factor in profile.risk_factors:
        rows.append(
            {
                "Enabled": factor.enabled,
                "Risk factor": factor.label,
                "Prevalence (%)": _percent(factor.prevalence.value),
                "Effect estimate": factor.effect_estimate.value,
                "Effect measure": None if factor.effect_measure is None else factor.effect_measure.value,
                "Source": factor.evidence_source,
                "Review status": REVIEW_LABELS[factor.review_status],
                "Status": risk_factor_status(factor),
                "Notes": factor.notes,
                "id": factor.risk_factor_id,
            }
        )
    return rows


def apply_risk_factor_rows(profile: PopulationProfile, rows: list[dict[str, Any]]) -> PopulationProfile:
    """Apply edited table rows, marking every changed field as user-defined.

    Blank numeric cells become *missing*; they are never converted to zero.
    Effect-measure types are stored exactly as selected.
    """
    existing = {factor.risk_factor_id: factor for factor in profile.risk_factors}
    updated: list[RiskFactor] = []
    used_ids: set[str] = set()
    for index, row in enumerate(rows):
        factor_id = _text(row.get("id"))
        factor = existing.get(factor_id) if factor_id else None
        if factor is None:
            label = _text(row.get("Risk factor"))
            if not label and _blank_number(row.get("Prevalence (%)")) and _blank_number(row.get("Effect estimate")):
                continue
            factor_id = _new_custom_id(label or f"risk_factor_{index + 1}", used_ids | set(existing))
            factor = RiskFactor(
                risk_factor_id=factor_id,
                label=label or "Custom risk factor",
                enabled=True,
                prevalence=user_value(None, "proportion"),
                effect_estimate=user_value(None, "ratio"),
                effect_measure=None,
                evidence_source="User-defined",
                user_modified_fields=("added",),
            )
        used_ids.add(factor.risk_factor_id)
        updated.append(_apply_row(factor, row))
    return replace(profile, risk_factors=tuple(updated))


def _apply_row(factor: RiskFactor, row: dict[str, Any]) -> RiskFactor:
    modified = list(factor.user_modified_fields)
    changes: dict[str, Any] = {}

    def mark(name: str) -> None:
        if name not in modified:
            modified.append(name)

    enabled = bool(row.get("Enabled")) if row.get("Enabled") is not None else factor.enabled
    if enabled != factor.enabled:
        changes["enabled"] = enabled
        mark("enabled")

    label = _text(row.get("Risk factor")) or factor.label
    if label != factor.label:
        changes["label"] = label
        mark("label")

    prevalence = _number(row.get("Prevalence (%)"))
    new_prevalence = None if prevalence is None else prevalence / 100.0
    if not _same_number(new_prevalence, factor.prevalence.value):
        changes["prevalence"] = user_value(new_prevalence, "proportion")

    effect = _number(row.get("Effect estimate"))
    if not _same_number(effect, factor.effect_estimate.value):
        changes["effect_estimate"] = user_value(effect, "ratio")

    measure_text = _text(row.get("Effect measure"))
    measure = EffectMeasure(measure_text) if measure_text else None
    if measure != factor.effect_measure:
        changes["effect_measure"] = measure
        mark("effectMeasure")

    source = _text(row.get("Source"))
    if source != factor.evidence_source:
        changes["evidence_source"] = source
        mark("source")

    notes = _text(row.get("Notes"))
    if notes != factor.notes:
        changes["notes"] = notes
        mark("notes")

    review = REVIEW_CODES.get(_text(row.get("Review status")), factor.review_status)
    if review != factor.review_status:
        changes["review_status"] = review
        mark("reviewStatus")
        for key in ("prevalence", "effect_estimate"):
            value = changes.get(key, getattr(factor, key))
            changes[key] = replace(value, review_status=review)

    if not changes:
        return factor
    changes["user_modified_fields"] = tuple(modified)
    return replace(factor, **changes)


def population_status_text(profile: PopulationProfile) -> str:
    size = profile.population_size
    return USER_DEFINED_MARK if size.is_user_override else "Demonstration default"


def _percent(value: float | None) -> float | None:
    return None if value is None else round(value * 100.0, 10)


def _number(value: Any) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, str):
        text = value.strip().rstrip("%").strip()
        if not text:
            return None
        try:
            value = float(text)
        except ValueError as exc:
            raise ValueError(f"Not a number: {value!r}") from exc
    number = float(value)
    return None if math.isnan(number) else number


def _blank_number(value: Any) -> bool:
    try:
        return _number(value) is None
    except ValueError:
        return False


def _same_number(a: float | None, b: float | None) -> bool:
    if a is None or b is None:
        return a is None and b is None
    return math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12)


def _text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return str(value).strip()


def _new_custom_id(label: str, used: set[str]) -> str:
    base = "custom_" + "".join(ch.lower() if ch.isalnum() else "_" for ch in label).strip("_")
    candidate, counter = base, 2
    while candidate in used:
        candidate = f"{base}_{counter}"
        counter += 1
    return candidate

