"""Pure helpers translating between population profiles and editable tables."""

from __future__ import annotations

import csv
from dataclasses import replace
import io
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
TRANSITION_LABELS = {
    "progression_to_disease": "Progression to disease",
    "infection": "Infection",
    "other": "Other",
}
TRANSITION_CODES = {label: code for code, label in TRANSITION_LABELS.items()}
RISK_FACTOR_COLUMNS = [
    "Enabled",
    "Risk factor",
    "Status",
    "Prevalence (%)",
    "Prevalence low (%)",
    "Prevalence high (%)",
    "Effect estimate",
    "Effect low",
    "Effect high",
    "Effect measure",
    "Affected transition",
    "Evidence year",
    "Applies to",
    "Source",
    "Review status",
    "Notes",
    "id",
]
CSV_COLUMNS = [column for column in RISK_FACTOR_COLUMNS if column != "Status"]


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
                "Prevalence low (%)": None if factor.prevalence_bounds is None else _percent(factor.prevalence_bounds[0]),
                "Prevalence high (%)": None if factor.prevalence_bounds is None else _percent(factor.prevalence_bounds[1]),
                "Effect estimate": factor.effect_estimate.value,
                "Effect low": None if factor.effect_bounds is None else factor.effect_bounds[0],
                "Effect high": None if factor.effect_bounds is None else factor.effect_bounds[1],
                "Effect measure": None if factor.effect_measure is None else factor.effect_measure.value,
                "Affected transition": TRANSITION_LABELS.get(factor.affected_transition, factor.affected_transition),
                "Evidence year": factor.evidence_year,
                "Applies to": factor.population_applicability,
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

    prevalence_bounds = _pair(row.get("Prevalence low (%)"), row.get("Prevalence high (%)"), scale=0.01, label="Prevalence")
    if not _same_pair(prevalence_bounds, factor.prevalence_bounds):
        changes["prevalence_bounds"] = prevalence_bounds
        mark("prevalenceBounds")
    effect_bounds = _pair(row.get("Effect low"), row.get("Effect high"), scale=1.0, label="Effect")
    if not _same_pair(effect_bounds, factor.effect_bounds):
        changes["effect_bounds"] = effect_bounds
        mark("effectBounds")

    transition_text = _text(row.get("Affected transition"))
    transition = TRANSITION_CODES.get(transition_text, transition_text or factor.affected_transition)
    if transition != factor.affected_transition:
        changes["affected_transition"] = transition
        mark("affectedTransition")

    year_value = _number(row.get("Evidence year"))
    evidence_year = None if year_value is None else int(year_value)
    if evidence_year != factor.evidence_year:
        changes["evidence_year"] = evidence_year
        mark("evidenceYear")

    applicability = _text(row.get("Applies to"))
    if applicability != factor.population_applicability:
        changes["population_applicability"] = applicability
        mark("applicability")

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


def risk_factor_editing_issues(profile: PopulationProfile) -> list[str]:
    """Engine-specific checks on edited risk factors (in addition to validate_profile)."""
    issues = []
    for factor in profile.risk_factors:
        if factor.engine_key is not None and factor.affected_transition != "progression_to_disease":
            issues.append(
                f"{factor.label}: the current engine applies this factor to progression to disease only; "
                "change 'Affected transition' back or add it as a separate custom factor."
            )
    return issues


def risk_factors_csv(profile: PopulationProfile) -> str:
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=CSV_COLUMNS, lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    for row in risk_factor_rows(profile):
        writer.writerow({key: "" if row.get(key) is None else row.get(key) for key in CSV_COLUMNS})
    return buffer.getvalue()


def risk_factor_template_csv() -> str:
    return ",".join(CSV_COLUMNS) + "\n"


def parse_risk_factor_csv(payload: bytes) -> tuple[list[dict[str, Any]], list[str]]:
    """Read an exported or template CSV into editor rows; returns (rows, errors)."""
    try:
        text = payload.decode("utf-8-sig")
    except UnicodeDecodeError:
        return [], ["The file is not UTF-8 encoded text."]
    reader = csv.DictReader(io.StringIO(text, newline=""))
    header = reader.fieldnames or []
    missing = [column for column in ("Risk factor", "Prevalence (%)", "Effect estimate", "Effect measure") if column not in header]
    if missing:
        return [], [f"Missing column(s): {', '.join(missing)}. Use the template or an exported table."]
    rows, errors = [], []
    for line, raw in enumerate(reader, start=2):
        row = {key: (value if value != "" else None) for key, value in raw.items() if key in CSV_COLUMNS}
        enabled = str(raw.get("Enabled", "True")).strip().lower()
        row["Enabled"] = enabled not in {"false", "0", "no"}
        measure = _text(row.get("Effect measure")).upper()
        if measure and measure not in EFFECT_MEASURE_OPTIONS:
            errors.append(f"Line {line}: effect measure must be RR, HR or OR (found {measure!r}).")
        row["Effect measure"] = measure or None
        for column in ("Prevalence (%)", "Prevalence low (%)", "Prevalence high (%)", "Effect estimate", "Effect low", "Effect high", "Evidence year"):
            try:
                _number(row.get(column))
            except ValueError:
                errors.append(f"Line {line}: {column} is not a number.")
        rows.append(row)
    return rows, errors


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


def _pair(low: Any, high: Any, *, scale: float, label: str) -> tuple[float, float] | None:
    low_value, high_value = _number(low), _number(high)
    if low_value is None and high_value is None:
        return None
    if low_value is None or high_value is None:
        raise ValueError(f"{label} bounds need both a low and a high value.")
    return (low_value * scale, high_value * scale)


def _same_pair(a: tuple[float, float] | None, b: tuple[float, float] | None) -> bool:
    if a is None or b is None:
        return a is None and b is None
    return _same_number(a[0], b[0]) and _same_number(a[1], b[1])


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

