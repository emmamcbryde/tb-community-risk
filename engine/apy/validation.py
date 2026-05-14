from __future__ import annotations

from math import isfinite
from typing import Any

from engine.apy.config import normalise_config, resolve_repo_path


REQUIRED_FIELDS = [
    "csvFile",
    "N",
    "nReps",
    "seed",
    "screenWindow",
    "followHorizon",
    "screenCoverage",
    "screeningStrategy",
    "ageDistributionFile",
    "ageDistributionTable",
    "age85PlusMax",
    "targetAgeOR",
    "riskPrev",
    "diseaseOR",
    "testType",
    "regimen",
]

TEST_TYPES = {"IGRA", "TST"}
REGIMENS = {"3HP", "4R", "3HR", "6H", "9H"}
SCREENING_STRATEGIES = {"random", "ltbi", "cure", "prevent"}
PARTIAL_SHORT_COURSE_MODES = {"threshold80", "linear", "none"}
RISK_PREV_FIELDS = ["contact", "MJ", "renal", "diabetes", "smoking", "cld", "alcohol"]
DISEASE_OR_FIELDS = ["contact", "MJ", "renal", "diabetes", "smoking", "cld", "alcohol"]


def collect_validation_issues(config: dict[str, Any]) -> dict[str, Any]:
    report = _init_report()
    if not isinstance(config, dict):
        _add_issue(
            report,
            "errors",
            "config",
            "Config",
            "invalid_type",
            "Config must be a dictionary.",
        )
        return _finalise_report(report)

    cfg = normalise_config(config)
    for field in REQUIRED_FIELDS:
        if field not in cfg:
            _add_issue(
                report,
                "errors",
                field,
                field,
                "required_missing",
                f"Required field is missing: {field}",
            )

    _validate_existing_file(report, cfg, "csvFile", "Main data file", required=True)
    _validate_existing_file(
        report, cfg, "ageDistributionFile", "Age distribution file", required=False
    )
    _validate_positive_scalar(report, cfg, "N", "Cohort size")
    _validate_positive_scalar(report, cfg, "nReps", "Number of replicates")
    _validate_nonnegative_scalar(report, cfg, "seed", "Random seed")
    _validate_positive_scalar(report, cfg, "screenWindow", "Screen window")
    _validate_positive_scalar(report, cfg, "followHorizon", "Follow-up horizon")
    _validate_fraction(report, cfg, "screenCoverage", "Screening coverage")
    _validate_positive_scalar(report, cfg, "age85PlusMax", "Age 85+ upper bound")
    _validate_positive_scalar(report, cfg, "targetAgeOR", "Age OR target")

    if _is_number(cfg.get("followHorizon")) and _is_number(cfg.get("screenWindow")):
        if float(cfg["followHorizon"]) <= float(cfg["screenWindow"]):
            _add_issue(
                report,
                "errors",
                "followHorizon",
                "Follow-up horizon",
                "invalid_range",
                "followHorizon must be greater than screenWindow.",
            )

    _validate_choice(report, cfg, "testType", "Test choice", TEST_TYPES)
    _validate_choice(report, cfg, "regimen", "Treatment choice", REGIMENS)
    _validate_choice(
        report, cfg, "screeningStrategy", "Screening strategy", SCREENING_STRATEGIES
    )
    _validate_optional_choice(
        report,
        cfg,
        "partialShortCourseMode",
        "Partial-course efficacy rule",
        PARTIAL_SHORT_COURSE_MODES,
    )

    for field, label in [
        ("pStartTPT", "Treatment start probability"),
        ("regimenPComplete", "Completion probability"),
        ("regimenADRstop", "ADR stop probability"),
        ("regimenEffFull", "Full-course efficacy"),
        ("partialDoseFractionADR", "Partial dose fraction after ADR stop"),
        ("partialDoseFractionOther", "Partial dose fraction after other stop"),
    ]:
        _validate_fraction(report, cfg, field, label, allow_empty=True)

    _validate_age_distribution_consistency(report, cfg)
    _validate_risk_prevalence_struct(report, cfg)
    _validate_optional_fraction(
        report, cfg.get("riskPrev"), "female", "Female prevalence", "riskPrev.female"
    )
    _validate_optional_fraction(
        report, cfg.get("riskPrev"), "BCG", "BCG prevalence", "riskPrev.BCG"
    )
    _validate_disease_or_struct(report, cfg)

    if _is_number(cfg.get("nReps")) and float(cfg["nReps"]) > 1000:
        _add_issue(
            report,
            "warnings",
            "nReps",
            "Number of replicates",
            "browser_large_nreps",
            "nReps is large for browser use and may lead to a slow UI response.",
        )

    return _finalise_report(report)


def validate_config(config: dict[str, Any]) -> dict[str, Any]:
    cfg = normalise_config(config)
    report = collect_validation_issues(cfg)
    if not report["isValid"]:
        fields = ", ".join(report["fatalFieldNames"])
        raise ValueError(f"Invalid APY config: {fields}")
    return cfg


def _validate_age_distribution_consistency(
    report: dict[str, Any], cfg: dict[str, Any]
) -> None:
    age_file = cfg.get("ageDistributionFile")
    age_table = cfg.get("ageDistributionTable")
    has_age_file = bool(str(age_file).strip()) if age_file is not None else False
    has_age_table = bool(age_table)
    if has_age_file and has_age_table:
        _add_issue(
            report,
            "errors",
            "ageDistributionTable",
            "Age distribution input",
            "conflict",
            "Provide either ageDistributionFile or ageDistributionTable, not both.",
        )


def _validate_risk_prevalence_struct(
    report: dict[str, Any], cfg: dict[str, Any]
) -> None:
    risk_prev = cfg.get("riskPrev")
    if not isinstance(risk_prev, dict):
        _add_issue(
            report,
            "errors",
            "riskPrev",
            "Risk-factor prevalences",
            "invalid_type",
            "riskPrev must be a dictionary.",
        )
        return

    labels = {
        "contact": "Close contact prevalence",
        "MJ": "Marijuana use prevalence",
        "renal": "Renal disease prevalence",
        "diabetes": "Diabetes prevalence",
        "smoking": "Smoking prevalence",
        "cld": "Chronic lung disease prevalence",
        "alcohol": "Alcohol / drugs prevalence",
    }
    for field in RISK_PREV_FIELDS:
        value = risk_prev.get(field)
        if value is None:
            continue
        values = _numeric_vector(value)
        dotted = f"riskPrev.{field}"
        if values is None:
            _add_issue(
                report,
                "errors",
                dotted,
                labels[field],
                "invalid_type",
                "Value must be numeric.",
            )
        elif len(values) not in {1, 3}:
            _add_issue(
                report,
                "errors",
                dotted,
                labels[field],
                "invalid_shape",
                "Value must be either a scalar or a 3-element vector.",
            )
        elif any(not isfinite(v) or v < 0 or v > 1 for v in values):
            _add_issue(
                report,
                "errors",
                dotted,
                labels[field],
                "invalid_range",
                "Value must be within [0,1].",
            )


def _validate_disease_or_struct(report: dict[str, Any], cfg: dict[str, Any]) -> None:
    disease_or = cfg.get("diseaseOR")
    if not isinstance(disease_or, dict):
        _add_issue(
            report,
            "errors",
            "diseaseOR",
            "Disease odds ratios",
            "invalid_type",
            "diseaseOR must be a dictionary.",
        )
        return

    for field in DISEASE_OR_FIELDS:
        _validate_optional_positive_scalar(
            report,
            disease_or,
            field,
            f"{field} progression OR",
            f"diseaseOR.{field}",
        )


def _validate_existing_file(
    report: dict[str, Any],
    parent: dict[str, Any],
    field: str,
    label: str,
    *,
    required: bool,
) -> None:
    value = parent.get(field)
    if value is None or str(value).strip() == "":
        if required:
            _add_issue(
                report,
                "errors",
                field,
                label,
                "required_empty",
                f"{label} must not be empty.",
            )
        return
    if not resolve_repo_path(str(value)).is_file():
        _add_issue(
            report,
            "errors",
            field,
            label,
            "file_missing",
            f"{label} not found: {value}",
        )


def _validate_positive_scalar(
    report: dict[str, Any], parent: dict[str, Any], field: str, label: str
) -> None:
    value = parent.get(field)
    if not _is_number(value) or float(value) <= 0:
        _add_issue(
            report,
            "errors",
            field,
            label,
            "invalid_range",
            f"{label} must be a positive scalar.",
        )


def _validate_optional_positive_scalar(
    report: dict[str, Any],
    parent: dict[str, Any],
    field: str,
    label: str,
    dotted: str,
) -> None:
    value = parent.get(field) if isinstance(parent, dict) else None
    if value is None:
        return
    if not _is_number(value) or float(value) <= 0:
        _add_issue(
            report,
            "errors",
            dotted,
            label,
            "invalid_range",
            f"{label} must be a positive scalar.",
        )


def _validate_nonnegative_scalar(
    report: dict[str, Any], parent: dict[str, Any], field: str, label: str
) -> None:
    value = parent.get(field)
    if not _is_number(value) or float(value) < 0:
        _add_issue(
            report,
            "errors",
            field,
            label,
            "invalid_range",
            f"{label} must be a non-negative scalar.",
        )


def _validate_fraction(
    report: dict[str, Any],
    parent: dict[str, Any],
    field: str,
    label: str,
    *,
    allow_empty: bool = False,
) -> None:
    value = parent.get(field)
    if value is None:
        if allow_empty:
            return
        _add_issue(
            report,
            "errors",
            field,
            label,
            "required_empty",
            f"{label} must not be empty.",
        )
        return
    if not _is_number(value) or not 0 <= float(value) <= 1:
        _add_issue(
            report,
            "errors",
            field,
            label,
            "invalid_range",
            f"{label} must be a scalar in [0,1].",
        )


def _validate_optional_fraction(
    report: dict[str, Any],
    parent: Any,
    field: str,
    label: str,
    dotted: str,
) -> None:
    value = parent.get(field) if isinstance(parent, dict) else None
    if value is None:
        return
    if not _is_number(value) or not 0 <= float(value) <= 1:
        _add_issue(
            report,
            "errors",
            dotted,
            label,
            "invalid_range",
            f"{label} must be a scalar in [0,1].",
        )


def _validate_choice(
    report: dict[str, Any],
    parent: dict[str, Any],
    field: str,
    label: str,
    allowed: set[str],
) -> None:
    value = parent.get(field)
    if value is None:
        _add_issue(
            report,
            "errors",
            field,
            label,
            "required_empty",
            f"{label} must not be empty.",
        )
        return
    if str(value) not in allowed:
        _add_issue(
            report,
            "errors",
            field,
            label,
            "invalid_choice",
            f"{label} must be one of: {', '.join(sorted(allowed))}.",
        )


def _validate_optional_choice(
    report: dict[str, Any],
    parent: dict[str, Any],
    field: str,
    label: str,
    allowed: set[str],
) -> None:
    value = parent.get(field)
    if value is None:
        return
    if str(value) not in allowed:
        _add_issue(
            report,
            "errors",
            field,
            label,
            "invalid_choice",
            f"{label} must be one of: {', '.join(sorted(allowed))}.",
        )


def _numeric_vector(value: Any) -> list[float] | None:
    if isinstance(value, (str, bytes, dict)):
        return None
    if isinstance(value, (list, tuple)):
        values = value
    else:
        values = [value]
    try:
        return [float(v) for v in values]
    except (TypeError, ValueError):
        return None


def _is_number(value: Any) -> bool:
    if isinstance(value, bool) or value is None:
        return False
    try:
        return isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _init_report() -> dict[str, Any]:
    return {
        "isValid": True,
        "hasWarnings": False,
        "errors": [],
        "warnings": [],
        "infos": [],
        "fatalFieldNames": [],
        "warningFieldNames": [],
    }


def _add_issue(
    report: dict[str, Any],
    bucket: str,
    field: str,
    label: str,
    code: str,
    message: str,
) -> None:
    severity = "error" if bucket == "errors" else "warning"
    report[bucket].append(
        {
            "fieldName": field,
            "fieldLabel": label,
            "severity": severity,
            "code": code,
            "message": message,
        }
    )


def _finalise_report(report: dict[str, Any]) -> dict[str, Any]:
    report["fatalFieldNames"] = [issue["fieldName"] for issue in report["errors"]]
    report["warningFieldNames"] = [
        issue["fieldName"] for issue in report["warnings"]
    ]
    report["isValid"] = not report["errors"]
    report["hasWarnings"] = bool(report["warnings"])
    return report
