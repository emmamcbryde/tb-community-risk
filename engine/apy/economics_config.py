from __future__ import annotations

from math import isfinite
from numbers import Number
from typing import Any


def build_default_economics_config() -> dict[str, Any]:
    return {
        "metadata": {
            "currencyCode": "",
            "priceYear": None,
            "locationLabel": "",
            "sourceNotes": "",
            "programCostBasis": "",
        },
        "costs": {
            "test": {
                "IGRA": None,
                "TST": None,
            },
            "regimen": {
                "x3HP": None,
                "x4R": None,
                "x3HR": None,
                "x6H": None,
                "x9H": None,
            },
            "falsePositiveIncrementalPerPerson": None,
            "activeTBDiseasePerCase": None,
            "programSetupTotal": None,
            "programRunningTotal": None,
        },
    }


def build_economics_preset_kwab150() -> dict[str, Any]:
    config = build_default_economics_config()
    config["metadata"].update(
        {
            "currencyCode": "AUD",
            "priceYear": 2019,
            "locationLabel": "Australia",
            "programCostBasis": "total",
            "sourceNotes": (
                "KWAB150 preset populated from local data/costs.csv: "
                "cscreenqft, cscreentst, ctreat*, and ctb mid values. "
                "Program setup/running and false-positive incremental costs "
                "are not specified in this preset."
            ),
        }
    )
    config["costs"]["test"].update(
        {
            "IGRA": 113.48,
            "TST": 116.07,
        }
    )
    config["costs"]["regimen"].update(
        {
            "x3HP": 165.5072,
            "x4R": 123.3172,
            "x3HR": 134.2272,
            "x6H": 187.7508,
            "x9H": 254.8544,
        }
    )
    config["costs"]["activeTBDiseasePerCase"] = 19079.6
    return config


def validate_economics_config(config: dict[str, Any]) -> dict[str, Any]:
    report = _init_report()

    if not isinstance(config, dict) or not config:
        _add_issue(report, "econConfig", "error", "Economics config must be a struct.")
        return _finalise_report(report)

    metadata = config.get("metadata")
    if not isinstance(metadata, dict):
        _add_issue(report, "metadata", "error", "metadata must be a struct.")
    else:
        _validate_text_field(report, metadata, "currencyCode", "metadata.currencyCode")
        _validate_optional_numeric_scalar(
            report, metadata, "priceYear", "metadata.priceYear"
        )
        _validate_text_field(report, metadata, "locationLabel", "metadata.locationLabel")
        _validate_text_field(report, metadata, "sourceNotes", "metadata.sourceNotes")
        _validate_text_field(
            report, metadata, "programCostBasis", "metadata.programCostBasis"
        )

    costs = config.get("costs")
    if not isinstance(costs, dict):
        _add_issue(report, "costs", "error", "costs must be a struct.")
    else:
        test_costs = costs.get("test")
        if not isinstance(test_costs, dict):
            _add_issue(report, "costs.test", "error", "costs.test must be a struct.")
        else:
            _validate_optional_cost(report, test_costs, "IGRA", "costs.test.IGRA")
            _validate_optional_cost(report, test_costs, "TST", "costs.test.TST")

        regimen_costs = costs.get("regimen")
        if not isinstance(regimen_costs, dict):
            _add_issue(
                report,
                "costs.regimen",
                "error",
                "costs.regimen must be a struct.",
            )
        else:
            for field in ["x3HP", "x4R", "x3HR", "x6H", "x9H"]:
                _validate_optional_cost(
                    report,
                    regimen_costs,
                    field,
                    f"costs.regimen.{field}",
                )

        _validate_optional_cost(
            report,
            costs,
            "falsePositiveIncrementalPerPerson",
            "costs.falsePositiveIncrementalPerPerson",
        )
        _validate_optional_cost(
            report, costs, "activeTBDiseasePerCase", "costs.activeTBDiseasePerCase"
        )
        _validate_optional_cost(
            report, costs, "programSetupTotal", "costs.programSetupTotal"
        )
        _validate_optional_cost(
            report, costs, "programRunningTotal", "costs.programRunningTotal"
        )

    return _finalise_report(report)


def _validate_text_field(
    report: dict[str, Any],
    parent: dict[str, Any],
    field: str,
    full_name: str,
) -> None:
    value = parent.get(field)
    if _is_blank(value):
        return
    if not isinstance(value, str):
        _add_issue(report, full_name, "error", f"{full_name} must be text.")


def _validate_optional_numeric_scalar(
    report: dict[str, Any],
    parent: dict[str, Any],
    field: str,
    full_name: str,
) -> None:
    value = parent.get(field)
    if _is_blank(value):
        return
    if not _is_finite_number(value):
        _add_issue(report, full_name, "error", f"{full_name} must be a finite scalar.")


def _validate_optional_cost(
    report: dict[str, Any],
    parent: dict[str, Any],
    field: str,
    full_name: str,
) -> None:
    value = parent.get(field)
    if _is_blank(value):
        return
    if not _is_finite_number(value) or float(value) < 0:
        _add_issue(
            report,
            full_name,
            "error",
            f"{full_name} must be empty or a non-negative scalar.",
        )


def _is_blank(value: Any) -> bool:
    return value is None or value == "" or (isinstance(value, list) and len(value) == 0)


def _is_finite_number(value: Any) -> bool:
    if isinstance(value, bool) or not isinstance(value, Number):
        return False
    return isfinite(float(value))


def _add_issue(
    report: dict[str, Any],
    field: str,
    severity: str,
    message: str,
) -> None:
    issue = {"field": field, "severity": severity, "message": message}
    if severity == "error":
        report["errors"].append(issue)
    else:
        report["warnings"].append(issue)


def _init_report() -> dict[str, Any]:
    return {
        "isValid": True,
        "hasWarnings": False,
        "errors": [],
        "warnings": [],
    }


def _finalise_report(report: dict[str, Any]) -> dict[str, Any]:
    report["isValid"] = not report["errors"]
    report["hasWarnings"] = bool(report["warnings"])
    return report
