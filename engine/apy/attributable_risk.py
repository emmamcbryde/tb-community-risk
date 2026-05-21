from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


KNOWN_MATLAB_ATTRIBUTABLE_METRICS = (
    "ORused",
    "PrevalencePopulation",
    "ExpectedPrevalenceAmongInfected",
    "ExpectedRisk20yAmongInfected_Exposed",
    "ExpectedRisk20yAmongInfected_Exposed_NoFactor",
    "AttributableRisk20y_AmongExposedInfected",
    "AttributableFraction20y_AmongExposedInfected",
    "PopulationAttributableFraction20y_ReactivationOnly",
    "ExpectedAttributableCases20y_Per1500",
    "ExpectedActiveCases20y_Per1500_DoNothing",
)

_SOURCE = "apy_python_attributable_risk_scaffold"


def run_attributable_risk(
    results: Mapping[str, Any],
    requested_metrics: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Return the Python APY attributable-risk scaffold contract.

    Reactivation attributable-risk calculations are intentionally not ported
    yet. This function reports the MATLAB table metrics as unsupported and
    validates only that the dynamic-comparison section is present.
    """
    missing_inputs: list[dict[str, str]] = []
    messages: list[str] = []

    if _dynamic_comparison(results) is None:
        missing_inputs.append(
            {
                "field": "technical.dynamicComparison",
                "message": (
                    "Attributable-risk scaffold expects "
                    "technical.dynamicComparison to be present in the results "
                    "bundle, even though dynamic active-TB comparison fields "
                    "are not attributable-risk metrics."
                ),
            }
        )
        messages.append(
            "Missing technical.dynamicComparison; no attributable-risk rows calculated."
        )

    metrics = _metrics_to_report(requested_metrics)
    calculated_rows, rows_source = _precomputed_rows(results)
    unsupported_metrics = [
        {
            "metric": metric,
            "reason": (
                "Reactivation attributable-risk calculation is not implemented "
                "in the Python APY port yet."
            ),
            "matlabSource": "run_tb_screening_reactivation_attributable_v9",
        }
        for metric in metrics
    ]

    if unsupported_metrics:
        messages.append(
            "Reactivation attributable-risk metrics are explicitly unsupported "
            "in this Python scaffold."
        )

    if calculated_rows:
        messages.append(
            f"Precomputed attributable-risk rows were passed through from {rows_source}."
        )

    return {
        "status": _status(missing_inputs, calculated_rows),
        "source": _SOURCE,
        "calculatedRows": calculated_rows,
        "missingInputs": missing_inputs,
        "unsupportedMetrics": unsupported_metrics,
        "messages": messages,
    }


def _dynamic_comparison(results: Mapping[str, Any]) -> Any:
    technical = results.get("technical")
    if not isinstance(technical, Mapping):
        return None
    return technical.get("dynamicComparison")


def _metrics_to_report(requested_metrics: Sequence[str] | None) -> list[str]:
    if requested_metrics is None:
        return list(KNOWN_MATLAB_ATTRIBUTABLE_METRICS)
    return [str(metric) for metric in requested_metrics]


def _precomputed_rows(results: Mapping[str, Any]) -> tuple[list[Any], str | None]:
    attributable_risk = results.get("attributableRisk")
    if not isinstance(attributable_risk, Mapping):
        return [], None

    for field in ("calculatedRows", "attributableRows"):
        rows = attributable_risk.get(field)
        if isinstance(rows, list) and rows:
            return [_json_like(row) for row in rows], f"attributableRisk.{field}"

    for field in ("calculatedRows", "attributableRows"):
        rows = attributable_risk.get(field)
        if isinstance(rows, list):
            return [], f"attributableRisk.{field}"

    return [], None


def _json_like(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_like(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_json_like(item) for item in value]
    if isinstance(value, list):
        return [_json_like(item) for item in value]
    return value


def _status(missing_inputs: list[dict[str, str]], calculated_rows: list[Any]) -> str:
    if missing_inputs:
        return "missing-input"
    if calculated_rows:
        return "partial"
    return "unsupported"
