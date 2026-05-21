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

    return {
        "status": "missing-input" if missing_inputs else "unsupported",
        "source": _SOURCE,
        "calculatedRows": [],
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
