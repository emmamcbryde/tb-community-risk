from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


METADATA_ALIASES = {
    "strategy": ("strategy", "strategyLabel", "strategy_label", "screeningStrategy"),
    "screeningTest": ("screeningTest", "screening_test", "test", "testType"),
    "regimen": ("regimen", "treatmentRegimen", "treatment_regimen"),
    "coverage": ("coverage", "screeningCoverage", "screening_coverage"),
    "population": ("population", "N", "n", "sampleSize", "sample_size"),
    "followHorizon": ("followHorizon", "follow_horizon", "timeHorizon", "time_horizon"),
    "seed": ("seed", "randomSeed", "random_seed"),
}


def extract_targeting_metadata(
    strategy: Any | None = None,
    config: Any | None = None,
) -> dict[str, Any]:
    """Extract stable targeting metadata from dict-like or attribute objects."""
    sources = (strategy, config)
    return {
        field: _json_scalar(_first_available(sources, aliases))
        for field, aliases in METADATA_ALIASES.items()
    }


def compare_summary_metrics(
    summary: Mapping[str, Any],
    metrics: Sequence[str],
    baseline: Mapping[str, Any] | None = None,
    *,
    include_relative: bool = False,
) -> list[dict[str, Any]]:
    """Build JSON-serialisable comparison rows for existing summary metrics."""
    if isinstance(metrics, str):
        raise TypeError("metrics must be a sequence of metric names, not a string")

    return [
        _comparison_row(summary, metric, baseline, include_relative=include_relative)
        for metric in metrics
    ]


def build_targeting_comparison(
    summary: Mapping[str, Any],
    metrics: Sequence[str],
    baseline: Mapping[str, Any] | None = None,
    *,
    strategy: Any | None = None,
    config: Any | None = None,
    include_relative: bool = False,
) -> dict[str, Any]:
    """Package metadata and summary comparisons without running a model."""
    return {
        "metadata": extract_targeting_metadata(strategy=strategy, config=config),
        "metricRows": compare_summary_metrics(
            summary,
            metrics,
            baseline=baseline,
            include_relative=include_relative,
        ),
    }


def _comparison_row(
    summary: Mapping[str, Any],
    metric: str,
    baseline: Mapping[str, Any] | None,
    *,
    include_relative: bool,
) -> dict[str, Any]:
    value_raw = _metric_value(summary, metric)
    value = _finite_number(value_raw)
    status = _metric_status(value_raw, value)

    baseline_raw = _metric_value(baseline, metric) if baseline is not None else None
    baseline_value = _finite_number(baseline_raw)
    baseline_status = (
        None if baseline is None else _metric_status(baseline_raw, baseline_value)
    )

    absolute = None
    relative = None
    if status == "ok" and baseline_status == "ok":
        absolute = value - baseline_value
        if include_relative and baseline_value != 0:
            relative = absolute / baseline_value

    return {
        "metric": str(metric),
        "status": status,
        "value": value if status == "ok" else None,
        "baselineStatus": baseline_status,
        "baseline": baseline_value if baseline_status == "ok" else None,
        "absoluteDifference": absolute,
        "relativeDifference": relative,
        "notes": _metric_notes(metric, status, baseline_status, include_relative),
    }


def _first_available(sources: tuple[Any | None, ...], aliases: Sequence[str]) -> Any:
    for source in sources:
        if source is None:
            continue
        for alias in aliases:
            value = _read_field(source, alias)
            if value is not _MISSING:
                return value
    return None


def _metric_value(summary: Mapping[str, Any] | None, metric: str) -> Any:
    if summary is None:
        return _MISSING
    value = _read_field(summary, metric)
    if value is _MISSING:
        return _MISSING
    if isinstance(value, Mapping):
        for field in ("Median", "median", "value", "Value"):
            nested = _read_field(value, field)
            if nested is not _MISSING:
                return nested
        return _MISSING
    return value


def _read_field(source: Any, field: str) -> Any:
    if isinstance(source, Mapping):
        return source[field] if field in source else _MISSING
    return getattr(source, field, _MISSING)


def _finite_number(value: Any) -> float | None:
    if isinstance(value, bool) or value is _MISSING:
        return None
    if isinstance(value, (int, float)):
        numeric = float(value)
        return numeric if math.isfinite(numeric) else None
    return None


def _metric_status(raw_value: Any, numeric_value: float | None) -> str:
    if raw_value is _MISSING:
        return "missing"
    if numeric_value is None:
        return "non_numeric"
    return "ok"


def _metric_notes(
    metric: str,
    status: str,
    baseline_status: str | None,
    include_relative: bool,
) -> str:
    if status == "missing":
        return f"Metric '{metric}' was not present in the supplied summary."
    if status == "non_numeric":
        return f"Metric '{metric}' was present but is not a finite numeric value."
    if baseline_status == "missing":
        return "Baseline metric was not present, so no difference was calculated."
    if baseline_status == "non_numeric":
        return "Baseline metric is not finite numeric, so no difference was calculated."
    if baseline_status == "ok" and not include_relative:
        return "Relative difference was not requested."
    return ""


def _json_scalar(value: Any) -> Any:
    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    return str(value)


class _Missing:
    pass


_MISSING = _Missing()
