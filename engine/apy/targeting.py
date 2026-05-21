from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from copy import deepcopy
from typing import Any


METADATA_ALIASES = {
    "strategy": (
        "name",
        "strategyName",
        "strategy_name",
        "strategy",
        "strategyLabel",
        "strategy_label",
        "screeningStrategy",
    ),
    "screeningTest": ("screeningTest", "screening_test", "test", "testType"),
    "regimen": ("regimen", "treatmentRegimen", "treatment_regimen"),
    "coverage": (
        "coverage",
        "screeningCoverage",
        "screening_coverage",
        "screenCoverage",
    ),
    "population": ("population", "N", "n", "sampleSize", "sample_size"),
    "followHorizon": ("followHorizon", "follow_horizon", "timeHorizon", "time_horizon"),
    "seed": ("seed", "randomSeed", "random_seed"),
}


def make_strategy_spec(
    name: str,
    config_overrides: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a minimal JSON-compatible APY strategy spec."""
    if not isinstance(name, str) or not name.strip():
        raise ValueError("strategy spec name must be a non-empty string")
    overrides = {} if config_overrides is None else dict(config_overrides)
    spec = {
        "name": name,
        "config_overrides": deepcopy(overrides),
    }
    _assert_json_compatible(spec)
    return spec


def apply_strategy_config_overrides(
    base_config: Mapping[str, Any],
    strategy_spec: Mapping[str, Any],
) -> dict[str, Any]:
    """Apply explicit strategy config overrides without mutating the base config."""
    if not isinstance(base_config, Mapping):
        raise TypeError("base_config must be a mapping")
    if not isinstance(strategy_spec, Mapping):
        raise TypeError("strategy_spec must be a mapping")
    if "config_overrides" not in strategy_spec:
        raise ValueError("strategy_spec must include explicit config_overrides")
    overrides = strategy_spec["config_overrides"]
    if not isinstance(overrides, Mapping):
        raise TypeError("strategy_spec config_overrides must be a mapping")

    config = deepcopy(dict(base_config))
    config.update(deepcopy(dict(overrides)))
    if "strategySpec" not in overrides and "strategy_spec" not in overrides:
        config["strategySpec"] = deepcopy(dict(strategy_spec))
    return config


def extract_targeting_metadata(
    strategy: Any | None = None,
    config: Any | None = None,
    result_bundle: Any | None = None,
) -> dict[str, Any]:
    """Extract stable targeting metadata from dict-like or attribute objects."""
    sources = _metadata_sources(
        strategy=strategy,
        config=config,
        result_bundle=result_bundle,
    )
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


def compare_targeting_result_bundles(
    baseline_result_bundle: Mapping[str, Any],
    comparator_result_bundles: Sequence[Mapping[str, Any]],
    metrics: Sequence[str],
) -> list[dict[str, Any]]:
    """Compare already-run result bundles using only existing summary sections."""
    if isinstance(comparator_result_bundles, (str, bytes, bytearray)):
        raise TypeError("comparator_result_bundles must be a sequence of bundles")

    baseline_summary = _bundle_summary_metrics(baseline_result_bundle)
    baseline_metadata = extract_targeting_metadata(
        result_bundle=baseline_result_bundle,
    )
    baseline_strategy = _strategy_name(baseline_metadata)

    rows: list[dict[str, Any]] = []
    for comparator_bundle in comparator_result_bundles:
        comparator_summary = _bundle_summary_metrics(comparator_bundle)
        comparator_metadata = extract_targeting_metadata(
            result_bundle=comparator_bundle,
        )
        comparator_strategy = _strategy_name(comparator_metadata)
        comparison_rows = compare_summary_metrics(
            comparator_summary,
            metrics,
            baseline=baseline_summary,
            include_relative=True,
        )
        for row in comparison_rows:
            rows.append(
                {
                    "baselineStrategy": baseline_strategy,
                    "comparatorStrategy": comparator_strategy,
                    "metric": row["metric"],
                    "status": row["status"],
                    "baselineStatus": row["baselineStatus"],
                    "comparatorStatus": row["status"],
                    "baselineValue": row["baseline"],
                    "comparatorValue": row["value"],
                    "absoluteDifference": row["absoluteDifference"],
                    "relativeDifference": row["relativeDifference"],
                    "notes": row["notes"],
                }
            )
    _assert_json_compatible(rows)
    return rows


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


def _bundle_summary_metrics(bundle: Any) -> dict[str, Any]:
    for source, field in (
        (bundle, "summary"),
        (bundle, "summaryRows"),
        (bundle, "keyMetricsRows"),
        (_read_optional(bundle, "results"), "summary"),
        (_read_optional(bundle, "headline"), "summaryRows"),
        (_read_optional(bundle, "headline"), "keyMetricsRows"),
    ):
        summary = _summary_mapping(_read_optional(source, field))
        if summary:
            return summary
    return {}


def _summary_mapping(summary: Any) -> dict[str, Any]:
    if summary is None:
        return {}

    rows = _summary_rows(summary)
    if rows is not None:
        mapped: dict[str, Any] = {}
        for row in rows:
            if not isinstance(row, Mapping):
                continue
            metric = _row_metric_name(row)
            if metric is not None:
                mapped[str(metric)] = row
        return mapped

    if isinstance(summary, Mapping):
        return dict(summary)

    return {}


def _summary_rows(summary: Any) -> list[Any] | None:
    if isinstance(summary, Sequence) and not isinstance(summary, (str, bytes, bytearray)):
        return list(summary)

    to_dict = getattr(summary, "to_dict", None)
    if callable(to_dict):
        try:
            records = to_dict(orient="records")
        except TypeError:
            return None
        if isinstance(records, list):
            return records

    return None


def _row_metric_name(row: Mapping[str, Any]) -> Any:
    for field in ("Metric", "metric", "name", "Name"):
        value = _read_field(row, field)
        if value is not _MISSING:
            return value
    return None


def _strategy_name(metadata: Mapping[str, Any]) -> Any:
    return _json_scalar(metadata.get("strategy"))


def _first_available(sources: tuple[Any | None, ...], aliases: Sequence[str]) -> Any:
    for source in sources:
        if source is None:
            continue
        for alias in aliases:
            value = _read_field(source, alias)
            if value is not _MISSING:
                return value
    return None


def _metadata_sources(
    *,
    strategy: Any | None,
    config: Any | None,
    result_bundle: Any | None,
) -> tuple[Any | None, ...]:
    bundle_sources = _bundle_metadata_sources(result_bundle)
    strategy_spec_sources = (
        _read_optional(strategy, "strategySpec"),
        _read_optional(strategy, "strategy_spec"),
        _read_optional(config, "strategySpec"),
        _read_optional(config, "strategy_spec"),
    )
    return (
        strategy,
        *strategy_spec_sources,
        config,
        *bundle_sources,
    )


def _bundle_metadata_sources(bundle: Any | None) -> tuple[Any | None, ...]:
    if bundle is None:
        return ()

    headline = _read_optional(bundle, "headline")
    technical = _read_optional(bundle, "technical")
    results = _read_optional(bundle, "results")
    result_strategy = _read_optional(results, "strategy")
    result_config = _read_optional(results, "interfaceConfig")
    return (
        _read_optional(bundle, "strategy"),
        _read_optional(bundle, "interfaceConfig"),
        _read_optional(bundle, "config"),
        _read_optional(headline, "strategy"),
        _read_optional(technical, "interfaceConfig"),
        result_strategy,
        result_config,
    )


def _read_optional(source: Any, field: str) -> Any | None:
    if source is None:
        return None
    value = _read_field(source, field)
    return None if value is _MISSING else value


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


def _assert_json_compatible(value: Any) -> None:
    if value is None or isinstance(value, (str, bool)):
        return
    if isinstance(value, int) and not isinstance(value, bool):
        return
    if isinstance(value, float):
        if math.isfinite(value):
            return
        raise ValueError("strategy spec must contain only finite numeric values")
    if isinstance(value, Mapping):
        for key, nested in value.items():
            if not isinstance(key, str):
                raise TypeError("strategy spec keys must be strings")
            _assert_json_compatible(nested)
        return
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        for nested in value:
            _assert_json_compatible(nested)
        return
    raise TypeError(
        "strategy spec contains non-JSON-compatible value: "
        f"{type(value).__name__}"
    )


class _Missing:
    pass


_MISSING = _Missing()
