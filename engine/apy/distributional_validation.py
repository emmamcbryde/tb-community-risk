from __future__ import annotations

from pathlib import Path
from typing import Any

import math
import pandas as pd

from engine.apy.config import build_default_config
from engine.apy.parity import DYNAMIC_COMPARISON_METRICS
from engine.apy.reference_loader import (
    list_reference_scenario_dirs,
    load_reference_dir,
    load_reference_scenario_suite,
)
from engine.apy.runner import run_scenario_with_do_nothing


DEFAULT_SUMMARY_METRICS = [
    "nScreened",
    "nInfected",
    "nLatentAtScreen",
    "nActiveAtScreen",
    "nTestPositive",
    "nTestPositiveNonActive",
    "nFalsePositiveTests",
    "nFalsePositiveTreated",
    "nTotalCoursesStarted",
    "nTotalCoursesCompleted",
    "nCuredInfection",
    "nPreventedActiveTB",
    "nActiveBy2y",
    "nActiveBy20y",
    "NNS_cureInfection",
    "NNS_preventActiveTB",
    "NNT_started_cureInfection",
    "NNT_started_preventActiveTB",
]


PORTABLE_REFERENCE_FIELDS = [
    "scenarioLabel",
    "N",
    "nReps",
    "seed",
    "screenWindow",
    "followHorizon",
    "screenCoverage",
    "screeningStrategy",
    "targetAgeOR",
    "testType",
    "regimen",
    "pStartTPT",
    "regimenPComplete",
    "regimenADRstop",
    "regimenEffFull",
    "partialShortCourseMode",
    "partialDoseFractionADR",
    "partialDoseFractionOther",
    "earlyLateRatio",
    "riskPrev",
    "diseaseOR",
]


def run_reference_distributional_validation(
    reference_dir: str | Path,
    config_overrides: dict[str, Any] | None = None,
    summary_metrics: list[str] | None = None,
    tolerance: dict[str, float] | None = None,
) -> dict[str, Any]:
    reference = load_reference_dir(reference_dir)
    config = portable_config_from_reference(reference["scenario_config"])
    if config_overrides:
        config.update(config_overrides)

    python_output = run_scenario_with_do_nothing(config)
    validation = build_distributional_validation_table(
        python_summary=python_output["results"]["summary"],
        matlab_summary=reference["summary"],
        python_dynamic_comparison=python_output["bundle"]["technical"][
            "dynamicComparison"
        ],
        matlab_dynamic_comparison=reference["dynamic_comparison"],
        summary_metrics=summary_metrics,
        tolerance=tolerance,
    )
    return {
        "config": config,
        "pythonOutput": python_output,
        "reference": reference,
        "validation": validation,
    }


def run_reference_suite_distributional_validation(
    reference_root: str | Path,
    suite_file: str | Path | None = None,
    scenario_ids: list[str] | None = None,
    config_overrides: dict[str, Any] | None = None,
    tolerance: dict[str, float] | None = None,
) -> dict[str, Any]:
    root = Path(reference_root)
    scenarios = _suite_scenarios(root, suite_file)
    requested_ids = _normalise_scenario_ids(scenario_ids)
    if requested_ids is not None:
        scenarios = [
            scenario for scenario in scenarios
            if scenario["scenario_id"] in requested_ids
        ]

    scenario_rows = []
    metric_frames = []
    scenario_results = {}
    for scenario in scenarios:
        scenario_id = str(scenario["scenario_id"])
        scenario_dir = root / scenario_id
        if not scenario_dir.is_dir():
            scenario_rows.append(
                {
                    "scenario_id": scenario_id,
                    "n_metrics": 0,
                    "n_pass": 0,
                    "n_fail": 0,
                    "pass_rate": math.nan,
                    "max_abs_relative_difference": math.nan,
                    "notes": f"Reference fixture directory missing: {scenario_dir}",
                }
            )
            continue

        validation_output = run_reference_distributional_validation(
            scenario_dir,
            config_overrides=config_overrides,
            tolerance=tolerance,
        )
        validation = validation_output["validation"].copy()
        validation.insert(0, "scenario_id", scenario_id)
        metric_frames.append(validation)
        scenario_results[scenario_id] = validation_output

        comparable_rel = validation["RelativeDifference"].abs().dropna()
        n_metrics = int(len(validation))
        n_pass = int(validation["Pass"].sum())
        scenario_rows.append(
            {
                "scenario_id": scenario_id,
                "n_metrics": n_metrics,
                "n_pass": n_pass,
                "n_fail": n_metrics - n_pass,
                "pass_rate": math.nan if n_metrics == 0 else n_pass / n_metrics,
                "max_abs_relative_difference": (
                    math.nan if comparable_rel.empty else float(comparable_rel.max())
                ),
                "notes": "Fixture validated against Python diagnostic tolerances.",
            }
        )

    metric_rows = (
        pd.concat(metric_frames, ignore_index=True)
        if metric_frames
        else pd.DataFrame(
            columns=[
                "scenario_id",
                "Metric",
                "PythonMedian",
                "MatlabMedian",
                "AbsoluteDifference",
                "RelativeDifference",
                "Tolerance",
                "Pass",
                "Notes",
            ]
        )
    )
    return {
        "scenarioRows": pd.DataFrame(scenario_rows),
        "metricRows": metric_rows,
        "scenarioResults": scenario_results,
    }


def portable_config_from_reference(reference_config: dict[str, Any]) -> dict[str, Any]:
    config = build_default_config()
    for field in PORTABLE_REFERENCE_FIELDS:
        if field in reference_config:
            config[field] = _matlab_empty_to_none(reference_config[field])
    return config


def build_distributional_validation_table(
    python_summary: pd.DataFrame,
    matlab_summary: pd.DataFrame,
    python_dynamic_comparison: dict[str, Any],
    matlab_dynamic_comparison: dict[str, Any],
    summary_metrics: list[str] | None = None,
    tolerance: dict[str, float] | None = None,
) -> pd.DataFrame:
    tol = {
        "absolute": 10.0,
        "relative": 0.50,
        "proportion_absolute": 0.15,
        "proportion_relative": 0.50,
    }
    if tolerance:
        tol.update(tolerance)

    rows = []
    summary_metrics = summary_metrics or DEFAULT_SUMMARY_METRICS
    python_summary_rows = _summary_rows_by_metric(python_summary)
    matlab_summary_rows = _summary_rows_by_metric(matlab_summary)
    for metric in summary_metrics:
        rows.append(
            _validation_row(
                metric=metric,
                python_value=_summary_median(python_summary_rows.get(metric)),
                matlab_value=_summary_median(matlab_summary_rows.get(metric)),
                tolerance=tol,
                notes=(
                    "Summary median comparison; strict equality is not expected "
                    "because MATLAB and NumPy random streams differ."
                ),
            )
        )

    for metric in DYNAMIC_COMPARISON_METRICS:
        rows.append(
            _validation_row(
                metric=metric,
                python_value=_dynamic_value(python_dynamic_comparison, metric),
                matlab_value=_dynamic_value(matlab_dynamic_comparison, metric),
                tolerance=tol,
                notes=(
                    "technical.dynamicComparison median comparison from "
                    "do-nothing derived rows."
                ),
            )
        )
    return pd.DataFrame(rows)


def _validation_row(
    metric: str,
    python_value,
    matlab_value,
    tolerance: dict[str, float],
    notes: str,
) -> dict[str, Any]:
    comparable = python_value is not None and matlab_value is not None
    if comparable:
        py = float(python_value)
        matlab = float(matlab_value)
        if math.isinf(py) or math.isinf(matlab):
            same_infinity = math.isinf(py) and math.isinf(matlab) and py == matlab
            abs_diff = 0.0 if same_infinity else math.inf
            rel_diff = 0.0 if same_infinity else math.inf
            allowed = 0.0
            passed = same_infinity
            tolerance_label = "same infinity"
        else:
            abs_diff = py - matlab
            rel_diff = math.nan if matlab == 0 else abs_diff / matlab
            allowed = _allowed_difference(metric, matlab, tolerance)
            passed = abs(abs_diff) <= allowed
            tolerance_label = f"abs <= {allowed:.6g}"
    else:
        abs_diff = math.nan
        rel_diff = math.nan
        passed = False
        tolerance_label = "not comparable"
        notes = "Metric missing from Python or MATLAB reference."

    return {
        "Metric": metric,
        "PythonMedian": python_value,
        "MatlabMedian": matlab_value,
        "AbsoluteDifference": abs_diff,
        "RelativeDifference": rel_diff,
        "Tolerance": tolerance_label,
        "Pass": bool(passed),
        "Notes": notes,
    }


def _allowed_difference(
    metric: str,
    matlab_value: float,
    tolerance: dict[str, float],
) -> float:
    if _is_proportion_metric(metric):
        return max(
            tolerance["proportion_absolute"],
            abs(matlab_value) * tolerance["proportion_relative"],
        )
    return max(tolerance["absolute"], abs(matlab_value) * tolerance["relative"])


def _is_proportion_metric(metric: str) -> bool:
    metric_lower = metric.lower()
    return (
        "prev" in metric_lower
        or "reduction" in metric_lower
        or metric_lower.startswith("rel")
    )


def _summary_rows_by_metric(summary: pd.DataFrame) -> dict[str, dict[str, Any]]:
    if summary is None or summary.empty or "Metric" not in summary.columns:
        return {}
    return {
        str(row["Metric"]): row
        for row in summary.to_dict(orient="records")
    }


def _summary_median(row: dict[str, Any] | None):
    if row is None or "Median" not in row:
        return None
    value = row["Median"]
    if pd.isna(value):
        return math.nan
    return float(value)


def _dynamic_value(dynamic_comparison: dict[str, Any], metric: str):
    if not dynamic_comparison:
        return None
    direct = dynamic_comparison.get(metric)
    if direct is not None:
        return float(direct)
    for row in dynamic_comparison.get("metricRows", []):
        if row.get("Metric") == metric and row.get("Median") is not None:
            return float(row["Median"])
    return None


def _matlab_empty_to_none(value):
    if value == []:
        return None
    if isinstance(value, dict):
        return {
            key: _matlab_empty_to_none(item)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [_matlab_empty_to_none(item) for item in value]
    return value


def _suite_scenarios(root: Path, suite_file: str | Path | None) -> list[dict[str, Any]]:
    if suite_file is not None:
        return load_reference_scenario_suite(suite_file)

    return [
        {
            "scenario_id": path.name,
            "description": "Discovered committed MATLAB reference fixture.",
            "config_overrides": {},
            "expected_focus": [],
            "notes": "Discovered from reference directory.",
        }
        for path in list_reference_scenario_dirs(root)
    ]


def _normalise_scenario_ids(scenario_ids: list[str] | None) -> set[str] | None:
    if not scenario_ids:
        return None
    out = set()
    for item in scenario_ids:
        for part in str(item).split(","):
            part = part.strip()
            if part:
                out.add(part)
    return out
