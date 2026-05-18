from __future__ import annotations

import copy
from collections.abc import Mapping
from numbers import Real
from typing import Any


_NOT_PORTED_MESSAGE = (
    "Only the tested minimal APY economics subset is implemented; the full "
    "APY health-economics model is not ported to Python."
)

_CONTRACT_VERSION = "apy_economics_minimal_v1"
_SOURCE = "apy_python_minimal_economics"
_CALCULATED_COMPONENTS = [
    "testingCost",
    "treatmentCost",
]
_UNSUPPORTED_COMPONENTS = [
    "baselineTBDiseaseCost",
    "interventionTBDiseaseCost",
    "tbDiseaseCostsAverted",
    "netCostVsBaseline",
    "costPerInfectionCured",
    "costPerTBCasePrevented",
]
_OPTIONAL_DIRECT_COMPONENTS = [
    "falsePositiveIncrementalCost",
    "programSetupCost",
    "programRunningCost",
    "totalProgramCost",
]
_TOTAL_PROGRAM_COST_COMPONENTS = (
    "testingCost",
    "treatmentCost",
    "falsePositiveIncrementalCost",
    "programSetupCost",
    "programRunningCost",
)
_UNSUPPORTED_COMPONENT_MESSAGE = (
    "Known MATLAB APY v9 economics component; not calculated by the current "
    "minimal Python economics port."
)


def calculate_economics(
    result_bundle: Mapping[str, Any],
    economics_config: Mapping[str, Any],
) -> dict[str, Any]:
    """Return the supported minimal APY economics payload."""
    _validate_result_bundle(result_bundle)
    _validate_economics_config(economics_config)

    results = result_bundle["results"]
    interface_config = results["interfaceConfig"]
    summary = results["summary"]

    test_type = _require_non_empty_string(interface_config, "testType")
    regimen = _require_non_empty_string(interface_config, "regimen")

    n_screened = _summary_median(summary, "nScreened")
    n_started = _summary_median(summary, "nTotalCoursesStarted")
    n_false_positive_treated = _optional_summary_median(
        summary,
        "nFalsePositiveTreated",
    )

    costs_config = economics_config["costs"]
    test_cost = _require_unit_cost(
        costs_config["test"],
        test_type,
        f"costs.test.{test_type}",
    )
    regimen_key = f"x{regimen}"
    regimen_cost = _require_unit_cost(
        costs_config["regimen"],
        regimen_key,
        f"costs.regimen.{regimen_key}",
    )
    false_positive_unit_cost = _optional_cost(
        costs_config,
        "falsePositiveIncrementalPerPerson",
        "costs.falsePositiveIncrementalPerPerson",
    )
    program_setup_total = _optional_cost(
        costs_config,
        "programSetupTotal",
        "costs.programSetupTotal",
    )
    program_running_total = _optional_cost(
        costs_config,
        "programRunningTotal",
        "costs.programRunningTotal",
    )

    testing_cost = n_screened * test_cost
    treatment_cost = n_started * regimen_cost
    costs = {
        "testingCost": testing_cost,
        "treatmentCost": treatment_cost,
    }
    quantities = {
        "nScreened": n_screened,
        "nTotalCoursesStarted": n_started,
    }
    unit_costs = {
        "testPerPerson": test_cost,
        "treatmentPerStarted": regimen_cost,
    }
    coverage = _coverage_metadata()
    summary_rows = [
        {
            "metric": "testingCost",
            "label": "Testing cost",
            "value": testing_cost,
        },
        {
            "metric": "treatmentCost",
            "label": "Treatment cost",
            "value": treatment_cost,
        },
    ]

    if n_false_positive_treated is not None:
        quantities["nFalsePositiveTreated"] = n_false_positive_treated
    if false_positive_unit_cost is not None:
        unit_costs["falsePositiveIncrementalPerPerson"] = false_positive_unit_cost
    _add_optional_direct_total_cost(
        costs,
        coverage,
        summary_rows,
        component="programSetupCost",
        label="Program setup cost",
        value=program_setup_total,
        missing_input="costs.programSetupTotal",
    )
    _add_optional_direct_total_cost(
        costs,
        coverage,
        summary_rows,
        component="programRunningCost",
        label="Program running cost",
        value=program_running_total,
        missing_input="costs.programRunningTotal",
    )

    if n_false_positive_treated is None:
        _mark_missing(
            coverage,
            "falsePositiveIncrementalCost",
            "summary.nFalsePositiveTreated",
        )
    if false_positive_unit_cost is None:
        _mark_missing(
            coverage,
            "falsePositiveIncrementalCost",
            "costs.falsePositiveIncrementalPerPerson",
        )
    if n_false_positive_treated is not None and false_positive_unit_cost is not None:
        false_positive_incremental_cost = (
            n_false_positive_treated * false_positive_unit_cost
        )
        costs["falsePositiveIncrementalCost"] = false_positive_incremental_cost
        _mark_calculated(coverage, "falsePositiveIncrementalCost")
        summary_rows.append(
            {
                "metric": "falsePositiveIncrementalCost",
                "label": "False-positive incremental cost",
                "value": false_positive_incremental_cost,
            }
        )

    if all(component in costs for component in _TOTAL_PROGRAM_COST_COMPONENTS):
        total_program_cost = (
            costs["testingCost"]
            + costs["treatmentCost"]
            + costs["falsePositiveIncrementalCost"]
            + costs["programSetupCost"]
            + costs["programRunningCost"]
        )
        costs["totalProgramCost"] = total_program_cost
        _mark_calculated(coverage, "totalProgramCost")
        summary_rows.append(
            {
                "metric": "totalProgramCost",
                "label": "Total program cost",
                "value": total_program_cost,
            }
        )
    else:
        _mark_not_calculated(coverage, "totalProgramCost")

    return {
        "available": True,
        "source": _SOURCE,
        "contractVersion": _CONTRACT_VERSION,
        "message": _NOT_PORTED_MESSAGE,
        "metadata": _json_like_copy(economics_config["metadata"]),
        "inputs": _json_like_copy(economics_config),
        "strategy": {
            "testType": test_type,
            "regimen": regimen,
        },
        "quantities": quantities,
        "unitCosts": unit_costs,
        "costs": costs,
        "costEffectiveness": {},
        "summaryRows": summary_rows,
        "coverage": coverage,
        "status": "partial",
    }


def _add_optional_direct_total_cost(
    costs: dict[str, float],
    coverage: dict[str, Any],
    summary_rows: list[dict[str, Any]],
    *,
    component: str,
    label: str,
    value: float | None,
    missing_input: str,
) -> None:
    if value is None:
        _mark_missing(coverage, component, missing_input)
        return

    costs[component] = value
    _mark_calculated(coverage, component)
    summary_rows.append(
        {
            "metric": component,
            "label": label,
            "value": value,
        }
    )


def _coverage_metadata() -> dict[str, Any]:
    unsupported = [
        {
            "component": component,
            "message": _UNSUPPORTED_COMPONENT_MESSAGE,
        }
        for component in _UNSUPPORTED_COMPONENTS
    ]
    return {
        "status": "partial",
        "calculatedComponents": list(_CALCULATED_COMPONENTS),
        "unsupportedComponents": unsupported,
        "notCalculated": (
            list(_OPTIONAL_DIRECT_COMPONENTS) + list(_UNSUPPORTED_COMPONENTS)
        ),
        "missingInputs": [],
        "messages": [
            _NOT_PORTED_MESSAGE,
            _UNSUPPORTED_COMPONENT_MESSAGE,
        ],
    }


def _mark_calculated(coverage: dict[str, Any], component: str) -> None:
    if component not in coverage["calculatedComponents"]:
        coverage["calculatedComponents"].append(component)
    if component in coverage["notCalculated"]:
        coverage["notCalculated"].remove(component)


def _mark_not_calculated(coverage: dict[str, Any], component: str) -> None:
    if component not in coverage["notCalculated"]:
        coverage["notCalculated"].append(component)


def _mark_missing(
    coverage: dict[str, Any],
    component: str,
    input_name: str,
) -> None:
    if input_name not in coverage["missingInputs"]:
        coverage["missingInputs"].append(input_name)
    _mark_not_calculated(coverage, component)


def _validate_result_bundle(result_bundle: Mapping[str, Any]) -> None:
    if not isinstance(result_bundle, Mapping):
        raise TypeError("result_bundle must be a mapping.")

    missing = [
        key for key in ("metadata", "results")
        if key not in result_bundle
    ]
    if missing:
        raise ValueError(
            "result_bundle is missing required top-level key(s): "
            + ", ".join(missing)
            + "."
        )

    for key in ("metadata", "results"):
        if not isinstance(result_bundle[key], Mapping):
            raise ValueError(f"result_bundle.{key} must be a mapping.")

    results = result_bundle["results"]
    missing_results = [
        key for key in ("interfaceConfig", "summary")
        if key not in results
    ]
    if missing_results:
        raise ValueError(
            "result_bundle.results is missing required key(s): "
            + ", ".join(missing_results)
            + "."
        )

    if not isinstance(results["interfaceConfig"], Mapping):
        raise ValueError("result_bundle.results.interfaceConfig must be a mapping.")
    if not isinstance(results["summary"], list):
        raise ValueError("result_bundle.results.summary must be a list.")


def _validate_economics_config(economics_config: Mapping[str, Any]) -> None:
    if not isinstance(economics_config, Mapping):
        raise TypeError("economics_config must be a mapping.")

    missing = [
        key for key in ("metadata", "costs")
        if key not in economics_config
    ]
    if missing:
        raise ValueError(
            "economics_config is missing required top-level key(s): "
            + ", ".join(missing)
            + "."
        )

    for key in ("metadata", "costs"):
        if not isinstance(economics_config[key], Mapping):
            raise ValueError(f"economics_config.{key} must be a mapping.")

    costs = economics_config["costs"]
    missing_costs = [
        key for key in ("test", "regimen")
        if key not in costs
    ]
    if missing_costs:
        raise ValueError(
            "economics_config.costs is missing required key(s): "
            + ", ".join(missing_costs)
            + "."
        )

    for key in ("test", "regimen"):
        if not isinstance(costs[key], Mapping):
            raise ValueError(f"economics_config.costs.{key} must be a mapping.")


def _require_non_empty_string(mapping: Mapping[str, Any], key: str) -> str:
    value = mapping.get(key)
    if not isinstance(value, str) or not value:
        raise ValueError(f"result_bundle.results.interfaceConfig.{key} is required.")
    return value


def _summary_median(summary: list[Any], metric: str) -> float:
    for row in summary:
        if not isinstance(row, Mapping):
            raise ValueError("result_bundle.results.summary rows must be mappings.")
        if row.get("Metric") == metric:
            if "Median" not in row:
                raise ValueError(f"summary metric {metric} is missing Median.")
            return _as_float(row["Median"], f"summary metric {metric} Median")

    raise ValueError(f"summary is missing required metric: {metric}.")


def _optional_summary_median(summary: list[Any], metric: str) -> float | None:
    for row in summary:
        if not isinstance(row, Mapping):
            raise ValueError("result_bundle.results.summary rows must be mappings.")
        if row.get("Metric") == metric:
            if "Median" not in row or row["Median"] is None:
                return None
            return _as_float(row["Median"], f"summary metric {metric} Median")

    return None


def _require_unit_cost(
    costs: Mapping[str, Any],
    key: str,
    path: str,
) -> float:
    if key not in costs or costs[key] is None:
        raise ValueError(f"{path} is required for minimal APY economics.")
    return _as_float(costs[key], path)


def _optional_cost(
    costs: Mapping[str, Any],
    key: str,
    path: str,
) -> float | None:
    if key not in costs or costs[key] is None:
        return None
    return _as_float(costs[key], path)


def _as_float(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{label} must be numeric.")
    return float(value)


def _json_like_copy(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_like_copy(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_like_copy(item) for item in value]
    if isinstance(value, tuple):
        return [_json_like_copy(item) for item in value]
    return copy.deepcopy(value)
