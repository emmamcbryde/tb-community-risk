from __future__ import annotations

import copy
from collections.abc import Mapping
from numbers import Real
from typing import Any


_NOT_PORTED_MESSAGE = (
    "Only the tested partial APY economics subset is implemented; the full "
    "APY health-economics model is not ported to Python."
)

_CONTRACT_VERSION = "apy_economics_partial_v1"
_SOURCE = "apy_python_partial_economics"
_LEGACY_CONTRACT_VERSION = "apy_economics_minimal_v1"
_LEGACY_SOURCE = "apy_python_minimal_economics"
_CALCULATED_COMPONENTS = [
    "testingCost",
    "treatmentCost",
]
_UNSUPPORTED_COMPONENTS = [
]
_OPTIONAL_DIRECT_COMPONENTS = [
    "falsePositiveIncrementalCost",
    "programSetupCost",
    "programRunningCost",
    "totalProgramCost",
]
_OPTIONAL_BENEFIT_COMPONENTS = [
    "baselineTBDiseaseCost",
    "interventionTBDiseaseCost",
    "tbDiseaseCostsAverted",
    "netCostVsBaseline",
    "costPerInfectionCured",
    "costPerTBCasePrevented",
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
    "partial Python economics port."
)


def calculate_economics(
    result_bundle: Mapping[str, Any],
    economics_config: Mapping[str, Any],
) -> dict[str, Any]:
    """Return the supported partial APY economics payload."""
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
    n_prevented_active_tb = _optional_summary_median(
        summary,
        "nPreventedActiveTB",
    )
    n_cured_infection = _optional_summary_median(
        summary,
        "nCuredInfection",
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
    active_tb_disease_per_case = _optional_cost(
        costs_config,
        "activeTBDiseasePerCase",
        "costs.activeTBDiseasePerCase",
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
        _implemented_direct_cost_row(
            "testingCost",
            "Testing cost",
            testing_cost,
        ),
        _implemented_direct_cost_row(
            "treatmentCost",
            "Treatment cost",
            treatment_cost,
        ),
    ]

    if n_false_positive_treated is not None:
        quantities["nFalsePositiveTreated"] = n_false_positive_treated
    if false_positive_unit_cost is not None:
        unit_costs["falsePositiveIncrementalPerPerson"] = false_positive_unit_cost
    if active_tb_disease_per_case is not None:
        unit_costs["activeTBDiseasePerCase"] = active_tb_disease_per_case
    if n_prevented_active_tb is not None:
        quantities["nPreventedActiveTB"] = n_prevented_active_tb
    if n_cured_infection is not None:
        quantities["nCuredInfection"] = n_cured_infection
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
            _implemented_direct_cost_row(
                "falsePositiveIncrementalCost",
                "False-positive incremental cost",
                false_positive_incremental_cost,
            )
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
            _implemented_direct_cost_row(
                "totalProgramCost",
                "Total program cost",
                total_program_cost,
                category="directCostTotal",
                included_in_total=False,
            )
        )
    else:
        _mark_not_calculated(coverage, "totalProgramCost")

    dynamic_comparison = _optional_dynamic_comparison(results)
    if dynamic_comparison is None:
        _add_coverage_message(
            coverage,
            "dynamicComparison is missing; baseline/intervention TB disease "
            "costs were not calculated.",
        )
        _mark_missing(
            coverage,
            "baselineTBDiseaseCost",
            "results.dynamicComparison",
        )
        _mark_missing(
            coverage,
            "interventionTBDiseaseCost",
            "results.dynamicComparison",
        )
    elif active_tb_disease_per_case is not None:
        baseline_cases = _optional_numeric_value(
            dynamic_comparison,
            "cumulative_baseline_active_tb_cases",
            "results.dynamicComparison.cumulative_baseline_active_tb_cases",
        )
        intervention_cases = _optional_numeric_value(
            dynamic_comparison,
            "cumulative_intervention_active_tb_cases",
            "results.dynamicComparison.cumulative_intervention_active_tb_cases",
        )
        if baseline_cases is not None:
            quantities["cumulativeBaselineActiveTBCases"] = baseline_cases
            baseline_tb_disease_cost = baseline_cases * active_tb_disease_per_case
            costs["baselineTBDiseaseCost"] = baseline_tb_disease_cost
            _mark_calculated(coverage, "baselineTBDiseaseCost")
            summary_rows.append(
                _implemented_direct_cost_row(
                    "baselineTBDiseaseCost",
                    "Baseline TB disease cost",
                    baseline_tb_disease_cost,
                    category="benefitCost",
                    included_in_total=False,
                )
            )
        else:
            _mark_missing(
                coverage,
                "baselineTBDiseaseCost",
                "results.dynamicComparison.cumulative_baseline_active_tb_cases",
            )
        if intervention_cases is not None:
            quantities["cumulativeInterventionActiveTBCases"] = intervention_cases
            intervention_tb_disease_cost = (
                intervention_cases * active_tb_disease_per_case
            )
            costs["interventionTBDiseaseCost"] = intervention_tb_disease_cost
            _mark_calculated(coverage, "interventionTBDiseaseCost")
            summary_rows.append(
                _implemented_direct_cost_row(
                    "interventionTBDiseaseCost",
                    "Intervention TB disease cost",
                    intervention_tb_disease_cost,
                    category="benefitCost",
                    included_in_total=False,
                )
            )
        else:
            _mark_missing(
                coverage,
                "interventionTBDiseaseCost",
                "results.dynamicComparison.cumulative_intervention_active_tb_cases",
            )
    else:
        _mark_missing(
            coverage,
            "baselineTBDiseaseCost",
            "costs.activeTBDiseasePerCase",
        )
        _mark_missing(
            coverage,
            "interventionTBDiseaseCost",
            "costs.activeTBDiseasePerCase",
        )

    if active_tb_disease_per_case is None:
        _mark_missing(
            coverage,
            "tbDiseaseCostsAverted",
            "costs.activeTBDiseasePerCase",
        )
    if n_prevented_active_tb is None:
        _mark_missing(
            coverage,
            "tbDiseaseCostsAverted",
            "summary.nPreventedActiveTB",
        )
    if active_tb_disease_per_case is not None and n_prevented_active_tb is not None:
        tb_disease_costs_averted = (
            n_prevented_active_tb * active_tb_disease_per_case
        )
        costs["tbDiseaseCostsAverted"] = tb_disease_costs_averted
        _mark_calculated(coverage, "tbDiseaseCostsAverted")
        summary_rows.append(
            _implemented_direct_cost_row(
                "tbDiseaseCostsAverted",
                "TB disease costs averted",
                tb_disease_costs_averted,
                category="benefitCost",
                included_in_total=False,
            )
        )

    if "totalProgramCost" in costs and "tbDiseaseCostsAverted" in costs:
        net_cost_vs_baseline = (
            costs["totalProgramCost"] - costs["tbDiseaseCostsAverted"]
        )
        costs["netCostVsBaseline"] = net_cost_vs_baseline
        _mark_calculated(coverage, "netCostVsBaseline")
        summary_rows.append(
            _implemented_direct_cost_row(
                "netCostVsBaseline",
                "Net cost vs baseline",
                net_cost_vs_baseline,
                category="netCost",
                included_in_total=False,
            )
        )
    else:
        _mark_not_calculated(coverage, "netCostVsBaseline")

    cost_effectiveness: dict[str, float] = {}
    _add_simple_ratio(
        cost_effectiveness,
        coverage,
        summary_rows,
        metric="costPerInfectionCured",
        label="Cost per infection cured",
        net_cost=costs.get("netCostVsBaseline"),
        denominator=n_cured_infection,
        denominator_input="summary.nCuredInfection",
    )
    _add_simple_ratio(
        cost_effectiveness,
        coverage,
        summary_rows,
        metric="costPerTBCasePrevented",
        label="Cost per TB case prevented",
        net_cost=costs.get("netCostVsBaseline"),
        denominator=n_prevented_active_tb,
        denominator_input="summary.nPreventedActiveTB",
    )

    summary_rows.append(_total_implemented_cost_row(summary_rows))

    return {
        "available": True,
        "source": _SOURCE,
        "contractVersion": _CONTRACT_VERSION,
        "message": _NOT_PORTED_MESSAGE,
        "metadata": _economics_metadata(economics_config["metadata"]),
        "inputs": _json_like_copy(economics_config),
        "strategy": {
            "testType": test_type,
            "regimen": regimen,
        },
        "quantities": quantities,
        "unitCosts": unit_costs,
        "costs": costs,
        "costEffectiveness": cost_effectiveness,
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
        _implemented_direct_cost_row(component, label, value)
    )


def _add_simple_ratio(
    cost_effectiveness: dict[str, float],
    coverage: dict[str, Any],
    summary_rows: list[dict[str, Any]],
    *,
    metric: str,
    label: str,
    net_cost: float | None,
    denominator: float | None,
    denominator_input: str,
) -> None:
    if net_cost is None:
        _add_coverage_message(
            coverage,
            f"{metric} was not calculated because netCostVsBaseline is unavailable.",
        )
        _mark_not_calculated(coverage, metric)
        return
    if denominator is None:
        _add_coverage_message(
            coverage,
            f"{metric} was not calculated because {denominator_input} is missing.",
        )
        _mark_missing(coverage, metric, denominator_input)
        return
    if denominator <= 0:
        _add_coverage_message(
            coverage,
            f"{metric} was not calculated because {denominator_input} is not positive.",
        )
        _mark_not_calculated(coverage, metric)
        return

    value = net_cost / denominator
    cost_effectiveness[metric] = value
    _mark_calculated(coverage, metric)
    summary_rows.append(
        _implemented_direct_cost_row(
            metric,
            label,
            value,
            category="costEffectiveness",
            included_in_total=False,
        )
    )


def _implemented_direct_cost_row(
    metric: str,
    label: str,
    value: float,
    *,
    category: str = "directCost",
    included_in_total: bool = True,
) -> dict[str, Any]:
    return {
        "metric": metric,
        "label": label,
        "value": value,
        "category": category,
        "status": "implemented",
        "includedInTotal": included_in_total,
    }


def _total_implemented_cost_row(summary_rows: list[dict[str, Any]]) -> dict[str, Any]:
    total = sum(
        row["value"]
        for row in summary_rows
        if row.get("includedInTotal") is True
    )
    return {
        "metric": "totalImplementedCost",
        "label": "Total implemented cost",
        "value": total,
        "category": "directCostTotal",
        "status": "implemented",
        "includedInTotal": False,
    }


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
            + list(_OPTIONAL_BENEFIT_COMPONENTS)
        ),
        "missingInputs": [],
        "messages": [
            _NOT_PORTED_MESSAGE,
            _UNSUPPORTED_COMPONENT_MESSAGE,
        ],
    }


def _economics_metadata(metadata: Mapping[str, Any]) -> dict[str, Any]:
    payload_metadata = _json_like_copy(metadata)
    payload_metadata["legacyIdentifiers"] = {
        "source": _LEGACY_SOURCE,
        "contractVersion": _LEGACY_CONTRACT_VERSION,
    }
    return payload_metadata


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


def _add_coverage_message(coverage: dict[str, Any], message: str) -> None:
    if message not in coverage["messages"]:
        coverage["messages"].append(message)


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
        raise ValueError(f"{path} is required for partial APY economics.")
    return _as_float(costs[key], path)


def _optional_cost(
    costs: Mapping[str, Any],
    key: str,
    path: str,
) -> float | None:
    if key not in costs or costs[key] is None:
        return None
    return _as_float(costs[key], path)


def _optional_dynamic_comparison(results: Mapping[str, Any]) -> Mapping[str, Any] | None:
    if "dynamicComparison" not in results or results["dynamicComparison"] is None:
        return None
    if not isinstance(results["dynamicComparison"], Mapping):
        raise ValueError("result_bundle.results.dynamicComparison must be a mapping.")
    return results["dynamicComparison"]


def _optional_numeric_value(
    mapping: Mapping[str, Any],
    key: str,
    path: str,
) -> float | None:
    if key not in mapping or mapping[key] is None:
        return None
    return _as_float(mapping[key], path)


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
