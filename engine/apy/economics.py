from __future__ import annotations

from copy import deepcopy
import math
from typing import Any

import pandas as pd

from engine.apy.costing import (
    TARGET_CURRENCY,
    TARGET_PRICE_YEAR,
    build_cost_item,
    normalise_cost_table,
    unresolved_cost_warnings,
    valid_converted_cost,
)
from engine.apy.runner import run_scenario_with_do_nothing
from engine.apy.scenario import DIRECT_EFFECTS_SCOPE_STATEMENT


ECONOMICS_CONTRACT_VERSION = "apy_v9_economics_results_v1"
ECONOMICS_SOURCE = "run_health_economics_v9_python_port"


def build_default_economics_config() -> dict[str, Any]:
    return {
        "metadata": {
            "currencyCode": TARGET_CURRENCY,
            "priceYear": TARGET_PRICE_YEAR,
            "targetCurrency": TARGET_CURRENCY,
            "targetPriceYear": TARGET_PRICE_YEAR,
            "perspective": "Australian health-system perspective",
            "locationLabel": "Australia",
            "sourceNotes": "",
            "programCostBasis": "",
            "scopeStatement": DIRECT_EFFECTS_SCOPE_STATEMENT,
        },
        "costNormalisation": {
            "contractVersion": "ltbi_cost_normalisation_v1",
            "formula": "converted_cost = original_cost * target_year_index / source_year_index",
            "defaultInflationIndexId": "AUD_HEALTH_SYSTEM_CPI",
            "indexTable": "data/inflation_indices.csv",
            "notes": (
                "Price-year conversion is separate from future discounting. "
                "Projected model costs remain in constant target-year prices."
            ),
        },
        "costItems": [],
        "discounting": {
            "selectedAnnualRate": 0.03,
            "availableAnnualRates": [0.03, 0.0],
            "primaryDisplayedRate": 0.03,
            "comparisonRate": 0.0,
            "notes": "Discounting applies to future costs and health outcomes only, not source-price conversion.",
        },
        "healthOutcome": {
            "primary": "DALYs averted",
            "unresolvedInputs": [],
            "notes": "Natural epidemiological outputs remain reported alongside DALYs.",
        },
        "threshold": {
            "metric": "cost per DALY averted",
            "concept": "1 x Australian GDP per capita per DALY averted",
            "value": None,
            "currency": TARGET_CURRENCY,
            "referenceYear": None,
            "source": "",
            "notes": "Illustrative benchmark only; not an official Australian funding threshold.",
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
            "priceYear": TARGET_PRICE_YEAR,
            "targetCurrency": TARGET_CURRENCY,
            "targetPriceYear": TARGET_PRICE_YEAR,
            "locationLabel": "Australia",
            "programCostBasis": "total",
            "sourceNotes": (
                "KWAB150 preset populated from local data/costs.csv: "
                "cscreenqft, cscreentst, ctreat*, and ctb mid values. "
                "The source price year is unresolved in the repository, so "
                "these costs are not treated as 2025-26 values until a verified "
                "source year and inflation-index values are supplied."
            ),
        }
    )
    config["costs"]["test"].update({"IGRA": 113.48, "TST": 116.07})
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
    config["costItems"] = normalise_cost_table(
        [
            build_cost_item(
                cost_item_id="test_igra",
                description="IGRA screening test per person",
                original_cost=113.48,
                original_currency="AUD",
                original_price_year=None,
                source="Local data/costs.csv row cscreenqft",
                page_table_item="cscreenqft",
                source_year_status="unknown",
                resource_category="screening_test",
                notes="Includes test kit and lab processing according to local row description.",
            ),
            build_cost_item(
                cost_item_id="test_tst",
                description="TST screening per person",
                original_cost=116.07,
                original_currency="AUD",
                original_price_year=None,
                source="Local data/costs.csv row cscreentst",
                page_table_item="cscreentst",
                source_year_status="unknown",
                resource_category="screening_test",
                notes="Local row description says TST includes PPD injection and reading visit.",
                resource_use={"returnVisitForReading": True},
            ),
            build_cost_item(
                cost_item_id="regimen_3hp",
                description="3HP preventive regimen per started course",
                original_cost=165.5072,
                original_currency="AUD",
                original_price_year=None,
                source="Local data/costs.csv row ctreat3HP",
                page_table_item="ctreat3HP",
                source_year_status="unknown",
                resource_category="preventive_regimen",
            ),
            build_cost_item(
                cost_item_id="regimen_4r",
                description="4R preventive regimen per started course",
                original_cost=123.3172,
                original_currency="AUD",
                original_price_year=None,
                source="Local data/costs.csv row ctreat4R",
                page_table_item="ctreat4R",
                source_year_status="unknown",
                resource_category="preventive_regimen",
            ),
            build_cost_item(
                cost_item_id="regimen_3hr",
                description="3HR preventive regimen per started course",
                original_cost=134.2272,
                original_currency="AUD",
                original_price_year=None,
                source="Local data/costs.csv row ctreat3HR",
                page_table_item="ctreat3HR",
                source_year_status="unknown",
                resource_category="preventive_regimen",
            ),
            build_cost_item(
                cost_item_id="regimen_6h",
                description="6H preventive regimen per started course",
                original_cost=187.7508,
                original_currency="AUD",
                original_price_year=None,
                source="Local data/costs.csv row ctreat6H",
                page_table_item="ctreat6H",
                source_year_status="unknown",
                resource_category="preventive_regimen",
            ),
            build_cost_item(
                cost_item_id="regimen_9h",
                description="9H preventive regimen per started course",
                original_cost=254.8544,
                original_currency="AUD",
                original_price_year=None,
                source="Local data/costs.csv row ctreat9H",
                page_table_item="ctreat9H",
                source_year_status="unknown",
                resource_category="preventive_regimen",
            ),
            build_cost_item(
                cost_item_id="active_tb_disease",
                description="Active TB disease management per case",
                original_cost=19079.6,
                original_currency="AUD",
                original_price_year=None,
                source="Local data/costs.csv row ctb",
                page_table_item="ctb",
                source_year_status="unknown",
                resource_category="tb_disease_care",
            ),
        ]
    )
    return config


def validate_economics_config(econ_config: dict[str, Any]) -> dict[str, Any]:
    errors: list[dict[str, str]] = []
    warnings: list[dict[str, str]] = []

    if not isinstance(econ_config, dict):
        errors.append(_issue("econConfig", "Economics config must be a dict."))
        return _report(errors, warnings)

    metadata = econ_config.get("metadata")
    if not isinstance(metadata, dict):
        errors.append(_issue("metadata", "metadata must be a dict."))
        metadata = {}
    else:
        for field in (
            "currencyCode",
            "targetCurrency",
            "targetPriceYear",
            "locationLabel",
            "sourceNotes",
            "programCostBasis",
            "perspective",
            "scopeStatement",
        ):
            _validate_text_field(errors, metadata, field, f"metadata.{field}")
        if metadata.get("priceYear") not in (None, "", []) and not isinstance(
            metadata.get("priceYear"), (str, int, float)
        ):
            errors.append(_issue("metadata.priceYear", "metadata.priceYear must be a year label."))

    costs = econ_config.get("costs")
    if not isinstance(costs, dict):
        errors.append(_issue("costs", "costs must be a dict."))
        costs = {}

    test_costs = costs.get("test")
    if not isinstance(test_costs, dict):
        errors.append(_issue("costs.test", "costs.test must be a dict."))
        test_costs = {}
    else:
        for field in ("IGRA", "TST"):
            _validate_optional_cost(errors, test_costs, field, f"costs.test.{field}")

    regimen_costs = costs.get("regimen")
    if not isinstance(regimen_costs, dict):
        errors.append(_issue("costs.regimen", "costs.regimen must be a dict."))
        regimen_costs = {}
    else:
        for field in ("x3HP", "x4R", "x3HR", "x6H", "x9H"):
            _validate_optional_cost(
                errors,
                regimen_costs,
                field,
                f"costs.regimen.{field}",
            )

    for field in (
        "falsePositiveIncrementalPerPerson",
        "activeTBDiseasePerCase",
        "programSetupTotal",
        "programRunningTotal",
    ):
        _validate_optional_cost(errors, costs, field, f"costs.{field}")

    cost_items = econ_config.get("costItems")
    if cost_items not in (None, []):
        if not isinstance(cost_items, list):
            errors.append(_issue("costItems", "costItems must be a list."))
        else:
            for idx, item in enumerate(cost_items):
                if not isinstance(item, dict):
                    errors.append(_issue(f"costItems.{idx}", "Each cost item must be a dict."))
                elif item.get("conversionStatus") != "valid":
                    warnings.append(
                        _issue(
                            f"costItems.{idx}",
                            f"{item.get('costItemId', idx)} conversion is unresolved.",
                            "warning",
                        )
                    )

    discounting = econ_config.get("discounting", {})
    if isinstance(discounting, dict):
        _validate_optional_cost(
            errors,
            discounting,
            "selectedAnnualRate",
            "discounting.selectedAnnualRate",
        )
        rate = _coerce_number(discounting.get("selectedAnnualRate"))
        if rate is not None and rate not in {0, 0.03}:
            warnings.append(
                _issue(
                    "discounting.selectedAnnualRate",
                    "Primary milestone rates are 0% and 3%; other numeric rates are retained for later review.",
                    "warning",
                )
            )

    threshold = econ_config.get("threshold", {})
    if isinstance(threshold, dict) and _coerce_number(threshold.get("value")) is None:
        warnings.append(
            _issue(
                "threshold.value",
                "GDP-per-capita willingness-to-pay value is unresolved; enter a verified value and provenance before using threshold-based conclusions.",
                "warning",
            )
        )

    return _report(errors, warnings)


def run_health_economics(
    results: dict[str, Any],
    econ_config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    if econ_config is None:
        econ_config = build_default_economics_config()
    else:
        econ_config = _normalise_empty_to_none(deepcopy(econ_config))

    validation_report = validate_economics_config(econ_config)
    if not validation_report["isValid"]:
        raise ValueError("Economics config contains fatal validation errors.")

    status = _empty_status(validation_report)
    test_type = get_config_text(results, "testType")
    regimen = get_config_text(results, "regimen")

    quantities = {
        "nScreened": summary_metric(results, "nScreened"),
        "nTotalCoursesStarted": summary_metric(results, "nTotalCoursesStarted"),
        "nFalsePositiveTreated": summary_metric(results, "nFalsePositiveTreated"),
        "nCuredInfection": summary_metric(results, "nCuredInfection"),
        "nPreventedActiveTB": summary_metric(results, "nPreventedActiveTB"),
    }
    baseline_cases, intervention_cases = baseline_intervention_cases(results)
    quantities["baselineActiveTBCases"] = baseline_cases
    quantities["interventionActiveTBCases"] = intervention_cases

    test_cost = selected_test_cost(econ_config, test_type, status)
    regimen_cost = selected_regimen_cost(econ_config, regimen, status)
    costs = econ_config.get("costs", {})
    normalised_items = normalise_cost_table(econ_config.get("costItems") or [])
    cost_item_lookup = {item.get("costItemId"): item for item in normalised_items}
    converted_test_cost = selected_converted_test_cost(cost_item_lookup, test_type)
    converted_regimen_cost = selected_converted_regimen_cost(cost_item_lookup, regimen)
    converted_active_tb_cost = valid_converted_cost(cost_item_lookup.get("active_tb_disease", {}))
    if normalised_items:
        test_cost = converted_test_cost
        regimen_cost = converted_regimen_cost
    active_tb_cost = converted_active_tb_cost if normalised_items else _optional_cost(costs, "activeTBDiseasePerCase")

    status["messages"].extend(unresolved_cost_warnings(normalised_items))

    unit_costs = {
        "testPerPerson": test_cost,
        "treatmentPerStarted": regimen_cost,
        "falsePositiveIncrementalPerPerson": _optional_cost(
            costs,
            "falsePositiveIncrementalPerPerson",
        ),
        "activeTBDiseasePerCase": active_tb_cost,
    }

    cost_values = {
        "testingCost": multiply_if_available(quantities["nScreened"], test_cost),
        "treatmentCost": multiply_if_available(
            quantities["nTotalCoursesStarted"],
            regimen_cost,
        ),
        "falsePositiveIncrementalCost": multiply_if_available(
            quantities["nFalsePositiveTreated"],
            unit_costs["falsePositiveIncrementalPerPerson"],
        ),
        "programSetupCost": _optional_cost(costs, "programSetupTotal"),
        "programRunningCost": _optional_cost(costs, "programRunningTotal"),
        "tbDiseaseCostsAverted": multiply_if_available(
            quantities["nPreventedActiveTB"],
            active_tb_cost,
        ),
        "baselineTBDiseaseCost": multiply_if_available(
            quantities["baselineActiveTBCases"],
            active_tb_cost,
        ),
        "interventionTBDiseaseCost": multiply_if_available(
            quantities["interventionActiveTBCases"],
            active_tb_cost,
        ),
    }

    _note_missing_if_empty(
        status,
        cost_values["testingCost"],
        "testingCost",
        "Testing cost not calculated because nScreened or selected test cost is missing.",
    )
    _note_missing_if_empty(
        status,
        cost_values["treatmentCost"],
        "treatmentCost",
        "Treatment cost not calculated because starts or selected regimen cost is missing.",
    )
    if _is_missing_number(unit_costs["falsePositiveIncrementalPerPerson"]):
        status["messages"].append(
            "No false-positive incremental unit cost supplied; no extra false-positive cost is included."
        )
    _note_missing_if_empty(
        status,
        cost_values["tbDiseaseCostsAverted"],
        "tbDiseaseCostsAverted",
        "TB disease costs averted not calculated because prevented cases or disease cost is missing.",
    )
    if _is_missing_number(cost_values["baselineTBDiseaseCost"]) or _is_missing_number(
        cost_values["interventionTBDiseaseCost"]
    ):
        status["partialCalculations"].append(
            "Baseline/intervention TB disease costs were not calculated because case counts are unavailable."
        )

    cost_values["totalProgramCost"] = sum_required_costs(
        [
            cost_values["testingCost"],
            cost_values["treatmentCost"],
            zero_if_empty(cost_values["falsePositiveIncrementalCost"]),
            cost_values["programSetupCost"],
            cost_values["programRunningCost"],
        ]
    )
    _note_missing_if_empty(
        status,
        cost_values["totalProgramCost"],
        "totalProgramCost",
        "Total program cost not calculated because one or more required cost components is missing.",
    )

    cost_values["netCostVsBaseline"] = subtract_if_available(
        cost_values["totalProgramCost"],
        cost_values["tbDiseaseCostsAverted"],
    )
    _note_missing_if_empty(
        status,
        cost_values["netCostVsBaseline"],
        "netCostVsBaseline",
        "Net cost versus baseline not calculated because total program cost or TB costs averted is missing.",
    )

    cost_effectiveness = {
        "costPerInfectionCured": divide_if_available(
            cost_values["netCostVsBaseline"],
            quantities["nCuredInfection"],
            "costPerInfectionCured",
            status,
        ),
        "costPerTBCasePrevented": divide_if_available(
            cost_values["netCostVsBaseline"],
            quantities["nPreventedActiveTB"],
            "costPerTBCasePrevented",
            status,
        ),
    }

    econ = {
        "available": True,
        "source": ECONOMICS_SOURCE,
        "contractVersion": ECONOMICS_CONTRACT_VERSION,
        "metadata": deepcopy(econ_config.get("metadata", {})),
        "inputs": deepcopy(econ_config),
        "costNormalisation": deepcopy(econ_config.get("costNormalisation", {})),
        "costItems": normalised_items,
        "discounting": deepcopy(econ_config.get("discounting", {})),
        "healthOutcome": deepcopy(econ_config.get("healthOutcome", {})),
        "threshold": deepcopy(econ_config.get("threshold", {})),
        "scopeStatement": econ_config.get("metadata", {}).get(
            "scopeStatement", DIRECT_EFFECTS_SCOPE_STATEMENT
        ),
        "strategy": {
            "testType": test_type,
            "regimen": regimen,
        },
        "quantities": quantities,
        "unitCosts": unit_costs,
        "costs": cost_values,
        "costEffectiveness": cost_effectiveness,
        "status": status,
    }
    econ["summaryRows"] = build_summary_rows(econ)
    econ["summaryTable"] = econ["summaryRows"]
    status["isComplete"] = not status["missingInputs"] and not status["notCalculated"]
    return econ


def run_health_economics_for_config(
    config: dict[str, Any],
    econ_config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    out = run_scenario_with_do_nothing(config)
    econ = run_health_economics(out, econ_config)
    econ["scenario"] = {
        "modelVersion": out["results"].get("modelVersion"),
        "backend": out["results"].get("backend"),
        "scenarioLabel": out["results"].get("interfaceConfig", {}).get("scenarioLabel"),
    }
    return econ


def summary_metric(results: dict[str, Any], metric_name: str) -> float | None:
    rows = _summary_rows(results)
    for row in rows:
        if str(row.get("Metric", row.get("metric", ""))) == metric_name:
            return _coerce_number(row.get("Median", row.get("median")))
    return None


def get_config_text(results: dict[str, Any], field_name: str) -> str:
    for config in _candidate_interface_configs(results):
        value = config.get(field_name)
        if value not in (None, "", []):
            return str(value)
    return ""


def selected_test_cost(
    econ_config: dict[str, Any],
    test_type: str,
    status: dict[str, Any],
) -> float | None:
    if not test_type:
        status["missingInputs"].append("results.interfaceConfig.testType")
        return None
    costs = econ_config.get("costs", {}).get("test", {})
    value = _coerce_number(costs.get(str(test_type).upper()))
    if value is None:
        status["missingInputs"].append(f"costs.test.{test_type}")
    return value


def selected_regimen_cost(
    econ_config: dict[str, Any],
    regimen: str,
    status: dict[str, Any],
) -> float | None:
    if not regimen:
        status["missingInputs"].append("results.interfaceConfig.regimen")
        return None
    field = _regimen_cost_field(regimen)
    if field is None:
        status["missingInputs"].append(f"Unknown regimen label: {regimen}")
        return None
    costs = econ_config.get("costs", {}).get("regimen", {})
    value = _coerce_number(costs.get(field))
    if value is None:
        status["missingInputs"].append(f"costs.regimen.{field}")
    return value


def selected_converted_test_cost(
    cost_items: dict[str, dict[str, Any]],
    test_type: str,
) -> float | None:
    key = {"IGRA": "test_igra", "TST": "test_tst"}.get(str(test_type).upper())
    if key is None:
        return None
    return valid_converted_cost(cost_items.get(key, {}))


def selected_converted_regimen_cost(
    cost_items: dict[str, dict[str, Any]],
    regimen: str,
) -> float | None:
    key = {
        "3HP": "regimen_3hp",
        "4R": "regimen_4r",
        "3HR": "regimen_3hr",
        "6H": "regimen_6h",
        "9H": "regimen_9h",
    }.get(str(regimen).upper())
    if key is None:
        return None
    return valid_converted_cost(cost_items.get(key, {}))


def baseline_intervention_cases(results: dict[str, Any]) -> tuple[float | None, float | None]:
    derived = _do_nothing_derived(results)
    if derived is not None:
        baseline = _median_from_rows_or_frame(derived, "nActiveBy20y_DoNothing")
        intervention = _median_from_rows_or_frame(derived, "nActiveBy20y_AfterStrategy")
        if baseline is not None or intervention is not None:
            return baseline, intervention

    dynamic = _dynamic_comparison(results)
    if dynamic:
        baseline = _coerce_number(dynamic.get("cumulative_baseline_active_tb_cases"))
        intervention = _coerce_number(
            dynamic.get("cumulative_intervention_active_tb_cases")
        )
        if baseline is not None or intervention is not None:
            return baseline, intervention

    raw = _raw_frame(results)
    if raw is not None and {"nActiveBy20y", "nPreventedActiveTB"}.issubset(raw.columns):
        baseline_values = pd.to_numeric(raw["nActiveBy20y"], errors="coerce")
        intervention_values = baseline_values - pd.to_numeric(
            raw["nPreventedActiveTB"],
            errors="coerce",
        )
        return _median_series(baseline_values), _median_series(intervention_values)

    return None, None


def multiply_if_available(a: Any, b: Any) -> float | None:
    a_num = _coerce_number(a)
    b_num = _coerce_number(b)
    if a_num is None or b_num is None:
        return None
    return a_num * b_num


def subtract_if_available(a: Any, b: Any) -> float | None:
    a_num = _coerce_number(a)
    b_num = _coerce_number(b)
    if a_num is None or b_num is None:
        return None
    return a_num - b_num


def divide_if_available(
    a: Any,
    b: Any,
    field_name: str,
    status: dict[str, Any],
) -> float | None:
    a_num = _coerce_number(a)
    b_num = _coerce_number(b)
    if a_num is None or b_num is None or b_num == 0:
        status["notCalculated"].append(field_name)
        status["messages"].append(
            f"{field_name} not calculated because the denominator is missing, zero, or NaN."
        )
        return None
    return a_num / b_num


def sum_required_costs(values: list[Any]) -> float | None:
    total = 0.0
    for value in values:
        value_num = _coerce_number(value)
        if value_num is None:
            return None
        total += value_num
    return total


def zero_if_empty(value: Any) -> float:
    value_num = _coerce_number(value)
    return 0.0 if value_num is None else value_num


def build_summary_rows(econ: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    metadata = econ.get("metadata", {})
    currency = metadata.get("currencyCode", "")
    price_year = metadata.get("priceYear")
    specs = [
        ("Costs", "testingCost", econ["costs"].get("testingCost")),
        ("Costs", "treatmentCost", econ["costs"].get("treatmentCost")),
        (
            "Costs",
            "falsePositiveIncrementalCost",
            econ["costs"].get("falsePositiveIncrementalCost"),
        ),
        ("Costs", "programSetupCost", econ["costs"].get("programSetupCost")),
        ("Costs", "programRunningCost", econ["costs"].get("programRunningCost")),
        ("Costs", "baselineTBDiseaseCost", econ["costs"].get("baselineTBDiseaseCost")),
        (
            "Costs",
            "interventionTBDiseaseCost",
            econ["costs"].get("interventionTBDiseaseCost"),
        ),
        ("Costs", "tbDiseaseCostsAverted", econ["costs"].get("tbDiseaseCostsAverted")),
        ("Costs", "totalProgramCost", econ["costs"].get("totalProgramCost")),
        ("Costs", "netCostVsBaseline", econ["costs"].get("netCostVsBaseline")),
        (
            "Cost-effectiveness",
            "costPerInfectionCured",
            econ["costEffectiveness"].get("costPerInfectionCured"),
        ),
        (
            "Cost-effectiveness",
            "costPerTBCasePrevented",
            econ["costEffectiveness"].get("costPerTBCasePrevented"),
        ),
    ]
    for section, metric, value in specs:
        rows.append(
            {
                "Section": section,
                "Metric": metric,
                "Value": _coerce_number(value),
                "Calculated": _coerce_number(value) is not None,
                "CurrencyCode": currency,
                "PriceYear": price_year,
                "Notes": "",
            }
        )
    return rows


def _issue(field: str, message: str, severity: str = "error") -> dict[str, str]:
    return {
        "field": field,
        "severity": severity,
        "message": message,
    }


def _report(
    errors: list[dict[str, str]],
    warnings: list[dict[str, str]],
) -> dict[str, Any]:
    return {
        "isValid": not errors,
        "hasWarnings": bool(warnings),
        "errors": errors,
        "warnings": warnings,
        "fatalFieldNames": [issue["field"] for issue in errors],
        "warningFieldNames": [issue["field"] for issue in warnings],
    }


def _validate_text_field(
    errors: list[dict[str, str]],
    parent: dict[str, Any],
    field: str,
    full_name: str,
) -> None:
    value = parent.get(field)
    if value in (None, "", []):
        return
    if not isinstance(value, str):
        errors.append(_issue(full_name, f"{full_name} must be text."))


def _validate_optional_numeric_scalar(
    errors: list[dict[str, str]],
    parent: dict[str, Any],
    field: str,
    full_name: str,
) -> None:
    value = parent.get(field)
    if value in (None, "", []):
        return
    if _coerce_number(value) is None:
        errors.append(_issue(full_name, f"{full_name} must be a finite scalar."))


def _validate_optional_cost(
    errors: list[dict[str, str]],
    parent: dict[str, Any],
    field: str,
    full_name: str,
) -> None:
    value = parent.get(field)
    if value in (None, "", []):
        return
    number = _coerce_number(value)
    if number is None or number < 0:
        errors.append(
            _issue(full_name, f"{full_name} must be empty or a non-negative scalar.")
        )


def _empty_status(validation_report: dict[str, Any]) -> dict[str, Any]:
    return {
        "isComplete": False,
        "missingInputs": [],
        "partialCalculations": [],
        "notCalculated": [],
        "messages": [],
        "validationReport": validation_report,
    }


def _optional_cost(parent: dict[str, Any], field: str) -> float | None:
    return _coerce_number(parent.get(field))


def _note_missing_if_empty(
    status: dict[str, Any],
    value: Any,
    name: str,
    message: str,
) -> None:
    if _is_missing_number(value):
        status["notCalculated"].append(name)
        status["messages"].append(message)


def _is_missing_number(value: Any) -> bool:
    return _coerce_number(value) is None


def _coerce_number(value: Any) -> float | None:
    if value in (None, "", []):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    if number.is_integer():
        return int(number)
    return number


def _normalise_empty_to_none(value: Any) -> Any:
    if value == []:
        return None
    if isinstance(value, dict):
        return {key: _normalise_empty_to_none(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_normalise_empty_to_none(item) for item in value]
    return value


def _summary_rows(results: dict[str, Any]) -> list[dict[str, Any]]:
    if "results" in results and isinstance(results["results"], dict):
        return _summary_rows(results["results"])
    summary = results.get("summary")
    if isinstance(summary, pd.DataFrame):
        return summary.to_dict(orient="records")
    if isinstance(summary, list):
        return summary
    headline = results.get("headline")
    if isinstance(headline, dict):
        rows = headline.get("summaryRows") or headline.get("keyMetricsRows")
        if isinstance(rows, list):
            return rows
    return []


def _candidate_interface_configs(results: dict[str, Any]) -> list[dict[str, Any]]:
    configs = []
    if "results" in results and isinstance(results["results"], dict):
        configs.extend(_candidate_interface_configs(results["results"]))
    for candidate in (
        results.get("interfaceConfig"),
        results.get("technical", {}).get("interfaceConfig")
        if isinstance(results.get("technical"), dict)
        else None,
    ):
        if isinstance(candidate, dict):
            configs.append(candidate)
    return configs


def _do_nothing_derived(results: dict[str, Any]) -> Any:
    do_nothing = results.get("doNothing")
    if isinstance(do_nothing, dict):
        return do_nothing.get("derived")
    return None


def _dynamic_comparison(results: dict[str, Any]) -> dict[str, Any]:
    if "bundle" in results and isinstance(results["bundle"], dict):
        dynamic = _dynamic_comparison(results["bundle"])
        if dynamic:
            return dynamic
    technical = results.get("technical")
    if isinstance(technical, dict) and isinstance(technical.get("dynamicComparison"), dict):
        return technical["dynamicComparison"]
    return {}


def _raw_frame(results: dict[str, Any]) -> pd.DataFrame | None:
    if "results" in results and isinstance(results["results"], dict):
        return _raw_frame(results["results"])
    raw = results.get("raw")
    if isinstance(raw, pd.DataFrame):
        return raw
    if isinstance(raw, list):
        return pd.DataFrame(raw)
    return None


def _median_from_rows_or_frame(rows_or_frame: Any, column: str) -> float | None:
    if isinstance(rows_or_frame, pd.DataFrame):
        if column not in rows_or_frame.columns:
            return None
        return _median_series(pd.to_numeric(rows_or_frame[column], errors="coerce"))
    if isinstance(rows_or_frame, list):
        values = [
            _coerce_number(row.get(column))
            for row in rows_or_frame
            if isinstance(row, dict)
        ]
        values = [value for value in values if value is not None]
        if not values:
            return None
        return _median_series(pd.Series(values, dtype="float64"))
    return None


def _median_series(values: pd.Series) -> float | None:
    clean = values.dropna()
    if clean.empty:
        return None
    value = float(clean.median())
    if value.is_integer():
        return int(value)
    return value


def _regimen_cost_field(regimen: str) -> str | None:
    normalised = str(regimen).strip().upper()
    mapping = {
        "3HP": "x3HP",
        "4R": "x4R",
        "3HR": "x3HR",
        "6H": "x6H",
        "9H": "x9H",
    }
    return mapping.get(normalised)
