from __future__ import annotations

from copy import deepcopy
import math
from typing import Any

from engine.apy.economics import summary_metric


def build_default_health_economic_assumptions() -> dict[str, Any]:
    return {
        "metadata": {
            "currencyCode": "AUD",
            "priceYear": 2019,
            "perspective": "Australian health-care system",
            "source": "Dale et al. KWAB150-informed assumptions",
            "discountRate": 0.03,
        },
        "daly": {
            "activeTBDisabilityWeight": 0.333,
            "activeTBDurationYears": 0.5,
            "includeYLL": True,
            "tbCaseFatalityRisk": 0.074,
            "yllPerTBDeath": 20.0,
            "notes": (
                "tbCaseFatalityRisk and yllPerTBDeath should be sensitivity-tested; "
                "default values are placeholders until APY-specific mortality/YLL "
                "assumptions are finalised."
            ),
        },
        "qaly": {
            "healthyUtility": 0.8733,
            "activeTBUtility": 0.8182,
            "activeTBDurationYears": 0.5,
            "qalyActiveTbDurationYears": 0.5,
            "saeUtility": 0.8685,
            "saeDurationYears": 7 / 365,
            "ltbiTreatmentUtility": 0.8733,
            "ltbiTreatmentDecrement": 0.0,
            "ltbiTreatmentDecrementSensitivity": 0.0133,
            "includeMortalityInQaly": True,
            "qalyMortalityMethod": "yll_times_healthy_utility",
            "qualityAdjustedLifeExpectancyPerTBDeath": None,
            "gbdAlignedQalyMorbiditySensitivity": True,
            "activeTbMorbidityQalyMethod": "dale_bauer",
        },
        "postTB": {
            "includePostTBSequelae": False,
            "postTBDalysPerCase": 5.8,
            "postTBQalysLostPerCase": 0.93,
            "totalTBQalysLostPerCase": 1.93,
            "postTBExcessMortalityMultiplier": 1.16,
            "postTBDurationYears": 10,
            "postTBDurationYearsSensitivityOptions": [10, 20, "lifetime"],
            "postTBLifetimeDurationYearsEquivalent": None,
            "postTBAnnualCareCost": 0.0,
            "pPTLD": None,
            "pPTLD_low": 0.20,
            "pPTLD_high": 0.68,
            "notes": (
                "Post-TB sequelae are handled as a scenario layer on top of acute APY "
                "health outcomes. APY-specific PTLD prevalence, duration, excess "
                "mortality, and annual care cost inputs are not yet available, so "
                "defaults should be treated as sensitivity assumptions only."
            ),
        },
        "thresholds": {
            "wtpLow": 45000,
            "wtpHigh": 75000,
        },
    }


def calculate_health_outcomes(
    results: dict[str, Any],
    assumptions: dict[str, Any] | None = None,
) -> dict[str, Any]:
    assumptions = _merge_assumptions(assumptions)
    notes: list[str] = []
    missing_inputs: list[str] = []

    tb_cases_prevented = _first_number(
        summary_metric(results, "nPreventedActiveTB"),
        _dynamic_value(results, "cumulative_cases_averted"),
    )
    n_courses_started = summary_metric(results, "nTotalCoursesStarted")
    n_adr_stop = summary_metric(results, "nADRstop")
    treatment_duration_years = _treatment_duration_years(results)

    if tb_cases_prevented is None:
        missing_inputs.append("nPreventedActiveTB")

    daly_assumptions = assumptions["daly"]
    qaly_assumptions = assumptions["qaly"]
    quality_adjusted_life_expectancy = _first_number(
        qaly_assumptions.get("qualityAdjustedLifeExpectancyPerTBDeath"),
        _multiply(
            daly_assumptions["yllPerTBDeath"],
            qaly_assumptions["healthyUtility"],
        ),
    )

    tb_deaths_prevented = _multiply(
        tb_cases_prevented,
        daly_assumptions["tbCaseFatalityRisk"],
    )
    yld_averted = _multiply(
        tb_cases_prevented,
        daly_assumptions["activeTBDisabilityWeight"],
        daly_assumptions["activeTBDurationYears"],
    )
    yll_averted = (
        _multiply(tb_deaths_prevented, daly_assumptions["yllPerTBDeath"])
        if daly_assumptions.get("includeYLL", True)
        else 0.0
    )
    dalys_averted = _add_if_available(yld_averted, yll_averted)

    active_tb_morbidity_qalys_gained_dale = _multiply(
        tb_cases_prevented,
        qaly_assumptions["healthyUtility"] - qaly_assumptions["activeTBUtility"],
        qaly_assumptions.get(
            "qalyActiveTbDurationYears",
            qaly_assumptions["activeTBDurationYears"],
        ),
    )
    active_tb_morbidity_qalys_gained_gbd = _multiply(
        tb_cases_prevented,
        daly_assumptions["activeTBDisabilityWeight"],
        daly_assumptions["activeTBDurationYears"],
    )
    mortality_qalys_gained = (
        _multiply(tb_deaths_prevented, quality_adjusted_life_expectancy)
        if qaly_assumptions.get("includeMortalityInQaly", True)
        else 0.0
    )
    qaly_loss_treatment = _multiply(
        n_courses_started,
        qaly_assumptions["ltbiTreatmentDecrement"],
        treatment_duration_years,
    )
    qaly_loss_sae = _multiply(
        n_adr_stop,
        qaly_assumptions["healthyUtility"] - qaly_assumptions["saeUtility"],
        qaly_assumptions["saeDurationYears"],
    )
    qaly_non_mortality_losses = _add_if_available(qaly_loss_treatment, qaly_loss_sae)
    morbidity_only_qalys_gained = _subtract_if_available(
        active_tb_morbidity_qalys_gained_dale,
        qaly_non_mortality_losses,
    )
    qalys_gained_dale_mortality_adjusted = _subtract_if_available(
        _add_if_available(active_tb_morbidity_qalys_gained_dale, mortality_qalys_gained),
        _add_if_available(qaly_loss_treatment, qaly_loss_sae),
    )
    qalys_gained_gbd_aligned_mortality_adjusted = _subtract_if_available(
        _add_if_available(active_tb_morbidity_qalys_gained_gbd, mortality_qalys_gained),
        qaly_non_mortality_losses,
    )
    qalys_gained = (
        qalys_gained_dale_mortality_adjusted
        if qaly_assumptions.get("includeMortalityInQaly", True)
        else morbidity_only_qalys_gained
    )

    if daly_assumptions.get("includeYLL", True):
        notes.append(
            "DALYs include YLD plus configurable YLL from assumed TB case fatality and YLL per TB death."
        )
    else:
        notes.append("DALYs use YLD only; YLL is excluded for this sensitivity analysis.")
    if qaly_assumptions.get("includeMortalityInQaly", True):
        notes.append(
            "Primary QALYs include mortality benefits using the same TB fatality and YLL assumptions as DALYs."
        )
    else:
        notes.append("Primary QALYs exclude mortality benefits for this sensitivity analysis.")
    notes.append(
        "Dale-compatible QALYs use the Dale/Bauer active TB utility decrement; GBD-aligned QALYs use the active TB disability weight as a morbidity sensitivity."
    )
    if qaly_assumptions["ltbiTreatmentDecrement"] == 0:
        notes.append("Base case assumes no LTBI treatment utility decrement.")
    notes.append(
        "Post-TB sequelae are available as an optional scenario layer and do not change the underlying APY epidemiological simulation."
    )

    return {
        "assumptions": assumptions,
        "healthOutcomes": {
            "tbCasesPrevented": tb_cases_prevented,
            "tbDeathsPrevented": tb_deaths_prevented,
            "yldAverted": yld_averted,
            "yllAverted": yll_averted,
            "dalysAverted": dalys_averted,
            "qualityAdjustedLifeExpectancyPerTBDeath": quality_adjusted_life_expectancy,
            "mortalityQalysGained": mortality_qalys_gained,
            "activeTbMorbidityQalysGained_Dale": active_tb_morbidity_qalys_gained_dale,
            "activeTbMorbidityQalysGained_GBD": active_tb_morbidity_qalys_gained_gbd,
            "morbidityOnlyQalysGained": morbidity_only_qalys_gained,
            "qalyLossActiveTBAverted": active_tb_morbidity_qalys_gained_dale,
            "qalyLossTreatment": qaly_loss_treatment,
            "qalyLossSAE": qaly_loss_sae,
            "treatmentDecrementQalyLoss": qaly_loss_treatment,
            "saeQalyLoss": qaly_loss_sae,
            "qalysGained_DaleMortalityAdjusted": qalys_gained_dale_mortality_adjusted,
            "qalysGained_GBDAlignedMortalityAdjusted": qalys_gained_gbd_aligned_mortality_adjusted,
            "qalysGained": qalys_gained,
        },
        "status": {
            "missingInputs": missing_inputs,
            "notes": notes,
        },
    }


def calculate_post_tb_sequelae(
    tb_cases_prevented: Any,
    acute_dalys_averted: Any,
    acute_qalys_gained: Any,
    acute_net_cost: Any,
    assumptions: dict[str, Any] | None = None,
) -> dict[str, Any]:
    assumptions = _merge_assumptions(assumptions)
    post_tb_assumptions = assumptions["postTB"]
    notes: list[str] = []

    include_post_tb = bool(post_tb_assumptions.get("includePostTBSequelae", False))
    tb_cases = _number_or_none(tb_cases_prevented)
    acute_dalys = _number_or_none(acute_dalys_averted)
    acute_qalys = _number_or_none(acute_qalys_gained)
    acute_cost = _number_or_none(acute_net_cost)
    post_tb_daly_per_case = _number_or_none(post_tb_assumptions.get("postTBDalysPerCase"))
    post_tb_qaly_per_case = _number_or_none(
        post_tb_assumptions.get("postTBQalysLostPerCase")
    )
    post_tb_annual_care_cost = _first_number(
        post_tb_assumptions.get("postTBAnnualCareCost"),
        0.0,
    )
    post_tb_duration_years = _duration_years_or_none(
        post_tb_assumptions.get("postTBDurationYears"),
        post_tb_assumptions.get("postTBLifetimeDurationYearsEquivalent"),
    )

    post_tb_dalys_averted = (
        _multiply(tb_cases, post_tb_daly_per_case) if include_post_tb else 0.0
    )
    post_tb_qalys_gained = (
        _multiply(tb_cases, post_tb_qaly_per_case) if include_post_tb else 0.0
    )
    total_dalys_including_post_tb = _add_optional_zero_safe(
        acute_dalys,
        post_tb_dalys_averted,
    )
    total_qalys_including_post_tb = _add_optional_zero_safe(
        acute_qalys,
        post_tb_qalys_gained,
    )
    post_tb_costs_averted = (
        _multiply(tb_cases, post_tb_annual_care_cost, post_tb_duration_years)
        if include_post_tb and post_tb_duration_years is not None
        else 0.0
    )
    net_cost_including_post_tb = _subtract_if_available(acute_cost, post_tb_costs_averted)

    if include_post_tb:
        notes.append(
            "Post-TB sequelae scenario applied as a tail on prevented incident TB cases; acute APY model outputs are unchanged."
        )
    else:
        notes.append(
            "Post-TB sequelae disabled; acute-only DALY/QALY and cost results are unchanged."
        )
    if (
        include_post_tb
        and str(post_tb_assumptions.get("postTBDurationYears")).strip().lower() == "lifetime"
        and post_tb_duration_years is None
    ):
        notes.append(
            "Lifetime post-TB duration was requested without a numeric equivalent, so post-TB annual care costs were not extrapolated."
        )
    if _first_number(post_tb_assumptions.get("pPTLD")) is None:
        notes.append(
            "PTLD prevalence remains unresolved for APY-specific use; the low/high evidence bounds are recorded for scenario review."
        )

    return {
        "assumptions": assumptions,
        "postTBScenarios": {
            "includePostTBSequelae": include_post_tb,
            "tbCasesPrevented": tb_cases,
            "acuteDALYsAverted": acute_dalys,
            "postTBDALYsAverted": post_tb_dalys_averted,
            "totalDALYsAvertedIncludingPostTB": total_dalys_including_post_tb,
            "acuteQALYsGained": acute_qalys,
            "postTBQALYsGained": post_tb_qalys_gained,
            "totalQALYsGainedIncludingPostTB": total_qalys_including_post_tb,
            "postTBCostsAverted": post_tb_costs_averted,
            "acuteNetCost": acute_cost,
            "netCostIncludingPostTB": net_cost_including_post_tb,
            "costPerDALYIncludingPostTB": _divide(
                net_cost_including_post_tb,
                total_dalys_including_post_tb,
            ),
            "costPerQALYIncludingPostTB": _divide(
                net_cost_including_post_tb,
                total_qalys_including_post_tb,
            ),
            "postTBDurationYearsApplied": post_tb_duration_years,
            "postTBExcessMortalityMultiplier": _number_or_none(
                post_tb_assumptions.get("postTBExcessMortalityMultiplier")
            ),
            "pPTLD": _number_or_none(post_tb_assumptions.get("pPTLD")),
            "pPTLD_low": _number_or_none(post_tb_assumptions.get("pPTLD_low")),
            "pPTLD_high": _number_or_none(post_tb_assumptions.get("pPTLD_high")),
        },
        "status": {
            "notes": notes,
        },
    }


def calculate_icers(
    economics: dict[str, Any],
    health_outcomes: dict[str, Any],
    assumptions: dict[str, Any] | None = None,
) -> dict[str, Any]:
    assumptions = _merge_assumptions(assumptions)
    costs = economics.get("costs", {})
    health = health_outcomes.get("healthOutcomes", health_outcomes)
    notes: list[str] = []
    missing_inputs: list[str] = []

    incremental_cost = _first_number(
        costs.get("netCostVsBaseline"),
        _subtract_if_available(
            costs.get("totalProgramCost"),
            costs.get("tbDiseaseCostsAverted"),
        ),
    )
    if incremental_cost is None:
        missing_inputs.append("incrementalCost")

    tb_cases_prevented = _number_or_none(health.get("tbCasesPrevented"))
    dalys_averted = _number_or_none(health.get("dalysAverted"))
    qalys_gained = _first_number(
        health.get("qalysGained_DaleMortalityAdjusted"),
        health.get("qalysGained"),
    )
    qalys_gained_gbd_aligned = _number_or_none(
        health.get("qalysGained_GBDAlignedMortalityAdjusted")
    )
    thresholds = assumptions["thresholds"]

    icers = {
        "costPerDALYAverted": _divide(incremental_cost, dalys_averted),
        "costPerQALYGained": _divide(incremental_cost, qalys_gained),
        "costPerQALYGained_GBDAligned": _divide(
            incremental_cost,
            qalys_gained_gbd_aligned,
        ),
        "costPerTBCasePrevented": _divide(incremental_cost, tb_cases_prevented),
    }
    nmb = {
        "netMonetaryBenefitQALY_low": _nmb(
            qalys_gained,
            thresholds["wtpLow"],
            incremental_cost,
        ),
        "netMonetaryBenefitQALY_high": _nmb(
            qalys_gained,
            thresholds["wtpHigh"],
            incremental_cost,
        ),
        "netMonetaryBenefitQALY_GBDAligned_low": _nmb(
            qalys_gained_gbd_aligned,
            thresholds["wtpLow"],
            incremental_cost,
        ),
        "netMonetaryBenefitQALY_GBDAligned_high": _nmb(
            qalys_gained_gbd_aligned,
            thresholds["wtpHigh"],
            incremental_cost,
        ),
        "netMonetaryBenefitDALY_low": _nmb(
            dalys_averted,
            thresholds["wtpLow"],
            incremental_cost,
        ),
        "netMonetaryBenefitDALY_high": _nmb(
            dalys_averted,
            thresholds["wtpHigh"],
            incremental_cost,
        ),
    }
    for metric, value in {
        "dalysAverted": dalys_averted,
        "qalysGained": qalys_gained,
        "tbCasesPrevented": tb_cases_prevented,
    }.items():
        if value in (None, 0):
            notes.append(f"{metric} is missing or zero; some ICERs may be unavailable.")

    result = {
        "inputs": {
            "incrementalCost": incremental_cost,
            "wtpLow": thresholds["wtpLow"],
            "wtpHigh": thresholds["wtpHigh"],
        },
        "healthOutcomes": health,
        "incrementalCosts": {
            "incrementalCost": incremental_cost,
            "source": "netCostVsBaseline or totalProgramCost minus tbDiseaseCostsAverted",
        },
        "icers": icers,
        "nmb": nmb,
        "status": {
            "missingInputs": missing_inputs,
            "notes": notes + health_outcomes.get("status", {}).get("notes", []),
        },
    }
    result["summaryRows"] = _summary_rows(result)
    return result


def _summary_rows(result: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for section, values in (
        ("Health outcomes", result["healthOutcomes"]),
        ("ICERs", result["icers"]),
        ("Net monetary benefit", result["nmb"]),
    ):
        for metric, value in values.items():
            if isinstance(value, (int, float)) or value is None:
                rows.append(
                    {
                        "Section": section,
                        "Metric": metric,
                        "Value": value,
                    }
                )
    rows.append(
        {
            "Section": "Costs",
            "Metric": "incrementalCost",
            "Value": result["incrementalCosts"]["incrementalCost"],
        }
    )
    return rows


def _merge_assumptions(assumptions: dict[str, Any] | None) -> dict[str, Any]:
    merged = build_default_health_economic_assumptions()
    if assumptions:
        _deep_update(merged, assumptions)
    return merged


def _deep_update(target: dict[str, Any], updates: dict[str, Any]) -> None:
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            _deep_update(target[key], value)
        else:
            target[key] = value


def _dynamic_value(results: dict[str, Any], field: str) -> float | None:
    for candidate in (
        results.get("bundle", {}),
        results,
    ):
        if not isinstance(candidate, dict):
            continue
        technical = candidate.get("technical")
        if isinstance(technical, dict):
            dynamic = technical.get("dynamicComparison")
            if isinstance(dynamic, dict):
                value = _number_or_none(dynamic.get(field))
                if value is not None:
                    return value
    return None


def _treatment_duration_years(results: dict[str, Any]) -> float | None:
    for candidate in (
        results.get("results", {}),
        results,
    ):
        if not isinstance(candidate, dict):
            continue
        strategy = candidate.get("strategy")
        if isinstance(strategy, dict):
            months = _number_or_none(strategy.get("regimenMonths"))
            if months is not None:
                return months / 12
    return None


def _first_number(*values: Any) -> float | None:
    for value in values:
        number = _number_or_none(value)
        if number is not None:
            return number
    return None


def _number_or_none(value: Any) -> float | None:
    if value in (None, "", []):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def _multiply(*values: Any) -> float | None:
    product = 1.0
    for value in values:
        number = _number_or_none(value)
        if number is None:
            return None
        product *= number
    return product


def _add_if_available(*values: Any) -> float | None:
    total = 0.0
    for value in values:
        number = _number_or_none(value)
        if number is None:
            return None
        total += number
    return total


def _add_optional_zero_safe(*values: Any) -> float | None:
    total = 0.0
    saw_number = False
    for value in values:
        number = _number_or_none(value)
        if number is None:
            continue
        saw_number = True
        total += number
    if not saw_number:
        return None
    return total


def _subtract_if_available(a: Any, b: Any) -> float | None:
    a_num = _number_or_none(a)
    b_num = _number_or_none(b)
    if a_num is None or b_num is None:
        return None
    return a_num - b_num


def _divide(a: Any, b: Any) -> float | None:
    a_num = _number_or_none(a)
    b_num = _number_or_none(b)
    if a_num is None or b_num is None or b_num == 0:
        return None
    return a_num / b_num


def _nmb(health_gain: Any, threshold: Any, incremental_cost: Any) -> float | None:
    benefit = _multiply(health_gain, threshold)
    return _subtract_if_available(benefit, incremental_cost)


def _duration_years_or_none(value: Any, lifetime_equivalent: Any = None) -> float | None:
    numeric = _number_or_none(value)
    if numeric is not None:
        return numeric
    if isinstance(value, str) and value.strip().lower() == "lifetime":
        return _number_or_none(lifetime_equivalent)
    return None
