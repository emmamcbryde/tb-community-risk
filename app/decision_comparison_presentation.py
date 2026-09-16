from __future__ import annotations

from typing import Any


STRATEGY_FIELDS = {
    "test",
    "testType",
    "regimen",
    "screeningStrategy",
    "screenCoverage",
    "screeningNumber",
    "screeningWindowYears",
    "followUpHorizonYears",
    "pStartTPT",
    "testSensitivity",
    "testSpecificity",
    "tstSensitivity",
    "tstSpecificityBCG",
    "tstSpecificityNoBCG",
    "economics",
}


def number(value: Any) -> float | None:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def arithmetic_icer(incremental_cost: Any, dalys_averted: Any) -> float | None:
    cost = number(incremental_cost)
    dalys = number(dalys_averted)
    if cost is None or dalys is None or abs(dalys) < 1e-12:
        return None
    return cost / dalys


def quadrant(dalys_averted: Any, incremental_cost: Any) -> str:
    x = number(dalys_averted)
    y = number(incremental_cost)
    if x is None or y is None:
        return "Unavailable"
    if abs(x) < 1e-12 and abs(y) < 1e-12:
        return "Origin"
    if x > 0 and y < 0:
        return "Lower right"
    if x > 0 and y > 0:
        return "Upper right"
    if x < 0 and y > 0:
        return "Upper left"
    if x < 0 and y < 0:
        return "Lower left"
    return "On axis"


def classify_incremental(incremental_cost: Any, dalys_averted: Any) -> str:
    q = quadrant(dalys_averted, incremental_cost)
    if q == "Lower right":
        return "Dominant"
    if q == "Upper right":
        return "Higher cost with health gain"
    if q == "Upper left":
        return "Dominated"
    if q == "Lower left":
        return "Lower cost with health loss"
    if q == "Origin":
        return "Comparator"
    return "ICER not calculable"


def format_money(value: Any) -> str:
    amount = number(value)
    if amount is None:
        return "N/A"
    sign = "-" if amount < 0 else ""
    return f"{sign}AUD {abs(amount):,.0f}"


def format_dalys(value: Any) -> str:
    amount = number(value)
    if amount is None:
        return "N/A"
    return f"{amount:,.1f}"


def format_icer(value: Any) -> str:
    amount = number(value)
    if amount is None:
        return "Not applicable"
    sign = "-" if amount < 0 else ""
    return f"{sign}AUD {abs(amount):,.0f} per DALY averted"


def strategy_display_name(
    *,
    test: str,
    regimen: str,
    coverage: float,
    prioritisation: str | None = None,
) -> str:
    percent = int(round(float(coverage) * 100))
    if prioritisation:
        return f"{test} + {regimen}, {percent}% coverage, {prioritisation}"
    return f"{test} + {regimen}, {percent}% coverage"


def common_comparator_status(comparison: dict[str, Any] | None, *, tolerance: float = 1e-6) -> dict[str, Any]:
    if not isinstance(comparison, dict):
        return {"valid": False, "message": "Run a strategy comparison first."}
    scenarios = comparison.get("scenarios") or []
    if len(scenarios) < 2:
        return {"valid": False, "message": "Compare at least two strategies."}

    for scenario in scenarios:
        changed = set((scenario.get("changedFields") or {}).keys())
        non_strategy = sorted(changed - STRATEGY_FIELDS)
        if non_strategy:
            return {
                "valid": False,
                "message": (
                    "A valid strategy comparison requires the same underlying population and natural-history "
                    f"assumptions. The following non-strategy inputs differ: {', '.join(non_strategy)}."
                ),
            }

    paired = comparison.get("pairedComparisons") or []
    if comparison.get("modelType") == "agent_based" or any(row.get("modelType") == "agent_based" for row in paired):
        for row in paired:
            if row.get("comparisonDesign") != "paired_shared_baseline" or row.get("pairingValid") is not True:
                return {
                    "valid": False,
                    "message": (
                        "A valid stochastic strategy comparison requires the same baseline cohort and untreated "
                        "natural history. Return to Set up or compare strategy-only changes."
                    ),
                }

    reference_metrics = scenarios[0].get("metrics") or {}
    for scenario in scenarios[1:]:
        metrics = scenario.get("metrics") or {}
        for metric in ["comparator_active_tb", "comparatorCost", "comparatorDALYs"]:
            left = number(reference_metrics.get(metric))
            right = number(metrics.get(metric))
            if left is None or right is None:
                continue
            if abs(left - right) > tolerance:
                return {
                    "valid": False,
                    "message": (
                        "The no-screening comparator differs between strategies. Use the same population, "
                        "natural-history basis, horizon, analysis method, seed and repetitions."
                    ),
                }

    return {"valid": True, "message": "Strategies share one business-as-usual comparator."}


def comparison_plane_rows(comparison: dict[str, Any]) -> list[dict[str, Any]]:
    rows = [
        {
            "Strategy": "Business as usual",
            "DALYs averted compared with business as usual": 0.0,
            "Incremental cost compared with business as usual (AUD)": 0.0,
            "Arithmetic ICER": "Not applicable",
            "Classification": "Comparator",
            "Quadrant": "Origin",
        }
    ]
    for scenario in comparison.get("scenarios") or []:
        metrics = scenario.get("metrics") or {}
        cost = number(metrics.get("incrementalCost"))
        dalys = number(metrics.get("dalysAverted"))
        rows.append(
            {
                "Strategy": scenario.get("label") or scenario.get("scenarioId"),
                "DALYs averted compared with business as usual": dalys,
                "Incremental cost compared with business as usual (AUD)": cost,
                "Arithmetic ICER": format_icer(arithmetic_icer(cost, dalys)),
                "Classification": classify_incremental(cost, dalys),
                "Quadrant": quadrant(dalys, cost),
            }
        )
    return rows


def comparison_table_rows(comparison: dict[str, Any]) -> list[dict[str, str]]:
    out = []
    for row in comparison_plane_rows(comparison):
        out.append(
            {
                "Strategy": str(row["Strategy"]),
                "Incremental cost": format_money(row["Incremental cost compared with business as usual (AUD)"]),
                "DALYs averted": format_dalys(row["DALYs averted compared with business as usual"]),
                "Arithmetic ICER": str(row["Arithmetic ICER"]),
                "Classification": str(row["Classification"]),
            }
        )
    return out


def head_to_head_row(comparison: dict[str, Any]) -> dict[str, str] | None:
    rows = comparison.get("scenarios") or []
    if len(rows) < 2:
        return None
    a = rows[0]
    b = rows[1]
    a_metrics = a.get("metrics") or {}
    b_metrics = b.get("metrics") or {}
    delta_cost = _subtract(b_metrics.get("incrementalCost"), a_metrics.get("incrementalCost"))
    delta_dalys = _subtract(b_metrics.get("dalysAverted"), a_metrics.get("dalysAverted"))
    if delta_cost is None or delta_dalys is None:
        return None
    if abs(delta_dalys) < 1e-12:
        if abs(delta_cost) < 1e-12:
            interpretation = "No material difference"
        elif delta_cost < 0:
            interpretation = "Same health effect; Strategy B is less costly"
        else:
            interpretation = "Same health effect; Strategy B is more costly"
        icer = "Not applicable"
    else:
        interpretation = classify_incremental(delta_cost, delta_dalys)
        icer = format_icer(arithmetic_icer(delta_cost, delta_dalys))
    return {
        "Comparison": f"{b.get('label', 'Strategy B')} compared with {a.get('label', 'Strategy A')}",
        "Additional incremental cost": format_money(delta_cost),
        "Additional DALYs averted": format_dalys(delta_dalys),
        "Pairwise ICER": icer,
        "Classification": interpretation,
    }


def _subtract(left: Any, right: Any) -> float | None:
    left_number = number(left)
    right_number = number(right)
    if left_number is None or right_number is None:
        return None
    return left_number - right_number
