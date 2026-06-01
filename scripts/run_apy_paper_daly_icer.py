from __future__ import annotations

import argparse
from copy import deepcopy
from pathlib import Path
import sys
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from engine.apy.config import build_default_config
from engine.apy.economics import (
    build_economics_preset_kwab150,
    run_health_economics,
)
from engine.apy.health_economics import (
    build_default_health_economic_assumptions,
    calculate_health_outcomes,
    calculate_icers,
)
from engine.apy.runner import run_scenario_with_do_nothing


DEFAULT_STRATEGIES = ["random", "ltbi", "cure", "prevent"]
OUTPUT_FILES = {
    "strategy": "apy_daly_icer_strategy_summary.csv",
    "health": "apy_daly_icer_health_outcomes.csv",
    "costs": "apy_daly_icer_costs.csv",
    "sensitivity": "apy_daly_icer_sensitivity.csv",
    "notes": "apy_daly_icer_notes.md",
}


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    n = 100 if args.quick else args.n
    n_reps = 5 if args.quick else args.n_reps
    strategies = parse_strategies(args.strategies)

    strategy_rows: list[dict[str, Any]] = []
    health_rows: list[dict[str, Any]] = []
    cost_rows: list[dict[str, Any]] = []
    sensitivity_rows: list[dict[str, Any]] = []

    for strategy in strategies:
        scenario = run_one_strategy(
            strategy=strategy,
            n=n,
            n_reps=n_reps,
            seed=args.seed,
            screen_coverage=args.screen_coverage,
            include_daly=args.include_daly,
            include_qaly=args.include_qaly,
        )
        strategy_rows.append(scenario["strategyRow"])
        health_rows.extend(scenario["healthRows"])
        cost_rows.extend(scenario["costRows"])
        sensitivity_rows.extend(scenario["sensitivityRows"])

    strategy_df = pd.DataFrame(strategy_rows)
    health_df = pd.DataFrame(health_rows)
    cost_df = pd.DataFrame(cost_rows)
    sensitivity_df = pd.DataFrame(sensitivity_rows)

    strategy_df.to_csv(output_dir / OUTPUT_FILES["strategy"], index=False)
    health_df.to_csv(output_dir / OUTPUT_FILES["health"], index=False)
    cost_df.to_csv(output_dir / OUTPUT_FILES["costs"], index=False)
    sensitivity_df.to_csv(output_dir / OUTPUT_FILES["sensitivity"], index=False)
    write_notes(output_dir / OUTPUT_FILES["notes"], strategy_df, sensitivity_df, args)

    print(f"Wrote APY DALY/QALY/ICER outputs to {output_dir}")
    print(strategy_df.head().to_string(index=False))
    return 0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run offline APY paper DALY/QALY/ICER economics outputs.",
    )
    parser.add_argument("--output-dir", default="outputs/apy_paper_daly_icer")
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--n", type=int, default=1500)
    parser.add_argument("--n-reps", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--screen-coverage", type=float, default=0.30)
    parser.add_argument("--strategies", default=",".join(DEFAULT_STRATEGIES))
    parser.add_argument("--include-qaly", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--include-daly", action=argparse.BooleanOptionalAction, default=True)
    return parser.parse_args(argv)


def parse_strategies(raw: str) -> list[str]:
    values = [item.strip().lower() for item in raw.split(",") if item.strip()]
    return values or DEFAULT_STRATEGIES


def run_one_strategy(
    strategy: str,
    n: int,
    n_reps: int,
    seed: int,
    screen_coverage: float,
    include_daly: bool = True,
    include_qaly: bool = True,
) -> dict[str, Any]:
    config = build_default_config()
    config.update(
        {
            "scenarioLabel": f"APY paper DALY ICER - {strategy}",
            "N": int(n),
            "nReps": int(n_reps),
            "seed": int(seed),
            "screenWindow": 2,
            "followHorizon": 20,
            "screenCoverage": float(screen_coverage),
            "screeningStrategy": strategy,
            "testType": "IGRA",
            "regimen": "3HP",
            "pStartTPT": 0.85,
            "regimenPComplete": 0.80,
            "regimenADRstop": 0.05,
            "regimenEffFull": 0.85,
        }
    )
    model_out = run_scenario_with_do_nothing(config)
    econ_config = build_economics_preset_kwab150()
    econ_config["costs"]["falsePositiveIncrementalPerPerson"] = 0.0
    econ_config["costs"]["programSetupTotal"] = 0.0
    econ_config["costs"]["programRunningTotal"] = 0.0
    economics = run_health_economics(model_out, econ_config)
    assumptions = build_default_health_economic_assumptions()
    health = calculate_health_outcomes(model_out, assumptions)
    icer = calculate_icers(economics, health, assumptions)

    scenario_id = f"apy_paper_{strategy}_igra_3hp"
    strategy_row = build_strategy_row(
        scenario_id,
        strategy,
        model_out,
        economics,
        health,
        icer,
        include_daly,
        include_qaly,
    )
    sensitivity_rows = build_sensitivity_rows(
        scenario_id,
        strategy,
        model_out,
        economics,
        include_daly,
        include_qaly,
    )
    return {
        "strategyRow": strategy_row,
        "healthRows": flatten_section_rows(scenario_id, strategy, "health", health["healthOutcomes"]),
        "costRows": flatten_section_rows(scenario_id, strategy, "costs", economics["costs"]),
        "sensitivityRows": sensitivity_rows,
    }


def build_strategy_row(
    scenario_id: str,
    strategy: str,
    model_out: dict[str, Any],
    economics: dict[str, Any],
    health: dict[str, Any],
    icer: dict[str, Any],
    include_daly: bool,
    include_qaly: bool,
) -> dict[str, Any]:
    results = model_out["results"]
    config = results["interfaceConfig"]
    dynamic = model_out["bundle"]["technical"]["dynamicComparison"]
    costs = economics["costs"]
    health_values = health["healthOutcomes"]
    icers = icer["icers"]
    nmb = icer["nmb"]
    row = {
        "scenario_id": scenario_id,
        "screeningStrategy": strategy,
        "testType": config["testType"],
        "regimen": config["regimen"],
        "N": config["N"],
        "nReps": config["nReps"],
        "screenCoverage": config["screenCoverage"],
        "screenWindow": config["screenWindow"],
        "followHorizon": config["followHorizon"],
        "nScreened_median": metric(results, "nScreened"),
        "nTotalCoursesStarted_median": metric(results, "nTotalCoursesStarted"),
        "nTotalCoursesCompleted_median": metric(results, "nTotalCoursesCompleted"),
        "nFalsePositiveTreated_median": metric(results, "nFalsePositiveTreated"),
        "nCuredInfection_median": metric(results, "nCuredInfection"),
        "nPreventedActiveTB_median": metric(results, "nPreventedActiveTB"),
        "nActiveBy20y_DoNothing_median": dynamic.get("cumulative_baseline_active_tb_cases"),
        "nActiveBy20y_AfterStrategy_median": dynamic.get("cumulative_intervention_active_tb_cases"),
        "relativeReduction20y_median": dynamic.get("relative_reduction_cumulative_active_tb_cases"),
        "testingCost": costs.get("testingCost"),
        "treatmentCost": costs.get("treatmentCost"),
        "tbDiseaseCostsAverted": costs.get("tbDiseaseCostsAverted"),
        "totalProgramCost": costs.get("totalProgramCost"),
        "netCostVsBaseline": costs.get("netCostVsBaseline"),
        "tbCasesPrevented": health_values.get("tbCasesPrevented") if include_daly or include_qaly else None,
        "tbDeathsPrevented": health_values.get("tbDeathsPrevented") if include_daly else None,
        "yldAverted": health_values.get("yldAverted") if include_daly else None,
        "yllAverted": health_values.get("yllAverted") if include_daly else None,
        "dalysAverted": health_values.get("dalysAverted") if include_daly else None,
        "qalysGained": health_values.get("qalysGained") if include_qaly else None,
        "costPerTBCasePrevented": icers.get("costPerTBCasePrevented"),
        "costPerDALYAverted": icers.get("costPerDALYAverted") if include_daly else None,
        "costPerQALYGained": icers.get("costPerQALYGained") if include_qaly else None,
        "nmbDALY_45000": nmb.get("netMonetaryBenefitDALY_low") if include_daly else None,
        "nmbDALY_75000": nmb.get("netMonetaryBenefitDALY_high") if include_daly else None,
        "nmbQALY_45000": nmb.get("netMonetaryBenefitQALY_low") if include_qaly else None,
        "nmbQALY_75000": nmb.get("netMonetaryBenefitQALY_high") if include_qaly else None,
    }
    row["dominance_status"] = dominance_status(row)
    row["notes"] = "Direct person-level prevented cases; no downstream transmission effects."
    return row


def build_sensitivity_rows(
    scenario_id: str,
    strategy: str,
    model_out: dict[str, Any],
    economics: dict[str, Any],
    include_daly: bool,
    include_qaly: bool,
) -> list[dict[str, Any]]:
    scenarios = [
        ("no_ltbi_treatment_utility_decrement", {}),
        (
            "ltbi_treatment_decrement_0_0133",
            {"qaly": {"ltbiTreatmentDecrement": 0.0133}},
        ),
        ("yld_only_no_yll", {"daly": {"includeYLL": False}}),
        ("conservative_yll_10", {"daly": {"yllPerTBDeath": 10.0}}),
        ("higher_yll_30", {"daly": {"yllPerTBDeath": 30.0}}),
        ("lower_tb_case_fatality", {"daly": {"tbCaseFatalityRisk": 0.037}}),
        ("higher_tb_case_fatality", {"daly": {"tbCaseFatalityRisk": 0.111}}),
    ]
    rows = []
    for sensitivity_name, overrides in scenarios:
        assumptions = build_default_health_economic_assumptions()
        deep_update(assumptions, overrides)
        health = calculate_health_outcomes(model_out, assumptions)
        icer = calculate_icers(economics, health, assumptions)
        rows.append(
            {
                "scenario_id": scenario_id,
                "screeningStrategy": strategy,
                "sensitivity": sensitivity_name,
                "tbCasesPrevented": health["healthOutcomes"].get("tbCasesPrevented"),
                "dalysAverted": health["healthOutcomes"].get("dalysAverted") if include_daly else None,
                "qalysGained": health["healthOutcomes"].get("qalysGained") if include_qaly else None,
                "costPerDALYAverted": icer["icers"].get("costPerDALYAverted") if include_daly else None,
                "costPerQALYGained": icer["icers"].get("costPerQALYGained") if include_qaly else None,
                "nmbDALY_45000": icer["nmb"].get("netMonetaryBenefitDALY_low") if include_daly else None,
                "nmbDALY_75000": icer["nmb"].get("netMonetaryBenefitDALY_high") if include_daly else None,
                "nmbQALY_45000": icer["nmb"].get("netMonetaryBenefitQALY_low") if include_qaly else None,
                "nmbQALY_75000": icer["nmb"].get("netMonetaryBenefitQALY_high") if include_qaly else None,
                "notes": "Sensitivity assumption; not a final policy setting.",
            }
        )
    return rows


def metric(results: dict[str, Any], metric_name: str) -> float | None:
    rows = results["summary"]
    matched = rows.loc[rows["Metric"] == metric_name, "Median"]
    if matched.empty:
        return None
    return matched.iloc[0]


def flatten_section_rows(
    scenario_id: str,
    strategy: str,
    section: str,
    values: dict[str, Any],
) -> list[dict[str, Any]]:
    return [
        {
            "scenario_id": scenario_id,
            "screeningStrategy": strategy,
            "section": section,
            "metric": key,
            "value": value,
        }
        for key, value in values.items()
        if isinstance(value, (int, float)) or value is None
    ]


def dominance_status(row: dict[str, Any]) -> str:
    net_cost = row.get("netCostVsBaseline")
    cases = row.get("tbCasesPrevented")
    if net_cost is None or cases is None:
        return "not assessed"
    if net_cost < 0 and cases > 0:
        return "dominant_cost_saving"
    if cases <= 0:
        return "no_health_gain"
    return "incremental_cost"


def deep_update(target: dict[str, Any], updates: dict[str, Any]) -> None:
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            deep_update(target[key], value)
        else:
            target[key] = value


def write_notes(
    path: Path,
    strategy_df: pd.DataFrame,
    sensitivity_df: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    best_daly = best_strategy(strategy_df, "costPerDALYAverted")
    best_qaly = best_strategy(strategy_df, "costPerQALYGained")
    lines = [
        "# APY Paper DALY/QALY/ICER Outputs",
        "",
        "Purpose: offline health-economic outputs for the APY paper/report.",
        "",
        "Model basis: Python APY closed-cohort, person-level prevention model using direct prevented TB cases. Downstream transmission effects are not included.",
        "",
        "Costing: Australian health-care system perspective, 2019 AUD, KWAB150/Dale-informed unit costs, 3% discounting assumption recorded for interpretation.",
        "",
        "Program setup, program running, and false-positive incremental costs are set to zero in this offline paper script because the KWAB150 preset does not specify them. This should be sensitivity-tested if programme implementation costs are added.",
        "",
        "DALY assumptions: active TB disability weight 0.333, active TB duration 0.5 years, configurable TB case fatality risk and YLL per TB death. YLL assumptions require explicit review.",
        "",
        "QALY assumptions: healthy utility 0.8733, active TB utility 0.8182, SAE utility 0.8685, and LTBI treatment decrement sensitivity based on Dale-informed assumptions.",
        "",
        "ICER formulas: incremental cost divided by DALYs averted, QALYs gained, or active TB cases prevented. Net monetary benefit equals health gain times WTP threshold minus incremental cost.",
        "",
        "Limitations: not an online tool output; not a full transmission model; not a final policy gate; mortality/YLL assumptions are sensitivity inputs.",
        "",
        f"Lowest cost per DALY averted among calculable base-case strategies: {best_daly}.",
        "",
        f"Lowest cost per QALY gained among calculable base-case strategies: {best_qaly}.",
        "",
        "The APY model supports sequencing and prioritisation, not exclusion from care.",
        "",
        f"Run settings: quick={args.quick}, N={100 if args.quick else args.n}, nReps={5 if args.quick else args.n_reps}, seed={args.seed}, screenCoverage={args.screen_coverage}.",
        "",
        f"Sensitivity rows generated: {len(sensitivity_df)}.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def best_strategy(df: pd.DataFrame, metric_name: str) -> str:
    if metric_name not in df.columns:
        return "not calculable"
    values = df.dropna(subset=[metric_name])
    if values.empty:
        return "not calculable"
    return str(values.loc[values[metric_name].idxmin(), "screeningStrategy"])


if __name__ == "__main__":
    raise SystemExit(main())
