from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from engine.apy.decision_analysis import run_scenario_comparison
from engine.apy.early_review import run_early_screening_review
from engine.apy.sensitivity import crossings_from_evaluated_grid, run_one_way_sensitivity
from tests.test_apy_event_ledger_economics import _reviewed_epi, _synthetic_econ


def main() -> None:
    base = _reviewed_epi({"N": 160, "screenCoverage": 0.3, "screeningStrategy": "prevent", "ltbiPrevalence": 0.08})
    prevalence = run_scenario_comparison(
        base,
        None,
        [
            {"scenarioId": "low", "changes": {"ltbiPrevalence": 0.02}},
            {"scenarioId": "reference", "changes": {}},
            {"scenarioId": "high", "changes": {"ltbiPrevalence": 0.16}},
        ],
        model_type="expected_value",
    )
    scenarios = prevalence["scenarios"]
    hazards = {(row["calibration"]["lambdaEarly"], row["calibration"]["lambdaLate"]) for row in scenarios}
    comparator_tb = [row["metrics"]["comparator_active_tb"] for row in scenarios]
    assert len(hazards) == 1, hazards
    assert comparator_tb[0] < comparator_tb[2], comparator_tb

    zero = run_scenario_comparison(
        base,
        _synthetic_econ({"threshold": 50000}),
        [{"scenarioId": "zero", "changes": {"ltbiPrevalence": 0.0, "testSpecificity": 0.8}}],
        model_type="expected_value",
    )
    zero_summary = zero["scenarioSummaries"][0]
    assert zero_summary["false_positive"] > 0
    assert zero_summary["tpt_started_total"] > 0
    assert zero_summary["active_tb_cases_prevented"] == 0

    paired = run_scenario_comparison(
        base,
        None,
        [
            {"scenarioId": "igra", "changes": {"test": "IGRA"}},
            {"scenarioId": "tst", "changes": {"test": "TST"}},
        ],
        model_type="agent_based",
        n_reps=3,
        master_seed=123,
    )
    assert paired["validation"]["isValid"]
    assert all(row["fingerprintMatch"] for row in paired["pairedReplicateComparisons"])

    spec = {
        "parameterId": "cost.test_igra.synthetic",
        "adapter": "testIGRACost",
        "label": "Synthetic IGRA target-year cost",
        "baseValue": 10,
        "lowValue": 1,
        "highValue": 50,
        "monetaryValueBasis": "target_year_converted_cost",
        "source": "Synthetic validation fixture",
        "reviewStatus": "configured_reviewed",
        "provisional": False,
    }
    sens = run_one_way_sensitivity(base, _synthetic_econ({"test_igra": 10}), [spec], ["incrementalCost"])
    row = sens["results"][0]
    assert row["highOutcome"] > row["lowOutcome"]

    early_low = run_early_screening_review(
        base,
        None,
        {
            "screenedToDate": 20,
            "testPositiveToDate": 0,
            "plannedTotalScreened": 60,
            "eligiblePopulation": 160,
            "prior": {"type": "beta", "mean": 0.08, "effectiveSampleSize": 20, "source": "Synthetic validation prior"},
            "prevalenceGrid": [0.02, 0.08, 0.14],
        },
    )
    early_high = run_early_screening_review(
        base,
        None,
        {
            "screenedToDate": 20,
            "testPositiveToDate": 8,
            "plannedTotalScreened": 60,
            "eligiblePopulation": 160,
            "prior": {"type": "beta", "mean": 0.08, "effectiveSampleSize": 20, "source": "Synthetic validation prior"},
            "prevalenceGrid": [0.02, 0.08, 0.14],
        },
    )
    assert early_low["posterior"]["summary"]["mean"] < early_low["prior"]["summary"]["mean"]
    assert early_high["posterior"]["summary"]["mean"] > early_high["prior"]["summary"]["mean"]

    gap_grid = [
        {"parameterValue": 0, "difference": -1, "complete": True},
        {"parameterValue": 1, "difference": None, "complete": False},
        {"parameterValue": 2, "difference": 1, "complete": True},
    ]
    assert crossings_from_evaluated_grid(gap_grid) == []
    print("APY M5 hardening validation passed")
    print(
        {
            "prevalenceHazards": list(hazards)[0],
            "comparatorTB": comparator_tb,
            "zeroFalsePositive": zero_summary["false_positive"],
            "pairedReplicates": len(paired["pairedReplicateComparisons"]),
        }
    )


if __name__ == "__main__":
    main()
