from __future__ import annotations

import py_compile
import unittest
from pathlib import Path

from app.icon_arrays import (
    ROUNDING_NOTE,
    build_100_person_visual_data,
    icon_grid_value,
)
from app.run_progress import (
    expected_value_status,
    finalising_status,
    replicate_status,
    status_from_event,
)
from app.epidemiology_inputs import apply_ltbi_state_assumption_update
from engine.apy.config import build_default_config
from engine.apy.runner import run_replicates


ROOT = Path(__file__).resolve().parents[1]


def _ledger(records: list[dict[str, object]]) -> dict[str, object]:
    return {"replicateTotals": records}


class IconArraySummaryTests(unittest.TestCase):
    def test_icon_array_scaling_to_per_100_and_rounding_rule(self) -> None:
        value = icon_grid_value(2.6)

        self.assertTrue(value["available"])
        self.assertEqual(value["filledIcons"], 3)
        self.assertEqual(value["displayPer100"], "2.6 per 100")
        self.assertEqual(value["roundingRule"], ROUNDING_NOTE)

    def test_comparator_intervention_visual_mapping_uses_eligible_denominator(self) -> None:
        rows = build_100_person_visual_data(
            _ledger(
                [
                    {"arm": "comparator", "eventName": "eligible_population", "value": 200},
                    {"arm": "intervention", "eventName": "eligible_population", "value": 200},
                    {"arm": "comparator", "eventName": "active_tb_cases", "value": 8},
                    {"arm": "intervention", "eventName": "active_tb_cases", "value": 5},
                    {"arm": "intervention", "eventName": "tpt_started_total", "value": 40},
                    {"arm": "intervention", "eventName": "tpt_completed_total", "value": 30},
                    {"arm": "intervention", "eventName": "tpt_started_false_positive", "value": 10},
                    {"arm": "intervention", "eventName": "tpt_adr_stop_total", "value": 4},
                    {"arm": "intervention", "eventName": "active_tb_cases_prevented", "value": 3},
                ]
            )
        )
        by_id = {row["outcomeId"]: row for row in rows}

        self.assertEqual(by_id["active_tb_cases"]["comparator"]["displayPer100"], "4.0 per 100")
        self.assertEqual(by_id["active_tb_cases"]["intervention"]["displayPer100"], "2.5 per 100")
        self.assertEqual(by_id["tpt_started"]["intervention"]["displayPer100"], "20.0 per 100")
        self.assertEqual(by_id["tpt_completed"]["intervention"]["displayPer100"], "15.0 per 100")
        self.assertEqual(by_id["false_positive_treatments"]["intervention"]["displayPer100"], "5.0 per 100")
        self.assertEqual(by_id["adverse_events"]["intervention"]["displayPer100"], "2.0 per 100")
        self.assertEqual(by_id["active_tb_cases_prevented"]["intervention"]["displayPer100"], "1.5 per 100")

    def test_intervention_only_outcome_comparator_zero_and_zero_values_do_not_crash(self) -> None:
        rows = build_100_person_visual_data(
            _ledger(
                [
                    {"arm": "comparator", "eventName": "eligible_population", "value": 100},
                    {"arm": "intervention", "eventName": "eligible_population", "value": 100},
                    {"arm": "intervention", "eventName": "active_tb_cases_prevented", "value": 0},
                ]
            ),
            outcomes=["active_tb_cases_prevented"],
        )

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["comparator"]["displayPer100"], "0.0 per 100")
        self.assertEqual(rows[0]["intervention"]["filledIcons"], 0)

    def test_unavailable_values_return_unavailable_without_crashing(self) -> None:
        self.assertEqual(build_100_person_visual_data({}), [])
        unavailable = icon_grid_value(None)
        self.assertFalse(unavailable["available"])
        self.assertEqual(unavailable["displayPer100"], "Unavailable")


class RunProgressTests(unittest.TestCase):
    def test_progress_status_helper_formats_replicates(self) -> None:
        status = replicate_status(1, 5)

        self.assertEqual(status.message, "Simulating cohort: replicate 1 of 5")
        self.assertGreater(status.progress, 0.0)
        self.assertLess(status.progress, 1.0)

    def test_expected_value_and_final_status_formatting(self) -> None:
        self.assertEqual(expected_value_status().message, "Running expected-outcomes analysis")
        self.assertEqual(finalising_status().progress, 1.0)

    def test_status_from_event_maps_replicate_completed(self) -> None:
        status = status_from_event(
            {"stage": "replicate_completed", "replicate": 3, "totalReplicates": 5}
        )

        self.assertEqual(status.message, "Simulating cohort: replicate 3 of 5")

    def test_real_runner_emits_replicate_progress_events(self) -> None:
        cfg = apply_ltbi_state_assumption_update(
            build_default_config(),
            baseline_recent_percent=20.0,
            transition_rate_per_year=0.2,
            source="Reviewed unit-test progress source",
            status="configured_reviewed",
            notes="",
        )
        cfg.update({"N": 12, "nReps": 2, "seed": 12, "screenCoverage": 0.1})
        events = []

        run_replicates(cfg, keep_example_cohort=False, progress_callback=events.append)

        replicate_events = [
            event for event in events
            if event.get("stage") == "replicate_completed"
        ]
        self.assertEqual([event["replicate"] for event in replicate_events], [1, 2])
        self.assertEqual(replicate_events[-1]["totalReplicates"], 2)


class PageIntegrationSmokeTests(unittest.TestCase):
    def test_results_decision_and_run_pages_compile(self) -> None:
        for relative in [
            "pages/2_Run_Model.py",
            "pages/3_Results.py",
            "pages/5_Decision_Analysis.py",
        ]:
            py_compile.compile(str(ROOT / relative), doraise=True)


if __name__ == "__main__":
    unittest.main()
