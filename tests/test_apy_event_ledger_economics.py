from __future__ import annotations

from io import BytesIO
import math
import unittest

from openpyxl import load_workbook

from app.results_workbook import build_results_workbook
from engine.apy.config import build_default_config
from engine.apy.economics import build_default_economics_config, run_health_economics
from engine.apy.event_ledger_economics import (
    HEALTH_ECONOMICS_CONTRACT_VERSION,
    classify_incremental_result,
    run_event_ledger_health_economics,
)
from engine.apy.expected_value import run_expected_value
from engine.apy.ltbi_state import apply_ltbi_state_edit, enable_development_compatibility_mode
from engine.apy.runner import run_replicates
from engine.apy.results_bundle import build_results_bundle


class ApyEventLedgerEconomicsTests(unittest.TestCase):
    def test_zero_intervention_has_zero_incremental_cost_and_dalys(self) -> None:
        result = run_expected_value(_reviewed_epi({"N": 100, "screenCoverage": 0}))
        econ = run_event_ledger_health_economics(result, _synthetic_econ({"program_setup": 0, "program_running": 0}))
        reps = econ["replicateResults"]
        primary = reps[reps["discountRate"] == 0.03].iloc[0]

        self.assertAlmostEqual(primary["incrementalCost"], 0)
        self.assertAlmostEqual(primary["dalysAverted"], 0)
        self.assertEqual(primary["classification"], "no material difference")

    def test_valid_synthetic_cost_and_daly_calculation(self) -> None:
        result = run_expected_value(_reviewed_epi({"N": 80, "screeningStrategy": "random"}))
        econ = run_event_ledger_health_economics(result, _synthetic_econ())
        annual = econ["annualByArm"]
        year0 = annual[
            (annual["arm"] == "intervention")
            & (annual["discountRate"] == 0.0)
            & (annual["modelYear"] == 0)
        ].iloc[0]

        self.assertEqual(econ["contractVersion"], HEALTH_ECONOMICS_CONTRACT_VERSION)
        self.assertAlmostEqual(year0["screeningTestCost"], year0["screened"] * 10)
        self.assertAlmostEqual(year0["programSetupCost"], 1000)
        self.assertEqual(econ["validation"]["isValid"], True)

    def test_unresolved_cost_prevents_complete_total(self) -> None:
        result = run_expected_value(_reviewed_epi({"N": 80}))
        econ_config = _synthetic_econ()
        for item in econ_config["costItems"]:
            if item["costItemId"] == "test_igra":
                item["originalPriceYear"] = None
                item["sourceYearStatus"] = "unknown"
        econ = run_event_ledger_health_economics(result, econ_config)

        self.assertIn("costItems.test_igra", {item["field"] for item in econ["unresolvedInputs"]})
        self.assertFalse(econ["status"]["isComplete"])
        self.assertIn("ICER", econ["status"]["notCalculated"])

    def test_active_tb_cost_reconciliation_and_no_double_subtraction(self) -> None:
        result = run_expected_value(_reviewed_epi({"N": 120, "screeningStrategy": "random"}))
        econ = run_event_ledger_health_economics(result, _synthetic_econ({"program_setup": 0, "program_running": 0}))
        annual = econ["annualByArm"]
        rows = econ["replicateResults"]
        r = rows[(rows["discountRate"] == 0.0)].iloc[0]
        tb = annual[annual["discountRate"] == 0.0].groupby("arm")["activeTBDiseaseCost"].sum()
        programme = annual[
            (annual["discountRate"] == 0.0)
            & (annual["arm"] == "intervention")
        ][["screeningTestCost", "tptRegimenCost", "falsePositiveIncrementalCost", "programSetupCost", "programRunningCost", "adrManagementCost"]].sum().sum()

        self.assertAlmostEqual(tb["comparator"] - tb["intervention"], (r["activeTBCasesPrevented"] * 1000))
        self.assertAlmostEqual(r["incrementalCost"], programme + tb["intervention"] - tb["comparator"])

    def test_programme_running_basis_and_setup_timing(self) -> None:
        result = run_expected_value(_reviewed_epi({"N": 100, "screeningWindowYears": 2}))
        econ = run_event_ledger_health_economics(
            result,
            _synthetic_econ({"program_running_basis": "annual_during_screening_window", "program_running": 50}),
        )
        annual = econ["annualByArm"]

        setup_rows = annual[(annual["programSetupCost"] > 0) & (annual["discountRate"] == 0.0)]
        running_years = annual[(annual["programRunningCost"] > 0) & (annual["discountRate"] == 0.0)]["modelYear"].unique()
        self.assertEqual(set(setup_rows["modelYear"]), {0})
        self.assertTrue(set(running_years).issubset({0, 1}))

        bad = _synthetic_econ({"program_running": 50, "program_running_basis": ""})
        unresolved = run_event_ledger_health_economics(result, bad)
        self.assertIn("costItems.program_running.resourceUse.costBasis", {item["field"] for item in unresolved["unresolvedInputs"]})

        total = run_event_ledger_health_economics(
            result,
            _synthetic_econ(
                {
                    "program_running": 120,
                    "program_running_basis": "total_over_screening_programme",
                    "program_running_allocation": "proportional_to_screening_volume",
                }
            ),
        )
        running = total["annualByArm"][
            (total["annualByArm"]["arm"] == "intervention")
            & (total["annualByArm"]["discountRate"] == 0.0)
        ]["programRunningCost"].sum()
        self.assertAlmostEqual(running, 120)

    def test_cost_target_year_must_match_analysis(self) -> None:
        result = run_expected_value(_reviewed_epi({"N": 80}))
        config = _synthetic_econ()
        for item in config["costItems"]:
            if item["costItemId"] == "test_igra":
                item["targetPriceYear"] = "2024-25"

        econ = run_event_ledger_health_economics(result, config)

        self.assertIn("costItems.test_igra.targetPriceYear", {item["field"] for item in econ["unresolvedInputs"]})
        self.assertIsNone(econ["costs"]["testingCost"])

    def test_false_positive_and_adr_costs_do_not_create_benefit(self) -> None:
        result = run_expected_value(_reviewed_epi({"N": 100, "testSpecificity": 0}))
        econ = run_event_ledger_health_economics(result, _synthetic_econ({"adr": 7}))
        annual = econ["annualByArm"]
        intervention = annual[(annual["arm"] == "intervention") & (annual["discountRate"] == 0.0)]

        self.assertGreater(intervention["falsePositiveIncrementalCost"].sum(), 0)
        self.assertAlmostEqual(
            intervention["adrManagementCost"].sum(),
            intervention["tpt_adr_stop_total"].sum() * 7,
        )
        self.assertLessEqual(intervention["active_tb_cases_prevented"].sum(), intervention["infection_effectively_treated_total"].sum())

    def test_discounting_known_year_and_no_quantity_change(self) -> None:
        result = run_expected_value(_reviewed_epi({"N": 80}))
        econ = run_event_ledger_health_economics(result, _synthetic_econ())
        annual = econ["annualByArm"]
        y10 = annual[(annual["modelYear"] == 10) & (annual["discountRate"] == 0.03)]
        if y10.empty:
            self.skipTest("Fixture has no year-10 economic events.")
        row = y10.iloc[0]
        self.assertAlmostEqual(row["discountFactor"], 1 / (1.03 ** 10))
        same = annual[
            (annual["modelYear"] == row["modelYear"])
            & (annual["arm"] == row["arm"])
            & (annual["discountRate"] == 0.0)
            & (annual["replicateId"] == row["replicateId"])
        ].iloc[0]
        self.assertEqual(row["active_tb_cases"], same["active_tb_cases"])

    def test_daly_equations_and_treatment_harms(self) -> None:
        result = run_expected_value(_reviewed_epi({"N": 80}))
        econ = run_event_ledger_health_economics(result, _synthetic_econ({"tpt_loss": 0.5}))
        annual = econ["annualByArm"]
        row = annual[(annual["active_tb_cases"] > 0) & (annual["discountRate"] == 0.0)].iloc[0]

        self.assertAlmostEqual(row["activeTBYLD"], row["active_tb_cases"] * 0.2 * 0.5)
        self.assertAlmostEqual(row["activeTBYLL"], row["active_tb_cases"] * 0.1 * 10)
        reps = econ["replicateResults"]
        r = reps[reps["discountRate"] == 0.0].iloc[0]
        self.assertAlmostEqual(r["dalysAverted"], r["comparatorDALYs"] - r["interventionDALYs"])

    def test_icer_ratio_of_means_not_mean_replicate_icer(self) -> None:
        result = run_replicates(_reviewed_epi({"N": 180, "nReps": 5, "seed": 5}), keep_example_cohort=False)
        econ = run_event_ledger_health_economics(result, _synthetic_econ())
        reps = econ["replicateResults"]
        reps = reps[(reps["discountRate"] == 0.03) & reps["replicateICER"].notna()]
        summary = econ["summaries"]
        primary = summary[
            (summary["discountRate"] == 0.03)
            & (summary["metric"] == "primaryICER_ratioOfMeans")
        ].iloc[0]["mean"]

        self.assertAlmostEqual(primary, reps["incrementalCost"].mean() / reps["dalysAverted"].mean())
        if len(reps) > 1 and not math.isclose(reps["replicateICER"].mean(), primary):
            self.assertNotAlmostEqual(primary, reps["replicateICER"].mean())

    def test_dominance_and_nmb(self) -> None:
        self.assertEqual(classify_incremental_result(-1, 1), "dominant")
        self.assertEqual(classify_incremental_result(1, -1), "dominated")
        self.assertEqual(classify_incremental_result(-1, -1), "cost saving with health loss")
        econ = run_event_ledger_health_economics(
            run_expected_value(_reviewed_epi({"N": 80})),
            _synthetic_econ({"threshold": 50000}),
        )
        reps = econ["replicateResults"]
        row = reps[reps["discountRate"] == 0.03].iloc[0]
        self.assertAlmostEqual(row["netMonetaryBenefit"], 50000 * row["dalysAverted"] - row["incrementalCost"])

        no_nmb = run_event_ledger_health_economics(run_expected_value(_reviewed_epi({"N": 80})), _synthetic_econ())
        self.assertIn("threshold.value", {item["field"] for item in no_nmb["unresolvedInputs"]})

    def test_provisional_epidemiology_propagates(self) -> None:
        result = run_expected_value(enable_development_compatibility_mode({"N": 80}))
        econ = run_event_ledger_health_economics(result, _synthetic_econ())

        self.assertTrue(econ["metadata"]["isProvisional"])
        self.assertFalse(econ["metadata"]["conclusionPermitted"])

    def test_existing_public_economics_routes_to_event_ledger(self) -> None:
        bundle = build_results_bundle(run_replicates(_reviewed_epi({"N": 80, "nReps": 2}), keep_example_cohort=False))
        econ = run_health_economics(bundle, _synthetic_econ())

        self.assertEqual(econ["contractVersion"], HEALTH_ECONOMICS_CONTRACT_VERSION)
        self.assertEqual(econ["source"], "event_ledger_health_economics_v2")

    def test_workbook_values_match_authoritative_economics(self) -> None:
        results = run_replicates(_reviewed_epi({"N": 80, "nReps": 2, "seed": 22}), keep_example_cohort=False)
        bundle = build_results_bundle(results)
        econ = run_event_ledger_health_economics(bundle, _synthetic_econ())
        payload = build_results_workbook(
            config=results["interfaceConfig"],
            bundle=bundle,
            backend_status={"name": "python_apy"},
            economics_results=econ,
            economics_config=_synthetic_econ(),
        )
        wb = load_workbook(BytesIO(payload), read_only=True, data_only=True)
        rows = list(wb["Economic_replicates"].iter_rows(values_only=True))
        headers = rows[0]
        first = dict(zip(headers, rows[1]))
        expected = econ["replicateResults"].iloc[0]

        self.assertIn("Economic_annual_by_arm", wb.sheetnames)
        self.assertEqual(first["incrementalCost"], expected["incrementalCost"])
        wb.close()


def _reviewed_epi(overrides: dict) -> dict:
    cfg = apply_ltbi_state_edit(
        build_default_config(),
        baseline_recent_proportion=0.35,
        transition_rate_per_year=0.2,
        source="Reviewed unit-test baseline recent LTBI source",
        status="configured_reviewed",
        notes="Synthetic fixture.",
    )
    cfg.update(overrides)
    return cfg


def _synthetic_econ(overrides: dict | None = None) -> dict:
    overrides = overrides or {}
    config = build_default_economics_config()
    values = {
        "test_igra": 10,
        "test_tst": 8,
        "regimen_3hp": 100,
        "regimen_4r": 90,
        "regimen_3hr": 95,
        "regimen_6h": 80,
        "regimen_9h": 70,
        "active_tb_disease": 1000,
        "false_positive_incremental": 5,
        "program_setup": overrides.get("program_setup", 1000),
        "program_running": overrides.get("program_running", 0),
        "tpt_adr_management": overrides.get("adr", 0),
    }
    for item in config["costItems"]:
        item_id = item["costItemId"]
        if item_id not in values:
            continue
        item["originalCost"] = values[item_id]
        item["originalPriceYear"] = "2025-26"
        item["sourceYearStatus"] = "explicit"
        item["sourceCitation"] = "Synthetic unit-test cost fixture"
        item["resourceUse"]["costBasis"] = item["resourceUse"].get("costBasis") or "per_event"
        if item_id == "program_running":
            item["resourceUse"]["costBasis"] = overrides.get("program_running_basis", "annual_during_screening_window")
            if "program_running_allocation" in overrides:
                item["resourceUse"]["allocationMethod"] = overrides["program_running_allocation"]
    config["dalyAssumptions"] = {
        "activeTBDisabilityWeight": _assumption(0.2),
        "activeTBDurationYears": _assumption(0.5),
        "tbCaseFatalityRisk": _assumption(0.1),
        "yllPerTBDeath": _assumption(10),
        "includeTPTHealthLoss": "tpt_loss" in overrides,
        "tptHealthLossExclusionStatus": "reviewed_exclusion" if "tpt_loss" not in overrides else "configured_reviewed",
        "tptHealthLossExclusionRationale": "Synthetic reviewed exclusion.",
        "dalyLossPerTPTStarted": _assumption(overrides.get("tpt_loss")),
        "dalyLossPerADRStop": _assumption(0.01),
    }
    if "threshold" in overrides:
        config["threshold"].update(
            {
                "value": overrides["threshold"],
                "currency": "AUD",
                "referenceYear": 2025,
                "source": "Synthetic unit-test threshold fixture",
                "status": "configured_reviewed",
            }
        )
    return config


def _assumption(value):
    return {
        "value": value,
        "source": "Synthetic unit-test DALY fixture" if value is not None else "",
        "status": "configured_reviewed" if value is not None else "unresolved",
        "notes": "Synthetic fixture.",
        "provisional": False if value is not None else True,
    }


if __name__ == "__main__":
    unittest.main()
