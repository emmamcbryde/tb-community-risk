from __future__ import annotations

from copy import deepcopy
import csv
import json
from pathlib import Path
import unittest

from app.health_economics_inputs import editable_assumption_rows
from engine.apy.event_ledger_economics import run_event_ledger_health_economics
from engine.apy.sa_health_reference_package import (
    PROGRAMME_NOT_COSTED_LABEL,
    SA_HEALTH_REFERENCE_PACKAGE_ID,
    build_sa_health_reference_package,
    build_same_ledger_economic_scenario_comparison,
    write_sa_health_reference_package,
    _replicate_row,
)

ROOT = Path(__file__).resolve().parents[1]


class SAHealthReferencePackageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.package = build_sa_health_reference_package(
            n_reps=50,
            include_technical_recent_remote=True,
        )

    def test_manifest_is_complete_and_provisional(self) -> None:
        manifest = self.package["manifest"]

        self.assertEqual(manifest["packageId"], SA_HEALTH_REFERENCE_PACKAGE_ID)
        self.assertEqual(manifest["modelType"], "agent_based")
        self.assertEqual(manifest["epidemiologicalAnchor"]["naturalHistorySemantics"], "matlab_v9_implicit_early_late")
        self.assertFalse(manifest["dynamicTransmissionIncluded"])
        self.assertEqual(manifest["eventLedgerContractVersion"], "ltbi_screening_event_ledger_v3")
        self.assertEqual(manifest["healthEconomicsContractVersion"], "ltbi_health_economics_results_v3")
        self.assertTrue(manifest["configurationHash"])
        self.assertTrue(manifest["economicsConfigurationHash"])
        self.assertTrue(manifest["evidenceRegistryHash"])
        self.assertFalse(manifest["readiness"]["overallClinicianReady"])
        self.assertIn("The inherited 10/770 active-TB calibration target remains unresolved and must not be described as validated incident progression from LTBI.", manifest["interpretationGuardrails"])

    def test_reference_config_matches_sa_health_decisions(self) -> None:
        config = self.package["config"]
        econ = self.package["economicsConfig"]

        self.assertEqual(config["N"], 1500)
        self.assertEqual(config["testType"], "IGRA")
        self.assertEqual(config["regimen"], "3HP")
        self.assertEqual(config["screeningStrategy"], "prevent")
        self.assertEqual(config["screenCoverage"], 0.30)
        self.assertEqual(config["screeningWindowYears"], 2)
        self.assertEqual(config["followUpHorizonYears"], 20)
        self.assertEqual(config["nReps"], 50)
        self.assertEqual(config["seed"], 1)
        self.assertEqual(config["naturalHistorySemantics"], "matlab_v9_implicit_early_late")
        self.assertEqual(econ["metadata"]["perspective"], "Australian health-care system")
        self.assertEqual(econ["metadata"]["targetCurrency"], "AUD")
        self.assertEqual(econ["metadata"]["targetPriceYear"], "2019")
        self.assertIsNone(econ["threshold"]["value"])

    def test_primary_anchor_uses_matlab_v9_compatible_hazards(self) -> None:
        calibration = self.package["epidemiology"]["calibration"]

        self.assertAlmostEqual(calibration["lambdaEarly"], 0.0037975986627, places=10)
        self.assertAlmostEqual(calibration["lambdaLate"], 0.0007595197325, places=10)
        self.assertNotAlmostEqual(calibration["lambdaEarly"], 0.0189879933116, places=5)

    def test_streamlit_cost_workspace_contains_requested_editable_costs(self) -> None:
        rows = editable_assumption_rows(economics_config=self.package["economicsConfig"])
        by_id = {row["assumptionId"]: row for row in rows}
        required = {
            "cost.test_igra",
            "cost.test_tst",
            "cost.regimen_3hp",
            "cost.regimen_4r",
            "cost.regimen_3hr",
            "cost.regimen_6h",
            "cost.regimen_9h",
            "cost.active_tb_disease",
            "cost.tpt_adr_management",
            "cost.false_positive_incremental",
            "cost.return_for_results",
            "cost.clinical_review",
            "cost.active_tb_exclusion_workup",
            "cost.program_setup",
            "cost.program_running",
            "cost.travel_outreach_staff_support",
        }

        self.assertTrue(required.issubset(by_id))
        for assumption_id in required:
            row = by_id[assumption_id]
            self.assertIn("currentValue", row)
            self.assertTrue(row.get("unit") or row.get("costBasis"))
            self.assertIn("originalCurrency", row)
            self.assertIn("originalPriceYear", row)
            self.assertIn("sourceCitation", row)
            self.assertIn("reviewStatus", row)
            self.assertIn("inclusionStatus", row)

    def test_report_tables_reconcile_to_authoritative_outputs(self) -> None:
        tables = self.package["tables"]
        primary = _replicate_row(self.package["economics"], "primary")
        executive = {row["metric"]: row for row in tables["executive_summary"]}

        self.assertAlmostEqual(executive["People screened"]["value"], 450.0)
        self.assertAlmostEqual(
            executive["Comparator active TB cases"]["value"] - executive["Intervention active TB cases"]["value"],
            executive["Active TB cases averted"]["value"],
        )
        self.assertAlmostEqual(executive["Incremental cost"]["value"], primary["incrementalCost"])
        self.assertAlmostEqual(
            executive["Net health-system cost or saving per active TB case averted"]["value"],
            primary["costPerActiveTBCasePrevented"],
        )
        annual = tables["annual_budget_impact"]
        self.assertAlmostEqual(
            annual[-1]["cumulativeIncrementalCostDiscounted"],
            primary["incrementalCost"],
        )

        category_total = next(row for row in tables["cost_categories"] if row["categoryId"] == "total")
        self.assertAlmostEqual(category_total["incrementalDiscountedCost"], primary["incrementalCost"])

    def test_unresolved_programme_costs_are_distinguishable_from_zero(self) -> None:
        categories = {
            row["categoryId"]: row
            for row in self.package["tables"]["cost_categories"]
        }

        for category_id in ["program_setup", "program_running", "travel_outreach_staff"]:
            self.assertEqual(categories[category_id]["interpretation"], PROGRAMME_NOT_COSTED_LABEL)
            self.assertIsNotNone(categories[category_id]["interventionDiscountedCost"])

    def test_missing_wtp_keeps_nmb_unavailable(self) -> None:
        economics = self.package["economics"]
        primary = economics["replicateResults"][economics["replicateResults"]["discountProfile"] == "primary"].iloc[0]
        unresolved = {item["field"] for item in economics["unresolvedInputs"]}

        self.assertIn("threshold.value", unresolved)
        self.assertIsNone(primary["netMonetaryBenefit"])
        self.assertNotIn("probabilityPositiveNMB_fixedParameterSimulation", set(economics["summaries"]["metric"]))

    def test_economic_scenarios_use_same_event_ledger_and_change_costs_only(self) -> None:
        comparison = build_same_ledger_economic_scenario_comparison(
            self.package["epidemiology"],
            self.package["economicsConfig"],
        )
        self.assertFalse(comparison["epidemiologyRerunRequired"])
        scenarios = {row["scenarioId"]: row for row in comparison["rows"]}
        primary = scenarios["primary_working_reference"]
        setup = scenarios["illustrative_500k_setup"]
        bundled = scenarios["pathway_components_bundled"]
        higher = scenarios["higher_burden_post_tb_daly"]

        self.assertTrue(all(row["usesSameEventLedger"] for row in scenarios.values()))
        self.assertAlmostEqual(setup["incrementalCost"] - primary["incrementalCost"], 500000.0)
        self.assertLess(bundled["interventionCost"], primary["interventionCost"])
        self.assertGreater(higher["dalysAverted"], primary["dalysAverted"])

    def test_distributional_agreement_with_stored_matlab_reference(self) -> None:
        rows = {row["metric"]: row for row in self.package["tables"]["matlab_reference_validation"]}

        for metric in [
            "nScreened",
            "nFalsePositiveTests",
            "nStartTPT",
            "nCompleteTPT",
            "nADRstop",
            "nCuredInfection",
            "nPreventedActiveTB",
            "nActiveBy20y",
        ]:
            self.assertTrue(rows[metric]["withinMatlabInterval"], metric)

    def test_technical_recent_remote_scenario_is_not_primary_anchor(self) -> None:
        technical = self.package["technicalRecentRemoteScenario"]

        self.assertFalse(technical["primarySAHealthAnchor"])
        self.assertFalse(technical["suitableAsPrimaryEstimate"])
        self.assertEqual(technical["packageId"], "sa_health_apy_technical_recent_remote_igra_3hp_prevent_30pct")
        self.assertNotEqual(
            technical["epidemiology"]["eventLedger"]["metadata"]["modelType"],
            self.package["eventLedger"]["metadata"]["modelType"],
        )

    def test_gross_and_net_ratio_tables_are_separate(self) -> None:
        gross = {row["ratioId"]: row for row in self.package["tables"]["gross_delivery_ratios"]}
        net = {row["ratioId"]: row for row in self.package["tables"]["net_health_system_ratios"]}

        self.assertGreater(gross["person_screened"]["value"], 0)
        self.assertIn("before active-TB care offsets", gross["person_screened"]["interpretation"])
        self.assertIn("after active-TB care offsets", net["person_screened"]["interpretation"])
        self.assertNotEqual(gross["person_screened"]["value"], net["person_screened"]["value"])

    def test_health_economics_page_uses_workspace_and_same_ledger_scenarios(self) -> None:
        page = (ROOT / "pages" / "4_Economics.py").read_text(encoding="utf-8")

        self.assertIn("Inputs required for this analysis", page)
        self.assertIn("Compare economic scenarios using current results", page)
        self.assertIn("Changing only economic assumptions does not rerun the epidemiological analysis", page)
        self.assertIn("frozen APY stochastic compatibility reference", page)
        self.assertIn("Gross delivery expenditure ratios are before active-TB care offsets", page)
        self.assertIn("build_same_ledger_economic_scenario_comparison", page)

    def test_economic_changes_do_not_alter_event_counts(self) -> None:
        ledger = self.package["epidemiology"]
        base_econ = self.package["economicsConfig"]
        edited_econ = deepcopy(base_econ)
        for item in edited_econ["costItems"]:
            if item["costItemId"] == "test_igra":
                item["originalCost"] = item["originalCost"] + 10
        base = run_event_ledger_health_economics(ledger, base_econ)
        edited = run_event_ledger_health_economics(ledger, edited_econ)

        self.assertEqual(
            base["replicateResults"].iloc[0]["activeTBCasesPrevented"],
            edited["replicateResults"].iloc[0]["activeTBCasesPrevented"],
        )
        self.assertNotEqual(
            base["replicateResults"].iloc[0]["incrementalCost"],
            edited["replicateResults"].iloc[0]["incrementalCost"],
        )

    def test_exported_and_reloaded_economics_config_reproduces_results(self) -> None:
        out_dir = ROOT / "outputs" / "test_tmp" / "sa_health_reference_package_test"
        result = write_sa_health_reference_package(out_dir, n_reps=10)
        econ_path = out_dir / "economics_config.json"
        loaded_econ = json.loads(econ_path.read_text(encoding="utf-8"))
        rerun = run_event_ledger_health_economics(result["package"]["epidemiology"], loaded_econ)
        original = result["package"]["economics"]["replicateResults"]
        reloaded = rerun["replicateResults"]

        self.assertAlmostEqual(
            original[original["discountProfile"] == "primary"].iloc[0]["incrementalCost"],
            reloaded[reloaded["discountProfile"] == "primary"].iloc[0]["incrementalCost"],
        )
        self.assertTrue((out_dir / "sa_health_apy_working_reference_outputs.xlsx").exists())
        self.assertTrue((out_dir / "technical_recent_remote_scenario.json").exists())
        self.assertTrue((out_dir / "figures" / "screening_treatment_cascade.svg").exists())
        with (out_dir / "tables" / "executive_summary.csv").open("r", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        self.assertTrue(any(row["metric"] == "Analysis status" and "provisional" in row["value"] for row in rows))


if __name__ == "__main__":
    unittest.main()
