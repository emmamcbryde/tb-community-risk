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
)

ROOT = Path(__file__).resolve().parents[1]


class SAHealthReferencePackageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.package = build_sa_health_reference_package()

    def test_manifest_is_complete_and_provisional(self) -> None:
        manifest = self.package["manifest"]

        self.assertEqual(manifest["packageId"], SA_HEALTH_REFERENCE_PACKAGE_ID)
        self.assertEqual(manifest["modelType"], "expected_value")
        self.assertFalse(manifest["dynamicTransmissionIncluded"])
        self.assertEqual(manifest["eventLedgerContractVersion"], "ltbi_screening_event_ledger_v3")
        self.assertEqual(manifest["healthEconomicsContractVersion"], "ltbi_health_economics_results_v3")
        self.assertTrue(manifest["configurationHash"])
        self.assertTrue(manifest["economicsConfigurationHash"])
        self.assertTrue(manifest["evidenceRegistryHash"])
        self.assertFalse(manifest["readiness"]["overallClinicianReady"])
        self.assertIn("Recent-LTBI proportion uses explicit 0% compatibility placeholder and is not an estimate.", manifest["interpretationGuardrails"])

    def test_reference_config_matches_sa_health_decisions(self) -> None:
        config = self.package["config"]
        econ = self.package["economicsConfig"]

        self.assertEqual(config["N"], 1500)
        self.assertEqual(config["testType"], "IGRA")
        self.assertEqual(config["regimen"], "3HP")
        self.assertEqual(config["screeningStrategy"], "prevent")
        self.assertEqual(config["screenCoverage"], 0.30)
        self.assertEqual(config["screeningWindowYears"], 3)
        self.assertEqual(config["followUpHorizonYears"], 20)
        self.assertEqual(econ["metadata"]["perspective"], "Australian health-care system")
        self.assertEqual(econ["metadata"]["targetCurrency"], "AUD")
        self.assertEqual(econ["metadata"]["targetPriceYear"], "2019")
        self.assertIsNone(econ["threshold"]["value"])

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
        econ_replicates = self.package["economics"]["replicateResults"]
        primary = econ_replicates[econ_replicates["discountProfile"] == "primary"].iloc[0]
        executive = {row["metric"]: row for row in tables["executive_summary"]}

        self.assertAlmostEqual(executive["People screened"]["value"], 450.0)
        self.assertAlmostEqual(
            executive["Comparator active TB cases"]["value"] - executive["Intervention active TB cases"]["value"],
            executive["Active TB cases averted"]["value"],
        )
        self.assertAlmostEqual(executive["Incremental cost"]["value"], primary["incrementalCost"])
        self.assertAlmostEqual(
            executive["Cost per active TB case averted"]["value"],
            primary["costPerActiveTBCasePrevented"],
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

    def test_health_economics_page_uses_workspace_and_same_ledger_scenarios(self) -> None:
        page = (ROOT / "pages" / "4_Economics.py").read_text(encoding="utf-8")

        self.assertIn("Inputs required for this analysis", page)
        self.assertIn("Compare economic scenarios using current results", page)
        self.assertIn("Changing only economic assumptions does not rerun the epidemiological analysis", page)
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
        result = write_sa_health_reference_package(out_dir)
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
        self.assertTrue((out_dir / "figures" / "screening_treatment_cascade.svg").exists())
        with (out_dir / "tables" / "executive_summary.csv").open("r", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        self.assertTrue(any(row["metric"] == "Analysis status" and "provisional" in row["value"] for row in rows))


if __name__ == "__main__":
    unittest.main()
