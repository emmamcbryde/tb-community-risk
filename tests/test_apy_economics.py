from __future__ import annotations

from copy import deepcopy
import unittest

from engine.apy.economics import (
    DALE_2019_PRESET_NAME,
    build_default_economics_config,
    build_economics_preset_dale2019_aud,
    build_economics_preset_kwab150,
    update_cost_item_original_values_from_legacy_fields,
    run_health_economics,
    run_health_economics_for_config,
    validate_economics_config,
)
from engine.apy.ltbi_state import enable_development_compatibility_mode
from engine.apy.runner import run_scenario_with_do_nothing


def _reviewed_daly_fixture() -> dict:
    def rec(value):
        return {
            "value": value,
            "source": "Synthetic legacy economics test fixture",
            "status": "configured_reviewed",
            "notes": "Synthetic fixture.",
            "provisional": False,
        }

    return {
        "activeTBDisabilityWeight": rec(0.2),
        "activeTBDurationYears": rec(0.5),
        "tbCaseFatalityRisk": rec(0.1),
        "yllPerTBDeath": rec(10),
        "includeTPTHealthLoss": False,
        "tptHealthLossExclusionStatus": "reviewed_exclusion",
        "tptHealthLossExclusionRationale": "Synthetic reviewed exclusion.",
        "includeADRHealthLoss": False,
        "adrHealthLossExclusionStatus": "reviewed_exclusion",
        "adrHealthLossExclusionRationale": "Synthetic reviewed exclusion.",
        "includePostTBSequelae": False,
        "postTBSequelaeStatus": "reviewed_exclusion",
        "postTBSequelaeExclusionRationale": "Synthetic reviewed exclusion.",
    }


class ApyEconomicsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        config = enable_development_compatibility_mode({"N": 100, "nReps": 5, "seed": 1})
        cls.python_output = run_scenario_with_do_nothing(config)
        cls.econ_config = build_economics_preset_kwab150()
        for item in cls.econ_config["costItems"]:
            item["originalPriceYear"] = "2025-26"
            item["sourceYearStatus"] = "explicit"
            item["conversionStatus"] = "not_converted"
            item["conversionApplied"] = False
            item["warnings"] = []
        cls.econ_config["costs"]["falsePositiveIncrementalPerPerson"] = 10.0
        cls.econ_config["costs"]["programSetupTotal"] = 1000.0
        cls.econ_config["costs"]["programRunningTotal"] = 2000.0
        cls.econ_config["dalyAssumptions"] = _reviewed_daly_fixture()
        cls.econ_config = update_cost_item_original_values_from_legacy_fields(cls.econ_config)
        for item in cls.econ_config["costItems"]:
            if item["originalCost"] not in (None, "", []):
                item["originalPriceYear"] = "2025-26"
                item["sourceYearStatus"] = "explicit"
            if item["costItemId"] == "program_running":
                item["resourceUse"]["costBasis"] = "annual_during_screening_window"
            if item["costItemId"] == "tpt_adr_management":
                item["originalCost"] = 0.0
                item["originalPriceYear"] = "2025-26"
                item["sourceYearStatus"] = "explicit"
        cls.econ = run_health_economics(cls.python_output, cls.econ_config)

    def test_default_economics_config_has_expected_blank_fields(self) -> None:
        config = build_default_economics_config()

        self.assertEqual(config["metadata"]["currencyCode"], "AUD")
        self.assertEqual(config["metadata"]["priceYear"], "2025-26")
        self.assertEqual(config["metadata"]["perspective"], "Australian health-system perspective")
        self.assertIsNone(config["costs"]["test"]["IGRA"])
        self.assertIsNone(config["costs"]["regimen"]["x3HP"])
        self.assertIsNone(config["costs"]["activeTBDiseasePerCase"])

    def test_kwab150_preset_has_expected_aud_target_year_metadata(self) -> None:
        config = build_economics_preset_kwab150()

        self.assertEqual(config["metadata"]["currencyCode"], "AUD")
        self.assertEqual(config["metadata"]["priceYear"], "2025-26")
        self.assertEqual(config["metadata"]["locationLabel"], "Australia")
        self.assertAlmostEqual(config["costs"]["test"]["IGRA"], 113.48)
        self.assertAlmostEqual(config["costs"]["test"]["TST"], 116.07)
        self.assertAlmostEqual(config["costs"]["regimen"]["x3HP"], 165.5072)
        self.assertAlmostEqual(config["costs"]["regimen"]["x4R"], 123.3172)
        self.assertAlmostEqual(config["costs"]["activeTBDiseasePerCase"], 19079.6)
        self.assertEqual(
            config["costItems"][0]["conversionStatus"],
            "not_converted",
        )

    def test_dale_2019_working_preset_has_expected_metadata_and_discounting(self) -> None:
        config = build_economics_preset_dale2019_aud()

        self.assertEqual(config["metadata"]["presetName"], DALE_2019_PRESET_NAME)
        self.assertEqual(config["metadata"]["currencyCode"], "AUD")
        self.assertEqual(config["metadata"]["priceYear"], "2019")
        self.assertEqual(config["metadata"]["targetCurrency"], "AUD")
        self.assertEqual(config["metadata"]["targetPriceYear"], "2019")
        self.assertEqual(config["metadata"]["perspective"], "Australian health-care system")
        self.assertTrue(config["metadata"]["workingDefault"])
        self.assertEqual(config["metadata"]["referenceStatus"], "working_default_not_final_apy_evidence")
        self.assertEqual(config["discounting"]["primaryDisplayedRate"], 0.03)
        self.assertEqual(config["discounting"]["comparisonRate"], 0.0)
        self.assertIsNone(config["threshold"]["value"])

    def test_dale_2019_working_preset_populates_exact_cost_values(self) -> None:
        config = build_economics_preset_dale2019_aud()
        by_id = {item["costItemId"]: item for item in config["costItems"]}

        expected = {
            "test_igra": 113.48,
            "test_tst": 116.07,
            "regimen_3hp": 165.5072,
            "regimen_4r": 123.3172,
            "regimen_3hr": 134.2272,
            "regimen_6h": 187.7508,
            "regimen_9h": 254.8544,
            "active_tb_disease": 19079.60,
        }
        for item_id, value in expected.items():
            with self.subTest(item_id=item_id):
                item = by_id[item_id]
                self.assertAlmostEqual(item["originalCost"], value)
                self.assertEqual(item["originalCurrency"], "AUD")
                self.assertEqual(item["originalPriceYear"], "2019")
                self.assertEqual(item["targetCurrency"], "AUD")
                self.assertEqual(item["targetPriceYear"], "2019")
                self.assertEqual(item["sourceYearStatus"], "explicit")

    def test_dale_2019_working_preset_converted_costs_equal_originals(self) -> None:
        from engine.apy.costing import normalise_cost_table

        config = build_economics_preset_dale2019_aud()
        by_id = {item["costItemId"]: item for item in normalise_cost_table(config["costItems"])}

        for item_id in ["test_igra", "test_tst", "regimen_3hp", "active_tb_disease", "program_setup"]:
            with self.subTest(item_id=item_id):
                item = by_id[item_id]
                self.assertEqual(item["conversionStatus"], "valid")
                self.assertEqual(item["inflationFactor"], 1.0)
                self.assertEqual(item["convertedTargetYearCost"], item["originalCost"])

    def test_dale_2019_adr_management_maps_selected_regimen(self) -> None:
        three_hp = {item["costItemId"]: item for item in build_economics_preset_dale2019_aud("3HP")["costItems"]}
        four_r = {item["costItemId"]: item for item in build_economics_preset_dale2019_aud("4R")["costItems"]}
        three_hr = {item["costItemId"]: item for item in build_economics_preset_dale2019_aud("3HR")["costItems"]}

        self.assertAlmostEqual(three_hp["tpt_adr_management"]["originalCost"], 39.4059)
        self.assertEqual(three_hp["tpt_adr_management"]["pageTableItem"], "csae3HP")
        self.assertAlmostEqual(four_r["tpt_adr_management"]["originalCost"], 23.141)
        self.assertEqual(four_r["tpt_adr_management"]["pageTableItem"], "csae4R")
        self.assertIsNone(three_hr["tpt_adr_management"]["originalCost"])
        self.assertEqual(three_hr["tpt_adr_management"]["sourceYearStatus"], "unknown")

    def test_dale_2019_working_preset_records_exclusions_and_dalys(self) -> None:
        config = build_economics_preset_dale2019_aud()
        rows = {row["assumptionId"]: row for row in config["assumptionEvidenceRegistry"]}

        self.assertEqual(rows["cost.false_positive_incremental"]["inclusionStatus"], "bundled")
        self.assertEqual(rows["cost.program_setup"]["inclusionStatus"], "excluded")
        self.assertIn("not evidence that implementation costs are truly zero", rows["cost.program_setup"]["notes"])
        self.assertEqual(config["dalyAssumptions"]["activeTBDisabilityWeight"]["value"], 0.333)
        self.assertEqual(config["dalyAssumptions"]["activeTBDurationYears"]["value"], 0.5)
        self.assertEqual(config["dalyAssumptions"]["tbCaseFatalityRisk"]["value"], 0.074)
        self.assertEqual(config["dalyAssumptions"]["yllPerTBDeath"]["value"], 20.0)
        self.assertEqual(config["dalyAssumptions"]["tptHealthLossExclusionStatus"], "reviewed_exclusion")
        self.assertTrue(config["dalyAssumptions"]["tptHealthLossExclusionRationale"])
        self.assertTrue(config["dalyAssumptions"]["adrHealthLossExclusionRationale"])
        self.assertTrue(config["dalyAssumptions"]["postTBSequelaeExclusionRationale"])

    def test_valid_default_economics_config_passes_validation(self) -> None:
        report = validate_economics_config(build_default_economics_config())

        self.assertTrue(report["isValid"])
        self.assertTrue(report["hasWarnings"])

    def test_negative_cost_fails_validation(self) -> None:
        config = build_default_economics_config()
        config["costs"]["test"]["IGRA"] = -1

        report = validate_economics_config(config)

        self.assertFalse(report["isValid"])
        self.assertEqual(report["errors"][0]["field"], "costs.test.IGRA")

    def test_invalid_metadata_type_fails_validation(self) -> None:
        config = build_default_economics_config()
        config["metadata"]["currencyCode"] = 123

        report = validate_economics_config(config)

        self.assertFalse(report["isValid"])
        self.assertEqual(report["errors"][0]["field"], "metadata.currencyCode")

    def test_run_health_economics_returns_summary_and_status(self) -> None:
        self.assertIn("summaryRows", self.econ)
        self.assertIn("status", self.econ)
        self.assertTrue(self.econ["status"]["validationReport"]["isValid"])

    def test_testing_cost_equals_screened_times_selected_test_cost(self) -> None:
        expected = (
            self.econ["quantities"]["nScreened"]
            * self.econ_config["costs"]["test"]["IGRA"]
        )

        self.assertAlmostEqual(self.econ["costs"]["testingCost"], expected)

    def test_treatment_cost_equals_started_times_selected_regimen_cost(self) -> None:
        expected = (
            self.econ["quantities"]["nTotalCoursesStarted"]
            * self.econ_config["costs"]["regimen"]["x3HP"]
        )

        self.assertAlmostEqual(self.econ["costs"]["treatmentCost"], expected)

    def test_tb_disease_costs_averted_uses_active_tb_cost(self) -> None:
        expected = (
            self.econ["quantities"]["nPreventedActiveTB"]
            * self.econ_config["costs"]["activeTBDiseasePerCase"]
        )

        self.assertAlmostEqual(self.econ["costs"]["tbDiseaseCostsAverted"], expected)

    def test_net_cost_vs_baseline_uses_program_cost_less_averted_cost(self) -> None:
        expected = (
            self.econ["costs"]["totalProgramCost"]
            - self.econ["costs"]["tbDiseaseCostsAverted"]
        )

        self.assertAlmostEqual(self.econ["costs"]["netCostVsBaseline"], expected)

    def test_missing_false_positive_incremental_cost_does_not_crash(self) -> None:
        config = build_economics_preset_kwab150()
        for item in config["costItems"]:
            item["originalPriceYear"] = "2025-26"
            item["sourceYearStatus"] = "explicit"
            item["conversionStatus"] = "not_converted"
            item["conversionApplied"] = False
            item["warnings"] = []
            if item["costItemId"] == "tpt_adr_management":
                item["originalCost"] = 0.0

        config = update_cost_item_original_values_from_legacy_fields(config)
        config["dalyAssumptions"] = _reviewed_daly_fixture()
        econ = run_health_economics(self.python_output, config)

        self.assertIsNone(econ["costs"]["falsePositiveIncrementalCost"])
        self.assertIn(
            "No false-positive incremental unit cost supplied",
            " ".join(econ["status"]["messages"]),
        )

    def test_missing_required_costs_produce_not_calculated_messages(self) -> None:
        econ = run_health_economics(self.python_output, build_default_economics_config())

        self.assertIn("testingCost", econ["status"]["notCalculated"])
        self.assertIn("treatmentCost", econ["status"]["notCalculated"])
        self.assertIn("tbDiseaseCostsAverted", econ["status"]["notCalculated"])
        self.assertFalse(econ["status"]["isComplete"])

    def test_run_health_economics_for_config_runs_without_matlab(self) -> None:
        config = enable_development_compatibility_mode({"N": 50, "nReps": 2, "seed": 3})

        econ = run_health_economics_for_config(config, self.econ_config)

        self.assertEqual(econ["source"], "event_ledger_health_economics_v2")
        self.assertIn("summaryRows", econ)


if __name__ == "__main__":
    unittest.main()
