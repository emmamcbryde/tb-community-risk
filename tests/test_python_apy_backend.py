from __future__ import annotations

import builtins
import json
import sys
import unittest

from adapters.paths import repo_root
from adapters.python_apy_backend import PythonApyBackend


class PythonApyBackendTests(unittest.TestCase):
    def setUp(self) -> None:
        self.backend = PythonApyBackend(repo_root())

    def test_status_reports_python_backend(self) -> None:
        status = self.backend.status()

        self.assertEqual(status["name"], "python_apy")
        self.assertTrue(status["started"])
        self.assertTrue(status["experimental"])
        self.assertFalse(status["matlabRequired"])

    def test_default_config_returns_dict(self) -> None:
        config = self.backend.default_config()

        self.assertIsInstance(config, dict)
        self.assertEqual(config["modelVersion"], "v9")

    def test_validate_default_config_passes(self) -> None:
        report = self.backend.validate_config(self.backend.default_config())

        self.assertTrue(report["isValid"])
        self.assertIn("warnings", report)

    def test_run_scenario_bundle_small_run(self) -> None:
        config = self.backend.default_config()
        config.update({"N": 100, "nReps": 5, "seed": 1})

        bundle = self.backend.run_scenario_bundle(config)

        self.assertIn("metadata", bundle)
        self.assertIn("headline", bundle)
        self.assertIn("technical", bundle)
        self.assertIn("dynamicComparison", bundle["technical"])
        self.assertEqual(bundle["metadata"]["backend"], "python")
        self.assertEqual(
            bundle["technical"]["dynamicComparison"]["source"],
            "doNothing.derived",
        )

    def test_default_economics_config_returns_blank_matlab_shape(self) -> None:
        config = self.backend.default_economics_config()

        self.assertEqual(config["metadata"]["currencyCode"], "")
        self.assertIsNone(config["metadata"]["priceYear"])
        self.assertIsNone(config["costs"]["test"]["IGRA"])
        self.assertIsNone(config["costs"]["regimen"]["x3HP"])
        self.assertIsNone(config["costs"]["falsePositiveIncrementalPerPerson"])

    def test_economics_preset_kwab150_returns_expected_values(self) -> None:
        config = self.backend.economics_preset_kwab150()

        self.assertEqual(config["metadata"]["currencyCode"], "AUD")
        self.assertEqual(config["metadata"]["priceYear"], 2019)
        self.assertEqual(config["metadata"]["locationLabel"], "Australia")
        self.assertEqual(config["metadata"]["programCostBasis"], "total")
        self.assertEqual(config["costs"]["test"]["IGRA"], 113.48)
        self.assertEqual(config["costs"]["test"]["TST"], 116.07)
        self.assertEqual(config["costs"]["regimen"]["x3HP"], 165.5072)
        self.assertEqual(config["costs"]["activeTBDiseasePerCase"], 19079.6)
        self.assertIsNone(config["costs"]["programSetupTotal"])

    def test_validate_economics_config_accepts_default_and_preset(self) -> None:
        default_report = self.backend.validate_economics_config(
            self.backend.default_economics_config()
        )
        preset_report = self.backend.validate_economics_config(
            self.backend.economics_preset_kwab150()
        )

        self.assertTrue(default_report["isValid"])
        self.assertFalse(default_report["hasWarnings"])
        self.assertEqual(default_report["errors"], [])
        self.assertTrue(preset_report["isValid"])
        self.assertEqual(preset_report["errors"], [])

    def test_validate_economics_config_reports_failures(self) -> None:
        config = self.backend.default_economics_config()
        config["metadata"]["priceYear"] = "2019"
        config["metadata"]["currencyCode"] = 123
        config["costs"]["test"]["IGRA"] = -1
        config["costs"]["regimen"] = "bad"

        report = self.backend.validate_economics_config(config)

        self.assertFalse(report["isValid"])
        self.assertFalse(report["hasWarnings"])
        self.assertEqual(
            [issue["field"] for issue in report["errors"]],
            [
                "metadata.currencyCode",
                "metadata.priceYear",
                "costs.test.IGRA",
                "costs.regimen",
            ],
        )

    def test_run_economics_returns_partial_json_payload(self) -> None:
        config = self.backend.default_config()
        config.update({"N": 50, "nReps": 1, "seed": 2})
        results = self.backend.run_scenario(config)

        payload = self.backend.run_economics(
            results,
            self.backend.economics_preset_kwab150(),
        )

        self.assertEqual(payload["status"], "partial")
        self.assertNotIn(
            "does not include economics execution",
            payload["message"].lower(),
        )
        self.assertIn("full", payload["message"].lower())
        self.assertIn("not ported", payload["message"].lower())
        self.assertEqual(payload["source"], "apy_python_partial_economics")
        self.assertEqual(payload["contractVersion"], "apy_economics_partial_v1")
        self.assertEqual(
            payload["metadata"]["legacyIdentifiers"],
            {
                "source": "apy_python_minimal_economics",
                "contractVersion": "apy_economics_minimal_v1",
            },
        )
        self.assertEqual(payload["strategy"], {"testType": "IGRA", "regimen": "3HP"})
        json.dumps(payload)

    def test_run_economics_preserves_partial_status_and_coverage_map(self) -> None:
        config = self.backend.default_config()
        config.update({"N": 50, "nReps": 1, "seed": 2})
        results = self.backend.run_scenario(config)

        payload = self.backend.run_economics(
            results,
            self.backend.economics_preset_kwab150(),
        )

        self.assertEqual(payload["status"], "partial")
        self.assertIn("coverage", payload)
        self.assertEqual(payload["coverage"]["status"], "partial")
        self.assertIn("testingCost", payload["coverage"]["calculatedComponents"])
        self.assertIn("treatmentCost", payload["coverage"]["calculatedComponents"])
        self.assertIn("notCalculated", payload["coverage"])
        self.assertIn("missingInputs", payload["coverage"])
        self.assertIn("messages", payload["coverage"])
        json.dumps(payload["coverage"])

    def test_run_economics_exposes_summary_rows_contract(self) -> None:
        config = self.backend.default_config()
        config.update({"N": 50, "nReps": 1, "seed": 2})
        results = self.backend.run_scenario(config)

        payload = self.backend.run_economics(
            results,
            self.backend.economics_preset_kwab150(),
        )

        rows = payload["summaryRows"]
        self.assertGreater(len(rows), 0)
        json.dumps(rows)

        for row in rows:
            self.assertEqual(
                set(row),
                {
                    "metric",
                    "label",
                    "value",
                    "category",
                    "status",
                    "includedInTotal",
                },
            )

        summary_by_metric = {row["metric"]: row for row in rows}
        self.assertIn("totalImplementedCost", summary_by_metric)
        self.assertAlmostEqual(
            summary_by_metric["totalImplementedCost"]["value"],
            sum(row["value"] for row in rows if row["includedInTotal"] is True),
        )

        self.assertEqual(payload["status"], "partial")
        coverage = payload["coverage"]
        json.dumps(coverage)
        self.assertTrue(
            coverage["unsupportedComponents"]
            or coverage["notCalculated"]
            or coverage["missingInputs"]
        )

    def test_run_economics_uses_app_bundle_dynamic_comparison_for_tb_disease_costs(self) -> None:
        config = self.backend.default_config()
        config.update({"N": 50, "nReps": 1, "seed": 2})
        bundle = self.backend.run_scenario_bundle(config)

        payload = self.backend.run_economics(
            bundle,
            self.backend.economics_preset_kwab150(),
        )

        summary_by_metric = {row["metric"]: row for row in payload["summaryRows"]}

        self.assertEqual(payload["status"], "partial")
        self.assertIn("baselineTBDiseaseCost", payload["costs"])
        self.assertIn("interventionTBDiseaseCost", payload["costs"])
        self.assertIn("baselineTBDiseaseCost", summary_by_metric)
        self.assertIn("interventionTBDiseaseCost", summary_by_metric)
        self.assertEqual(
            summary_by_metric["baselineTBDiseaseCost"]["category"],
            "benefitCost",
        )
        self.assertEqual(
            summary_by_metric["interventionTBDiseaseCost"]["category"],
            "benefitCost",
        )
        self.assertFalse(summary_by_metric["baselineTBDiseaseCost"]["includedInTotal"])
        self.assertFalse(
            summary_by_metric["interventionTBDiseaseCost"]["includedInTotal"]
        )
        self.assertNotIn("results.dynamicComparison", payload["coverage"]["missingInputs"])
        json.dumps(payload)

    def test_run_economics_for_config_runs_python_scenario_then_economics(self) -> None:
        config = self.backend.default_config()
        config.update({"N": 50, "nReps": 1, "seed": 2})

        payload = self.backend.run_economics_for_config(
            config,
            self.backend.economics_preset_kwab150(),
        )

        self.assertEqual(payload["status"], "partial")
        self.assertIn("not ported", payload["message"].lower())
        self.assertGreaterEqual(payload["quantities"]["nScreened"], 0.0)
        self.assertGreaterEqual(payload["costs"]["testingCost"], 0.0)
        json.dumps(payload)

    def test_python_backend_does_not_import_matlab_engine(self) -> None:
        sys.modules.pop("matlab.engine", None)
        original_import = builtins.__import__

        def fail_on_matlab_engine_import(name, *args, **kwargs):
            if name == "matlab.engine":
                raise AssertionError("Python backend attempted to import matlab.engine")
            return original_import(name, *args, **kwargs)

        config = self.backend.default_config()
        config.update({"N": 50, "nReps": 1, "seed": 2})
        economics_config = self.backend.economics_preset_kwab150()

        try:
            builtins.__import__ = fail_on_matlab_engine_import
            payload = self.backend.run_economics_for_config(config, economics_config)
        finally:
            builtins.__import__ = original_import

        self.assertNotIn("matlab.engine", sys.modules)
        self.assertEqual(payload["status"], "partial")


if __name__ == "__main__":
    unittest.main()
