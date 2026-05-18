from __future__ import annotations

import sys
import unittest

from adapters.paths import repo_root
from adapters.python_apy_backend import (
    PYTHON_ECONOMICS_UNSUPPORTED,
    PythonApyBackend,
)


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

    def test_full_economics_execution_methods_are_clear(self) -> None:
        methods = [
            lambda: self.backend.run_economics({}, {}),
            lambda: self.backend.run_economics_for_config({}, {}),
        ]
        for method in methods:
            with self.subTest(method=method):
                with self.assertRaisesRegex(
                    NotImplementedError,
                    "does not yet include economics execution",
                ):
                    method()

    def test_python_backend_does_not_import_matlab_engine(self) -> None:
        sys.modules.pop("matlab.engine", None)
        config = self.backend.default_config()
        config.update({"N": 50, "nReps": 1, "seed": 2})

        self.backend.run_scenario_bundle(config)

        self.assertNotIn("matlab.engine", sys.modules)
        self.assertIn("economics", PYTHON_ECONOMICS_UNSUPPORTED)


if __name__ == "__main__":
    unittest.main()
