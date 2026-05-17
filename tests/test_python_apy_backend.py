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

    def test_unsupported_economics_methods_are_clear(self) -> None:
        methods = [
            lambda: self.backend.default_economics_config(),
            lambda: self.backend.economics_preset_kwab150(),
            lambda: self.backend.run_economics({}, {}),
            lambda: self.backend.run_economics_for_config({}, {}),
        ]
        for method in methods:
            with self.subTest(method=method):
                with self.assertRaisesRegex(
                    NotImplementedError,
                    "does not yet include economics",
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
