from __future__ import annotations

import unittest

from engine.apy.config import build_default_config, normalise_config, strip_empty_fields
from engine.apy.validation import collect_validation_issues, validate_config


class ApyConfigValidationTests(unittest.TestCase):
    def test_default_config_has_expected_key_values(self) -> None:
        cfg = build_default_config()

        self.assertEqual(cfg["configVersion"], "apy_v9_config_v1")
        self.assertEqual(cfg["modelVersion"], "v9")
        self.assertEqual(cfg["scenarioLabel"], "APY default scenario")
        self.assertTrue(cfg["useDefaults"])
        self.assertEqual(cfg["N"], 1500)
        self.assertEqual(cfg["nReps"], 2000)
        self.assertEqual(cfg["seed"], 1)
        self.assertEqual(cfg["screenWindow"], 2)
        self.assertEqual(cfg["followHorizon"], 20)
        self.assertEqual(cfg["screenCoverage"], 0.30)
        self.assertEqual(cfg["screeningStrategy"], "prevent")
        self.assertEqual(cfg["targetAgeOR"], 7.54)
        self.assertEqual(cfg["testType"], "IGRA")
        self.assertEqual(cfg["regimen"], "3HP")
        self.assertIsNone(cfg["riskPrev"]["smoking"])
        self.assertEqual(cfg["diseaseOR"]["renal"], 3.6)
        self.assertIsNone(cfg["pStartTPT"])
        self.assertIsNone(cfg["earlyLateRatio"])

    def test_normalised_config_fills_missing_fields_from_defaults(self) -> None:
        cfg = normalise_config({"scenarioLabel": "Custom", "riskPrev": {"smoking": 0.2}})

        self.assertEqual(cfg["scenarioLabel"], "Custom")
        self.assertEqual(cfg["N"], 1500)
        self.assertEqual(cfg["riskPrev"]["smoking"], 0.2)
        self.assertIsNone(cfg["riskPrev"]["diabetes"])

    def test_strip_empty_fields_removes_empty_nested_values(self) -> None:
        stripped = strip_empty_fields(
            {"a": None, "b": "", "c": [], "d": {"x": None, "y": 1}}
        )

        self.assertEqual(stripped, {"d": {"y": 1}})

    def test_valid_default_config_passes_validation(self) -> None:
        cfg = validate_config(build_default_config())
        report = collect_validation_issues(cfg)

        self.assertTrue(report["isValid"])

    def test_invalid_n_fails(self) -> None:
        report = collect_validation_issues({"N": 0})

        self.assertFalse(report["isValid"])
        self.assertIn("N", report["fatalFieldNames"])

    def test_invalid_screen_coverage_fails(self) -> None:
        report = collect_validation_issues({"screenCoverage": 1.5})

        self.assertFalse(report["isValid"])
        self.assertIn("screenCoverage", report["fatalFieldNames"])

    def test_invalid_test_type_fails(self) -> None:
        report = collect_validation_issues({"testType": "XPERT"})

        self.assertFalse(report["isValid"])
        self.assertIn("testType", report["fatalFieldNames"])

    def test_follow_horizon_must_exceed_screen_window(self) -> None:
        report = collect_validation_issues({"followHorizon": 2, "screenWindow": 2})

        self.assertFalse(report["isValid"])
        self.assertIn("followHorizon", report["fatalFieldNames"])

    def test_invalid_risk_prevalence_fails(self) -> None:
        report = collect_validation_issues({"riskPrev": {"smoking": [0.1, 0.2, 2.0]}})

        self.assertFalse(report["isValid"])
        self.assertIn("riskPrev.smoking", report["fatalFieldNames"])

    def test_invalid_disease_or_fails(self) -> None:
        report = collect_validation_issues({"diseaseOR": {"renal": 0}})

        self.assertFalse(report["isValid"])
        self.assertIn("diseaseOR.renal", report["fatalFieldNames"])

    def test_high_nreps_emits_warning(self) -> None:
        report = collect_validation_issues({"nReps": 2001})

        self.assertTrue(report["isValid"])
        self.assertTrue(report["hasWarnings"])
        self.assertIn("nReps", report["warningFieldNames"])

    def test_validate_config_raises_for_invalid_config(self) -> None:
        with self.assertRaisesRegex(ValueError, "Invalid APY config"):
            validate_config({"testType": "invalid"})


if __name__ == "__main__":
    unittest.main()
