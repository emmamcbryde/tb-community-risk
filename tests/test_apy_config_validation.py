from __future__ import annotations

import unittest

from engine.apy.config import build_default_config, normalise_config, strip_empty_fields
from engine.apy.ltbi_state import (
    COMPATIBILITY_MODE_WARNING,
    LTBI_STATE_MODEL,
    apply_ltbi_state_edit,
    enable_development_compatibility_mode,
    is_clinician_ready_ltbi_state,
    require_numeric_ltbi_state_assumptions,
    resolve_ltbi_state_assumptions,
)
from engine.apy.validation import collect_validation_issues, validate_config


class ApyConfigValidationTests(unittest.TestCase):
    def test_default_config_has_expected_key_values(self) -> None:
        cfg = build_default_config()

        self.assertEqual(cfg["configVersion"], "apy_v9_config_v1")
        self.assertEqual(cfg["modelVersion"], "v9")
        self.assertEqual(cfg["scenarioLabel"], "APY demonstration LTBI screening scenario")
        self.assertEqual(cfg["populationPresetId"], "apy_demonstration")
        self.assertTrue(cfg["useDefaults"])
        self.assertEqual(cfg["N"], 1500)
        self.assertEqual(cfg["nReps"], 2000)
        self.assertEqual(cfg["seed"], 1)
        self.assertEqual(cfg["screenWindow"], 3)
        self.assertEqual(cfg["followHorizon"], 20)
        self.assertEqual(cfg["screenCoverage"], 0.30)
        self.assertEqual(cfg["screeningStrategy"], "prevent")
        self.assertEqual(cfg["targetAgeOR"], 7.54)
        self.assertEqual(cfg["testType"], "IGRA")
        self.assertEqual(cfg["regimen"], "3HP")
        self.assertIsNone(cfg["riskPrev"]["smoking"])
        self.assertEqual(cfg["diseaseOR"]["renal"], 3.6)
        self.assertEqual(cfg["testSensitivity"], 0.95)
        self.assertEqual(cfg["testSpecificity"], 0.98)
        self.assertEqual(cfg["tstSensitivity"], 0.80)
        self.assertIsNone(cfg["pStartTPT"])
        self.assertIsNone(cfg["earlyLateRatio"])
        self.assertIsNone(cfg["ltbiStateAssumptions"]["baselineRecentLTBIProportion"])
        self.assertFalse(cfg["ltbiStateAssumptions"]["developmentCompatibilityMode"])

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
        self.assertTrue(report["hasWarnings"])
        self.assertIn(
            "ltbiStateAssumptions.baselineRecentLTBIProportion",
            report["warningFieldNames"],
        )

    def test_nested_ltbi_state_fields_are_authoritative(self) -> None:
        cfg = normalise_config(
            {
                "ltbiStateAssumptions": {
                    "baselineRecentLTBIProportion": 0.42,
                    "recentToRemoteTransitionRatePerYear": 0.25,
                    "source": "test source",
                    "status": "configured_reviewed",
                    "developmentCompatibilityMode": False,
                }
            }
        )
        state = resolve_ltbi_state_assumptions(cfg)

        self.assertEqual(cfg["baselineRecentLTBIProportion"], 0.42)
        self.assertEqual(cfg["recentToRemoteTransitionRatePerYear"], 0.25)
        self.assertEqual(state["baselineRecentLTBIProportion"], 0.42)
        self.assertEqual(state["source"], "test source")
        self.assertFalse(state["provisional"])

    def test_legacy_only_ltbi_state_fields_are_migrated_as_provisional(self) -> None:
        cfg = normalise_config(
            {
                "baselineRecentLTBIProportion": 0.3,
                "recentToRemoteTransitionRatePerYear": 0.4,
            }
        )

        self.assertEqual(cfg["ltbiStateAssumptions"]["baselineRecentLTBIProportion"], 0.3)
        self.assertEqual(cfg["ltbiStateAssumptions"]["recentToRemoteTransitionRatePerYear"], 0.4)
        self.assertEqual(cfg["ltbiStateAssumptions"]["status"], "migrated_legacy_unverified")
        self.assertTrue(cfg["ltbiStateAssumptions"]["provisional"])

    def test_conflicting_ltbi_state_fields_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "Conflicting LTBI-state configuration"):
            normalise_config(
                {
                    "baselineRecentLTBIProportion": 0.1,
                    "ltbiStateAssumptions": {
                        "baselineRecentLTBIProportion": 0.2,
                        "recentToRemoteTransitionRatePerYear": 0.2,
                    },
                }
            )

    def test_unresolved_ltbi_state_is_flagged_and_clinician_ready_rejects(self) -> None:
        cfg = build_default_config()
        report = collect_validation_issues(cfg)
        clinician_report = collect_validation_issues(cfg, clinician_ready=True)

        self.assertTrue(report["isValid"])
        self.assertIn(
            "ltbiStateAssumptions.baselineRecentLTBIProportion",
            report["warningFieldNames"],
        )
        self.assertFalse(clinician_report["isValid"])
        self.assertIn(
            "ltbiStateAssumptions.baselineRecentLTBIProportion",
            clinician_report["fatalFieldNames"],
        )
        self.assertFalse(is_clinician_ready_ltbi_state(cfg))

    def test_development_compatibility_mode_is_visibly_provisional(self) -> None:
        state = resolve_ltbi_state_assumptions(
            enable_development_compatibility_mode(build_default_config())
        )

        self.assertEqual(state["baselineRecentLTBIProportion"], 0.0)
        self.assertTrue(state["developmentCompatibilityMode"])
        self.assertTrue(state["provisional"])
        self.assertIn(COMPATIBILITY_MODE_WARNING, state["warning"])

    def test_unresolved_default_run_requires_explicit_compatibility_opt_in(self) -> None:
        state = resolve_ltbi_state_assumptions(build_default_config())

        self.assertIsNone(state["baselineRecentLTBIProportion"])
        self.assertFalse(state["developmentCompatibilityMode"])
        with self.assertRaisesRegex(ValueError, "developmentCompatibilityMode=true"):
            require_numeric_ltbi_state_assumptions(build_default_config())

    def test_legacy_migration_fails_clinician_ready_validation(self) -> None:
        cfg = normalise_config(
            {
                "baselineRecentLTBIProportion": 0.3,
                "recentToRemoteTransitionRatePerYear": 0.2,
            }
        )
        report = collect_validation_issues(cfg, clinician_ready=True)

        self.assertFalse(report["isValid"])
        self.assertFalse(is_clinician_ready_ltbi_state(cfg))
        self.assertIn(
            "ltbiStateAssumptions.baselineRecentLTBIProportion",
            report["fatalFieldNames"],
        )

    def test_provisional_status_remains_non_ready(self) -> None:
        cfg = apply_ltbi_state_edit(
            build_default_config(),
            baseline_recent_proportion=0.2,
            transition_rate_per_year=0.2,
            source="test source",
            status="provisional",
            notes="not reviewed",
        )
        state = resolve_ltbi_state_assumptions(cfg)
        report = collect_validation_issues(cfg, clinician_ready=True)

        self.assertTrue(state["provisional"])
        self.assertFalse(report["isValid"])

    def test_transition_provenance_alone_does_not_validate_recent_fraction(self) -> None:
        cfg = normalise_config(
            {
                "ltbiStateAssumptions": {
                    "baselineRecentLTBIProportion": 0.2,
                    "recentToRemoteTransitionRatePerYear": 0.2,
                    "transitionModelSource": "Reviewed transition source",
                    "transitionModelStatus": "configured_reviewed",
                    "status": "configured_reviewed",
                    "baselineRecentLTBIProportionStatus": "configured_reviewed",
                    "baselineRecentLTBIProportionSource": "",
                    "developmentCompatibilityMode": False,
                }
            }
        )

        self.assertFalse(is_clinician_ready_ltbi_state(cfg))
        self.assertFalse(collect_validation_issues(cfg, clinician_ready=True)["isValid"])

    def test_explicit_reviewed_value_can_pass_clinician_ready_validation(self) -> None:
        cfg = apply_ltbi_state_edit(
            build_default_config(),
            baseline_recent_proportion=0.24,
            transition_rate_per_year=0.2,
            source="Scientific model-derived fixture source",
            status="configured_reviewed",
            notes="Reviewed for clinician-ready validation test.",
        )
        report = collect_validation_issues(cfg, clinician_ready=True)

        self.assertTrue(report["isValid"])
        self.assertTrue(is_clinician_ready_ltbi_state(cfg))

    def test_markov_terminology_and_implied_mean_residence_time(self) -> None:
        state = resolve_ltbi_state_assumptions(
            normalise_config(
                {
                    "ltbiStateAssumptions": {
                        "baselineRecentLTBIProportion": 0.25,
                        "recentToRemoteTransitionRatePerYear": 0.2,
                        "status": "configured_reviewed",
                        "source": "test",
                        "developmentCompatibilityMode": False,
                    }
                }
            )
        )

        self.assertEqual(state["transitionModel"], LTBI_STATE_MODEL)
        self.assertAlmostEqual(state["impliedMeanResidenceTimeYears"], 5.0)
        self.assertIn("mean residence time", state["stateDefinition"])
        self.assertNotIn("infection acquired within the last 5 years", state["stateDefinition"])

    def test_ltbi_state_save_load_round_trip_preserves_provenance(self) -> None:
        cfg = apply_ltbi_state_edit(
            build_default_config(),
            baseline_recent_proportion=0.27,
            transition_rate_per_year=0.18,
            source="round-trip source",
            status="configured",
            notes="round-trip notes",
        )
        loaded = normalise_config(cfg)

        self.assertEqual(loaded["ltbiStateAssumptions"]["baselineRecentLTBIProportion"], 0.27)
        self.assertEqual(loaded["ltbiStateAssumptions"]["source"], "round-trip source")
        self.assertEqual(loaded["ltbiStateAssumptions"]["status"], "configured")
        self.assertEqual(loaded["ltbiStateAssumptions"]["notes"], "round-trip notes")

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

    def test_follow_horizon_must_be_at_least_screen_window(self) -> None:
        report = collect_validation_issues({"followHorizon": 2, "screenWindow": 2})

        self.assertTrue(report["isValid"])

        report = collect_validation_issues({"followHorizon": 1.9, "screenWindow": 2})

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
