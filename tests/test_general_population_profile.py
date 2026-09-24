from __future__ import annotations

from dataclasses import replace
import json
import unittest

from app.general.profile_editing import apply_risk_factor_rows, risk_factor_rows
from app.general.terminology import USER_DEFINED_MARK
from engine.profiles.demonstration import (
    DEMONSTRATION_PROFILE_LABEL,
    GENERAL_DEFAULT_POPULATION,
    GENERAL_DEFAULT_STOCHASTIC_SIMULATIONS,
    build_demonstration_profile,
    demonstration_profile_without_risk_factors,
)
from engine.profiles.engine_mapping import (
    ProfileMappingError,
    build_engine_config,
    default_analysis,
    effect_interpretation_rows,
)
from engine.profiles.population_profile import (
    PROFILE_SCHEMA_VERSION,
    EffectMeasure,
    PopulationProfile,
    ProfileValidationError,
    Provenance,
    ReviewStatus,
    ValueState,
    bundled_value,
    unresolved_inputs,
    user_override_fields,
    user_value,
    validate_profile,
    with_population_size,
)


class PopulationProfileContractTests(unittest.TestCase):
    def setUp(self) -> None:
        self.profile = build_demonstration_profile()

    def test_general_defaults(self) -> None:
        self.assertEqual(GENERAL_DEFAULT_POPULATION, 10_000)
        self.assertEqual(self.profile.population_size.value, 10_000)
        self.assertEqual(GENERAL_DEFAULT_STOCHASTIC_SIMULATIONS, 1_000)
        self.assertEqual(default_analysis()["nReps"], 1_000)
        self.assertEqual(default_analysis()["analysisMethod"], "agent_based")
        self.assertEqual(self.profile.name, DEMONSTRATION_PROFILE_LABEL)
        self.assertTrue(self.profile.demonstration)
        self.assertEqual(validate_profile(self.profile), [])
        self.assertFalse(self.profile.has_user_overrides)

    def test_json_round_trip(self) -> None:
        text = self.profile.to_json()
        restored = PopulationProfile.from_json(text)
        self.assertEqual(restored, self.profile)
        self.assertEqual(restored.to_json(), text)
        self.assertEqual(json.loads(text)["schemaVersion"], PROFILE_SCHEMA_VERSION)

    def test_hash_is_stable_and_content_sensitive(self) -> None:
        first = self.profile.profile_hash()
        self.assertEqual(first, build_demonstration_profile().profile_hash())
        self.assertEqual(first, PopulationProfile.from_json(self.profile.to_json()).profile_hash())
        self.assertRegex(first, r"^[0-9a-f]{64}$")
        changed = with_population_size(self.profile, 9_999)
        self.assertNotEqual(first, changed.profile_hash())

    def test_missing_zero_and_not_reviewed_remain_distinct(self) -> None:
        zero = user_value(0.0, "proportion")
        missing = user_value(None, "proportion")
        not_reviewed = bundled_value(0.2, "proportion", "Bundled")
        self.assertEqual(zero.state, ValueState.VALUE)
        self.assertEqual(zero.value, 0.0)
        self.assertEqual(missing.state, ValueState.MISSING)
        self.assertIsNone(missing.value)
        self.assertIn("zero", zero.status_codes())
        self.assertIn("missing", missing.status_codes())
        self.assertNotIn("zero", missing.status_codes())
        self.assertIn("not_reviewed", not_reviewed.status_codes())
        self.assertEqual(len({zero.status_codes(), missing.status_codes(), not_reviewed.status_codes()}), 3)

        factors = list(self.profile.risk_factors)
        factors[0] = replace(factors[0], prevalence=zero)
        factors[1] = replace(factors[1], prevalence=missing, enabled=False)
        profile = replace(self.profile, risk_factors=tuple(factors))
        restored = PopulationProfile.from_json(profile.to_json())
        self.assertEqual(restored.risk_factors[0].prevalence.value, 0.0)
        self.assertEqual(restored.risk_factors[0].prevalence.state, ValueState.VALUE)
        self.assertIsNone(restored.risk_factors[1].prevalence.value)
        self.assertEqual(restored.risk_factors[1].prevalence.state, ValueState.MISSING)
        self.assertEqual(restored.risk_factors[2].prevalence.review_status, ReviewStatus.NOT_REVIEWED)

    def test_state_and_value_must_agree(self) -> None:
        payload = self.profile.to_dict()
        payload["ltbiPrevalence"]["state"] = "missing"
        profile = PopulationProfile.from_dict(payload)
        self.assertTrue(any(issue["field"] == "ltbiPrevalence" for issue in validate_profile(profile)))
        payload["ltbiPrevalence"]["value"] = "0.1"
        with self.assertRaises(ProfileValidationError):
            PopulationProfile.from_dict(payload)

    def test_effect_measure_types_survive_round_trip_without_conversion(self) -> None:
        factors = list(self.profile.risk_factors)
        for factor, measure in zip(factors[:3], (EffectMeasure.RR, EffectMeasure.HR, EffectMeasure.OR)):
            factors[factors.index(factor)] = replace(factor, effect_measure=measure, effect_estimate=user_value(2.5, "ratio"))
        profile = replace(self.profile, risk_factors=tuple(factors))
        restored = PopulationProfile.from_json(profile.to_json())
        self.assertEqual([f.effect_measure for f in restored.risk_factors[:3]], [EffectMeasure.RR, EffectMeasure.HR, EffectMeasure.OR])
        self.assertEqual([f.effect_estimate.value for f in restored.risk_factors[:3]], [2.5, 2.5, 2.5])
        notes = {row["Declared measure"]: row["Note"] for row in effect_interpretation_rows(restored)}
        self.assertIn("without conversion", notes["OR"])
        self.assertIn("without conversion", notes["RR"])
        config = build_engine_config(restored)
        engine_keys = [f.engine_key for f in restored.risk_factors[:3]]
        self.assertEqual([config["diseaseOR"][key] for key in engine_keys], [2.5, 2.5, 2.5])

    def test_rejects_unknown_effect_measure_and_schema(self) -> None:
        payload = self.profile.to_dict()
        payload["riskFactors"][0]["effectMeasure"] = "IRR"
        with self.assertRaises(ProfileValidationError):
            PopulationProfile.from_dict(payload)
        payload = self.profile.to_dict()
        payload["schemaVersion"] = "population_profile_v0"
        with self.assertRaises(ProfileValidationError):
            PopulationProfile.from_dict(payload)

    def test_profile_without_risk_factors_is_valid_and_maps(self) -> None:
        profile = demonstration_profile_without_risk_factors()
        self.assertEqual(profile.risk_factors, ())
        self.assertEqual(validate_profile(profile), [])
        self.assertEqual(PopulationProfile.from_json(profile.to_json()), profile)
        config = build_engine_config(profile)
        self.assertTrue(all(value == 0.0 for value in config["riskPrev"].values() if value is not None))
        from engine.apy.validation import collect_validation_issues

        self.assertTrue(collect_validation_issues(config)["isValid"])

    def test_unchanged_demonstration_profile_preserves_engine_inputs(self) -> None:
        from engine.apy.working_defaults import build_unified_working_default_preset

        reference = build_unified_working_default_preset()["config"]
        config = build_engine_config(self.profile)
        for key in ("riskPrev", "diseaseOR", "ltbiPrevalence", "testType", "regimen", "screenCoverage", "screeningStrategy"):
            self.assertEqual(config[key], reference[key], key)
        self.assertEqual(config["N"], 10_000)
        self.assertEqual(config["nReps"], 1_000)
        self.assertEqual(config["generalProfileLink"]["profileHash"], self.profile.profile_hash())

    def test_enabled_risk_factor_with_missing_value_blocks_mapping(self) -> None:
        factors = list(self.profile.risk_factors)
        factors[0] = replace(factors[0], prevalence=user_value(None, "proportion"))
        profile = replace(self.profile, risk_factors=tuple(factors))
        self.assertIn("blocking", {issue["severity"] for issue in validate_profile(profile)})
        with self.assertRaises(ProfileMappingError):
            build_engine_config(profile)

    def test_unresolved_inputs_listed_for_demonstration(self) -> None:
        items = {row["Input"] for row in unresolved_inputs(self.profile)}
        self.assertIn("LTBI prevalence", items)
        self.assertIn("TB incidence", items)
        self.assertIn("Age distribution", items)


class RiskFactorEditingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.profile = build_demonstration_profile()

    def test_edits_are_marked_user_defined(self) -> None:
        rows = risk_factor_rows(self.profile)
        rows[0]["Prevalence (%)"] = 25.0
        rows[1]["Effect measure"] = "HR"
        rows[2]["Prevalence (%)"] = None
        rows[3]["Prevalence (%)"] = 0.0
        edited = apply_risk_factor_rows(self.profile, rows)
        factors = edited.risk_factors
        self.assertEqual(factors[0].prevalence.value, 0.25)
        self.assertEqual(factors[0].prevalence.provenance, Provenance.USER_DEFINED)
        self.assertEqual(factors[1].effect_measure, EffectMeasure.HR)
        self.assertEqual(factors[1].effect_estimate.value, self.profile.risk_factors[1].effect_estimate.value)
        self.assertEqual(factors[2].prevalence.state, ValueState.MISSING)
        self.assertTrue(factors[3].prevalence.is_zero)
        self.assertFalse(factors[4].is_user_override)
        statuses = [row["Status"] for row in risk_factor_rows(edited)]
        self.assertTrue(all(status.startswith(USER_DEFINED_MARK) for status in statuses[:4]))
        self.assertFalse(statuses[4].startswith(USER_DEFINED_MARK))
        self.assertEqual(len(user_override_fields(edited)), 4)

    def test_unchanged_rows_leave_profile_identical(self) -> None:
        self.assertEqual(apply_risk_factor_rows(self.profile, risk_factor_rows(self.profile)), self.profile)

    def test_custom_and_deleted_rows(self) -> None:
        rows = risk_factor_rows(self.profile)[:2]
        rows.append({"Enabled": True, "Risk factor": "Silicosis", "Prevalence (%)": 3.0, "Effect estimate": 3.1,
                     "Effect measure": "RR", "id": None})
        edited = apply_risk_factor_rows(self.profile, rows)
        self.assertEqual(len(edited.risk_factors), 3)
        custom = edited.risk_factors[2]
        self.assertIsNone(custom.engine_key)
        self.assertEqual(custom.effect_measure, EffectMeasure.RR)
        self.assertTrue(custom.is_user_override)
        restored = PopulationProfile.from_json(edited.to_json())
        self.assertEqual(restored, edited)
        config = build_engine_config(edited)
        self.assertEqual(config["riskPrev"]["diabetes"], 0.0)


if __name__ == "__main__":
    unittest.main()
