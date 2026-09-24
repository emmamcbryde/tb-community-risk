from __future__ import annotations

from dataclasses import replace
import json
import unittest

from app.general.profile_editing import (
    CSV_COLUMNS,
    apply_risk_factor_rows,
    parse_risk_factor_csv,
    risk_factor_editing_issues,
    risk_factor_rows,
    risk_factor_template_csv,
    risk_factors_csv,
)
from engine.profiles.country import (
    KEEP_CURRENT,
    USE_NEW,
    ConflictRequiresChoice,
    apply_snapshot_country,
    preview_country_application,
)
from engine.profiles.demonstration import build_demonstration_profile
from engine.profiles.effect_measures import (
    ConversionPolicy,
    convert_effect,
    crosswalk_rows,
    effect_warnings,
)
from engine.profiles.engine_mapping import build_engine_config, epidemiological_config_hash
from engine.profiles.local_incidence import apply_local_incidence, parse_local_incidence
from engine.profiles.population_profile import (
    EffectMeasure,
    PopulationProfile,
    Provenance,
    user_override_fields,
    user_value,
    validate_profile,
)
from engine.who_incidence.snapshot import FIXTURE_DIR, find_manifests, load_snapshot


LOCAL = (
    b"location,iso3,year,measure,incidence_per_100k,lower,upper,population,source,notes\n"
    b"North District,,2022,estimated_incidence,50,40,60,250000,District survey 2023,\n"
    b"North District,,2023,estimated_incidence,48,39,58,252000,District survey 2023,\n"
)


class RiskFactorTableTests(unittest.TestCase):
    def setUp(self) -> None:
        self.profile = build_demonstration_profile()

    def test_bounds_and_evidence_fields_round_trip(self) -> None:
        rows = risk_factor_rows(self.profile)
        rows[0].update({"Prevalence low (%)": 10.0, "Prevalence high (%)": 16.0, "Effect low": 3.0, "Effect high": 8.0, "Evidence year": 2021, "Applies to": "Adults 15+"})
        edited = apply_risk_factor_rows(self.profile, rows)
        factor = edited.risk_factors[0]
        self.assertEqual(factor.prevalence_bounds, (0.10, 0.16))
        self.assertEqual(factor.effect_bounds, (3.0, 8.0))
        self.assertEqual(factor.evidence_year, 2021)
        self.assertTrue(factor.is_user_override)
        self.assertEqual(factor.prevalence.provenance, Provenance.BUNDLED)
        self.assertEqual(validate_profile(edited), [])
        self.assertEqual(PopulationProfile.from_json(edited.to_json()), edited)

    def test_invalid_bounds_detected(self) -> None:
        rows = risk_factor_rows(self.profile)
        rows[0].update({"Effect low": 6.0, "Effect high": 8.0})
        edited = apply_risk_factor_rows(self.profile, rows)
        self.assertTrue(any("outside its bounds" in issue["message"] for issue in validate_profile(edited)))
        rows = risk_factor_rows(self.profile)
        rows[0].update({"Effect low": 3.0})
        with self.assertRaises(ValueError):
            apply_risk_factor_rows(self.profile, rows)

    def test_engine_factor_transition_is_fixed(self) -> None:
        rows = risk_factor_rows(self.profile)
        rows[1]["Affected transition"] = "Infection"
        edited = apply_risk_factor_rows(self.profile, rows)
        self.assertTrue(risk_factor_editing_issues(edited))

    def test_csv_export_template_and_import(self) -> None:
        exported = risk_factors_csv(self.profile)
        self.assertEqual(exported.splitlines()[0].split(","), CSV_COLUMNS)
        rows, errors = parse_risk_factor_csv(exported.encode("utf-8"))
        self.assertEqual(errors, [])
        self.assertEqual(apply_risk_factor_rows(self.profile, rows), self.profile)
        self.assertEqual(risk_factor_template_csv().strip().split(","), CSV_COLUMNS)
        bad = exported.replace(",OR,", ",IRR,", 1).encode("utf-8")
        _, errors = parse_risk_factor_csv(bad)
        self.assertTrue(any("RR, HR or OR" in error for error in errors))

    def test_run_without_stratification(self) -> None:
        profile = replace(self.profile, risk_factors=())
        self.assertEqual(validate_profile(profile), [])
        config = build_engine_config(profile)
        self.assertTrue(all(value == 0.0 for value in config["riskPrev"].values() if value is not None))


class EffectMeasureTests(unittest.TestCase):
    def test_warnings_do_not_change_numbers(self) -> None:
        profile = build_demonstration_profile()
        before = epidemiological_config_hash(build_engine_config(profile))
        codes = {warning.code for warning in effect_warnings(profile)}
        self.assertIn("or_as_hazard_multiplier", codes)
        self.assertIn("multiplied_factors", codes)
        self.assertIn("extreme_combined_multiplier", codes)
        self.assertEqual(before, epidemiological_config_hash(build_engine_config(profile)))

    def test_missing_measure_warned_and_hr_not_flagged(self) -> None:
        profile = build_demonstration_profile()
        factors = list(profile.risk_factors)
        factors[0] = replace(factors[0], effect_measure=None, effect_estimate=user_value(2.0, "ratio"))
        factors[1] = replace(factors[1], effect_measure=EffectMeasure.HR)
        profile = replace(profile, risk_factors=tuple(factors[:2]))
        codes = [(w.code, w.risk_factor) for w in effect_warnings(profile)]
        self.assertIn(("missing_measure_type", factors[0].label), codes)
        self.assertNotIn(("or_as_hazard_multiplier", factors[1].label), codes)

    def test_conversion_policy_inactive(self) -> None:
        self.assertEqual(convert_effect(3.0, EffectMeasure.OR), 3.0)
        with self.assertRaises(NotImplementedError):
            convert_effect(3.0, EffectMeasure.OR, ConversionPolicy.OR_TO_RR_BASELINE_RISK)

    def test_crosswalk_covers_every_factor(self) -> None:
        rows = crosswalk_rows(build_demonstration_profile())
        self.assertEqual(len(rows), 7)
        self.assertTrue(all(row["Cap applied"] == "No" for row in rows))
        self.assertTrue(all("without conversion" in row["Engine interpretation"] for row in rows))


class CountryAndLocalDataTests(unittest.TestCase):
    def setUp(self) -> None:
        self.snapshot = load_snapshot(find_manifests(FIXTURE_DIR)[0])
        self.demo = build_demonstration_profile()

    def test_country_application_changes_only_supported_fields(self) -> None:
        changes = preview_country_application(self.demo, self.snapshot, "AUS")
        self.assertEqual({change.key for change in changes}, {"location", "incidence", "national_population"})
        applied = apply_snapshot_country(self.demo, self.snapshot, "AUS")
        for field in ("population_size", "age_distribution", "ltbi_prevalence", "risk_factors"):
            self.assertEqual(getattr(applied, field), getattr(self.demo, field), field)
        self.assertEqual(applied.incidence.provenance, Provenance.WHO_SNAPSHOT)
        self.assertEqual(dict(applied.incidence.source_detail)["snapshotId"], self.snapshot.snapshot_id)
        self.assertEqual(applied.location.national_population_year, 2024)
        self.assertEqual(epidemiological_config_hash(build_engine_config(applied)), epidemiological_config_hash(build_engine_config(self.demo)))

    def test_local_upload_validation_provenance_and_hash(self) -> None:
        upload = parse_local_incidence(LOCAL, filename="north.csv")
        self.assertTrue(upload.is_valid, upload.errors)
        self.assertEqual(upload.content_hash, parse_local_incidence(LOCAL, filename="north.csv").content_hash)
        profile = apply_local_incidence(self.demo, upload)
        self.assertEqual(profile.incidence.provenance, Provenance.LOCAL_UPLOAD)
        self.assertIsNone(profile.incidence.snapshot_id)
        self.assertEqual(profile.location.kind.value, "subnational")
        self.assertEqual(PopulationProfile.from_json(profile.to_json()), profile)
        detail = dict(profile.incidence.source_detail)
        self.assertEqual(detail["fileSha256"], upload.file_sha256)
        self.assertIn("Incidence series (local file)", user_override_fields(profile))

    def test_local_upload_rejections(self) -> None:
        cases = {
            "notification data": LOCAL.replace(b"estimated_incidence,50", b"notification,50"),
            "lower bound exceeds": LOCAL.replace(b"50,40,60", b"50,55,60"),
            "both lower and upper": LOCAL.replace(b"50,40,60", b"50,,60"),
            "more than one location": LOCAL.replace(b"North District,,2023", b"South,,2023"),
            "source is required": LOCAL.replace(b"District survey 2023,\nNorth District,,2023", b",\nNorth District,,2023"),
            "duplicate year": LOCAL.replace(b",2023,", b",2022,"),
        }
        for expected, payload in cases.items():
            with self.subTest(expected=expected):
                result = parse_local_incidence(payload)
                self.assertFalse(result.is_valid)
                self.assertTrue(any(expected in error for error in result.errors), result.errors)
        missing_bounds = parse_local_incidence(LOCAL.replace(b"50,40,60", b"50,,").replace(b"48,39,58", b"48,,"))
        self.assertTrue(missing_bounds.is_valid)
        self.assertTrue(any("propagated trend uncertainty will be unavailable" in w for w in missing_bounds.warnings))
        who = parse_local_incidence(LOCAL.replace(b"District survey 2023", b"WHO copy"))
        self.assertTrue(any("still treated as a local file" in w for w in who.warnings))

    def test_conflicts_require_explicit_choice(self) -> None:
        local = apply_local_incidence(self.demo, parse_local_incidence(LOCAL, filename="north.csv"))
        with self.assertRaises(ConflictRequiresChoice):
            apply_snapshot_country(local, self.snapshot, "AUS")
        kept = apply_snapshot_country(local, self.snapshot, "AUS", resolutions={"incidence": KEEP_CURRENT, "location": KEEP_CURRENT})
        self.assertEqual(kept.incidence, local.incidence)
        replaced = apply_snapshot_country(local, self.snapshot, "AUS", resolutions={"incidence": USE_NEW, "location": USE_NEW})
        self.assertEqual(replaced.incidence.provenance, Provenance.WHO_SNAPSHOT)
        who = apply_snapshot_country(self.demo, self.snapshot, "AUS")
        with self.assertRaises(ConflictRequiresChoice):
            apply_local_incidence(who, parse_local_incidence(LOCAL, filename="north.csv"))

    def test_v1_profile_migrates(self) -> None:
        applied = apply_snapshot_country(self.demo, self.snapshot, "ZAF")
        payload = json.loads(applied.to_json())
        payload["schemaVersion"] = "population_profile_v1"
        payload["incidence"]["provenance"] = "bundled"
        payload["incidence"].pop("dataHash")
        migrated = PopulationProfile.from_dict(payload)
        self.assertEqual(migrated.schema_version, "population_profile_v2")
        self.assertEqual(migrated.incidence.provenance, Provenance.WHO_SNAPSHOT)


if __name__ == "__main__":
    unittest.main()
