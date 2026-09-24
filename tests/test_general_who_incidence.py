from __future__ import annotations

import csv
import io
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from engine.profiles.country import apply_snapshot_country, apply_user_incidence
from engine.profiles.demonstration import build_demonstration_profile
from engine.profiles.population_profile import LocationKind, PopulationProfile, Provenance, validate_profile
from engine.who_incidence.adapters import WhoBurdenCsvAdapter, parse_user_incidence_csv, transform
from engine.who_incidence.schema import (
    RECORD_COLUMNS,
    check_vintage_compatibility,
    parse_rows,
    validate_columns,
)
from engine.who_incidence.snapshot import (
    BUNDLED_MANIFEST,
    SnapshotError,
    load_bundled_snapshot,
    load_snapshot,
    write_snapshot,
)
from engine.who_incidence.trend import (
    INCIDENCE_TO_INFECTION_POLICY,
    TrendMethod,
    TrendSettings,
    covid_disruption_flags,
    estimate_trend,
)


def _row(iso3="AUS", year=2020, est=6.0, lo=5.0, hi=7.0, **extra):
    row = {column: "" for column in RECORD_COLUMNS}
    row.update(
        {"iso3": iso3, "country": "Example", "year": str(year), "incidence_per_100k": str(est),
         "incidence_per_100k_lo": str(lo), "incidence_per_100k_hi": str(hi)}
    )
    row.update(extra)
    return row


class BundledSnapshotTests(unittest.TestCase):
    def test_fixture_loads_offline_with_manifest_provenance(self) -> None:
        snapshot = load_bundled_snapshot()
        manifest = snapshot.manifest
        self.assertEqual(manifest["snapshotKind"], "test_fixture")
        self.assertFalse(snapshot.is_complete_dataset)
        self.assertEqual(manifest["sourceReportYear"], 2025)
        self.assertEqual(manifest["upstream"]["commit"], "666088cac1e20dd1e9e52016ea857a9c48e0ba8d")
        for key in ("repository", "files"):
            self.assertTrue(manifest["upstream"][key])
        for key in ("extractionDate", "transformationVersion", "schemaVersion", "licence", "citation"):
            self.assertTrue(manifest[key])
        self.assertGreaterEqual(len(snapshot.countries()), 5)
        self.assertTrue(snapshot.validation["isValid"])
        series = snapshot.series_for("ZAF")
        self.assertEqual((series[0].year, series[-1].year), (2000, 2024))
        self.assertTrue(all(p.lower <= p.estimate <= p.upper for p in series))

    def test_checksum_mismatch_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp)
            manifest = json.loads(BUNDLED_MANIFEST.read_text(encoding="utf-8"))
            shutil.copy(BUNDLED_MANIFEST, target / BUNDLED_MANIFEST.name)
            data = (BUNDLED_MANIFEST.parent / manifest["dataFile"]["filename"]).read_text(encoding="utf-8")
            (target / manifest["dataFile"]["filename"]).write_text(data.replace("6.82", "6.83", 1), encoding="utf-8", newline="\n")
            with self.assertRaises(SnapshotError):
                load_snapshot(target / BUNDLED_MANIFEST.name)

    def test_changed_schema_version_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp)
            manifest = json.loads(BUNDLED_MANIFEST.read_text(encoding="utf-8"))
            shutil.copy(BUNDLED_MANIFEST.parent / manifest["dataFile"]["filename"], target)
            manifest["schemaVersion"] = "who_incidence_snapshot_v2"
            (target / "m.json").write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaises(SnapshotError):
                load_snapshot(target / "m.json")

    def test_country_selection_attaches_incidence_but_not_other_evidence(self) -> None:
        snapshot = load_bundled_snapshot()
        demo = build_demonstration_profile()
        profile = apply_snapshot_country(demo, snapshot, "PHL")
        self.assertEqual(profile.location.iso3, "PHL")
        self.assertEqual(profile.location.kind, LocationKind.COUNTRY)
        self.assertEqual(profile.incidence.snapshot_id, snapshot.snapshot_id)
        self.assertEqual(profile.incidence.data_year_range, (2000, 2024))
        self.assertEqual(profile.ltbi_prevalence, demo.ltbi_prevalence)
        self.assertEqual(profile.risk_factors, demo.risk_factors)
        self.assertEqual(validate_profile(profile), [])
        self.assertEqual(PopulationProfile.from_json(profile.to_json()), profile)
        again = apply_snapshot_country(profile, snapshot, "ZAF")
        self.assertEqual(again.profile_id, "zaf-demonstration-working-defaults")
        self.assertTrue(again.name.startswith("South Africa: "))


class IncidenceValidationTests(unittest.TestCase):
    def test_valid_rows(self) -> None:
        records, report = parse_rows([_row(year=2020), _row(year=2021)])
        self.assertTrue(report.is_valid, report.to_dict())
        self.assertEqual(len(records), 2)

    def test_malformed_bounds(self) -> None:
        _, report = parse_rows([_row(est=6.0, lo=6.5, hi=7.0), _row(year=2021, est=6.0, lo=5.0, hi=5.5)])
        self.assertFalse(report.is_valid)
        self.assertIn("lower_exceeds_estimate", report.codes())
        self.assertIn("upper_below_estimate", report.codes())

    def test_duplicate_rows(self) -> None:
        _, report = parse_rows([_row(year=2020), _row(year=2020)])
        self.assertIn("duplicate_row", report.codes())
        self.assertFalse(report.is_valid)

    def test_invalid_iso3_non_numeric_negative(self) -> None:
        _, report = parse_rows([_row(iso3="AU"), _row(year=2021, est="six"), _row(year=2022, est=-1, lo=-2, hi=1)])
        self.assertTrue({"invalid_iso3", "non_numeric", "negative_value"} <= report.codes())

    def test_missing_years_and_discontinuities_warn(self) -> None:
        _, report = parse_rows([_row(year=2018), _row(year=2020), _row(year=2021, est=30, lo=25, hi=35)])
        self.assertTrue(report.is_valid)
        self.assertIn("missing_years", report.codes())
        self.assertIn("implausible_discontinuity", report.codes())

    def test_schema_change_detected(self) -> None:
        self.assertFalse(validate_columns(["iso3", "year", "e_inc_100k"]).is_valid)
        row = _row()
        row.pop("incidence_per_100k_hi")
        _, report = parse_rows([row])
        self.assertIn("schema_changed", report.codes())

    def test_incompatible_vintages(self) -> None:
        base = {"sourceDataset": "WHO", "schemaVersion": "who_incidence_snapshot_v1", "transformationVersion": "t1"}
        report = check_vintage_compatibility([{**base, "sourceReportYear": 2024}, {**base, "sourceReportYear": 2025}])
        self.assertIn("incompatible_vintage", report.codes())
        self.assertTrue(check_vintage_compatibility([{**base, "sourceReportYear": 2025}] * 2).is_valid)

    def test_who_public_csv_adapter_and_snapshot_writer(self) -> None:
        header = ["country", "iso2", "iso3", "g_whoregion", "year", "e_pop_num", "e_inc_100k", "e_inc_100k_lo",
                  "e_inc_100k_hi", "e_inc_num", "e_inc_num_lo", "e_inc_num_hi"]
        buffer = io.StringIO()
        writer = csv.writer(buffer)
        writer.writerow(header)
        writer.writerow(["Exampleland", "EX", "EXA", "AFR", 2023, 1000000, 100, 80, 120, 1000, 800, 1200])
        writer.writerow(["Exampleland", "EX", "EXA", "AFR", 2024, 1010000, 95, 76, 115, 960, 770, 1160])
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "burden.csv"
            source.write_text(buffer.getvalue(), encoding="utf-8")
            records, report = transform(WhoBurdenCsvAdapter(source))
            self.assertTrue(report.is_valid)
            manifest = write_snapshot(
                records,
                out_dir=Path(tmp) / "out",
                data_filename="s.csv",
                manifest_filename="s_manifest.json",
                manifest_fields={
                    "snapshotId": "s", "snapshotKind": "test_fixture", "isCompleteDataset": False,
                    "sourceDataset": "WHO", "sourceReportYear": 2025, "upstream": {}, "extractionDate": "2026-01-01",
                    "licence": {"status": "to_be_confirmed"}, "citation": "WHO",
                },
            )
            loaded = load_snapshot(Path(tmp) / "out" / "s_manifest.json")
            self.assertEqual(manifest["dataFile"]["rows"], 2)
            self.assertEqual(loaded.series_for("EXA")[1].estimate, 95.0)

    def test_user_subnational_upload(self) -> None:
        text = "country,year,incidence_per_100k,incidence_per_100k_lo,incidence_per_100k_hi\nNorth District,2022,50,40,60\nNorth District,2023,48,39,58\n"
        records, report = parse_user_incidence_csv(text)
        self.assertTrue(report.is_valid, report.to_dict())
        profile = apply_user_incidence(build_demonstration_profile(), records, source_label="north.csv")
        self.assertEqual(profile.location.kind, LocationKind.SUBNATIONAL)
        self.assertEqual(profile.incidence.provenance, Provenance.USER_DEFINED)
        self.assertEqual(PopulationProfile.from_json(profile.to_json()), profile)
        _, bad = parse_user_incidence_csv(text.replace("40,60", "55,60", 1))
        self.assertFalse(bad.is_valid)


class TrendInterfaceTests(unittest.TestCase):
    def test_observed_values_preserved_and_unimplemented_methods_refuse(self) -> None:
        series = load_bundled_snapshot().series_for("IDN")
        result = estimate_trend(series, TrendSettings())
        self.assertEqual([p.observed for p in result.points], [p.estimate for p in series])
        self.assertTrue(all(p.fitted is None for p in result.points))
        self.assertEqual(result.quantity, "estimated_tb_disease_incidence")
        for method in (TrendMethod.LOG_LINEAR_RECENT, TrendMethod.PENALISED_SPLINE, TrendMethod.STATE_SPACE):
            with self.assertRaises(NotImplementedError):
                estimate_trend(series, TrendSettings(method=method))
        override = estimate_trend(series, TrendSettings(method=TrendMethod.USER_OVERRIDE, user_annual_percent_change=-2.0))
        self.assertEqual(override.annual_percent_change, -2.0)
        self.assertEqual(covid_disruption_flags([2019, 2020, 2023]), [False, True, False])
        self.assertIn("not used as the slope of infection pressure", INCIDENCE_TO_INFECTION_POLICY)


if __name__ == "__main__":
    unittest.main()
