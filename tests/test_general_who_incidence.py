from __future__ import annotations

import csv
import io
import json
from pathlib import Path
import shutil
import socket
import tempfile
import unittest
from unittest.mock import patch

from engine.who_incidence import snapshot as snapshot_module
from engine.who_incidence.contract import COLUMNS, SNAPSHOT_CONTRACT_VERSION, rounding_half_unit
from engine.who_incidence.snapshot import (
    FIXTURE_DIR,
    PRODUCTION_DIR,
    SnapshotError,
    SnapshotUnavailable,
    _as_who_csv,
    find_manifests,
    load_bundled_snapshot,
    load_snapshot,
)
from engine.who_incidence.who_import import (
    ImportFailed,
    build_manifest,
    check_dictionary,
    cross_check,
    parse_who_estimates,
    snapshot_bytes,
    write_snapshot_files,
)


ROOT = Path(__file__).resolve().parents[1]
FIXTURE_MANIFEST = find_manifests(FIXTURE_DIR)[0]


def _no_network(*args, **kwargs):
    raise AssertionError("network access attempted")


def fixture_who_csv() -> bytes:
    return _as_who_csv(list(load_snapshot(FIXTURE_MANIFEST).rows))


def mutate(payload: bytes, code: str, at_year: int, **changes: str) -> bytes:
    rows = list(csv.DictReader(io.StringIO(payload.decode("utf-8"))))
    for row in rows:
        if row["iso3"] == code and row["year"] == str(at_year):
            row.update(changes)
    return write_rows(rows)


def write_rows(rows, fieldnames=None) -> bytes:
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames or list(rows[0].keys()), lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


class InstalledSnapshotTests(unittest.TestCase):
    def test_production_snapshot_complete_and_valid(self) -> None:
        manifests = find_manifests(PRODUCTION_DIR)
        if not manifests:
            self.skipTest("No production snapshot installed.")
        with patch.object(socket.socket, "connect", _no_network):
            snapshot = load_snapshot(manifests[0])
        manifest = snapshot.manifest
        self.assertEqual(manifest["contractVersion"], SNAPSHOT_CONTRACT_VERSION)
        self.assertEqual(manifest["snapshotKind"], "production")
        self.assertTrue(manifest["isCompleteDataset"])
        self.assertEqual(manifest["reportYear"], 2025)
        coverage = manifest["coverage"]
        self.assertEqual(coverage["countriesAndAreas"], len(snapshot.countries()))
        self.assertEqual(coverage["yearRange"], [2000, 2024])
        self.assertGreater(coverage["countriesAndAreas"], 200)
        self.assertEqual(manifest["crossCheck"]["disagreements"], 0)
        self.assertTrue(manifest["crossCheck"]["commit"].startswith("666088c"))
        self.assertIn("not reviewed or endorsed", manifest["terms"]["noEndorsement"])
        for key in ("sourceFiles", "importerVersion", "transformation", "citation", "fields", "knownExclusions", "analyticalHash"):
            self.assertTrue(manifest[key], key)
        self.assertEqual({item["role"] for item in manifest["sourceFiles"]}, {"estimates", "dictionary"})

    def test_default_loader_prefers_production(self) -> None:
        path = snapshot_module.default_manifest_path()
        expected = find_manifests(PRODUCTION_DIR) or find_manifests(FIXTURE_DIR)
        self.assertEqual(path, expected[0])

    def test_published_precision_is_retained(self) -> None:
        payload = (FIXTURE_MANIFEST.parent / json.loads(FIXTURE_MANIFEST.read_text())["dataFile"]["filename"]).read_bytes()
        rows = list(csv.DictReader(io.StringIO(payload.decode("utf-8"))))
        aus2000 = next(row for row in rows if row["iso3"] == "AUS" and row["year"] == "2000")
        self.assertEqual((aus2000["incidence_per_100k"], aus2000["incident_cases"], aus2000["iso_numeric"]), ("6.8", "1300", "036"))

    def test_fixture_edge_cases(self) -> None:
        snapshot = load_snapshot(FIXTURE_MANIFEST)
        self.assertEqual(snapshot.series_for("PRK"), ())
        self.assertIsNone(snapshot.country_summary("PRK")["latest"])
        self.assertEqual(snapshot.country_summary("ANT")["years"], [2000, 2009])
        self.assertIn({"iso3": "PRK", "reason": "No incidence estimates published in this report round."}, snapshot.manifest["knownExclusions"])

    def test_missing_snapshot_gives_maintainer_message(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(snapshot_module, "PRODUCTION_DIR", Path(tmp) / "none"), patch.object(snapshot_module, "FIXTURE_DIR", Path(tmp) / "none2"):
                with self.assertRaises(SnapshotUnavailable) as caught:
                    load_bundled_snapshot()
        self.assertIn("scripts/import_who_incidence.py", str(caught.exception))
        self.assertIn("docs/who_data_update_guide.md", str(caught.exception))

    def test_corruption_detected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp)
            manifest = json.loads(FIXTURE_MANIFEST.read_text(encoding="utf-8"))
            shutil.copy(FIXTURE_MANIFEST, target / FIXTURE_MANIFEST.name)
            data = (FIXTURE_DIR / manifest["dataFile"]["filename"]).read_text(encoding="utf-8")
            (target / manifest["dataFile"]["filename"]).write_text(data.replace("6.8", "6.9", 1), encoding="utf-8", newline="\n")
            with self.assertRaises(SnapshotError) as caught:
                load_snapshot(target / FIXTURE_MANIFEST.name)
            self.assertIn("Checksum mismatch", str(caught.exception))

    def test_invalid_manifest_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp)
            manifest = json.loads(FIXTURE_MANIFEST.read_text(encoding="utf-8"))
            shutil.copy(FIXTURE_DIR / manifest["dataFile"]["filename"], target)
            for broken in ({**manifest, "contractVersion": "who_incidence_snapshot_v1"}, {k: v for k, v in manifest.items() if k != "terms"}):
                (target / "m_manifest.json").write_text(json.dumps(broken), encoding="utf-8")
                with self.assertRaises(SnapshotError):
                    load_snapshot(target / "m_manifest.json")


class ImporterValidationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.payload = fixture_who_csv()

    def test_clean_import_and_determinism(self) -> None:
        rows1, report1 = parse_who_estimates(self.payload, report_year=2025)
        rows2, _ = parse_who_estimates(self.payload, report_year=2025)
        self.assertTrue(report1.is_valid, report1.summary_markdown())
        self.assertEqual(snapshot_bytes(rows1), snapshot_bytes(rows2))
        kwargs = dict(report_year=2025, source_files=[], crosscheck=None, importer_commit="x", kind="test_fixture")
        m1 = build_manifest(rows1, report1, data_bytes=snapshot_bytes(rows1), access_date="2026-01-01", **kwargs)
        m2 = build_manifest(rows2, report1, data_bytes=snapshot_bytes(rows2), access_date="2026-02-02", **kwargs)
        self.assertEqual(m1["analyticalHash"], m2["analyticalHash"])
        self.assertEqual(m1["snapshotId"], m2["snapshotId"])
        self.assertNotEqual(m1["volatile"], m2["volatile"])

    def test_malformed_inputs_fail_clearly(self) -> None:
        cases = {
            "duplicate_row": self.payload + self.payload.split(b"\n", 2)[1] + b"\n",
            "invalid_iso3": mutate(self.payload, "AUS", 2000, iso3="AU1"),
            "lower_exceeds_estimate": mutate(self.payload, "AUS", 2001, e_inc_100k_lo="99"),
            "upper_below_estimate": mutate(self.payload, "AUS", 2001, e_inc_100k_hi="1"),
            "non_numeric": mutate(self.payload, "AUS", 2002, e_inc_100k="six"),
            "negative_value": mutate(self.payload, "AUS", 2002, e_inc_num="-5"),
            "zero_population": mutate(self.payload, "AUS", 2003, e_pop_num="0"),
            "rate_count_inconsistent": mutate(self.payload, "ZAF", 2010, e_inc_num="10"),
            "mixed_vintage": mutate(self.payload, "AUS", 2024, year="2025"),
            "incomplete_triple": mutate(self.payload, "AUS", 2004, e_inc_100k_hi=""),
            "inconsistent_identity": mutate(self.payload, "AUS", 2005, g_whoregion="EUR"),
            "invalid_region": mutate(self.payload, "AUS", 2006, g_whoregion="XXX"),
        }
        for code, payload in cases.items():
            with self.subTest(code=code):
                _, report = parse_who_estimates(payload, report_year=2025)
                self.assertFalse(report.is_valid)
                self.assertIn(code, report.codes("fatal"))

    def test_duplicate_iso3_mapping(self) -> None:
        rows = list(csv.DictReader(io.StringIO(self.payload.decode("utf-8"))))
        for row in rows:
            if row["iso3"] == "BRA":
                row["country"] = "Australia"
        _, report = parse_who_estimates(write_rows(rows), report_year=2025)
        self.assertIn("duplicate_iso3_mapping", report.codes("fatal"))

    def test_changed_schema_encoding_and_empty(self) -> None:
        rows = list(csv.DictReader(io.StringIO(self.payload.decode("utf-8"))))
        fields = [f for f in rows[0] if f != "e_inc_100k_hi"]
        _, report = parse_who_estimates(write_rows(rows, fields), report_year=2025)
        self.assertIn("schema_changed", report.codes("fatal"))
        _, report = parse_who_estimates(b"\xff\xfe\x00c\x00o", report_year=2025)
        self.assertIn("corrupt_encoding", report.codes("fatal"))
        _, report = parse_who_estimates(self.payload.split(b"\n", 1)[0] + b"\n", report_year=2025)
        self.assertIn("empty_file", report.codes("fatal"))

    def test_unexpected_year_range(self) -> None:
        _, report = parse_who_estimates(self.payload, report_year=2026)
        self.assertIn("unexpected_year_range", report.codes("fatal"))

    def test_categories_distinguish_missingness(self) -> None:
        _, report = parse_who_estimates(self.payload, report_year=2025)
        self.assertIn("no_estimates", report.codes("expected_missingness"))
        self.assertIn("short_series", report.codes("incomplete_series"))
        self.assertEqual(report.fatal, [])
        summary = report.summary_markdown()
        self.assertIn("Result: PASSED", summary)
        self.assertIn("Expected missingness (1)", summary)

    def test_crosscheck_tolerance(self) -> None:
        rows, report = parse_who_estimates(self.payload, report_year=2025)
        reference = {}
        for row in rows:
            values = {c: row.number(c) for c in ("incidence_per_100k", "incidence_per_100k_lo", "incidence_per_100k_hi", "incident_cases", "incident_cases_lo", "incident_cases_hi", "population")}
            reference[(row.iso3, row.year)] = values
        result = cross_check(rows, reference, report)
        self.assertEqual(result["disagreements"], 0)
        key = next(k for k in reference if k[0] == "ZAF")
        reference[key] = {**reference[key], "incidence_per_100k": reference[key]["incidence_per_100k"] * 1.2}
        _, fresh = parse_who_estimates(self.payload, report_year=2025)
        self.assertEqual(cross_check(rows, reference, fresh)["disagreements"], 1)
        self.assertIn("crosscheck_disagreement", fresh.codes("fatal"))
        self.assertEqual(rounding_half_unit("6.8"), 0.05)
        self.assertEqual(rounding_half_unit("1300"), 50)

    def test_dictionary_check(self) -> None:
        _, report = parse_who_estimates(self.payload, report_year=2025)
        check_dictionary(b"variable_name,definition\niso3,code\n", report)
        self.assertIn("dictionary_missing_variable", report.codes("fatal"))

    def test_writer_refuses_invalid_report(self) -> None:
        rows, report = parse_who_estimates(mutate(self.payload, "AUS", 2001, e_inc_100k_lo="99"), report_year=2025)
        with tempfile.TemporaryDirectory() as tmp, self.assertRaises(ImportFailed):
            write_snapshot_files(Path(tmp), snapshot_bytes(rows), {"snapshotId": "x"}, report)

    def test_snapshot_columns_match_contract(self) -> None:
        payload = (FIXTURE_DIR / json.loads(FIXTURE_MANIFEST.read_text())["dataFile"]["filename"]).read_bytes()
        self.assertEqual(tuple(payload.split(b"\n", 1)[0].decode().split(",")), COLUMNS)


class ImporterIsOfflineTests(unittest.TestCase):
    def test_import_and_load_without_network_or_gtbreport(self) -> None:
        with patch.object(socket.socket, "connect", _no_network), patch("socket.create_connection", _no_network):
            rows, report = parse_who_estimates(fixture_who_csv(), report_year=2025)
            self.assertTrue(report.is_valid)
            load_snapshot(FIXTURE_MANIFEST)
        for module in ("engine/who_incidence/snapshot.py", "engine/who_incidence/who_import.py", "engine/who_incidence/contract.py"):
            source = (ROOT / module).read_text(encoding="utf-8")
            self.assertNotRegex(source, r"import (requests|urllib|http\.client)")
            self.assertNotIn("../gtbreport2025", source)


if __name__ == "__main__":
    unittest.main()
