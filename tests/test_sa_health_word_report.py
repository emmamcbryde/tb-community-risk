from __future__ import annotations

import csv
import json
from pathlib import Path
import shutil
import unittest
import zipfile
import xml.etree.ElementTree as ET

from engine.apy.sa_health_reference_package import write_sa_health_reference_package
from engine.apy.sa_health_report import (
    build_sa_health_word_report,
    load_report_data,
    report_ready_tables,
    validate_report_data,
)


ROOT = Path(__file__).resolve().parents[1]


class SAHealthWordReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.package_dir = ROOT / "outputs" / "test_tmp" / "sa_health_word_report_package"
        cls.report_dir = ROOT / "reports" / "test_tmp" / "sa_health_word_report"
        shutil.rmtree(cls.package_dir, ignore_errors=True)
        shutil.rmtree(cls.report_dir, ignore_errors=True)
        write_sa_health_reference_package(cls.package_dir, n_reps=10)
        cls.result = build_sa_health_word_report(cls.package_dir, cls.report_dir)

    def test_report_outputs_are_created_and_parseable(self) -> None:
        docx = Path(self.result["docxPath"])
        pdf = Path(self.result["pdfPath"])
        manifest = self.report_dir / "report_manifest.json"
        technical_archive = self.report_dir / "technical_event_ledger_and_economics_csv.zip"

        self.assertTrue(docx.exists())
        self.assertTrue(pdf.exists())
        self.assertTrue(manifest.exists())
        self.assertTrue(technical_archive.exists())
        payload = json.loads(manifest.read_text(encoding="utf-8"))
        self.assertEqual(payload["reportPackageVersion"], "sa_health_apy_health_economic_report_v1")
        self.assertEqual(payload["nReps"], 10)
        self.assertEqual(payload["pageCount"], len(payload["previewPages"]))
        self.assertEqual(payload["technicalCsvArchive"], str(technical_archive))

        with zipfile.ZipFile(technical_archive) as archive:
            names = set(archive.namelist())
        self.assertIn("event_ledger_totals.csv", names)
        self.assertIn("event_ledger_annual.csv", names)
        self.assertIn("economic_replicates.csv", names)
        self.assertIn("economic_annual_by_arm.csv", names)

    def test_word_report_contains_hyperlink_and_no_raw_markdown_placeholders(self) -> None:
        with zipfile.ZipFile(self.result["docxPath"]) as docx:
            document = docx.read("word/document.xml").decode("utf-8")
            rels = docx.read("word/_rels/document.xml.rels").decode("utf-8")

        self.assertIn("Targeted latent tuberculosis infection screening", document)
        self.assertIn("https://tb-community-risk.streamlit.app", rels)
        for forbidden in ["TODO", "FIXME", "INSERT", "**"]:
            self.assertNotIn(forbidden, document)

    def test_word_xml_is_well_formed_and_drawing_ids_are_unique(self) -> None:
        with zipfile.ZipFile(self.result["docxPath"]) as docx:
            document = docx.read("word/document.xml")
            styles = docx.read("word/styles.xml").decode("utf-8")

        root = ET.fromstring(document)
        namespaces = {
            "wp": "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing",
        }
        drawing_ids = [
            item.attrib["id"]
            for item in root.findall(".//wp:docPr", namespaces)
        ]
        self.assertGreater(len(drawing_ids), 1)
        self.assertEqual(len(drawing_ids), len(set(drawing_ids)))
        self.assertIn('w:type="character" w:styleId="Hyperlink"', styles)
        self.assertIn('w:tblW w:w="9638"', document.decode("utf-8"))

    def test_report_tables_reconcile_to_authoritative_package(self) -> None:
        data = load_report_data(self.package_dir)
        validate_report_data(data)
        tables = report_ready_tables(data)

        primary = data["scenarioById"]["primary_working_reference"]
        annual = data["annualBudget"][-1]
        self.assertAlmostEqual(
            float(annual["cumulativeIncrementalCostDiscounted"]),
            float(primary["incrementalCost"]),
        )
        headline = {row["Measure"]: row for row in tables["headline_results"]}
        self.assertIn("net saving", headline["Net included health-system result"]["Value"])
        self.assertEqual(headline["Willingness-to-pay threshold"]["Value"], "Not supplied")

    def test_report_ready_csvs_preserve_provisional_and_not_locally_costed_labels(self) -> None:
        path = self.report_dir / "tables" / "economic_inputs.csv"
        with path.open("r", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        by_item = {row["Cost item"]: row for row in rows}

        programme_rows = [
            row for row in rows
            if any(
                label in row["Cost item"]
                for label in ["Programme setup", "Programme running", "Travel, outreach and staff-support"]
            )
        ]
        self.assertTrue(programme_rows)
        for row in programme_rows:
            self.assertIn("Not yet locally costed", row["Caveat"])
        self.assertIn("Reviewed", by_item["IGRA screening test per person"]["Status"])


if __name__ == "__main__":
    unittest.main()
