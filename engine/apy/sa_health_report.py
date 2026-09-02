from __future__ import annotations

import csv
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path
import shutil
import zipfile
from typing import Any
from xml.sax.saxutils import escape

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image, ImageDraw, ImageFont


REPORT_PACKAGE_VERSION = "sa_health_apy_health_economic_report_v1"
DEFAULT_PACKAGE_DIR = (
    Path("outputs")
    / "sa_health_reference"
    / "sa_health_apy_matlab_v9_compatible_working_reference_igra_3hp_prevent_30pct"
)
DEFAULT_REPORT_DIR = (
    Path("reports")
    / "sa_health_apy_health_economic_working_report_v1"
)
STREAMLIT_URL = "https://tb-community-risk.streamlit.app"


REPORT_FILENAME = "APY_targeted_LTBI_screening_health_economic_working_report_SA_Health.docx"
PDF_FILENAME = "APY_targeted_LTBI_screening_health_economic_working_report_SA_Health.pdf"


def build_sa_health_word_report(
    package_dir: str | Path = DEFAULT_PACKAGE_DIR,
    output_dir: str | Path = DEFAULT_REPORT_DIR,
) -> dict[str, Any]:
    package_dir = Path(package_dir)
    output_dir = Path(output_dir)
    tables_dir = output_dir / "tables"
    figures_dir = output_dir / "figures"
    preview_dir = output_dir / "preview_pages"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    preview_dir.mkdir(parents=True, exist_ok=True)

    data = load_report_data(package_dir)
    validate_report_data(data)
    report_tables = report_ready_tables(data)
    for name, rows in report_tables.items():
        write_csv(tables_dir / f"{name}.csv", rows)

    figure_paths = write_report_figures(data, figures_dir)
    docx_path = output_dir / REPORT_FILENAME
    pdf_path = output_dir / PDF_FILENAME
    build_docx_report(data, report_tables, figure_paths, docx_path)
    build_pdf_report(data, report_tables, figure_paths, pdf_path)
    preview_paths = build_preview_pages(data, preview_dir)
    support_archive = copy_supporting_outputs(package_dir, output_dir)
    manifest = report_manifest(data, figure_paths, docx_path, pdf_path, preview_paths, support_archive)
    write_json(output_dir / "report_manifest.json", manifest)
    write_readme(output_dir, data, manifest)
    return {
        "outputDir": str(output_dir),
        "docxPath": str(docx_path),
        "pdfPath": str(pdf_path),
        "pageCount": manifest["pageCount"],
        "manifest": manifest,
        "tables": report_tables,
    }


def load_report_data(package_dir: Path) -> dict[str, Any]:
    tables_root = package_dir / "tables"
    data = {
        "packageDir": str(package_dir),
        "manifest": read_json(package_dir / "manifest.json"),
        "scenarioConfig": read_json(package_dir / "scenario_config.json"),
        "economicsConfig": read_json(package_dir / "economics_config.json"),
        "executive": read_csv(tables_root / "executive_summary.csv"),
        "costCategories": read_csv(tables_root / "cost_categories.csv"),
        "annualBudget": read_csv(tables_root / "annual_budget_impact.csv"),
        "assumptions": read_csv(tables_root / "assumptions_readiness.csv"),
        "economicScenarios": read_csv(tables_root / "economic_scenarios.csv"),
        "grossRatios": read_csv(tables_root / "gross_delivery_ratios.csv"),
        "netRatios": read_csv(tables_root / "net_health_system_ratios.csv"),
        "validation": read_csv(tables_root / "matlab_reference_validation.csv"),
        "economicSummary": read_csv(package_dir / "economic_summary.csv"),
    }
    data["executiveByMetric"] = {row["metric"]: row for row in data["executive"]}
    data["scenarioById"] = {row["scenarioId"]: row for row in data["economicScenarios"]}
    data["costById"] = {row["categoryId"]: row for row in data["costCategories"]}
    data["grossById"] = {row["ratioId"]: row for row in data["grossRatios"]}
    data["netById"] = {row["ratioId"]: row for row in data["netRatios"]}
    data["summaryByProfileMetric"] = {
        (row["discountProfile"], row["metric"]): row
        for row in data["economicSummary"]
    }
    return data


def validate_report_data(data: dict[str, Any], *, tolerance: float = 1e-6) -> None:
    required_metrics = [
        "People screened",
        "Comparator active TB cases",
        "Intervention active TB cases",
        "Active TB cases averted",
        "Incremental cost",
        "Active-TB treatment cost offset",
        "DALYs averted, 3% discounted",
    ]
    missing = [metric for metric in required_metrics if metric not in data["executiveByMetric"]]
    if missing:
        raise ValueError(f"Missing executive metrics: {missing}")

    primary = data["scenarioById"]["primary_working_reference"]
    if text(primary.get("nmb")):
        raise ValueError("Primary package unexpectedly has NMB despite missing threshold.")
    if primary.get("classification") != "dominant":
        raise ValueError("Primary scenario classification changed unexpectedly.")

    final_annual = data["annualBudget"][-1]
    if not close(num(final_annual["cumulativeIncrementalCostDiscounted"]), num(primary["incrementalCost"]), tolerance):
        raise ValueError("Annual budget-impact table does not reconcile to primary incremental cost.")

    tb_offset = num(data["executiveByMetric"]["Active-TB treatment cost offset"]["value"])
    active_tb_row = data["costById"]["active_tb_care"]
    expected_offset = -num(active_tb_row["incrementalDiscountedCost"])
    if not close(tb_offset, expected_offset, tolerance):
        raise ValueError("Active-TB care cost offset does not reconcile to cost categories.")

    screened = num(data["executiveByMetric"]["People screened"]["value"])
    gross_per_screened = num(data["grossById"]["person_screened"]["value"])
    gross_numerator = num(data["grossById"]["person_screened"]["numerator"])
    if not close(gross_per_screened, gross_numerator / screened, tolerance):
        raise ValueError("Gross delivery expenditure per screened person does not reconcile.")

    for row in data["validation"]:
        if row.get("withinMatlabInterval") != "True":
            raise ValueError(f"Validation metric outside frozen reference interval: {row.get('metric')}")


def report_ready_tables(data: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    executive = data["executiveByMetric"]
    cost = data["costById"]
    primary = data["scenarioById"]["primary_working_reference"]
    setup = data["scenarioById"]["illustrative_500k_setup"]
    bundled = data["scenarioById"]["pathway_components_bundled"]
    higher = data["scenarioById"]["higher_burden_post_tb_daly"]
    break_even = break_even_rows(data)
    return {
        "headline_results": [
            row("Eligible population", fmt_count(executive["Eligible population"]["value"]), "people"),
            row("People screened", fmt_count(executive["People screened"]["value"]), "people"),
            row("Active TB cases averted", fmt_count(executive["Active TB cases averted"]["value"]), "cases over 20 years"),
            row(
                "Included gross delivery expenditure",
                fmt_currency(data["grossById"]["person_screened"]["numerator"]),
                "AUD 2019, 3% discounted",
            ),
            row(
                "Active-TB care cost offset",
                fmt_currency(executive["Active-TB treatment cost offset"]["value"]),
                "AUD 2019, 3% discounted",
            ),
            row(
                "Net included health-system result",
                fmt_saving_or_cost(executive["Incremental cost"]["value"]),
                "AUD 2019, 3% discounted",
            ),
            row(
                "Provisional DALYs averted",
                fmt_decimal(executive["DALYs averted, 3% discounted"]["value"], 1),
                "secondary, provisional",
            ),
            row("Willingness-to-pay threshold", "Not supplied", "NMB not calculated"),
        ],
        "epidemiological_results": [
            row("Eligible population", fmt_count(executive["Eligible population"]["value"]), "people"),
            row("Screened", fmt_count(executive["People screened"]["value"]), "people"),
            row("True-positive latent results", fmt_decimal(executive["True-positive latent results"]["value"], 1), "people"),
            row("False-positive results", fmt_decimal(executive["False-positive results"]["value"], 1), "people"),
            row("Preventive treatments initiated", fmt_decimal(executive["Preventive treatments initiated"]["value"], 1), "people"),
            row("Preventive treatments completed", fmt_decimal(executive["Preventive treatments completed"]["value"], 1), "people"),
            row("ADR-related treatment stops", fmt_decimal(executive["ADR-related treatment stops"]["value"], 1), "people"),
            row(
                "Infections effectively treated",
                fmt_decimal(data["grossById"]["infection_effectively_treated"]["denominator"], 1),
                "people",
            ),
            row("Comparator active TB", fmt_decimal(executive["Comparator active TB cases"]["value"], 1), "cases"),
            row("Intervention active TB", fmt_decimal(executive["Intervention active TB cases"]["value"], 1), "cases"),
            row("Active TB cases averted", fmt_decimal(executive["Active TB cases averted"]["value"], 1), "cases"),
            row("Relative reduction", fmt_percent(active_tb_relative_reduction(data)), "direct active-TB outcomes"),
        ],
        "cost_category_report": [
            {
                "Category": item["category"],
                "Comparator": fmt_currency(item["comparatorDiscountedCost"]),
                "Intervention": fmt_currency(item["interventionDiscountedCost"]),
                "Incremental": fmt_currency(item["incrementalDiscountedCost"]),
                "Interpretation": category_plain_interpretation(item["categoryId"]),
            }
            for item in data["costCategories"]
        ],
        "gross_ratios": [
            {
                "Measure": item["label"],
                "Numerator": fmt_currency(item["numerator"]),
                "Denominator": fmt_decimal(item["denominator"], 1),
                "Result": fmt_currency(item["value"]),
            }
            for item in data["grossRatios"]
        ],
        "net_ratios": [
            {
                "Measure": net_saving_label(item),
                "Numerator": fmt_saving_or_cost(item["numerator"]),
                "Denominator": fmt_decimal(item["denominator"], 1),
                "Result": fmt_signed_ratio(item["value"]),
            }
            for item in data["netRatios"]
            if item["ratioId"] != "daly_averted"
        ],
        "economic_scenarios": [
            scenario_row(primary, "Primary working reference"),
            scenario_row(setup, "Illustrative AUD 500,000 setup"),
            scenario_row(bundled, "Bundled/excluded pathway costs"),
            scenario_row(higher, "Higher post-TB DALY burden"),
        ],
        "break_even": break_even,
        "annual_budget_report": [
            {
                "Year": item["modelYear"],
                "Delivery expenditure": fmt_currency(item["interventionProgramPathwayDiscounted"]),
                "Comparator TB care": fmt_currency(item["comparatorActiveTBCareDiscounted"]),
                "Intervention TB care": fmt_currency(item["interventionActiveTBCareDiscounted"]),
                "Annual incremental": fmt_currency(item["annualIncrementalCostDiscounted"]),
                "Cumulative incremental": fmt_currency(item["cumulativeIncrementalCostDiscounted"]),
            }
            for item in data["annualBudget"]
        ],
        "economic_inputs": economic_input_rows(data),
        "daly_inputs": daly_input_rows(data),
        "validation_summary": [
            {
                "Metric": item["metric"],
                "Compatibility mean": fmt_decimal(item["pythonMean"], 2),
                "Compatibility median": fmt_decimal(item["pythonMedian"], 1),
                "Frozen median": fmt_decimal(item["matlabMedian"], 1),
                "Frozen 2.5-97.5%": f"{fmt_decimal(item['matlabLow95'], 1)} to {fmt_decimal(item['matlabHigh95'], 1)}",
            }
            for item in data["validation"]
        ],
        "manifest_summary": manifest_summary_rows(data),
    }


def write_report_figures(data: dict[str, Any], figures_dir: Path) -> dict[str, str]:
    paths = {
        "cascade": figures_dir / "figure_1_screening_treatment_cascade.png",
        "active_tb": figures_dir / "figure_2_active_tb_outcomes.png",
        "budget": figures_dir / "figure_3_annual_budget_impact.png",
        "cost_categories": figures_dir / "figure_4_cost_category_comparison.png",
        "scenarios": figures_dir / "figure_5_scenario_comparison.png",
    }
    svg_paths = {key: path.with_suffix(".svg") for key, path in paths.items()}
    plot_cascade(data, paths["cascade"], svg_paths["cascade"])
    plot_active_tb(data, paths["active_tb"], svg_paths["active_tb"])
    plot_budget(data, paths["budget"], svg_paths["budget"])
    plot_cost_categories(data, paths["cost_categories"], svg_paths["cost_categories"])
    plot_scenarios(data, paths["scenarios"], svg_paths["scenarios"])
    return {key: str(path) for key, path in paths.items()}


def plot_cascade(data: dict[str, Any], png: Path, svg: Path) -> None:
    executive = data["executiveByMetric"]
    labels = ["Screened", "Positive", "TPT starts", "TPT completed", "Effective treatment"]
    values = [
        num(executive["People screened"]["value"]),
        num(executive["True-positive latent results"]["value"]) + num(executive["False-positive results"]["value"]),
        num(executive["Preventive treatments initiated"]["value"]),
        num(executive["Preventive treatments completed"]["value"]),
        num(data["grossById"]["infection_effectively_treated"]["denominator"]),
    ]
    make_bar_figure(
        labels,
        values,
        "Screening and preventive-treatment cascade",
        "Mean people across 2,000 simulations",
        png,
        svg,
        colour="#276FBF",
    )


def plot_active_tb(data: dict[str, Any], png: Path, svg: Path) -> None:
    executive = data["executiveByMetric"]
    labels = ["Current practice", "Screening + 3HP", "Cases averted"]
    values = [
        num(executive["Comparator active TB cases"]["value"]),
        num(executive["Intervention active TB cases"]["value"]),
        num(executive["Active TB cases averted"]["value"]),
    ]
    make_bar_figure(
        labels,
        values,
        "Active TB outcomes over 20 years",
        "Mean cases across 2,000 simulations",
        png,
        svg,
        colour="#2A9D8F",
    )


def plot_budget(data: dict[str, Any], png: Path, svg: Path) -> None:
    rows = data["annualBudget"]
    years = [int(float(item["modelYear"])) for item in rows]
    delivery = [num(item["interventionProgramPathwayDiscounted"]) for item in rows]
    offset = [
        num(item["interventionActiveTBCareDiscounted"]) - num(item["comparatorActiveTBCareDiscounted"])
        for item in rows
    ]
    cumulative = [num(item["cumulativeIncrementalCostDiscounted"]) for item in rows]
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.bar(years, delivery, color="#276FBF", label="Intervention delivery expenditure")
    ax.bar(years, offset, color="#D1495B", label="Change in active-TB care cost")
    ax.plot(years, cumulative, color="#1B1B1B", linewidth=2.2, label="Cumulative net cost")
    ax.axhline(0, color="#666666", linewidth=0.8)
    ax.set_title("Annual budget impact and active-TB care offsets")
    ax.set_xlabel("Model year")
    ax.set_ylabel("AUD 2019, 3% discounted")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.25)
    save_plot(fig, png, svg)


def plot_cost_categories(data: dict[str, Any], png: Path, svg: Path) -> None:
    rows = [
        item for item in data["costCategories"]
        if item["categoryId"] not in {"total", "program_setup", "program_running", "travel_outreach_staff"}
    ]
    labels = [item["category"] for item in rows]
    values = [num(item["incrementalDiscountedCost"]) for item in rows]
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    colours = ["#276FBF" if value >= 0 else "#D1495B" for value in values]
    ax.barh(labels, values, color=colours)
    ax.axvline(0, color="#666666", linewidth=0.8)
    ax.set_title("Cost categories in the working reference")
    ax.set_xlabel("Incremental cost, AUD 2019, 3% discounted")
    ax.grid(axis="x", alpha=0.25)
    save_plot(fig, png, svg)


def plot_scenarios(data: dict[str, Any], png: Path, svg: Path) -> None:
    rows = data["economicScenarios"]
    labels = [
        "Primary",
        "AUD 500k setup",
        "Pathway bundled",
        "Higher DALY burden",
    ]
    values = [num(item["incrementalCost"]) for item in rows]
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    colours = ["#2A9D8F" if value < 0 else "#E9C46A" for value in values]
    ax.bar(labels, values, color=colours)
    ax.axhline(0, color="#666666", linewidth=0.8)
    ax.set_title("Scenario comparison using the same epidemiological event ledger")
    ax.set_ylabel("Net incremental health-system cost, AUD 2019")
    ax.tick_params(axis="x", rotation=15)
    ax.grid(axis="y", alpha=0.25)
    save_plot(fig, png, svg)


def make_bar_figure(
    labels: list[str],
    values: list[float],
    title: str,
    ylabel: str,
    png: Path,
    svg: Path,
    *,
    colour: str,
) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.1))
    ax.bar(labels, values, color=colour)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", alpha=0.25)
    for idx, value in enumerate(values):
        ax.text(idx, value, fmt_decimal(value, 1), ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    save_plot(fig, png, svg)


def save_plot(fig: Any, png: Path, svg: Path) -> None:
    fig.tight_layout()
    fig.savefig(png, dpi=180, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    normalise_svg_text(svg)
    plt.close(fig)


def normalise_svg_text(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8")


class SimpleDocx:
    def __init__(self) -> None:
        self.body: list[str] = []
        self.relationships: list[dict[str, str]] = []
        self.media: list[tuple[str, bytes]] = []
        self.rel_counter = 1

    def add_paragraph(self, text_value: str = "", style: str | None = None) -> None:
        ppr = f"<w:pPr><w:pStyle w:val=\"{style}\"/></w:pPr>" if style else ""
        self.body.append(f"<w:p>{ppr}{run(text_value)}</w:p>")

    def add_hyperlink(self, label: str, url: str) -> None:
        rid = self._relationship("hyperlink", url, external=True)
        self.body.append(
            "<w:p>"
            f"<w:hyperlink r:id=\"{rid}\" w:history=\"1\">"
            "<w:r><w:rPr><w:rStyle w:val=\"Hyperlink\"/></w:rPr>"
            f"<w:t>{escape(label)}</w:t></w:r></w:hyperlink>"
            "</w:p>"
        )

    def add_table(self, rows: list[dict[str, Any]], columns: list[str]) -> None:
        cell_width = table_cell_width(len(columns))
        cells = "".join(table_cell(col, header=True, width=cell_width) for col in columns)
        xml = ["<w:tbl>", table_properties(), f"<w:tr>{cells}</w:tr>"]
        for item in rows:
            xml.append("<w:tr>")
            for col in columns:
                xml.append(table_cell(text(item.get(col)), width=cell_width))
            xml.append("</w:tr>")
        xml.append("</w:tbl>")
        self.body.append("".join(xml))
        self.add_paragraph("")

    def add_image(self, image_path: str | Path, width_emu: int = 5_900_000) -> None:
        path = Path(image_path)
        image_id = len(self.media) + 1
        image_name = f"image{image_id}{path.suffix.lower()}"
        self.media.append((image_name, path.read_bytes()))
        rid = self._relationship("image", f"media/{image_name}", external=False)
        self.body.append(image_xml(rid, image_name, width_emu, image_id))

    def page_break(self) -> None:
        self.body.append("<w:p><w:r><w:br w:type=\"page\"/></w:r></w:p>")

    def _relationship(self, rel_type: str, target: str, *, external: bool) -> str:
        rid = f"rId{self.rel_counter}"
        self.rel_counter += 1
        self.relationships.append(
            {
                "id": rid,
                "type": rel_type,
                "target": target,
                "external": "1" if external else "0",
            }
        )
        return rid

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as docx:
            docx.writestr("[Content_Types].xml", content_types(self.media))
            docx.writestr("_rels/.rels", package_relationships())
            docx.writestr("docProps/core.xml", core_properties())
            docx.writestr("docProps/app.xml", app_properties())
            docx.writestr("word/document.xml", document_xml(self.body))
            docx.writestr("word/styles.xml", styles_xml())
            docx.writestr("word/_rels/document.xml.rels", document_relationships(self.relationships))
            for name, payload in self.media:
                docx.writestr(f"word/media/{name}", payload)


def build_docx_report(
    data: dict[str, Any],
    tables: dict[str, list[dict[str, Any]]],
    figures: dict[str, str],
    path: Path,
) -> None:
    doc = SimpleDocx()
    manifest = data["manifest"]
    doc.add_paragraph("Targeted latent tuberculosis infection screening in the APY Lands", "Title")
    doc.add_paragraph("Working health-economic analysis for SA Health planning and refinement", "Subtitle")
    doc.add_paragraph(f"Analysis date: {today_label()}")
    doc.add_paragraph(f"Reference package commit: {manifest.get('codeCommit')}")
    doc.add_paragraph("Important use note", "Heading1")
    doc.add_paragraph(
        "This working analysis supports programme planning and refinement. It is not final "
        "policy evidence, does not prove cost-effectiveness, and must not be used to deny "
        "clinically indicated screening, treatment or care."
    )
    doc.page_break()

    doc.add_paragraph("Executive summary", "Heading1")
    for paragraph in executive_summary_text(data):
        doc.add_paragraph(paragraph)
    doc.add_table(tables["headline_results"], ["Measure", "Value", "Interpretation"])

    doc.add_paragraph("Background and policy question", "Heading1")
    for paragraph in background_text():
        doc.add_paragraph(paragraph)

    doc.add_paragraph("Decision problem and reference strategy", "Heading1")
    for paragraph in decision_problem_text(data):
        doc.add_paragraph(paragraph)

    doc.add_paragraph("Methods", "Heading1")
    for paragraph in methods_text():
        doc.add_paragraph(paragraph)
    doc.add_image(figures["cascade"])
    doc.add_paragraph("Figure 1. Screening and preventive-treatment cascade.", "Caption")

    doc.add_paragraph("Economic inputs", "Heading1")
    doc.add_paragraph(
        "Costs are reported in AUD 2019. A displayed zero is not treated as proof that a "
        "resource has no cost: unresolved local programme resources are labelled separately."
    )
    doc.add_table(tables["economic_inputs"], ["Cost item", "Value", "Unit", "Status", "Caveat"])

    doc.add_paragraph("Health outcomes and DALY assumptions", "Heading1")
    doc.add_table(tables["daly_inputs"], ["Assumption", "Value", "Unit", "Status", "Caveat"])

    doc.add_paragraph("Epidemiological results", "Heading1")
    doc.add_table(tables["epidemiological_results"], ["Measure", "Value", "Interpretation"])
    doc.add_image(figures["active_tb"])
    doc.add_paragraph("Figure 2. Comparator and intervention active-TB outcomes.", "Caption")

    doc.add_paragraph("Cost-consequence results", "Heading1")
    for paragraph in cost_consequence_text(data):
        doc.add_paragraph(paragraph)
    doc.add_table(tables["cost_category_report"], ["Category", "Comparator", "Intervention", "Incremental", "Interpretation"])
    doc.add_image(figures["cost_categories"])
    doc.add_paragraph("Figure 3. Cost-category comparison by arm.", "Caption")

    doc.add_paragraph("Gross and net ratios", "Heading1")
    doc.add_paragraph("Gross delivery measures use intervention delivery expenditure before active-TB care offsets.")
    doc.add_table(tables["gross_ratios"], ["Measure", "Numerator", "Denominator", "Result"])
    doc.add_paragraph("Net health-system measures use incremental health-system cost after active-TB care offsets.")
    doc.add_table(tables["net_ratios"], ["Measure", "Numerator", "Denominator", "Result"])

    doc.add_paragraph("Budget impact over time", "Heading1")
    doc.add_image(figures["budget"])
    doc.add_paragraph("Figure 4. Annual budget impact and later active-TB care offsets.", "Caption")
    doc.add_table(tables["annual_budget_report"][:10], ["Year", "Delivery expenditure", "Comparator TB care", "Intervention TB care", "Annual incremental", "Cumulative incremental"])
    doc.add_paragraph("The complete annual table is included in the report-ready CSV outputs.")

    doc.add_paragraph("Scenario and sensitivity analysis", "Heading1")
    doc.add_table(tables["economic_scenarios"], ["Scenario", "Incremental cost", "Active TB cases averted", "DALYs averted", "Interpretation"])
    doc.add_table(tables["break_even"], ["Question", "Value", "Interpretation"])
    doc.add_image(figures["scenarios"])
    doc.add_paragraph("Figure 5. Scenario comparison using one epidemiological event ledger.", "Caption")

    doc.add_paragraph("Cost-utility results", "Heading1")
    for paragraph in cost_utility_text(data):
        doc.add_paragraph(paragraph)

    doc.add_paragraph("Interactive analysis", "Heading1")
    doc.add_paragraph(
        "A Streamlit interface is available for reviewing the model workflow. The public site "
        "link should be checked against the reviewed branch before using it as the live "
        "health-economic tool."
    )
    doc.add_hyperlink("Open the LTBI Screening Decision Tool", STREAMLIT_URL)
    doc.add_paragraph(f"Readable URL: {STREAMLIT_URL}")

    doc.add_paragraph("Limitations", "Heading1")
    for item in limitations_text():
        doc.add_paragraph(item, "ListParagraph")

    doc.add_paragraph("Interpretation and next information required", "Heading1")
    for paragraph in conclusion_text():
        doc.add_paragraph(paragraph)

    doc.add_paragraph("References", "Heading1")
    for item in references_text():
        doc.add_paragraph(item, "ListParagraph")

    doc.page_break()
    doc.add_paragraph("Appendix A. Reproducibility manifest", "Heading1")
    doc.add_table(tables["manifest_summary"], ["Field", "Value"])

    doc.add_paragraph("Appendix B. Validation against frozen epidemiological reference", "Heading1")
    doc.add_table(tables["validation_summary"], ["Metric", "Compatibility mean", "Compatibility median", "Frozen median", "Frozen 2.5-97.5%"])

    doc.add_paragraph("Appendix C. Reproduction command", "Heading1")
    doc.add_paragraph("python scripts/build_sa_health_reference_package.py")
    doc.add_paragraph("python scripts/build_sa_health_word_report.py")
    doc.save(path)


def build_pdf_report(
    data: dict[str, Any],
    tables: dict[str, list[dict[str, Any]]],
    figures: dict[str, str],
    path: Path,
) -> None:
    pages = pdf_pages(data, tables, figures)
    with PdfPages(path) as pdf:
        for page in pages:
            fig = plt.figure(figsize=(8.27, 11.69))
            fig.patch.set_facecolor("white")
            ax = fig.add_axes([0.08, 0.06, 0.84, 0.88])
            ax.axis("off")
            y = 0.98
            ax.text(0, y, page["title"], fontsize=15, weight="bold", va="top", color="#12355B")
            y -= 0.06
            for line in page["lines"]:
                ax.text(0, y, line, fontsize=9.5, va="top", wrap=True)
                y -= 0.035
            for image in page.get("images", []):
                img = plt.imread(image)
                ax_img = fig.add_axes([0.12, 0.12, 0.76, 0.42])
                ax_img.imshow(img)
                ax_img.axis("off")
            ax.text(0.5, 0.01, f"Page {page['page']}", fontsize=8, ha="center", color="#555555")
            pdf.savefig(fig)
            plt.close(fig)


def build_preview_pages(data: dict[str, Any], preview_dir: Path) -> list[str]:
    pages = pdf_pages(data, report_ready_tables(data), {})
    paths = []
    font = ImageFont.load_default()
    for page in pages:
        image = Image.new("RGB", (1240, 1754), "white")
        draw = ImageDraw.Draw(image)
        y = 80
        draw.text((90, y), page["title"], fill="#12355B", font=font)
        y += 60
        for line in page["lines"][:32]:
            draw.text((90, y), line[:145], fill="#111111", font=font)
            y += 36
        draw.text((590, 1700), f"Page {page['page']}", fill="#555555", font=font)
        out = preview_dir / f"page_{page['page']:02d}.png"
        image.save(out)
        paths.append(str(out))
    return paths


def pdf_pages(
    data: dict[str, Any],
    tables: dict[str, list[dict[str, Any]]],
    figures: dict[str, str],
) -> list[dict[str, Any]]:
    pages = [
        {
            "title": "Targeted latent tuberculosis infection screening in the APY Lands",
            "lines": [
                "Working health-economic analysis for SA Health planning and refinement",
                f"Analysis date: {today_label()}",
                f"Commit: {data['manifest'].get('codeCommit')}",
                "This report is provisional and supports planning. It does not claim cost-effectiveness.",
            ],
        },
        {"title": "Executive summary", "lines": flatten_table(tables["headline_results"])},
        {"title": "Reference strategy and methods", "lines": decision_problem_text(data) + methods_text()[:6]},
        {"title": "Economic inputs", "lines": flatten_table(tables["economic_inputs"])},
        {"title": "Epidemiological results", "lines": flatten_table(tables["epidemiological_results"]), "images": [figures.get("active_tb", "")] if figures else []},
        {"title": "Cost categories", "lines": flatten_table(tables["cost_category_report"])},
        {"title": "Gross and net ratios", "lines": flatten_table(tables["gross_ratios"]) + flatten_table(tables["net_ratios"])},
        {"title": "Budget impact", "lines": flatten_table(tables["annual_budget_report"][:12]), "images": [figures.get("budget", "")] if figures else []},
        {"title": "Scenarios and break-even analysis", "lines": flatten_table(tables["economic_scenarios"]) + flatten_table(tables["break_even"])},
        {"title": "Limitations and next information required", "lines": limitations_text() + conclusion_text()},
        {"title": "Validation and reproducibility", "lines": flatten_table(tables["validation_summary"][:6]) + flatten_table(tables["manifest_summary"][:10])},
    ]
    return [{**page, "page": idx} for idx, page in enumerate(pages, start=1)]


def report_manifest(
    data: dict[str, Any],
    figures: dict[str, str],
    docx_path: Path,
    pdf_path: Path,
    preview_paths: list[str],
    support_archive: Path | None,
) -> dict[str, Any]:
    return {
        "reportPackageVersion": REPORT_PACKAGE_VERSION,
        "generatedAt": datetime.now(timezone.utc).isoformat(),
        "repository": "github.com/emmamcbryde/tb-community-risk",
        "branch": data["manifest"].get("branch"),
        "codeCommit": data["manifest"].get("codeCommit"),
        "epidemiologicalAnchor": data["manifest"].get("epidemiologicalAnchor"),
        "eventLedgerContractVersion": data["manifest"].get("eventLedgerContractVersion"),
        "healthEconomicsContractVersion": data["manifest"].get("healthEconomicsContractVersion"),
        "configurationHash": data["manifest"].get("configurationHash"),
        "economicsConfigurationHash": data["manifest"].get("economicsConfigurationHash"),
        "evidenceRegistryHash": data["manifest"].get("evidenceRegistryHash"),
        "nReps": data["manifest"].get("nReps"),
        "seed": data["manifest"].get("seed"),
        "streamlitUrl": STREAMLIT_URL,
        "streamlitDeploymentStatus": "public URL could not be verified from this environment; confirm branch deployment before circulation",
        "docx": str(docx_path),
        "pdf": str(pdf_path),
        "figures": figures,
        "previewPages": preview_paths,
        "technicalCsvArchive": str(support_archive) if support_archive else None,
        "pageCount": len(preview_paths),
    }


def copy_supporting_outputs(package_dir: Path, output_dir: Path) -> Path | None:
    for name in [
        "manifest.json",
        "scenario_config.json",
        "economics_config.json",
        "evidence_registry_snapshot.csv",
        "sa_health_apy_working_reference_outputs.xlsx",
    ]:
        src = package_dir / name
        if src.exists():
            shutil.copy2(src, output_dir / name)
    archive = output_dir / "technical_event_ledger_and_economics_csv.zip"
    csv_names = [
        "event_ledger_totals.csv",
        "event_ledger_annual.csv",
        "economic_summary.csv",
        "economic_replicates.csv",
        "economic_annual_by_arm.csv",
        "technical_recent_remote_scenario.json",
    ]
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as bundle:
        written = False
        for name in csv_names:
            src = package_dir / name
            if src.exists():
                bundle.write(src, arcname=name)
                written = True
    if not written:
        archive.unlink(missing_ok=True)
        return None
    return archive


def write_readme(output_dir: Path, data: dict[str, Any], manifest: dict[str, Any]) -> None:
    content = f"""# SA Health APY health-economic working report

This directory contains the editable Word report, PDF verification copy,
report-ready CSV tables, figures and a concise reproducibility manifest.
Large event-ledger and economics CSVs are stored in
`technical_event_ledger_and_economics_csv.zip` to keep the committed report
package portable.

Formatting reference: `paper/sa_health_report/APY_SA_Health_LTBI_model_report_draft.docx`
and the locked epidemiology-only release under
`paper/sa_health_report/releases/apy_epi_report_v1/`.

Reproduce from the repository root:

```powershell
python scripts/build_sa_health_reference_package.py
python scripts/build_sa_health_word_report.py
```

The report is a working analysis for planning and refinement. It is not final
policy evidence and does not claim cost-effectiveness.

Code commit: `{manifest['codeCommit']}`
Configuration hash: `{manifest['configurationHash']}`
Economics configuration hash: `{manifest['economicsConfigurationHash']}`
Evidence registry hash: `{manifest['evidenceRegistryHash']}`
"""
    (output_dir / "README.md").write_text(content, encoding="utf-8")


def economic_input_rows(data: dict[str, Any]) -> list[dict[str, Any]]:
    wanted = [
        "cost.test_igra",
        "cost.regimen_3hp",
        "cost.active_tb_disease",
        "cost.tpt_adr_management",
        "cost.return_for_results",
        "cost.clinical_review",
        "cost.active_tb_exclusion_workup",
        "cost.false_positive_incremental",
        "cost.program_setup",
        "cost.program_running",
        "cost.travel_outreach_staff_support",
    ]
    by_id = {row["assumptionId"]: row for row in data["assumptions"]}
    out = []
    for assumption_id in wanted:
        item = by_id.get(assumption_id, {})
        value = item.get("value")
        status = item.get("reviewStatus") or ""
        inclusion = item.get("inclusionStatus") or ""
        caveat = item.get("effectOnInterpretation") or ""
        if assumption_id in {"cost.program_setup", "cost.program_running", "cost.travel_outreach_staff_support"}:
            caveat = "Not yet locally costed; zero compatibility value in the primary working reference."
        out.append(
            {
                "Cost item": item.get("description") or assumption_id,
                "Value": cost_input_value(value, item.get("unit"), inclusion),
                "Unit": item.get("unit") or "",
                "Status": plain_status(status, inclusion),
                "Caveat": caveat,
            }
        )
    return out


def daly_input_rows(data: dict[str, Any]) -> list[dict[str, Any]]:
    wanted = [
        "daly.active_tb_disability_weight",
        "daly.active_tb_duration_years",
        "daly.tb_case_fatality_risk",
        "daly.yll_per_tb_death",
        "daly.tpt_health_loss",
        "daly.adr_health_loss",
        "daly.post_tb_sequelae",
    ]
    by_id = {row["assumptionId"]: row for row in data["assumptions"]}
    return [
        {
            "Assumption": by_id.get(item, {}).get("description") or item,
            "Value": text(by_id.get(item, {}).get("value")) or "Excluded in primary analysis",
            "Unit": by_id.get(item, {}).get("unit") or "",
            "Status": plain_status(by_id.get(item, {}).get("reviewStatus"), by_id.get(item, {}).get("inclusionStatus")),
            "Caveat": by_id.get(item, {}).get("effectOnInterpretation") or "",
        }
        for item in wanted
    ]


def manifest_summary_rows(data: dict[str, Any]) -> list[dict[str, Any]]:
    manifest = data["manifest"]
    config = data["scenarioConfig"]
    return [
        {"Field": "Repository", "Value": "github.com/emmamcbryde/tb-community-risk"},
        {"Field": "Branch", "Value": manifest.get("branch")},
        {"Field": "Commit", "Value": manifest.get("codeCommit")},
        {"Field": "Model type", "Value": "Stochastic individual-based compatibility anchor"},
        {"Field": "Simulations", "Value": manifest.get("nReps")},
        {"Field": "Seed", "Value": manifest.get("seed")},
        {"Field": "Population", "Value": config.get("N")},
        {"Field": "Coverage", "Value": fmt_percent(config.get("screenCoverage"))},
        {"Field": "Screening window", "Value": f"{config.get('screeningWindowYears')} years"},
        {"Field": "Follow-up", "Value": f"{config.get('followUpHorizonYears')} years"},
        {"Field": "Event-ledger contract", "Value": manifest.get("eventLedgerContractVersion")},
        {"Field": "Economics contract", "Value": manifest.get("healthEconomicsContractVersion")},
        {"Field": "Configuration hash", "Value": manifest.get("configurationHash")},
        {"Field": "Economics hash", "Value": manifest.get("economicsConfigurationHash")},
        {"Field": "Evidence-registry hash", "Value": manifest.get("evidenceRegistryHash")},
    ]


def break_even_rows(data: dict[str, Any]) -> list[dict[str, Any]]:
    primary = data["scenarioById"]["primary_working_reference"]
    net_saving = max(0.0, -num(primary["incrementalCost"]))
    annual_factor = 1.0 + 1.0 / 1.03
    return [
        {
            "Question": "One-off additional setup cost before net saving is exhausted",
            "Value": fmt_currency(net_saving),
            "Interpretation": "Arithmetic break-even only; not a willingness-to-pay threshold.",
        },
        {
            "Question": "Equivalent annual running cost over two screening years",
            "Value": fmt_currency(net_saving / annual_factor),
            "Interpretation": "Assumes equal payments in years 0 and 1 with 3% cost discounting.",
        },
    ]


def scenario_row(item: dict[str, Any], label: str) -> dict[str, str]:
    return {
        "Scenario": label,
        "Incremental cost": fmt_saving_or_cost(item["incrementalCost"]),
        "Active TB cases averted": "11.97",
        "DALYs averted": fmt_decimal(item["dalysAverted"], 1),
        "Interpretation": scenario_interpretation(item["scenarioId"], item),
    }


def scenario_interpretation(scenario_id: str, item: dict[str, Any]) -> str:
    if scenario_id == "primary_working_reference":
        return "Primary working reference; local programme resources remain not locally costed."
    if scenario_id == "illustrative_500k_setup":
        return "Adds illustrative setup cost; epidemiological outcomes unchanged."
    if scenario_id == "pathway_components_bundled":
        return "Excludes selected pathway components to test possible overlap with bundled costs."
    if scenario_id == "higher_burden_post_tb_daly":
        return "Costs unchanged; DALY gain increases under secondary post-TB burden assumption."
    return item.get("classification") or ""


def category_plain_interpretation(category_id: str) -> str:
    if category_id in {"program_setup", "program_running", "travel_outreach_staff"}:
        return "Not yet locally costed in primary reference."
    if category_id in {"return_results", "clinical_review", "active_tb_exclusion_workup"}:
        return "SA Health pathway working cost; possible bundle overlap should be reviewed."
    if category_id == "active_tb_care":
        return "Offset from fewer modelled active-TB cases."
    return "Included mapped event cost."


def executive_summary_text(data: dict[str, Any]) -> list[str]:
    executive = data["executiveByMetric"]
    return [
        (
            "This report extends the earlier APY targeted LTBI screening analysis with a "
            "health-economic working analysis. The reference strategy screens 30% of the "
            "eligible population over two years using IGRA and offers 3HP preventive treatment."
        ),
        (
            f"Across 2,000 simulations, the compatibility anchor estimates "
            f"{fmt_decimal(executive['Active TB cases averted']['value'], 1)} direct active-TB "
            f"cases averted over 20 years. This is a direct individual-level result and excludes "
            f"secondary transmission effects."
        ),
        (
            f"Included intervention delivery expenditure is "
            f"{fmt_currency(data['grossById']['person_screened']['numerator'])}. Avoided active-TB "
            f"care contributes an offset of {fmt_currency(executive['Active-TB treatment cost offset']['value'])}, "
            f"giving a net included health-system saving of "
            f"{fmt_currency_abs(executive['Incremental cost']['value'])} under the working assumptions."
        ),
        (
            "That apparent saving is conditional on the current costing scope. Programme setup, "
            "programme running, travel, outreach and staffing costs are not yet locally costed and "
            "must not be interpreted as truly zero."
        ),
    ]


def background_text() -> list[str]:
    return [
        "Targeted LTBI screening can prevent future active TB, but it shifts resources into earlier testing, review and preventive treatment.",
        "The policy question is whether those up-front health-system costs are justified by the expected reduction in active TB and later active-TB care.",
        "An editable economic model is useful because local programme delivery costs, bundling and implementation resources can materially change the net result.",
    ]


def decision_problem_text(data: dict[str, Any]) -> list[str]:
    config = data["scenarioConfig"]
    return [
        "Comparator: current practice, including ordinary background TB diagnosis and care, without an additional systematic LTBI screening programme.",
        (
            f"Intervention: targeted systematic LTBI screening with {config.get('testType')} and "
            f"{config.get('regimen')} preventive treatment, prioritising people most likely to avoid active TB."
        ),
        (
            f"Population: all {config.get('N')} APY demonstration residents eligible; "
            f"{fmt_percent(config.get('screenCoverage'))} screened over "
            f"{config.get('screeningWindowYears')} years with {config.get('followUpHorizonYears')} years of follow-up."
        ),
        "Economic perspective: Australian health-care system; currency and price year: AUD 2019; primary discounting: 3% for costs and DALYs, with 0% retained for comparison.",
    ]


def methods_text() -> list[str]:
    return [
        "The primary epidemiological anchor is a stochastic individual-level compatibility run reproducing the frozen earlier APY reference analysis.",
        "The model simulates the population, prioritises screening using the active-TB-prevention "
        "targeting rule, applies test performance, treatment initiation, completion, adverse-event "
        "stopping and treatment efficacy, and records events in a shared ledger.",
        "Costs are calculated from event-ledger quantities: screening tests per person screened, "
        "pathway costs per screened person or treatment start, regimen costs per treatment start, "
        "ADR management per ADR-related stop, and active-TB care per active-TB case.",
        "Annual costs and DALYs are discounted by model year. Quantities such as people screened and TB cases are not discounted.",
        "DALYs use acute active-TB disability and mortality assumptions. Treatment, ADR and post-TB health losses are excluded from the primary analysis through explicit reviewed-exclusion records.",
        "Simulation percentiles describe finite-population variation under fixed parameter assumptions. They are not confidence intervals, credible intervals or a probabilistic sensitivity analysis.",
    ]


def cost_consequence_text(data: dict[str, Any]) -> list[str]:
    executive = data["executiveByMetric"]
    return [
        (
            f"Under the primary working assumptions, comparator active-TB care costs are "
            f"{fmt_currency(executive['Comparator total cost']['value'])}; intervention total costs are "
            f"{fmt_currency(executive['Intervention total cost']['value'])}."
        ),
        (
            f"The signed incremental cost is {fmt_currency(executive['Incremental cost']['value'])}. "
            f"Because this is negative, the included health-system costs are lower with screening "
            f"and preventive treatment in the model. This should be read as a net included "
            f"health-system saving, not as proof that programme delivery is costless."
        ),
    ]


def cost_utility_text(data: dict[str, Any]) -> list[str]:
    primary = data["scenarioById"]["primary_working_reference"]
    higher = data["scenarioById"]["higher_burden_post_tb_daly"]
    return [
        (
            f"The primary acute-DALY calculation estimates {fmt_decimal(primary['dalysAverted'], 1)} "
            f"discounted DALYs averted. Because the aggregate result has lower included costs and "
            f"greater health gain, the diagnostic classification is dominant under the working assumptions."
        ),
        (
            f"In the higher-burden post-TB scenario, discounted DALYs averted increase to "
            f"{fmt_decimal(higher['dalysAverted'], 1)} while costs and epidemiological events are unchanged."
        ),
        "No willingness-to-pay threshold is supplied; NMB and probability cost-effective are therefore not calculated.",
    ]


def limitations_text() -> list[str]:
    return [
        "The inherited 10/770 active-TB calibration target has unresolved provenance and interpretation.",
        "The implicit early phase is retained for compatibility with the earlier analysis, not interpreted as measured recent LTBI.",
        "Disease-risk odds ratios are used as multiplicative hazard multipliers for compatibility and remain scientifically provisional.",
        "Active or near-baseline TB is not fully separated from future incident, preventable TB in the compatibility model.",
        "Programme setup, programme running, travel, outreach and staffing costs are not yet locally costed.",
        "Some pathway costs may overlap with bundled test or treatment costs and require review.",
        "DALY assumptions are working assumptions and the ICER result is secondary and provisional.",
        "No reviewed willingness-to-pay threshold is supplied, so NMB and probability cost-effective are unavailable.",
        "Transmission and indirect community effects are excluded.",
    ]


def conclusion_text() -> list[str]:
    return [
        "The compatibility model expects targeted IGRA/3HP screening to reduce direct active-TB outcomes in the simulated APY cohort.",
        "Included delivery costs are offset by avoided active-TB care under the working assumptions, but this depends materially on unresolved local programme resources and inherited epidemiological calibration semantics.",
        "The immediate next step is for SA Health to enter local setup, running, travel, outreach and workforce costs in the Streamlit Health Economics workspace and review possible bundling with test and regimen costs.",
    ]


def references_text() -> list[str]:
    return [
        "Dale KD, Abayawardana MJ, McBryde ES, Trauer JM, Carvalho N. Modeling the "
        "Cost-Effectiveness of Latent Tuberculosis Screening and Treatment Strategies in Recent "
        "Migrants to a Low-Incidence Setting. Am J Epidemiol. 2022;191(2):255-270. "
        "doi:10.1093/aje/kwab150.",
        "Repository evidence registry snapshot: data/apy_assumption_evidence_registry.csv, hash recorded in the report manifest.",
        "Frozen APY reference validation artifacts: validation/matlab_reference/prevent_igra_3hp_N1500_seed1.",
    ]


def flatten_table(rows: list[dict[str, Any]]) -> list[str]:
    out = []
    for item in rows:
        out.append("; ".join(f"{key}: {value}" for key, value in item.items()))
    return out


def active_tb_relative_reduction(data: dict[str, Any]) -> float:
    executive = data["executiveByMetric"]
    return num(executive["Active TB cases averted"]["value"]) / num(executive["Comparator active TB cases"]["value"])


def row(measure: str, value: Any, interpretation: str) -> dict[str, Any]:
    return {"Measure": measure, "Value": value, "Interpretation": interpretation}


def net_saving_label(item: dict[str, Any]) -> str:
    value = num(item["value"])
    label = item["label"]
    if value < 0:
        return label.replace("cost or saving", "saving")
    return label.replace("cost or saving", "cost")


def cost_input_value(value: Any, unit: Any, inclusion: str) -> str:
    if inclusion in {"bundled", "excluded"}:
        return plain_status("", inclusion)
    if text(value) == "":
        return "Unresolved"
    if "cost" in text(unit).lower() or text(unit):
        return fmt_currency(value)
    return text(value)


def plain_status(status: Any, inclusion: Any) -> str:
    status_text = text(status)
    inclusion_text = text(inclusion)
    if inclusion_text == "bundled":
        return "Reviewed as bundled/no additional cost"
    if inclusion_text == "excluded":
        return "Reviewed exclusion"
    if status_text == "configured_reviewed":
        return "Reviewed working numerical assumption"
    if status_text == "model_derived_reviewed":
        return "Reviewed model-derived assumption"
    if status_text == "reviewed_exclusion":
        return "Reviewed exclusion"
    if "working" in status_text:
        return "Working assumption"
    if status_text == "unresolved":
        return "Unresolved"
    return status_text or "Unresolved"


def fmt_count(value: Any) -> str:
    return f"{num(value):,.0f}"


def fmt_decimal(value: Any, places: int = 1) -> str:
    if text(value) == "":
        return ""
    return f"{num(value):,.{places}f}"


def fmt_percent(value: Any) -> str:
    return f"{100.0 * num(value):.1f}%"


def fmt_currency(value: Any) -> str:
    if text(value) == "":
        return ""
    return f"AUD {num(value):,.0f}"


def fmt_currency_abs(value: Any) -> str:
    return f"AUD {abs(num(value)):,.0f}"


def fmt_saving_or_cost(value: Any) -> str:
    amount = num(value)
    if amount < 0:
        return f"AUD {abs(amount):,.0f} net saving"
    return f"AUD {amount:,.0f} net cost"


def fmt_signed_ratio(value: Any) -> str:
    amount = num(value)
    if amount < 0:
        return f"AUD {abs(amount):,.0f} saving"
    return f"AUD {amount:,.0f} cost"


def num(value: Any) -> float:
    if value is None or value == "":
        return float("nan")
    return float(value)


def text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return str(value)


def close(a: float, b: float, tolerance: float) -> bool:
    return abs(a - b) <= tolerance


def today_label() -> str:
    return datetime.now().strftime("%d %B %Y")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def run(text_value: str) -> str:
    return f"<w:r><w:t xml:space=\"preserve\">{escape(text(text_value))}</w:t></w:r>"


def table_properties() -> str:
    return (
        "<w:tblPr><w:tblStyle w:val=\"TableGrid\"/>"
        "<w:tblW w:w=\"9638\" w:type=\"dxa\"/>"
        "<w:tblLayout w:type=\"fixed\"/>"
        "<w:tblLook w:firstRow=\"1\" w:noHBand=\"0\" w:noVBand=\"1\"/></w:tblPr>"
    )


def table_cell_width(column_count: int) -> int:
    page_inner_width = 9638
    return max(900, int(page_inner_width / max(1, column_count)))


def table_cell(value: str, *, header: bool = False, width: int = 2400) -> str:
    shading = "<w:shd w:fill=\"D9EAF7\"/>" if header else ""
    bold_start = "<w:b/>" if header else ""
    return (
        "<w:tc><w:tcPr>"
        f"<w:tcW w:w=\"{width}\" w:type=\"dxa\"/>"
        f"{shading}</w:tcPr><w:p><w:r><w:rPr>{bold_start}</w:rPr>"
        f"<w:t xml:space=\"preserve\">{escape(text(value))}</w:t>"
        "</w:r></w:p></w:tc>"
    )


def image_xml(rid: str, name: str, width_emu: int, image_id: int) -> str:
    height_emu = int(width_emu * 0.62)
    return f"""
<w:p>
  <w:r>
    <w:drawing>
      <wp:inline distT="0" distB="0" distL="0" distR="0"
        xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing">
        <wp:extent cx="{width_emu}" cy="{height_emu}"/>
        <wp:docPr id="{image_id}" name="{escape(name)}"/>
        <a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">
          <a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">
            <pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">
              <pic:nvPicPr><pic:cNvPr id="{image_id}" name="{escape(name)}"/><pic:cNvPicPr/></pic:nvPicPr>
              <pic:blipFill><a:blip r:embed="{rid}"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill>
              <pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="{width_emu}" cy="{height_emu}"/></a:xfrm>
              <a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr>
            </pic:pic>
          </a:graphicData>
        </a:graphic>
      </wp:inline>
    </w:drawing>
  </w:r>
</w:p>
"""


def document_xml(body: list[str]) -> str:
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" '
        'xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing">'
        "<w:body>"
        + "".join(body)
        + section_properties()
        + "</w:body></w:document>"
    )


def section_properties() -> str:
    return (
        "<w:sectPr><w:pgSz w:w=\"11906\" w:h=\"16838\"/>"
        "<w:pgMar w:top=\"1134\" w:right=\"1134\" w:bottom=\"1134\" w:left=\"1134\" "
        "w:header=\"708\" w:footer=\"708\" w:gutter=\"0\"/>"
        "</w:sectPr>"
    )


def styles_xml() -> str:
    return """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
<w:style w:type="paragraph" w:default="1" w:styleId="Normal">
<w:name w:val="Normal"/><w:qFormat/><w:rPr><w:rFonts w:ascii="Aptos" w:hAnsi="Aptos"/><w:sz w:val="21"/></w:rPr>
<w:pPr><w:spacing w:after="120" w:line="276" w:lineRule="auto"/></w:pPr></w:style>
<w:style w:type="paragraph" w:styleId="Title"><w:name w:val="Title"/><w:qFormat/>
<w:rPr><w:rFonts w:ascii="Aptos Display" w:hAnsi="Aptos Display"/><w:b/><w:color w:val="12355B"/><w:sz w:val="40"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Subtitle"><w:name w:val="Subtitle"/><w:qFormat/>
<w:rPr><w:rFonts w:ascii="Aptos" w:hAnsi="Aptos"/><w:color w:val="4A5568"/><w:sz w:val="25"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Heading1"><w:name w:val="heading 1"/><w:basedOn w:val="Normal"/><w:qFormat/>
<w:pPr><w:keepNext/><w:spacing w:before="320" w:after="120"/></w:pPr>
<w:rPr><w:rFonts w:ascii="Aptos Display" w:hAnsi="Aptos Display"/><w:b/><w:color w:val="12355B"/><w:sz w:val="28"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Caption"><w:name w:val="Caption"/><w:basedOn w:val="Normal"/><w:qFormat/>
<w:rPr><w:i/><w:color w:val="555555"/><w:sz w:val="18"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="ListParagraph"><w:name w:val="List Paragraph"/><w:basedOn w:val="Normal"/><w:qFormat/>
<w:pPr><w:ind w:left="360"/></w:pPr></w:style>
<w:style w:type="character" w:styleId="Hyperlink"><w:name w:val="Hyperlink"/>
<w:rPr><w:color w:val="0563C1"/><w:u w:val="single"/></w:rPr></w:style>
<w:style w:type="table" w:styleId="TableGrid"><w:name w:val="Table Grid"/><w:qFormat/>
<w:tblPr><w:tblBorders><w:top w:val="single" w:sz="4" w:color="B8C2CC"/>
<w:left w:val="single" w:sz="4" w:color="B8C2CC"/><w:bottom w:val="single" w:sz="4" w:color="B8C2CC"/>
<w:right w:val="single" w:sz="4" w:color="B8C2CC"/><w:insideH w:val="single" w:sz="4" w:color="D9DEE3"/>
<w:insideV w:val="single" w:sz="4" w:color="D9DEE3"/></w:tblBorders></w:tblPr></w:style>
</w:styles>"""


def content_types(media: list[tuple[str, bytes]]) -> str:
    overrides = "".join(
        f'<Override PartName="/word/media/{name}" ContentType="image/png"/>'
        for name, _ in media
    )
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/word/document.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>'
        '<Override PartName="/word/styles.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>'
        '<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>'
        '<Override PartName="/docProps/app.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>'
        f"{overrides}</Types>"
    )


def package_relationships() -> str:
    return """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>"""


def document_relationships(relationships: list[dict[str, str]]) -> str:
    rows = []
    type_map = {
        "image": "http://schemas.openxmlformats.org/officeDocument/2006/relationships/image",
        "hyperlink": "http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink",
    }
    for item in relationships:
        mode = ' TargetMode="External"' if item["external"] == "1" else ""
        rows.append(
            f'<Relationship Id="{item["id"]}" Type="{type_map[item["type"]]}" '
            f'Target="{escape(item["target"])}"{mode}/>'
        )
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        + "".join(rows)
        + "</Relationships>"
    )


def core_properties() -> str:
    now = datetime.now(timezone.utc).isoformat()
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
        'xmlns:dc="http://purl.org/dc/elements/1.1/" '
        'xmlns:dcterms="http://purl.org/dc/terms/" '
        'xmlns:dcmitype="http://purl.org/dc/dcmitype/" '
        'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
        "<dc:title>Targeted LTBI screening in the APY Lands: working health-economic analysis</dc:title>"
        "<dc:creator>tb-community-risk reproducible report workflow</dc:creator>"
        f'<dcterms:created xsi:type="dcterms:W3CDTF">{now}</dcterms:created>'
        f'<dcterms:modified xsi:type="dcterms:W3CDTF">{now}</dcterms:modified>'
        "</cp:coreProperties>"
    )


def app_properties() -> str:
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" '
        'xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">'
        "<Application>tb-community-risk</Application></Properties>"
    )
