from __future__ import annotations

from copy import deepcopy
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Any

import pandas as pd

from engine.apy.evidence import assess_apy_reference_readiness, load_apy_evidence_registry
from engine.apy.event_ledger import EVENT_LEDGER_CONTRACT_VERSION
from engine.apy.event_ledger_economics import (
    HEALTH_ECONOMICS_CONTRACT_VERSION,
    PROGRAM_COMPONENTS,
    run_event_ledger_health_economics,
)
from engine.apy.distributional_validation import portable_config_from_reference
from engine.apy.expected_value import run_expected_value
from engine.apy.ltbi_state import enable_development_compatibility_mode, resolve_ltbi_state_assumptions
from engine.apy.reference_loader import load_reference_scenario_config, load_reference_summary
from engine.apy.runner import run_scenario_with_do_nothing
from engine.apy.working_defaults import (
    HIGHER_BURDEN_DALY_OUTCOME_PRESET,
    UNIFIED_WORKING_DEFAULT_PRESET_ID,
    UNIFIED_WORKING_DEFAULT_PRESET_VERSION,
    apply_outcome_preset,
    apply_program_setup_preset,
    build_unified_working_default_preset,
)


SA_HEALTH_REFERENCE_PACKAGE_VERSION = "sa_health_apy_reference_package_v1"
SA_HEALTH_REFERENCE_PACKAGE_ID = "sa_health_apy_matlab_v9_compatible_working_reference_igra_3hp_prevent_30pct"
SA_HEALTH_TECHNICAL_RECENT_REMOTE_PACKAGE_ID = "sa_health_apy_technical_recent_remote_igra_3hp_prevent_30pct"
SA_HEALTH_REFERENCE_OUTPUT_DIR = (
    Path("outputs")
    / "sa_health_reference"
    / SA_HEALTH_REFERENCE_PACKAGE_ID
)
MATLAB_REFERENCE_DIR = Path("validation") / "matlab_reference" / "prevent_igra_3hp_N1500_seed1"
REFERENCE_LABEL = "SA Health APY MATLAB-v9-compatible working-reference analysis"
TECHNICAL_RECENT_REMOTE_LABEL = "SA Health APY technical explicit recent/remote calibration scenario"
PROVISIONAL_LABEL = "working-default analysis for planning - provisional"
PROGRAMME_NOT_COSTED_LABEL = "not yet locally costed; zero compatibility value used"
MATLAB_COMPATIBILITY_SEMANTICS = "matlab_v9_implicit_early_late"


def build_sa_health_reference_package(
    *,
    n_reps: int = 2000,
    seed: int = 1,
    include_technical_recent_remote: bool = False,
) -> dict[str, Any]:
    """Build the MATLAB-v9-compatible SA Health APY working-reference package.

    This is an export layer over the Python stochastic event ledger and
    health-economics engine using MATLAB v9's implicit early/late natural
    history semantics for compatibility with the earlier plain SA Health
    report. It does not treat the inherited calibration target or implicit
    early phase as reviewed APY evidence.
    """
    preset = build_unified_working_default_preset()
    config = _matlab_compatible_reference_config(
        preset["config"],
        n_reps=n_reps,
        seed=seed,
    )
    economics_config = _reference_economics_config(preset["economicsConfig"])
    epi_bundle = run_scenario_with_do_nothing(config)
    epi = epi_bundle["results"]
    ledger = epi["eventLedger"]
    economics = run_event_ledger_health_economics(epi, economics_config)
    readiness = assess_apy_reference_readiness(
        config,
        economics_config,
        economics_config.get("assumptionEvidenceRegistry"),
    )

    scenarios = _economic_scenarios(epi, economics_config)
    technical_recent_remote = (
        build_sa_health_technical_recent_remote_package()
        if include_technical_recent_remote
        else technical_recent_remote_metadata()
    )
    validation = matlab_reference_validation_rows(epi)
    tables = {
        "executive_summary": executive_summary_rows(config, ledger, economics, readiness),
        "cost_categories": cost_category_rows(economics, scenario_label="Primary working reference"),
        "annual_budget_impact": annual_budget_impact_rows(economics, scenario_label="Primary working reference"),
        "gross_delivery_ratios": gross_delivery_ratio_rows(ledger, economics),
        "net_health_system_ratios": net_health_system_ratio_rows(ledger, economics),
        "assumptions_readiness": assumptions_readiness_rows(config, economics_config, readiness),
        "economic_scenarios": economic_scenario_rows(scenarios),
        "matlab_reference_validation": validation,
    }
    figures = {
        "screening_treatment_cascade.svg": cascade_svg(tables["executive_summary"]),
        "annual_budget_impact.svg": annual_budget_svg(tables["annual_budget_impact"]),
        "cost_category_comparison.svg": cost_category_svg(tables["cost_categories"]),
    }
    manifest = manifest_payload(
        config=config,
        economics_config=economics_config,
        economics=economics,
        readiness=readiness,
        scenarios=scenarios,
        tables=tables,
        figures=figures,
    )
    return {
        "packageVersion": SA_HEALTH_REFERENCE_PACKAGE_VERSION,
        "packageId": SA_HEALTH_REFERENCE_PACKAGE_ID,
        "config": config,
        "economicsConfig": economics_config,
        "epidemiology": epi,
        "eventLedger": ledger,
        "economics": economics,
        "readiness": readiness,
        "economicScenarios": scenarios,
        "technicalRecentRemoteScenario": technical_recent_remote,
        "matlabReferenceValidation": validation,
        "tables": tables,
        "figures": figures,
        "manifest": manifest,
    }


def build_sa_health_technical_recent_remote_package() -> dict[str, Any]:
    """Build the superseded deterministic recent/remote technical scenario."""
    preset = build_unified_working_default_preset()
    config = _technical_recent_remote_config(preset["config"])
    economics_config = _reference_economics_config(preset["economicsConfig"])
    epi = run_expected_value(config)
    economics = run_event_ledger_health_economics(epi, economics_config)
    return {
        "packageId": SA_HEALTH_TECHNICAL_RECENT_REMOTE_PACKAGE_ID,
        "label": TECHNICAL_RECENT_REMOTE_LABEL,
        "primarySAHealthAnchor": False,
        "suitableAsPrimaryEstimate": False,
        "reason": (
            "This scenario uses the explicit recent/remote Markov structure with "
            "the 0% recent-LTBI compatibility placeholder. It is preserved for "
            "technical comparison and is not the primary SA Health epidemiological anchor."
        ),
        "config": config,
        "epidemiology": epi,
        "economics": economics,
    }


def technical_recent_remote_metadata() -> dict[str, Any]:
    return {
        "packageId": SA_HEALTH_TECHNICAL_RECENT_REMOTE_PACKAGE_ID,
        "label": TECHNICAL_RECENT_REMOTE_LABEL,
        "primarySAHealthAnchor": False,
        "suitableAsPrimaryEstimate": False,
        "reason": (
            "The deterministic explicit recent/remote package generated before "
            "the epidemiological-anchor correction is preserved under the prior "
            "output directory and can be rebuilt explicitly for technical review."
        ),
    }


def build_same_ledger_economic_scenario_comparison(
    results_or_bundle: dict[str, Any],
    economics_config: dict[str, Any],
) -> dict[str, Any]:
    """Compare economic scenarios using one existing APY event ledger."""
    scenarios = _economic_scenarios(results_or_bundle, economics_config)
    return {
        "comparisonType": "same_event_ledger_economic_scenarios",
        "epidemiologyRerunRequired": False,
        "scenarios": scenarios,
        "rows": economic_scenario_rows(scenarios),
        "notes": (
            "Economic assumptions are varied against the current event ledger. "
            "Epidemiological event counts are not recalculated."
        ),
    }


def write_sa_health_reference_package(
    output_dir: str | Path | None = None,
    *,
    n_reps: int = 2000,
    seed: int = 1,
) -> dict[str, Any]:
    from app.results_workbook import build_results_workbook

    package = build_sa_health_reference_package(n_reps=n_reps, seed=seed)
    root = Path(output_dir) if output_dir is not None else SA_HEALTH_REFERENCE_OUTPUT_DIR
    tables_dir = root / "tables"
    figures_dir = root / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    _write_json(root / "scenario_config.json", package["config"])
    _write_json(root / "economics_config.json", package["economicsConfig"])
    _write_json(root / "manifest.json", package["manifest"])
    _write_json(root / "technical_recent_remote_scenario.json", package["technicalRecentRemoteScenario"])
    _write_csv(root / "evidence_registry_snapshot.csv", load_apy_evidence_registry())
    _write_csv(root / "event_ledger_totals.csv", _table_rows(package["eventLedger"]["replicateTotals"]))
    _write_csv(root / "event_ledger_annual.csv", _table_rows(package["eventLedger"]["annualEvents"]))
    _write_csv(root / "economic_replicates.csv", _table_rows(package["economics"]["replicateResults"]))
    _write_csv(root / "economic_annual_by_arm.csv", _table_rows(package["economics"]["annualByArm"]))
    _write_csv(root / "economic_summary.csv", _table_rows(package["economics"]["summaries"]))
    for name, rows in package["tables"].items():
        _write_csv(tables_dir / f"{name}.csv", rows)
    for name, svg in package["figures"].items():
        (figures_dir / name).write_text(svg, encoding="utf-8")

    workbook_bundle = {
        "metadata": {
            "scenarioLabel": REFERENCE_LABEL,
            "modelVersion": package["epidemiology"].get("modelVersion"),
            "backend": package["epidemiology"].get("backend"),
        },
        "headline": {},
        "technical": {
            "eventLedger": {
                "metadata": package["eventLedger"].get("metadata"),
                "definitions": package["eventLedger"].get("definitions"),
                "replicateTotals": package["eventLedger"].get("replicateTotals"),
                "annualEvents": [],
                "validation": package["eventLedger"].get("validation"),
            }
        },
        "downloads": {
            "eventLedgerTotals": _table_rows(package["eventLedger"]["replicateTotals"]),
            "eventLedgerAnnual": [],
            "eventDefinitions": _table_rows(package["eventLedger"]["definitions"]),
        },
    }
    workbook = build_results_workbook(
        config=package["config"],
        bundle=workbook_bundle,
        backend_status={"name": "python_apy", "experimental": False},
        economics_results=package["economics"],
        economics_config=package["economicsConfig"],
    )
    (root / "sa_health_apy_working_reference_outputs.xlsx").write_bytes(workbook)
    (root / "README.md").write_text(reproduction_readme(package), encoding="utf-8")
    return {"outputDir": str(root), "manifest": package["manifest"], "package": package}


def _matlab_compatible_reference_config(
    config: dict[str, Any],
    *,
    n_reps: int,
    seed: int,
) -> dict[str, Any]:
    reference_cfg = load_reference_scenario_config(MATLAB_REFERENCE_DIR / "scenario_config.json")
    cfg = portable_config_from_reference(reference_cfg)
    cfg.update(
        {
            "scenarioLabel": REFERENCE_LABEL,
            "analysisMethod": "agent_based",
            "N": 1500,
            "nReps": int(n_reps),
            "seed": int(seed),
            "screenCoverage": 0.30,
            "screeningStrategy": "prevent",
            "testType": "IGRA",
            "regimen": "3HP",
            "screeningWindowYears": 2,
            "screenWindow": 2,
            "followUpHorizonYears": 20,
            "followHorizon": 20,
            "naturalHistorySemantics": MATLAB_COMPATIBILITY_SEMANTICS,
            "eligible": {"proportion": 1.0, "number": 1500, "description": "All APY demonstration residents eligible"},
        }
    )
    return cfg


def _reference_config(config: dict[str, Any]) -> dict[str, Any]:
    return _technical_recent_remote_config(config)


def _technical_recent_remote_config(config: dict[str, Any]) -> dict[str, Any]:
    cfg = deepcopy(config)
    cfg.update(
        {
            "scenarioLabel": TECHNICAL_RECENT_REMOTE_LABEL,
            "analysisMethod": "expected_value",
            "N": 1500,
            "screenCoverage": 0.30,
            "screeningStrategy": "prevent",
            "testType": "IGRA",
            "regimen": "3HP",
            "screeningWindowYears": 3,
            "screenWindow": 3,
            "followUpHorizonYears": 20,
            "followHorizon": 20,
            "eligible": {"proportion": 1.0, "number": 1500, "description": "All APY demonstration residents eligible"},
        }
    )
    return enable_development_compatibility_mode(cfg)


def _reference_economics_config(economics_config: dict[str, Any]) -> dict[str, Any]:
    econ = deepcopy(economics_config)
    econ["scenarioConfig"] = {
        "testType": "IGRA",
        "regimen": "3HP",
        "screeningStrategy": "prevent",
        "screenCoverage": 0.30,
        "screeningWindowYears": 2,
        "followUpHorizonYears": 20,
    }
    econ.setdefault("metadata", {}).update(
        {
            "scenarioLabel": REFERENCE_LABEL,
            "perspective": "Australian health-care system",
            "targetCurrency": "AUD",
            "targetPriceYear": "2019",
            "currencyCode": "AUD",
            "priceYear": "2019",
            "referencePackageId": SA_HEALTH_REFERENCE_PACKAGE_ID,
            "referencePackageVersion": SA_HEALTH_REFERENCE_PACKAGE_VERSION,
            "interpretation": PROVISIONAL_LABEL,
        }
    )
    econ.setdefault("discounting", {})["profiles"] = {
        "primary": {"costRate": 0.03, "healthRate": 0.03},
        "comparison": {"costRate": 0.0, "healthRate": 0.0},
    }
    econ.setdefault("threshold", {}).update(
        {"value": None, "currency": "AUD", "referenceYear": None, "source": "", "status": "unresolved"}
    )
    return econ


def _economic_scenarios(epi: dict[str, Any], base_econ: dict[str, Any]) -> list[dict[str, Any]]:
    definitions = [
        ("primary_working_reference", "Primary working reference", base_econ),
        (
            "illustrative_500k_setup",
            "Illustrative AUD 500,000 setup scenario",
            _illustrative_setup_config(base_econ),
        ),
        (
            "pathway_components_bundled",
            "Sensitivity: return/results, clinical review and work-up bundled/excluded",
            _bundle_pathway_components(base_econ),
        ),
        (
            "higher_burden_post_tb_daly",
            "Sensitivity: higher TB burden with post-TB DALYs",
            apply_outcome_preset(base_econ, HIGHER_BURDEN_DALY_OUTCOME_PRESET),
        ),
    ]
    out = []
    for scenario_id, label, econ_config in definitions:
        result = run_event_ledger_health_economics(epi, econ_config)
        out.append(
            {
                "scenarioId": scenario_id,
                "label": label,
                "economicsConfig": econ_config,
                "economics": result,
                "usesSameEventLedger": True,
            }
        )
    return out


def _illustrative_setup_config(economics_config: dict[str, Any]) -> dict[str, Any]:
    try:
        return apply_program_setup_preset(economics_config, "new_programme_implementation")
    except ValueError:
        econ = deepcopy(economics_config)
        for item in econ.get("costItems") or []:
            if item.get("costItemId") == "program_setup":
                item.update(
                    {
                        "originalCost": 500000.0,
                        "originalCurrency": "AUD",
                        "originalPriceYear": "2019",
                        "targetCurrency": "AUD",
                        "targetPriceYear": "2019",
                        "sourceCitation": "Illustrative SA Health setup-cost scenario",
                        "pageTableItem": "user-defined illustrative setup scenario",
                        "sourceYearStatus": "explicit",
                        "notes": "Illustrative only; not evidence-backed and not the primary reference calculation.",
                    }
                )
                item.setdefault("resourceUse", {})["costBasis"] = "total_once_at_program_start"
        econ.setdefault("metadata", {})["programSetupPreset"] = "illustrative_500000_setup"
        return econ


def _bundle_pathway_components(economics_config: dict[str, Any]) -> dict[str, Any]:
    econ = deepcopy(economics_config)
    notes = {
        "return_for_results": "Sensitivity: return/results cost treated as bundled into the selected test cost.",
        "clinical_review": "Sensitivity: clinical review cost treated as bundled into the selected regimen cost.",
        "active_tb_exclusion_workup": "Sensitivity: active-TB exclusion/work-up cost treated as bundled into the positive-test pathway.",
        "travel_outreach_staff_support": "Sensitivity: travel/outreach/staff support remains not locally costed.",
    }
    for item in econ.get("costItems") or []:
        if item.get("costItemId") in notes:
            item["originalCost"] = 0.0
            item["notes"] = notes[item["costItemId"]]
            item["sourceCitation"] = "SA Health pathway-cost bundling sensitivity"
            item["pageTableItem"] = "bundled/excluded sensitivity"
            item["originalPriceYear"] = "2019"
            item["targetPriceYear"] = "2019"
            item["originalCurrency"] = "AUD"
            item["targetCurrency"] = "AUD"
    for row in econ.get("assumptionEvidenceRegistry") or []:
        mapping = {
            "cost.return_for_results": ("bundled", "cost.test_igra"),
            "cost.clinical_review": ("bundled", "cost.regimen_3hp"),
            "cost.active_tb_exclusion_workup": ("bundled", "cost.regimen_3hp"),
            "cost.travel_outreach_staff_support": ("excluded", ""),
        }
        if row.get("assumptionId") in mapping:
            status, bundled_into = mapping[row["assumptionId"]]
            row.update(
                {
                    "currentValue": 0.0,
                    "reviewStatus": "reviewed_exclusion",
                    "provisional": False,
                    "inclusionStatus": status,
                    "bundledIntoAssumptionId": bundled_into,
                    "notes": notes[row["assumptionId"].split(".", 1)[1]],
                    "unresolvedReason": "",
                }
            )
    econ.setdefault("metadata", {})["economicScenario"] = "pathway_components_bundled"
    return econ


def executive_summary_rows(
    config: dict[str, Any],
    ledger: dict[str, Any],
    economics: dict[str, Any],
    readiness: dict[str, Any],
) -> list[dict[str, Any]]:
    totals = _totals_by_arm(ledger)
    comp = totals.get("comparator", {})
    inter = totals.get("intervention", {})
    primary = _replicate_row(economics, "primary")
    comparison = _replicate_row(economics, "comparison")
    gross_delivery = _gross_delivery_expenditure(primary)
    ltbi = resolve_ltbi_state_assumptions(config)
    rows = [
        _summary("Population", "Eligible population", inter.get("eligible_population"), "people", "All 1,500 APY demonstration residents are eligible."),
        _summary("Population", "People screened", inter.get("screened"), "people", "30% coverage over two years."),
        _summary("Screening cascade", "True-positive latent results", inter.get("true_positive_latent"), "people", ""),
        _summary("Screening cascade", "False-positive results", inter.get("false_positive"), "people", ""),
        _summary("Treatment cascade", "Preventive treatments initiated", inter.get("tpt_started_total"), "people", ""),
        _summary("Treatment cascade", "Preventive treatments completed", inter.get("tpt_completed_total"), "people", ""),
        _summary("Treatment cascade", "ADR-related treatment stops", inter.get("tpt_adr_stop_total"), "people", ""),
        _summary("TB outcomes", "Comparator active TB cases", comp.get("active_tb_cases"), "cases", "Current practice without additional systematic LTBI screening."),
        _summary("TB outcomes", "Intervention active TB cases", inter.get("active_tb_cases"), "cases", "Targeted IGRA + 3HP screening and treatment."),
        _summary("TB outcomes", "Active TB cases averted", inter.get("active_tb_cases_prevented"), "cases", "Direct effects only; no transmission-mediated effects."),
        _summary("Costs, 3% discounted", "Comparator total cost", primary.get("comparatorCost"), "AUD 2019", ""),
        _summary("Costs, 3% discounted", "Intervention total cost", primary.get("interventionCost"), "AUD 2019", ""),
        _summary("Costs, 3% discounted", "Incremental cost", primary.get("incrementalCost"), "AUD 2019", "Intervention minus comparator."),
        _summary("Costs, 3% discounted", "Active-TB treatment cost offset", _tb_cost_offset(primary), "AUD 2019", "Comparator active-TB care cost minus intervention active-TB care cost."),
        _summary(
            "Gross delivery measures",
            "Intervention delivery expenditure per person screened",
            _divide(gross_delivery, inter.get("screened")),
            "AUD 2019 per person screened",
            "Programme/pathway expenditure before active-TB care offsets.",
        ),
        _summary(
            "Gross delivery measures",
            "Intervention delivery expenditure per TPT start",
            _divide(gross_delivery, inter.get("tpt_started_total")),
            "AUD 2019 per treatment start",
            "Programme/pathway expenditure before active-TB care offsets.",
        ),
        _summary(
            "Gross delivery measures",
            "Intervention delivery expenditure per TPT completion",
            _divide(gross_delivery, inter.get("tpt_completed_total")),
            "AUD 2019 per completed course",
            "Programme/pathway expenditure before active-TB care offsets.",
        ),
        _summary(
            "Gross delivery measures",
            "Intervention delivery expenditure per active TB case averted",
            _divide(gross_delivery, inter.get("active_tb_cases_prevented")),
            "AUD 2019 per case averted",
            "Programme/pathway expenditure before active-TB care offsets.",
        ),
        _summary(
            "Net health-system measures",
            "Net health-system cost or saving per person screened",
            _divide(primary.get("incrementalCost"), inter.get("screened")),
            "AUD 2019 per person screened",
            "Net incremental cost after active-TB care offsets.",
        ),
        _summary(
            "Net health-system measures",
            "Net health-system cost or saving per TPT completion",
            _divide(primary.get("incrementalCost"), inter.get("tpt_completed_total")),
            "AUD 2019 per completed course",
            "Net incremental cost after active-TB care offsets.",
        ),
        _summary(
            "Net health-system measures",
            "Net health-system cost or saving per active TB case averted",
            primary.get("costPerActiveTBCasePrevented"),
            "AUD 2019 per case averted",
            "Net incremental cost after active-TB care offsets.",
        ),
        _summary("Secondary DALY/ICER", "DALYs averted, 3% discounted", primary.get("dalysAverted"), "DALYs", "Provisional secondary output; primary calculation excludes post-TB sequelae."),
        _summary("Secondary DALY/ICER", "ICER classification, 3% discounted", primary.get("classification"), "", "No cost-effectiveness conclusion is permitted while evidence remains provisional and threshold is unresolved."),
        _summary("Secondary DALY/ICER", "Numerical ICER, 3% discounted", primary.get("replicateICER"), "AUD 2019 per DALY averted", "Blank when the quadrant is dominant/dominated or otherwise not interpretable."),
        _summary("Secondary DALY/ICER", "NMB", primary.get("netMonetaryBenefit"), "AUD 2019", "Unavailable because no reviewed WTP threshold is supplied."),
        _summary("Comparison, 0% undiscounted", "Incremental cost, 0% discounting", comparison.get("incrementalCost"), "AUD 2019", ""),
        _summary("Evidence/readiness", "Epidemiological anchor", "MATLAB-v9-compatible stochastic APY reference", "", "Compatible with the frozen earlier plain SA Health APY v9 report."),
        _summary("Evidence/readiness", "Analysis status", PROVISIONAL_LABEL, "", "Inherited calibration semantics and evidence gaps remain provisional."),
        _summary(
            "Evidence/readiness",
            "Baseline recent-LTBI proportion used",
            ltbi.get("baselineRecentLTBIProportion"),
            "proportion",
            "Compatibility metadata only; the primary anchor uses MATLAB v9 "
            "implicit early/late semantics rather than measured recent-LTBI composition.",
        ),
        _summary("Evidence/readiness", "Overall clinician-ready", readiness.get("overallClinicianReady"), "boolean", "False until unresolved evidence is reviewed."),
    ]
    return rows


def cost_category_rows(economics: dict[str, Any], *, scenario_label: str) -> list[dict[str, Any]]:
    primary = _replicate_row(economics, "primary")
    categories = [
        ("screening_test", "Screening test", "screeningTestCost"),
        ("return_results", "Return/results", "returnForResultsCost"),
        ("clinical_review", "Clinical review", "clinicalReviewCost"),
        ("active_tb_exclusion_workup", "Active-TB exclusion/work-up", "activeTBExclusionWorkupCost"),
        ("preventive_regimen", "Preventive regimen", "tptRegimenCost"),
        ("adr_management", "ADR management", "adrManagementCost"),
        ("program_setup", "Programme setup", "programSetupCost"),
        ("program_running", "Programme running", "programRunningCost"),
        ("travel_outreach_staff", "Travel/outreach/staff", "travelOutreachStaffSupportCost"),
        ("active_tb_care", "Active-TB care", "activeTBDiseaseCost"),
        ("total", "Total", None),
    ]
    rows = []
    for category_id, label, component in categories:
        if component is None:
            comp_value = primary.get("comparatorCost")
            inter_value = primary.get("interventionCost")
            comp_undiscounted = primary.get("comparator_totalArmCostUndiscounted")
            inter_undiscounted = primary.get("intervention_totalArmCostUndiscounted")
        else:
            comp_value = primary.get(f"comparator_{component}Discounted")
            inter_value = primary.get(f"intervention_{component}Discounted")
            comp_undiscounted = primary.get(f"comparator_{component}Undiscounted")
            inter_undiscounted = primary.get(f"intervention_{component}Undiscounted")
        rows.append(
            {
                "scenario": scenario_label,
                "categoryId": category_id,
                "category": label,
                "comparatorDiscountedCost": comp_value,
                "interventionDiscountedCost": inter_value,
                "incrementalDiscountedCost": _subtract(inter_value, comp_value),
                "comparatorUndiscountedCost": comp_undiscounted,
                "interventionUndiscountedCost": inter_undiscounted,
                "incrementalUndiscountedCost": _subtract(inter_undiscounted, comp_undiscounted),
                "interpretation": _category_interpretation(category_id),
            }
        )
    return rows


def gross_delivery_ratio_rows(ledger: dict[str, Any], economics: dict[str, Any]) -> list[dict[str, Any]]:
    totals = _totals_by_arm(ledger)
    inter = totals.get("intervention", {})
    primary = _replicate_row(economics, "primary")
    gross = _gross_delivery_expenditure(primary)
    denominators = [
        ("person_screened", "Intervention delivery expenditure per person screened", inter.get("screened")),
        ("tpt_start", "Intervention delivery expenditure per TPT start", inter.get("tpt_started_total")),
        ("tpt_completion", "Intervention delivery expenditure per TPT completion", inter.get("tpt_completed_total")),
        ("infection_effectively_treated", "Intervention delivery expenditure per infection effectively treated", inter.get("infection_effectively_treated_total")),
        ("active_tb_case_averted", "Intervention delivery expenditure per active TB case averted", inter.get("active_tb_cases_prevented")),
    ]
    return [
        {
            "ratioId": ratio_id,
            "label": label,
            "numerator": gross,
            "denominator": denominator,
            "value": _divide(gross, denominator),
            "unit": "AUD 2019",
            "interpretation": "Gross programme/pathway expenditure before active-TB care offsets.",
        }
        for ratio_id, label, denominator in denominators
    ]


def net_health_system_ratio_rows(ledger: dict[str, Any], economics: dict[str, Any]) -> list[dict[str, Any]]:
    totals = _totals_by_arm(ledger)
    inter = totals.get("intervention", {})
    primary = _replicate_row(economics, "primary")
    net = primary.get("incrementalCost")
    denominators = [
        ("person_screened", "Net health-system cost or saving per person screened", inter.get("screened")),
        ("tpt_completion", "Net health-system cost or saving per TPT completion", inter.get("tpt_completed_total")),
        ("active_tb_case_averted", "Net health-system cost or saving per active TB case averted", inter.get("active_tb_cases_prevented")),
        ("daly_averted", "Net incremental cost per DALY averted", primary.get("dalysAverted")),
    ]
    return [
        {
            "ratioId": ratio_id,
            "label": label,
            "numerator": net,
            "denominator": denominator,
            "value": _divide(net, denominator),
            "unit": "AUD 2019",
            "classification": primary.get("classification"),
            "interpretation": "Net incremental health-system cost after active-TB care offsets.",
        }
        for ratio_id, label, denominator in denominators
    ]


def annual_budget_impact_rows(economics: dict[str, Any], *, scenario_label: str) -> list[dict[str, Any]]:
    annual = economics["annualByArm"]
    if not isinstance(annual, pd.DataFrame):
        annual = pd.DataFrame(annual)
    primary = annual[annual["discountProfile"] == "primary"].copy()
    rows = []
    cumulative_discounted = 0.0
    cumulative_undiscounted = 0.0
    for year in sorted(primary["modelYear"].dropna().unique()):
        by_year = primary[primary["modelYear"] == year]
        comp = by_year[by_year["arm"] == "comparator"]
        inter = by_year[by_year["arm"] == "intervention"]
        comp_total = _strict_numeric_sum(comp["totalDiscountedCost"])
        inter_total = _strict_numeric_sum(inter["totalDiscountedCost"])
        comp_total_u = _strict_numeric_sum(comp["totalUndiscountedCost"])
        inter_total_u = _strict_numeric_sum(inter["totalUndiscountedCost"])
        annual_inc = _subtract(inter_total, comp_total)
        annual_inc_u = _subtract(inter_total_u, comp_total_u)
        if annual_inc is not None:
            cumulative_discounted += annual_inc
        if annual_inc_u is not None:
            cumulative_undiscounted += annual_inc_u
        rows.append(
            {
                "scenario": scenario_label,
                "modelYear": int(year),
                "comparatorActiveTBCareDiscounted": _strict_numeric_sum(comp["activeTBDiseaseCost"] * comp["costDiscountFactor"]),
                "interventionProgramPathwayDiscounted": _annual_programme_cost(inter, discounted=True),
                "interventionActiveTBCareDiscounted": _strict_numeric_sum(inter["activeTBDiseaseCost"] * inter["costDiscountFactor"]),
                "comparatorTotalCostDiscounted": comp_total,
                "interventionTotalCostDiscounted": inter_total,
                "annualIncrementalCostDiscounted": annual_inc,
                "cumulativeIncrementalCostDiscounted": cumulative_discounted,
                "comparatorTotalCostUndiscounted": comp_total_u,
                "interventionTotalCostUndiscounted": inter_total_u,
                "annualIncrementalCostUndiscounted": annual_inc_u,
                "cumulativeIncrementalCostUndiscounted": cumulative_undiscounted,
            }
        )
    return rows


def _annual_programme_cost(frame: pd.DataFrame, *, discounted: bool) -> float | None:
    if frame.empty:
        return 0.0
    values = sum(frame[component] for component in PROGRAM_COMPONENTS if component in frame)
    if discounted:
        values = values * frame["costDiscountFactor"]
    return _strict_numeric_sum(values)


def assumptions_readiness_rows(config: dict[str, Any], economics_config: dict[str, Any], readiness: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for row in economics_config.get("assumptionEvidenceRegistry") or []:
        assumption_id = row.get("assumptionId")
        if not _is_report_relevant_assumption(assumption_id):
            continue
        rows.append(
            {
                "assumptionId": assumption_id,
                "description": row.get("description"),
                "value": row.get("currentValue"),
                "unit": row.get("unit"),
                "sourceProvenance": row.get("sourceCitation") or row.get("sourceLocation"),
                "reviewStatus": row.get("reviewStatus"),
                "provisional": row.get("provisional"),
                "inclusionStatus": row.get("inclusionStatus"),
                "effectOnInterpretation": _assumption_effect(row),
            }
        )
    ltbi = resolve_ltbi_state_assumptions(config)
    rows.append(
        {
            "assumptionId": "epi.baseline_recent_ltbi_proportion",
            "description": "Baseline recent-LTBI proportion",
            "value": ltbi.get("baselineRecentLTBIProportion"),
            "unit": "proportion",
            "sourceProvenance": ltbi.get("baselineRecentLTBIProportionSource") or ltbi.get("source"),
            "reviewStatus": ltbi.get("status"),
            "provisional": ltbi.get("provisional"),
            "inclusionStatus": "unresolved compatibility metadata",
            "effectOnInterpretation": "Blocks clinician-ready interpretation; primary anchor uses MATLAB v9 implicit early/late semantics rather than a reviewed recent-LTBI fraction.",
        }
    )
    rows.append(
        {
            "assumptionId": "threshold.gdp_per_capita",
            "description": "Willingness-to-pay threshold",
            "value": (economics_config.get("threshold") or {}).get("value"),
            "unit": "AUD 2019 per DALY averted",
            "sourceProvenance": (economics_config.get("threshold") or {}).get("source"),
            "reviewStatus": (economics_config.get("threshold") or {}).get("status", "unresolved"),
            "provisional": True,
            "inclusionStatus": "unresolved",
            "effectOnInterpretation": "NMB and probability cost-effective unavailable.",
        }
    )
    return rows


def economic_scenario_rows(scenarios: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for scenario in scenarios:
        primary = _replicate_row(scenario["economics"], "primary")
        rows.append(
            {
                "scenarioId": scenario["scenarioId"],
                "label": scenario["label"],
                "usesSameEventLedger": scenario["usesSameEventLedger"],
                "comparatorCost": primary.get("comparatorCost"),
                "interventionCost": primary.get("interventionCost"),
                "incrementalCost": primary.get("incrementalCost"),
                "dalysAverted": primary.get("dalysAverted"),
                "costPerActiveTBCasePrevented": primary.get("costPerActiveTBCasePrevented"),
                "classification": primary.get("classification"),
                "nmb": primary.get("netMonetaryBenefit"),
                "provisional": True,
            }
        )
    return rows


def matlab_reference_validation_rows(epi: dict[str, Any]) -> list[dict[str, Any]]:
    reference = load_reference_summary(MATLAB_REFERENCE_DIR / "matlab_summary.csv")
    reference_rows = {str(row["Metric"]): row for row in reference.to_dict(orient="records")}
    raw = epi.get("raw")
    if not isinstance(raw, pd.DataFrame):
        raw = pd.DataFrame(raw)
    metrics = [
        "nScreened",
        "nTestPositive",
        "nTestPositiveNonActive",
        "nFalsePositiveTests",
        "nStartTPT",
        "nCompleteTPT",
        "nADRstop",
        "nCuredInfection",
        "nPreventedActiveTB",
        "nActiveBy20y",
    ]
    rows = []
    for metric in metrics:
        ref = reference_rows.get(metric, {})
        py_median = _number_or_none(raw[metric].median()) if metric in raw else None
        low = _number_or_none(ref.get("Low95"))
        high = _number_or_none(ref.get("High95"))
        rows.append(
            {
                "metric": metric,
                "pythonMean": _number_or_none(raw[metric].mean()) if metric in raw else None,
                "pythonMedian": py_median,
                "matlabMedian": _number_or_none(ref.get("Median")),
                "matlabLow95": low,
                "matlabHigh95": high,
                "withinMatlabInterval": None if py_median is None or low is None or high is None else low <= py_median <= high,
                "notes": "Distributional validation against frozen MATLAB APY v9 reference; row-for-row RNG identity is not expected.",
            }
        )
    return rows


def manifest_payload(
    *,
    config: dict[str, Any],
    economics_config: dict[str, Any],
    economics: dict[str, Any],
    readiness: dict[str, Any],
    scenarios: list[dict[str, Any]],
    tables: dict[str, list[dict[str, Any]]],
    figures: dict[str, str],
) -> dict[str, Any]:
    return {
        "packageId": SA_HEALTH_REFERENCE_PACKAGE_ID,
        "packageVersion": SA_HEALTH_REFERENCE_PACKAGE_VERSION,
        "label": REFERENCE_LABEL,
        "analysisStatus": PROVISIONAL_LABEL,
        "codeCommit": _git_value(["git", "rev-parse", "HEAD"]),
        "branch": _git_value(["git", "branch", "--show-current"]),
        "modelType": "agent_based",
        "nReps": config.get("nReps"),
        "seed": config.get("seed"),
        "epidemiologicalAnchor": {
            "type": "matlab_v9_compatible_python_stochastic_event_ledger",
            "storedReferenceDirectory": str(MATLAB_REFERENCE_DIR),
            "naturalHistorySemantics": MATLAB_COMPATIBILITY_SEMANTICS,
            "interpretation": "Compatibility with frozen MATLAB APY v9 reference; not evidence that the inherited calibration target is scientifically validated.",
        },
        "technicalRecentRemoteScenario": {
            "packageId": SA_HEALTH_TECHNICAL_RECENT_REMOTE_PACKAGE_ID,
            "label": TECHNICAL_RECENT_REMOTE_LABEL,
            "suitableAsPrimaryEstimate": False,
        },
        "dynamicTransmissionIncluded": False,
        "eventLedgerContractVersion": EVENT_LEDGER_CONTRACT_VERSION,
        "healthEconomicsContractVersion": HEALTH_ECONOMICS_CONTRACT_VERSION,
        "workingDefaultPresetId": UNIFIED_WORKING_DEFAULT_PRESET_ID,
        "workingDefaultPresetVersion": UNIFIED_WORKING_DEFAULT_PRESET_VERSION,
        "configurationHash": _hash_payload(config),
        "economicsConfigurationHash": _hash_payload(economics_config),
        "evidenceRegistryHash": _hash_payload(load_apy_evidence_registry()),
        "economicScenarioIds": [scenario["scenarioId"] for scenario in scenarios],
        "readiness": {
            "overallClinicianReady": readiness.get("overallClinicianReady"),
            "epidemiologyReady": readiness.get("epidemiologyReady"),
            "costReady": readiness.get("costReady"),
            "dalyReady": readiness.get("dalyReady"),
            "thresholdReady": readiness.get("thresholdReady"),
        },
        "unresolvedInputs": economics.get("unresolvedInputs") or [],
        "warnings": economics.get("warnings") or [],
        "generatedFiles": {
            "scenarioConfig": "scenario_config.json",
            "economicsConfig": "economics_config.json",
            "manifest": "manifest.json",
            "workbook": "sa_health_apy_working_reference_outputs.xlsx",
            "tables": [f"tables/{name}.csv" for name in tables],
            "figures": [f"figures/{name}" for name in figures],
        },
        "interpretationGuardrails": [
            "Working-default analysis for planning; provisional.",
            "Primary epidemiology uses MATLAB-v9-compatible implicit early/late semantics for compatibility with the earlier plain SA Health report.",
            "The inherited 10/770 active-TB calibration target remains unresolved and must not be described as validated incident progression from LTBI.",
            "The implicit early phase is not measured APY recent-LTBI composition.",
            "Multiplicative OR-as-hazard and joint-risk assumptions are preserved for compatibility and remain scientifically provisional.",
            "Programme setup, running, travel, outreach and staff-support costs are not yet locally costed.",
            "NMB and probability cost-effective are unavailable without a reviewed threshold.",
            "No dynamic transmission effects are included.",
        ],
    }


def reproduction_readme(package: dict[str, Any]) -> str:
    return f"""# SA Health APY Working-Reference Package

Package: `{SA_HEALTH_REFERENCE_PACKAGE_ID}`

Run from the repository root:

```powershell
python scripts/build_sa_health_reference_package.py
```

This package is generated from the non-dynamic Python APY stochastic event
ledger and the event-ledger health-economics engine. The epidemiological
anchor uses MATLAB v9 compatible implicit early/late natural-history
semantics for compatibility with the frozen earlier plain SA Health report.
It excludes the dynamic community transmission model.

Interpretation: {PROVISIONAL_LABEL}. The inherited `10/770` active-TB target
is treated only as a MATLAB software calibration target for active TB observed
over the historical two-year screening window. Its scientific provenance and
interpretation remain unresolved. The implicit early phase is not measured
recent LTBI. Programme setup, programme running, travel, outreach and staff
support remain not locally costed. No willingness-to-pay threshold is supplied,
so NMB and probability cost-effective are unavailable.

Configuration hash: `{package["manifest"]["configurationHash"]}`

Economics configuration hash: `{package["manifest"]["economicsConfigurationHash"]}`

Evidence registry hash: `{package["manifest"]["evidenceRegistryHash"]}`
"""


def cascade_svg(rows: list[dict[str, Any]]) -> str:
    values = {
        row["metric"]: row.get("value")
        for row in rows
        if row.get("section") in {"Population", "Screening cascade", "Treatment cascade", "TB outcomes"}
    }
    labels = [
        ("People screened", values.get("People screened")),
        ("True-positive latent results", values.get("True-positive latent results")),
        ("False-positive results", values.get("False-positive results")),
        ("Preventive treatments initiated", values.get("Preventive treatments initiated")),
        ("Preventive treatments completed", values.get("Preventive treatments completed")),
        ("Active TB cases averted", values.get("Active TB cases averted")),
    ]
    return _bar_svg("Screening and treatment cascade", labels, width=820)


def annual_budget_svg(rows: list[dict[str, Any]]) -> str:
    labels = [(str(row["modelYear"]), row.get("annualIncrementalCostDiscounted")) for row in rows]
    return _bar_svg("Annual incremental budget impact, 3% discounted", labels, width=900)


def cost_category_svg(rows: list[dict[str, Any]]) -> str:
    labels = [
        (row["category"], row.get("incrementalDiscountedCost"))
        for row in rows
        if row["categoryId"] not in {"total"}
    ]
    return _bar_svg("Incremental cost by category, 3% discounted", labels, width=940)


def _bar_svg(title: str, labels: list[tuple[str, Any]], *, width: int) -> str:
    bar_h = 22
    gap = 14
    left = 260
    height = 70 + len(labels) * (bar_h + gap)
    numbers = [abs(float(value)) for _, value in labels if _number_or_none(value) is not None]
    max_value = max(numbers) if numbers else 1.0
    body = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="20" y="30" font-family="Arial" font-size="20" font-weight="700" fill="#1f2933">{_escape(title)}</text>',
    ]
    for i, (label, value) in enumerate(labels):
        y = 58 + i * (bar_h + gap)
        num = _number_or_none(value)
        bar_w = 0 if num is None else max(2, int((abs(num) / max_value) * (width - left - 80)))
        color = "#2f6f9f" if (num or 0) >= 0 else "#b65c3a"
        body.append(f'<text x="20" y="{y + 16}" font-family="Arial" font-size="13" fill="#1f2933">{_escape(label)}</text>')
        body.append(f'<rect x="{left}" y="{y}" width="{bar_w}" height="{bar_h}" fill="{color}"/>')
        text = "" if num is None else f"{num:,.1f}"
        body.append(f'<text x="{left + bar_w + 8}" y="{y + 16}" font-family="Arial" font-size="12" fill="#1f2933">{_escape(text)}</text>')
    body.append("</svg>")
    return "\n".join(body)


def _totals_by_arm(ledger: dict[str, Any]) -> dict[str, dict[str, Any]]:
    totals = ledger.get("replicateTotals")
    if not isinstance(totals, pd.DataFrame):
        totals = pd.DataFrame(totals)
    if {"arm", "eventName", "value"}.issubset(totals.columns):
        wide = totals.pivot_table(
            index=["arm"],
            columns="eventName",
            values="value",
            aggfunc="mean",
        ).reset_index()
        totals = wide
    elif "arm" in totals.columns:
        numeric = [col for col in totals.columns if col != "arm" and pd.api.types.is_numeric_dtype(totals[col])]
        totals = totals.groupby("arm", dropna=False)[numeric].mean().reset_index()
    out = {}
    for _, row in totals.iterrows():
        out[str(row.get("arm"))] = row.to_dict()
    return out


def _replicate_row(economics: dict[str, Any], profile: str) -> dict[str, Any]:
    reps = economics["replicateResults"]
    if not isinstance(reps, pd.DataFrame):
        reps = pd.DataFrame(reps)
    subset = reps[reps["discountProfile"] == profile]
    if subset.empty:
        return {}
    if len(subset) == 1:
        return subset.iloc[0].to_dict()
    if "economicPairComplete" in subset:
        complete = subset[subset["economicPairComplete"].astype(bool)].copy()
    else:
        complete = subset.copy()
    if complete.empty:
        return {"classification": "incomplete / not calculated"}
    out: dict[str, Any] = {"discountProfile": profile}
    for col in complete.columns:
        if col in {"replicateId", "pairedReplicateId", "replicateSeed"}:
            continue
        numeric = pd.to_numeric(complete[col], errors="coerce").dropna()
        if not numeric.empty:
            out[col] = _number_or_none(numeric.mean())
    out["incrementalCost"] = _subtract(out.get("interventionCost"), out.get("comparatorCost"))
    out["dalysAverted"] = _subtract(out.get("comparatorDALYs"), out.get("interventionDALYs"))
    out["costPerActiveTBCasePrevented"] = _divide(out.get("incrementalCost"), out.get("activeTBCasesPrevented"))
    out["costPerInfectionEffectivelyTreated"] = _divide(out.get("incrementalCost"), out.get("infectionsEffectivelyTreated"))
    summaries = economics.get("summaries")
    if not isinstance(summaries, pd.DataFrame):
        summaries = pd.DataFrame(summaries)
    primary = summaries[
        (summaries.get("discountProfile") == profile)
        & (summaries.get("metric") == "primaryICER_ratioOfMeans")
    ] if not summaries.empty else pd.DataFrame()
    if not primary.empty:
        row = primary.iloc[0].to_dict()
        out["classification"] = row.get("classification")
        out["replicateICER"] = row.get("mean")
    return out


def _summary(section: str, metric: str, value: Any, unit: str, interpretation: str) -> dict[str, Any]:
    return {
        "section": section,
        "metric": metric,
        "value": value,
        "unit": unit,
        "interpretation": interpretation,
        "status": PROVISIONAL_LABEL if section in {"Secondary DALY/ICER", "Evidence/readiness"} else "working default",
    }


def _category_interpretation(category_id: str) -> str:
    if category_id in {"program_setup", "program_running", "travel_outreach_staff"}:
        return PROGRAMME_NOT_COSTED_LABEL
    if category_id in {"return_results", "clinical_review", "active_tb_exclusion_workup"}:
        return "SA Health working pathway component shown separately; verify overlap with bundled test/regimen costs."
    if category_id == "active_tb_care":
        return "Ordinary active-TB care included independently in both arms."
    return "Included in authoritative economics when mapped event quantity is present."


def _gross_delivery_expenditure(row: dict[str, Any]) -> float | None:
    if not row:
        return None
    total = 0.0
    for component in PROGRAM_COMPONENTS:
        value = _number_or_none(row.get(f"intervention_{component}Discounted"))
        if value is None:
            return None
        total += value
    return total


def _assumption_effect(row: dict[str, Any]) -> str:
    assumption_id = str(row.get("assumptionId") or "")
    if assumption_id in {"cost.program_setup", "cost.program_running", "cost.travel_outreach_staff_support"}:
        return PROGRAMME_NOT_COSTED_LABEL
    if assumption_id == "threshold.gdp_per_capita":
        return "Blocks NMB only; does not block cost-consequence or ICER calculation."
    if assumption_id.startswith("daly."):
        return "Supports secondary provisional DALY/ICER outputs."
    return "Supports current cost-consequence calculation."


def _is_report_relevant_assumption(assumption_id: Any) -> bool:
    return str(assumption_id or "").startswith(("cost.", "daly.", "threshold."))


def _tb_cost_offset(row: dict[str, Any]) -> float | None:
    return _subtract(row.get("comparator_activeTBDiseaseCostDiscounted"), row.get("intervention_activeTBDiseaseCostDiscounted"))


def _strict_numeric_sum(values: Any) -> float | None:
    series = pd.to_numeric(pd.Series(values), errors="coerce")
    if series.isna().any():
        return None
    return float(series.sum())


def _subtract(a: Any, b: Any) -> float | None:
    a_num = _number_or_none(a)
    b_num = _number_or_none(b)
    if a_num is None or b_num is None:
        return None
    return a_num - b_num


def _divide(a: Any, b: Any) -> float | None:
    a_num = _number_or_none(a)
    b_num = _number_or_none(b)
    if a_num is None or b_num is None or abs(b_num) < 1e-12:
        return None
    return a_num / b_num


def _number_or_none(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, str) and value == "":
        return None
    if isinstance(value, list) and value == []:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if pd.notna(number) else None


def _table_rows(value: Any) -> list[dict[str, Any]]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if hasattr(value, "to_dict"):
        return value.to_dict(orient="records")
    return []


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(_json_safe(payload), indent=2, sort_keys=True), encoding="utf-8")


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    headers = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(_json_safe(row) for row in rows)


def _json_safe(value: Any) -> Any:
    if isinstance(value, pd.DataFrame):
        return value.to_dict(orient="records")
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_json_safe(v) for v in value]
    if isinstance(value, float) and (not pd.notna(value) or not math.isfinite(value)):
        return None
    if hasattr(value, "item"):
        try:
            return _json_safe(value.item())
        except Exception:
            pass
    return value


def _hash_payload(payload: Any) -> str:
    text = json.dumps(_json_safe(payload), sort_keys=True, default=str, separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _git_value(args: list[str]) -> str:
    try:
        return subprocess.check_output(args, cwd=Path(__file__).resolve().parents[2], text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return ""


def _escape(value: Any) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )
