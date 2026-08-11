from __future__ import annotations

import csv
from copy import deepcopy
import math
from pathlib import Path
from typing import Any

import pandas as pd

from engine.apy.costing import (
    TARGET_CURRENCY,
    TARGET_PRICE_YEAR,
    build_cost_item,
    normalise_cost_table,
    unresolved_cost_warnings,
    valid_converted_cost,
)
from engine.apy.runner import run_scenario_with_do_nothing
from engine.apy.scenario import DIRECT_EFFECTS_SCOPE_STATEMENT
from engine.apy.event_ledger_economics import (
    HEALTH_ECONOMICS_CONTRACT_VERSION,
    default_daly_assumptions,
    run_event_ledger_health_economics,
)


ECONOMICS_CONTRACT_VERSION = "apy_v9_economics_results_v1"
ECONOMICS_SOURCE = "run_health_economics_v9_python_port"
COMPATIBILITY_COST_NOTES = (
    "Compatibility mirror only. Authoritative Milestone 1 calculations use "
    "costItems and valid convertedTargetYearCost values."
)
LEGACY_COST_ITEM_MAP = {
    ("costs", "test", "IGRA"): "test_igra",
    ("costs", "test", "TST"): "test_tst",
    ("costs", "regimen", "x3HP"): "regimen_3hp",
    ("costs", "regimen", "x4R"): "regimen_4r",
    ("costs", "regimen", "x3HR"): "regimen_3hr",
    ("costs", "regimen", "x6H"): "regimen_6h",
    ("costs", "regimen", "x9H"): "regimen_9h",
    ("costs", "activeTBDiseasePerCase"): "active_tb_disease",
    ("costs", "falsePositiveIncrementalPerPerson"): "false_positive_incremental",
    ("costs", "programSetupTotal"): "program_setup",
    ("costs", "programRunningTotal"): "program_running",
}

DALE_2019_PRESET_NAME = "Dale 2019 AUD working defaults"
DALE_2019_REFERENCE_STATUS = "working_default_not_final_apy_evidence"
DALE_2019_CITATION = (
    "Dale KD, Abayawardana MJ, McBryde ES, Trauer JM, Carvalho N. "
    "Modeling the Cost-Effectiveness of Latent Tuberculosis Screening and "
    "Treatment Strategies in Recent Migrants to a Low-Incidence Setting. "
    "Am J Epidemiol. 2022;191(2):255-270. doi:10.1093/aje/kwab150."
)
DALE_2019_SOURCE_LOCATION = "Dale et al. 2022 Web Appendix 3 and Web Tables 16-24; data/costs.csv"
DALE_2019_COST_VALUES = {
    "test_igra": (113.48, "cscreenqft"),
    "test_tst": (116.07, "cscreentst"),
    "regimen_3hp": (165.5072, "ctreat3HP"),
    "regimen_4r": (123.3172, "ctreat4R"),
    "regimen_3hr": (134.2272, "ctreat3HR"),
    "regimen_6h": (187.7508, "ctreat6H"),
    "regimen_9h": (254.8544, "ctreat9H"),
    "active_tb_disease": (19079.60, "ctb"),
}
DALE_2019_ADR_COST_VALUES = {
    "3HP": (39.4059, "csae3HP"),
    "4R": (23.141, "csae4R"),
    "6H": (71.42, "csae6H"),
    "9H": (71.42, "csae9H"),
}
DALE_2019_PROGRAM_EXCLUSION_NOTE = (
    "Reviewed exclusion for this Dale-reproduction working default. The Dale "
    "model did not include the costs and logistics of programme implementation "
    "and expansion; this is not evidence that implementation costs are truly zero."
)
DALE_2019_FALSE_POSITIVE_BUNDLING_NOTE = (
    "Reviewed bundling decision: additional false-positive-specific costs are "
    "represented through screening-test cost, preventive-treatment course cost "
    "and applicable ADR-management cost. No additional false-positive-only cost "
    "is added."
)


def build_default_economics_config() -> dict[str, Any]:
    return {
        "metadata": {
            "currencyCode": TARGET_CURRENCY,
            "priceYear": TARGET_PRICE_YEAR,
            "targetCurrency": TARGET_CURRENCY,
            "targetPriceYear": TARGET_PRICE_YEAR,
            "perspective": "Australian health-system perspective",
            "locationLabel": "Australia",
            "sourceNotes": "",
            "programCostBasis": "",
            "scopeStatement": DIRECT_EFFECTS_SCOPE_STATEMENT,
        },
        "costNormalisation": {
            "contractVersion": "ltbi_cost_normalisation_v1",
            "formula": "converted_cost = original_cost * target_year_index / source_year_index",
            "defaultInflationIndexId": "AUD_HEALTH_SYSTEM_CPI",
            "indexTable": "data/inflation_indices.csv",
            "notes": (
                "Price-year conversion is separate from future discounting. "
                "Projected model costs remain in constant target-year prices."
            ),
        },
        "costItems": default_cost_items(),
        "costsCompatibilityNotes": COMPATIBILITY_COST_NOTES,
        "discounting": {
            "selectedAnnualRate": 0.03,
            "availableAnnualRates": [0.03, 0.0],
            "primaryDisplayedRate": 0.03,
            "comparisonRate": 0.0,
            "notes": "Discounting applies to future costs and health outcomes only, not source-price conversion.",
        },
        "healthOutcome": {
            "primary": "DALYs averted",
            "unresolvedInputs": [],
            "notes": "Natural epidemiological outputs remain reported alongside DALYs.",
        },
        "dalyAssumptions": default_daly_assumptions(),
        "threshold": {
            "metric": "cost per DALY averted",
            "concept": "1 x Australian GDP per capita per DALY averted",
            "value": None,
            "currency": TARGET_CURRENCY,
            "referenceYear": None,
            "source": "",
            "notes": "Illustrative benchmark only; not an official Australian funding threshold.",
        },
        "costs": {
            "test": {
                "IGRA": None,
                "TST": None,
            },
            "regimen": {
                "x3HP": None,
                "x4R": None,
                "x3HR": None,
                "x6H": None,
                "x9H": None,
            },
            "falsePositiveIncrementalPerPerson": None,
            "activeTBDiseasePerCase": None,
            "programSetupTotal": None,
            "programRunningTotal": None,
        },
    }


def build_economics_preset_kwab150() -> dict[str, Any]:
    config = build_default_economics_config()
    config["metadata"].update(
        {
            "currencyCode": "AUD",
            "priceYear": TARGET_PRICE_YEAR,
            "targetCurrency": TARGET_CURRENCY,
            "targetPriceYear": TARGET_PRICE_YEAR,
            "locationLabel": "Australia",
            "programCostBasis": "total",
            "sourceNotes": (
                "KWAB150 preset populated from local data/costs.csv: "
                "cscreenqft, cscreentst, ctreat*, and ctb mid values. "
                "The source price year is unresolved in the repository, so "
                "these costs are not treated as 2025-26 values until a verified "
                "source year and inflation-index values are supplied."
            ),
        }
    )
    config["costs"]["test"].update({"IGRA": 113.48, "TST": 116.07})
    config["costs"]["regimen"].update(
        {
            "x3HP": 165.5072,
            "x4R": 123.3172,
            "x3HR": 134.2272,
            "x6H": 187.7508,
            "x9H": 254.8544,
        }
    )
    config["costs"]["activeTBDiseasePerCase"] = 19079.6
    config["costItems"] = default_cost_items(
        {
            "test_igra": {
                "originalCost": 113.48,
                "sourceCitation": "Local data/costs.csv row cscreenqft",
                "pageTableItem": "cscreenqft",
                "notes": "Includes test kit and lab processing according to local row description.",
            },
            "test_tst": {
                "originalCost": 116.07,
                "sourceCitation": "Local data/costs.csv row cscreentst",
                "pageTableItem": "cscreentst",
                "notes": "Local row description says TST includes PPD injection and reading visit.",
            },
            "regimen_3hp": {
                "originalCost": 165.5072,
                "sourceCitation": "Local data/costs.csv row ctreat3HP",
                "pageTableItem": "ctreat3HP",
            },
            "regimen_4r": {
                "originalCost": 123.3172,
                "sourceCitation": "Local data/costs.csv row ctreat4R",
                "pageTableItem": "ctreat4R",
            },
            "regimen_3hr": {
                "originalCost": 134.2272,
                "sourceCitation": "Local data/costs.csv row ctreat3HR",
                "pageTableItem": "ctreat3HR",
            },
            "regimen_6h": {
                "originalCost": 187.7508,
                "sourceCitation": "Local data/costs.csv row ctreat6H",
                "pageTableItem": "ctreat6H",
            },
            "regimen_9h": {
                "originalCost": 254.8544,
                "sourceCitation": "Local data/costs.csv row ctreat9H",
                "pageTableItem": "ctreat9H",
            },
            "active_tb_disease": {
                "originalCost": 19079.6,
                "sourceCitation": "Local data/costs.csv row ctb",
                "pageTableItem": "ctb",
            },
        }
    )
    return sync_legacy_cost_fields_from_cost_items(config)


def build_economics_preset_dale2019_aud(regimen: str | None = "3HP") -> dict[str, Any]:
    """Build the user-selectable Dale/KWAB150 2019 AUD working preset.

    The preset is a transparent working default for reproducing a first APY
    ICER analysis. It does not make the overall APY evidence base
    clinician-ready.
    """
    selected_regimen = _normalise_regimen(regimen)
    config = build_default_economics_config()
    config["metadata"].update(
        {
            "presetName": DALE_2019_PRESET_NAME,
            "currencyCode": "AUD",
            "priceYear": "2019",
            "targetCurrency": "AUD",
            "targetPriceYear": "2019",
            "perspective": "Australian health-care system",
            "locationLabel": "Australia",
            "programCostBasis": "Dale reproduction working default; programme implementation costs excluded.",
            "sourceCitation": DALE_2019_CITATION,
            "sourceNotes": (
                "Cost rows use the existing KWAB150 model inputs from data/costs.csv. "
                "Dale et al. specify 2019 AUD and refer to Web Appendix 3 and Web "
                "Tables 16-24 for detailed costs; the main article does not print "
                "each unit cost."
            ),
            "workingDefault": True,
            "referenceStatus": DALE_2019_REFERENCE_STATUS,
            "provisional": True,
            "provisionalReason": "Working defaults do not resolve the full APY evidence-readiness gate.",
        }
    )
    config["discounting"].update(
        {
            "selectedAnnualRate": 0.03,
            "availableAnnualRates": [0.03, 0.0],
            "primaryDisplayedRate": 0.03,
            "comparisonRate": 0.0,
        }
    )
    config["threshold"].update(
        {
            "value": None,
            "currency": "AUD",
            "referenceYear": None,
            "source": "",
            "status": "unresolved",
            "notes": (
                "GDP-per-capita DALY threshold remains unresolved. Dale QALY "
                "thresholds are not silently reused as the authoritative DALY threshold."
            ),
        }
    )
    config["dalyAssumptions"] = _dale_2019_daly_assumptions()
    config["costItems"] = default_cost_items(_dale_2019_cost_overrides(selected_regimen))
    config["assumptionEvidenceRegistry"] = _dale_2019_registry(selected_regimen)
    config["assumptionEvidenceValidation"] = {
        "workingDefault": True,
        "referenceStatus": DALE_2019_REFERENCE_STATUS,
        "overallClinicianReady": False,
        "thresholdReady": False,
        "notes": "ICER may be calculable, but NMB and overall APY clinician readiness remain unavailable.",
    }
    config = sync_legacy_cost_fields_from_cost_items(config)
    config["costs"]["falsePositiveIncrementalPerPerson"] = 0.0
    config["costs"]["programSetupTotal"] = 0.0
    config["costs"]["programRunningTotal"] = 0.0
    return config


def default_cost_items(overrides: dict[str, dict[str, Any]] | None = None) -> list[dict[str, Any]]:
    overrides = overrides or {}
    specs = [
        ("test_igra", "IGRA screening test per person", "screening_test", {"costBasis": "per_person_screened"}),
        (
            "test_tst",
            "TST screening per person",
            "screening_test",
            {"costBasis": "per_person_screened", "returnVisitForReading": True},
        ),
        ("regimen_3hp", "3HP preventive regimen per started course", "preventive_regimen", {"costBasis": "per_started_course"}),
        ("regimen_4r", "4R preventive regimen per started course", "preventive_regimen", {"costBasis": "per_started_course"}),
        ("regimen_3hr", "3HR preventive regimen per started course", "preventive_regimen", {"costBasis": "per_started_course"}),
        ("regimen_6h", "6H preventive regimen per started course", "preventive_regimen", {"costBasis": "per_started_course"}),
        ("regimen_9h", "9H preventive regimen per started course", "preventive_regimen", {"costBasis": "per_started_course"}),
        ("active_tb_disease", "Active TB disease management per case", "tb_disease_care", {"costBasis": "per_active_tb_case"}),
        (
            "false_positive_incremental",
            "False-positive incremental resource use or cost per treated false-positive person",
            "false_positive_care",
            {"costBasis": "per_false_positive_tpt_started"},
        ),
        ("program_setup", "Programme setup total cost", "program_setup", {"costBasis": "total_once_at_program_start"}),
        ("program_running", "Programme running total cost", "program_running", {"costBasis": "annual_during_screening_window"}),
        (
            "tpt_adr_management",
            "TPT adverse-event management per ADR-related stop",
            "adverse_event_management",
            {"costBasis": "per_adr_stop"},
        ),
    ]
    items = []
    for cost_item_id, description, category, resource_use in specs:
        override = overrides.get(cost_item_id, {})
        items.append(
            build_cost_item(
                cost_item_id=cost_item_id,
                description=str(override.get("description", description)),
                original_cost=override.get("originalCost"),
                original_currency=str(override.get("originalCurrency", TARGET_CURRENCY)),
                original_price_year=override.get("originalPriceYear"),
                source=str(override.get("sourceCitation", "Not specified in repository.")),
                page_table_item=str(override.get("pageTableItem", "")),
                source_year_status=str(override.get("sourceYearStatus", "unknown")),
                resource_category=str(override.get("resourceCategory", category)),
                target_currency=str(override.get("targetCurrency", TARGET_CURRENCY)),
                target_price_year=str(override.get("targetPriceYear", TARGET_PRICE_YEAR)),
                inflation_index_id=str(override.get("inflationIndexId", "AUD_HEALTH_SYSTEM_CPI")),
                notes=str(
                    override.get(
                        "notes",
                        "Unresolved until a source value, source year and provenance are supplied.",
                    )
                ),
                resource_use=override.get("resourceUse", resource_use),
            )
        )
    return items


def _dale_2019_cost_overrides(regimen: str) -> dict[str, dict[str, Any]]:
    base = {
        item_id: _dale_2019_cost_override(item_id, value, variable)
        for item_id, (value, variable) in DALE_2019_COST_VALUES.items()
    }
    base["false_positive_incremental"] = {
        "originalCost": 0.0,
        "originalCurrency": "AUD",
        "originalPriceYear": "2019",
        "targetCurrency": "AUD",
        "targetPriceYear": "2019",
        "sourceCitation": DALE_2019_CITATION,
        "pageTableItem": "reviewed false-positive bundling decision",
        "sourceYearStatus": "explicit",
        "notes": DALE_2019_FALSE_POSITIVE_BUNDLING_NOTE,
        "resourceUse": {"costBasis": "per_false_positive_tpt_started"},
    }
    base["program_setup"] = {
        "originalCost": 0.0,
        "originalCurrency": "AUD",
        "originalPriceYear": "2019",
        "targetCurrency": "AUD",
        "targetPriceYear": "2019",
        "sourceCitation": DALE_2019_CITATION,
        "pageTableItem": "reviewed programme setup exclusion",
        "sourceYearStatus": "explicit",
        "notes": DALE_2019_PROGRAM_EXCLUSION_NOTE,
        "resourceUse": {"costBasis": "total_once_at_program_start"},
    }
    base["program_running"] = {
        "originalCost": 0.0,
        "originalCurrency": "AUD",
        "originalPriceYear": "2019",
        "targetCurrency": "AUD",
        "targetPriceYear": "2019",
        "sourceCitation": DALE_2019_CITATION,
        "pageTableItem": "reviewed programme running exclusion",
        "sourceYearStatus": "explicit",
        "notes": DALE_2019_PROGRAM_EXCLUSION_NOTE,
        "resourceUse": {"costBasis": "annual_during_screening_window"},
    }
    adr = DALE_2019_ADR_COST_VALUES.get(regimen)
    if adr:
        value, variable = adr
        base["tpt_adr_management"] = {
            "originalCost": value,
            "originalCurrency": "AUD",
            "originalPriceYear": "2019",
            "targetCurrency": "AUD",
            "targetPriceYear": "2019",
            "sourceCitation": DALE_2019_CITATION,
            "pageTableItem": variable,
            "sourceYearStatus": "explicit",
            "notes": (
                f"Dale/KWAB150 ADR-management cost mapped from {variable}. "
                "The APY event is an ADR-related treatment stop and may not "
                "exactly reproduce the severe-adverse-event definition in Dale et al."
            ),
            "resourceUse": {"costBasis": "per_adr_stop"},
        }
    else:
        base["tpt_adr_management"] = {
            "originalCost": None,
            "originalCurrency": "AUD",
            "originalPriceYear": None,
            "targetCurrency": "AUD",
            "targetPriceYear": "2019",
            "sourceCitation": "",
            "pageTableItem": "no Dale/KWAB150 3HR ADR-management mapping supplied",
            "sourceYearStatus": "unknown",
            "notes": "No 3HR ADR-management cost is supplied by this working default.",
            "resourceUse": {"costBasis": "per_adr_stop"},
        }
    return base


def _dale_2019_cost_override(item_id: str, value: float, variable: str) -> dict[str, Any]:
    notes = (
        f"Existing KWAB150 model input {variable} from data/costs.csv; "
        "treated as 2019 AUD only when the Dale 2019 AUD working-default preset is explicitly selected."
    )
    if item_id == "test_tst":
        notes += " TST row includes the reading visit according to the local row description."
    return {
        "originalCost": value,
        "originalCurrency": "AUD",
        "originalPriceYear": "2019",
        "targetCurrency": "AUD",
        "targetPriceYear": "2019",
        "sourceCitation": DALE_2019_CITATION,
        "pageTableItem": variable,
        "sourceYearStatus": "explicit",
        "notes": notes,
    }


def _dale_2019_daly_assumptions() -> dict[str, Any]:
    def rec(value: float, unit: str, source: str, notes: str) -> dict[str, Any]:
        return {
            "value": value,
            "unit": unit,
            "source": source,
            "status": "configured_reviewed",
            "notes": notes,
            "provisional": False,
        }

    return {
        "activeTBDisabilityWeight": rec(
            0.333,
            "disability weight",
            DALE_2019_CITATION,
            "Dale et al. discuss the GBD disability weight of 0.333 for active TB.",
        ),
        "activeTBDurationYears": rec(
            0.5,
            "years",
            DALE_2019_CITATION,
            "Dale et al. discuss active-TB disability experienced for six months.",
        ),
        "tbCaseFatalityRisk": rec(
            0.074,
            "probability",
            "User-approved APY working assumption",
            "Scalar model assumption requiring sensitivity analysis; not a Dale-derived age-specific mortality estimate.",
        ),
        "yllPerTBDeath": rec(
            20.0,
            "years",
            "User-approved APY working assumption",
            "Scalar model assumption requiring sensitivity analysis; not a Dale-derived age-specific mortality estimate.",
        ),
        "includeTPTHealthLoss": False,
        "tptHealthLossExclusionStatus": "reviewed_exclusion",
        "tptHealthLossExclusionRationale": (
            "Primary working analysis follows the Dale base case of no general "
            "LTBI-treatment health decrement beyond adverse events."
        ),
        "includeADRHealthLoss": False,
        "adrHealthLossExclusionStatus": "reviewed_exclusion",
        "adrHealthLossExclusionRationale": (
            "The APY ledger currently records ADR-related treatment stopping, "
            "which is not an unambiguous match to the severe-adverse-event utility "
            "state in Dale et al. ADR health loss will be explored when an "
            "appropriate severity mapping is available."
        ),
        "includePostTBSequelae": False,
        "postTBSequelaeStatus": "reviewed_exclusion",
        "postTBSequelaeExclusionRationale": (
            "Primary analysis estimates acute active-TB disability and mortality. "
            "Post-TB sequelae will be considered as a separate sensitivity scenario."
        ),
        "method": "scalar_expected_daly_per_active_tb_case",
        "notes": "Dale 2019 AUD working defaults; DALY assumptions remain working assumptions, not final APY evidence.",
    }


def _dale_2019_registry(regimen: str) -> list[dict[str, Any]]:
    rows = _load_evidence_registry_rows()
    by_id = {row.get("assumptionId"): row for row in rows}
    for item_id, (value, variable) in DALE_2019_COST_VALUES.items():
        _update_registry_cost_row(
            by_id,
            f"cost.{item_id}",
            value=value,
            variable=variable,
            status="configured_reviewed",
            provisional=False,
            inclusion_status="included",
            notes=(
                f"Reviewed working-default value from existing KWAB150 input {variable}; "
                "Dale et al. specify 2019 AUD and Web Appendix cost tables. "
                "This remains a working default, not final APY evidence."
            ),
        )
    adr = DALE_2019_ADR_COST_VALUES.get(regimen)
    if adr:
        value, variable = adr
        _update_registry_cost_row(
            by_id,
            "cost.tpt_adr_management",
            value=value,
            variable=variable,
            status="configured_reviewed",
            provisional=False,
            inclusion_status="included",
            notes=(
                f"Regimen-specific Dale/KWAB150 ADR-management mapping from {variable}. "
                "APY ADR stops may not exactly match Dale severe-adverse-event definitions."
            ),
        )
    else:
        _update_registry_cost_row(
            by_id,
            "cost.tpt_adr_management",
            value=None,
            variable="",
            status="unresolved",
            provisional=True,
            inclusion_status="included",
            source_citation="",
            notes="No Dale/KWAB150 ADR-management cost is supplied for 3HR.",
            unresolved_reason="3HR ADR-management cost remains unresolved.",
        )
    _update_registry_cost_row(
        by_id,
        "cost.false_positive_incremental",
        value="",
        variable="reviewed false-positive bundling decision",
        status="reviewed_exclusion",
        provisional=False,
        inclusion_status="bundled",
        bundled_into=f"cost.regimen_{regimen.lower()}",
        notes=DALE_2019_FALSE_POSITIVE_BUNDLING_NOTE,
    )
    _update_registry_cost_row(
        by_id,
        "cost.program_setup",
        value=0.0,
        variable="reviewed programme setup exclusion",
        status="reviewed_exclusion",
        provisional=False,
        inclusion_status="excluded",
        notes=DALE_2019_PROGRAM_EXCLUSION_NOTE
        + " This working default is likely optimistic where a new programme requires additional setup, staffing, engagement or delivery expenditure.",
    )
    _update_registry_cost_row(
        by_id,
        "cost.program_running",
        value=0.0,
        variable="reviewed programme running exclusion",
        status="reviewed_exclusion",
        provisional=False,
        inclusion_status="excluded",
        notes=DALE_2019_PROGRAM_EXCLUSION_NOTE
        + " This working default is likely optimistic where a new programme requires additional setup, staffing, engagement or delivery expenditure.",
    )
    _update_registry_daly_rows(by_id)
    row = by_id.get("threshold.gdp_per_capita")
    if row:
        row.update(
            {
                "currentValue": "",
                "targetCurrency": "AUD",
                "targetPriceYear": "2019",
                "originalPriceYear": "",
                "sourceCitation": "",
                "sourceLocation": "Not supplied in Dale 2019 AUD working defaults",
                "reviewStatus": "unresolved",
                "provisional": True,
                "notes": "No verified GDP-per-capita DALY threshold is supplied by this working default.",
                "unresolvedReason": "Threshold value, currency/year and source remain unresolved; NMB unavailable.",
            }
        )
    return [by_id.get(row.get("assumptionId"), row) for row in rows]


def _update_registry_cost_row(
    by_id: dict[str, dict[str, Any]],
    assumption_id: str,
    *,
    value: Any,
    variable: str,
    status: str,
    provisional: bool,
    inclusion_status: str,
    notes: str,
    source_citation: str | None = None,
    bundled_into: str = "",
    unresolved_reason: str = "",
) -> None:
    row = by_id.get(assumption_id)
    if not row:
        return
    row.update(
        {
            "currentValue": value,
            "originalCurrency": "AUD",
            "targetCurrency": "AUD",
            "originalPriceYear": "2019",
            "targetPriceYear": "2019",
            "sourceCitation": DALE_2019_CITATION if source_citation is None else source_citation,
            "sourceLocation": DALE_2019_SOURCE_LOCATION,
            "pageTableItem": variable,
            "repositoryPath": "data/costs.csv",
            "repositoryVariable": variable,
            "inflationIndexId": "AUD_HEALTH_SYSTEM_CPI",
            "inclusionStatus": inclusion_status,
            "reviewStatus": status,
            "provisional": provisional,
            "bundledIntoAssumptionId": bundled_into,
            "notes": notes,
            "unresolvedReason": unresolved_reason,
        }
    )


def _update_registry_daly_rows(by_id: dict[str, dict[str, Any]]) -> None:
    daly_rows = {
        "daly.active_tb_disability_weight": (
            0.333,
            "disability weight",
            DALE_2019_CITATION,
            "Dale et al. discussion of the GBD disability weight of 0.333.",
            "configured_reviewed",
            "included",
        ),
        "daly.active_tb_duration": (
            0.5,
            "years",
            DALE_2019_CITATION,
            "Dale et al. discussion of six months active-TB disability duration.",
            "configured_reviewed",
            "included",
        ),
        "daly.tb_case_fatality_risk": (
            0.074,
            "probability",
            "User-approved APY working assumption",
            "Scalar working assumption requiring sensitivity analysis.",
            "configured_reviewed",
            "included",
        ),
        "daly.yll_per_tb_death": (
            20.0,
            "years",
            "User-approved APY working assumption",
            "Scalar working assumption requiring sensitivity analysis.",
            "configured_reviewed",
            "included",
        ),
        "daly.tpt_health_loss": (
            "",
            "DALYs per TPT start",
            DALE_2019_CITATION,
            "Primary working analysis follows the Dale base case of no general LTBI-treatment health decrement beyond adverse events.",
            "reviewed_exclusion",
            "excluded",
        ),
        "daly.adr_health_loss": (
            "",
            "DALYs per ADR stop",
            DALE_2019_CITATION,
            "APY ADR-related stopping is not an unambiguous match to Dale severe-adverse-event utility state.",
            "reviewed_exclusion",
            "excluded",
        ),
        "daly.post_tb_sequelae": (
            "",
            "DALYs",
            "User-approved APY working assumption",
            "Primary analysis estimates acute active-TB disability and mortality; post-TB sequelae reserved for sensitivity scenarios.",
            "reviewed_exclusion",
            "excluded",
        ),
    }
    for assumption_id, (value, unit, source, notes, status, inclusion) in daly_rows.items():
        row = by_id.get(assumption_id)
        if not row:
            continue
        row.update(
            {
                "currentValue": value,
                "unit": unit,
                "sourceCitation": source,
                "sourceLocation": "Dale 2019 AUD working-default preset" if source == DALE_2019_CITATION else "APY working assumption",
                "reviewStatus": status,
                "provisional": False,
                "inclusionStatus": inclusion,
                "notes": notes,
                "unresolvedReason": "",
            }
        )


def _load_evidence_registry_rows() -> list[dict[str, Any]]:
    path = Path(__file__).resolve().parents[2] / "data" / "apy_assumption_evidence_registry.csv"
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def _normalise_regimen(regimen: str | None) -> str:
    text = str(regimen or "3HP").strip().upper().removeprefix("X")
    return text if text in {"3HP", "4R", "3HR", "6H", "9H"} else "3HP"


def sync_legacy_cost_fields_from_cost_items(econ_config: dict[str, Any]) -> dict[str, Any]:
    config = deepcopy(econ_config)
    costs = config.setdefault("costs", {})
    costs.setdefault("test", {})
    costs.setdefault("regimen", {})
    by_id = {
        item.get("costItemId"): item
        for item in ensure_authoritative_cost_items(config.get("costItems") or [])
    }
    for path, cost_item_id in LEGACY_COST_ITEM_MAP.items():
        _set_nested(config, path, by_id.get(cost_item_id, {}).get("originalCost"))
    config["costsCompatibilityNotes"] = COMPATIBILITY_COST_NOTES
    return config


def update_cost_item_original_values_from_legacy_fields(econ_config: dict[str, Any]) -> dict[str, Any]:
    config = deepcopy(econ_config)
    items = ensure_authoritative_cost_items(config.get("costItems") or [])
    by_id = {item.get("costItemId"): item for item in items}
    for path, cost_item_id in LEGACY_COST_ITEM_MAP.items():
        value = _get_nested(config, path)
        if value in (None, "", []):
            continue
        item = by_id.get(cost_item_id)
        if item is None:
            continue
        item["originalCost"] = value
        item["conversionStatus"] = "not_converted"
        item["conversionApplied"] = False
        item["costRecordType"] = "source"
        item["warnings"] = []
    config["costItems"] = list(by_id.values())
    return sync_legacy_cost_fields_from_cost_items(config)


def ensure_authoritative_cost_items(cost_items: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_id = {item.get("costItemId"): deepcopy(item) for item in cost_items if isinstance(item, dict)}
    for item in default_cost_items():
        by_id.setdefault(item["costItemId"], item)
    ordered_ids = [item["costItemId"] for item in default_cost_items()]
    return [by_id[item_id] for item_id in ordered_ids] + [
        item for item_id, item in by_id.items() if item_id not in ordered_ids
    ]


def validate_economics_config(econ_config: dict[str, Any]) -> dict[str, Any]:
    errors: list[dict[str, str]] = []
    warnings: list[dict[str, str]] = []

    if not isinstance(econ_config, dict):
        errors.append(_issue("econConfig", "Economics config must be a dict."))
        return _report(errors, warnings)

    metadata = econ_config.get("metadata")
    if not isinstance(metadata, dict):
        errors.append(_issue("metadata", "metadata must be a dict."))
        metadata = {}
    else:
        for field in (
            "currencyCode",
            "targetCurrency",
            "targetPriceYear",
            "locationLabel",
            "sourceNotes",
            "programCostBasis",
            "perspective",
            "scopeStatement",
        ):
            _validate_text_field(errors, metadata, field, f"metadata.{field}")
        if metadata.get("priceYear") not in (None, "", []) and not isinstance(
            metadata.get("priceYear"), (str, int, float)
        ):
            errors.append(_issue("metadata.priceYear", "metadata.priceYear must be a year label."))

    costs = econ_config.get("costs")
    if not isinstance(costs, dict):
        errors.append(_issue("costs", "costs must be a dict."))
        costs = {}

    test_costs = costs.get("test")
    if not isinstance(test_costs, dict):
        errors.append(_issue("costs.test", "costs.test must be a dict."))
        test_costs = {}
    else:
        for field in ("IGRA", "TST"):
            _validate_optional_cost(errors, test_costs, field, f"costs.test.{field}")

    regimen_costs = costs.get("regimen")
    if not isinstance(regimen_costs, dict):
        errors.append(_issue("costs.regimen", "costs.regimen must be a dict."))
        regimen_costs = {}
    else:
        for field in ("x3HP", "x4R", "x3HR", "x6H", "x9H"):
            _validate_optional_cost(
                errors,
                regimen_costs,
                field,
                f"costs.regimen.{field}",
            )

    for field in (
        "falsePositiveIncrementalPerPerson",
        "activeTBDiseasePerCase",
        "programSetupTotal",
        "programRunningTotal",
    ):
        _validate_optional_cost(errors, costs, field, f"costs.{field}")

    cost_items = econ_config.get("costItems")
    if cost_items not in (None, []):
        if not isinstance(cost_items, list):
            errors.append(_issue("costItems", "costItems must be a list."))
        else:
            for idx, item in enumerate(cost_items):
                if not isinstance(item, dict):
                    errors.append(_issue(f"costItems.{idx}", "Each cost item must be a dict."))
                elif item.get("conversionStatus") != "valid":
                    warnings.append(
                        _issue(
                            f"costItems.{idx}",
                            f"{item.get('costItemId', idx)} conversion is unresolved.",
                            "warning",
                        )
                    )

    discounting = econ_config.get("discounting", {})
    if isinstance(discounting, dict):
        _validate_optional_cost(
            errors,
            discounting,
            "selectedAnnualRate",
            "discounting.selectedAnnualRate",
        )
        rate = _coerce_number(discounting.get("selectedAnnualRate"))
        if rate is not None and rate not in {0, 0.03}:
            warnings.append(
                _issue(
                    "discounting.selectedAnnualRate",
                    "Primary milestone rates are 0% and 3%; other numeric rates are retained for later review.",
                    "warning",
                )
            )

    threshold = econ_config.get("threshold", {})
    if isinstance(threshold, dict) and _coerce_number(threshold.get("value")) is None:
        warnings.append(
            _issue(
                "threshold.value",
                "GDP-per-capita willingness-to-pay value is unresolved; enter a verified value and provenance before using threshold-based conclusions.",
                "warning",
            )
        )

    return _report(errors, warnings)


def run_health_economics(
    results: dict[str, Any],
    econ_config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    if econ_config is None:
        econ_config = build_default_economics_config()
    else:
        econ_config = _normalise_empty_to_none(deepcopy(econ_config))
    econ_config = sync_legacy_cost_fields_from_cost_items(
        {
            **econ_config,
            "costItems": ensure_authoritative_cost_items(econ_config.get("costItems") or []),
        }
    )

    if _has_event_ledger(results):
        routed_config = deepcopy(econ_config)
        configs = _candidate_interface_configs(results)
        if configs:
            routed_config["scenarioConfig"] = configs[0]
        return run_event_ledger_health_economics(results, routed_config)

    validation_report = validate_economics_config(econ_config)
    if not validation_report["isValid"]:
        raise ValueError("Economics config contains fatal validation errors.")

    status = _empty_status(validation_report)
    test_type = get_config_text(results, "testType")
    regimen = get_config_text(results, "regimen")

    quantities = {
        "nScreened": summary_metric(results, "nScreened"),
        "nTotalCoursesStarted": summary_metric(results, "nTotalCoursesStarted"),
        "nFalsePositiveTreated": summary_metric(results, "nFalsePositiveTreated"),
        "nCuredInfection": summary_metric(results, "nCuredInfection"),
        "nPreventedActiveTB": summary_metric(results, "nPreventedActiveTB"),
    }
    baseline_cases, intervention_cases = baseline_intervention_cases(results)
    quantities["baselineActiveTBCases"] = baseline_cases
    quantities["interventionActiveTBCases"] = intervention_cases

    normalised_items = normalise_cost_table(econ_config.get("costItems") or [])
    cost_item_lookup = {item.get("costItemId"): item for item in normalised_items}
    test_cost = selected_converted_test_cost(cost_item_lookup, test_type, status)
    regimen_cost = selected_converted_regimen_cost(cost_item_lookup, regimen, status)
    false_positive_cost = cost_item_value(cost_item_lookup, "false_positive_incremental")
    setup_cost = cost_item_value(cost_item_lookup, "program_setup")
    running_cost = cost_item_value(cost_item_lookup, "program_running")
    active_tb_cost = cost_item_value(cost_item_lookup, "active_tb_disease")

    status["messages"].extend(unresolved_cost_warnings(normalised_items))

    unit_costs = {
        "testPerPerson": test_cost,
        "treatmentPerStarted": regimen_cost,
        "falsePositiveIncrementalPerPerson": false_positive_cost,
        "activeTBDiseasePerCase": active_tb_cost,
    }

    cost_values = {
        "testingCost": multiply_if_available(quantities["nScreened"], test_cost),
        "treatmentCost": multiply_if_available(
            quantities["nTotalCoursesStarted"],
            regimen_cost,
        ),
        "falsePositiveIncrementalCost": multiply_if_available(
            quantities["nFalsePositiveTreated"],
            unit_costs["falsePositiveIncrementalPerPerson"],
        ),
        "programSetupCost": setup_cost,
        "programRunningCost": running_cost,
        "tbDiseaseCostsAverted": multiply_if_available(
            quantities["nPreventedActiveTB"],
            active_tb_cost,
        ),
        "baselineTBDiseaseCost": multiply_if_available(
            quantities["baselineActiveTBCases"],
            active_tb_cost,
        ),
        "interventionTBDiseaseCost": multiply_if_available(
            quantities["interventionActiveTBCases"],
            active_tb_cost,
        ),
    }

    _note_missing_if_empty(
        status,
        cost_values["testingCost"],
        "testingCost",
        "Testing cost not calculated because nScreened or selected test cost is missing.",
    )
    _note_missing_if_empty(
        status,
        cost_values["treatmentCost"],
        "treatmentCost",
        "Treatment cost not calculated because starts or selected regimen cost is missing.",
    )
    if _is_missing_number(unit_costs["falsePositiveIncrementalPerPerson"]):
        status["messages"].append(
            "No false-positive incremental unit cost supplied; no extra false-positive cost is included."
        )
    _note_missing_if_empty(
        status,
        cost_values["tbDiseaseCostsAverted"],
        "tbDiseaseCostsAverted",
        "TB disease costs averted not calculated because prevented cases or disease cost is missing.",
    )
    if _is_missing_number(cost_values["baselineTBDiseaseCost"]) or _is_missing_number(
        cost_values["interventionTBDiseaseCost"]
    ):
        status["partialCalculations"].append(
            "Baseline/intervention TB disease costs were not calculated because case counts are unavailable."
        )

    cost_values["totalProgramCost"] = sum_required_costs(
        [
            cost_values["testingCost"],
            cost_values["treatmentCost"],
            zero_if_empty(cost_values["falsePositiveIncrementalCost"]),
            cost_values["programSetupCost"],
            cost_values["programRunningCost"],
        ]
    )
    _note_missing_if_empty(
        status,
        cost_values["totalProgramCost"],
        "totalProgramCost",
        "Total program cost not calculated because one or more required cost components is missing.",
    )

    cost_values["netCostVsBaseline"] = subtract_if_available(
        cost_values["totalProgramCost"],
        cost_values["tbDiseaseCostsAverted"],
    )
    _note_missing_if_empty(
        status,
        cost_values["netCostVsBaseline"],
        "netCostVsBaseline",
        "Net cost versus baseline not calculated because total program cost or TB costs averted is missing.",
    )

    cost_effectiveness = {
        "costPerInfectionCured": divide_if_available(
            cost_values["netCostVsBaseline"],
            quantities["nCuredInfection"],
            "costPerInfectionCured",
            status,
        ),
        "costPerTBCasePrevented": divide_if_available(
            cost_values["netCostVsBaseline"],
            quantities["nPreventedActiveTB"],
            "costPerTBCasePrevented",
            status,
        ),
    }

    econ = {
        "available": True,
        "source": ECONOMICS_SOURCE,
        "contractVersion": ECONOMICS_CONTRACT_VERSION,
        "metadata": deepcopy(econ_config.get("metadata", {})),
        "inputs": deepcopy(econ_config),
        "costNormalisation": deepcopy(econ_config.get("costNormalisation", {})),
        "costItems": normalised_items,
        "discounting": deepcopy(econ_config.get("discounting", {})),
        "healthOutcome": deepcopy(econ_config.get("healthOutcome", {})),
        "threshold": deepcopy(econ_config.get("threshold", {})),
        "scopeStatement": econ_config.get("metadata", {}).get(
            "scopeStatement", DIRECT_EFFECTS_SCOPE_STATEMENT
        ),
        "strategy": {
            "testType": test_type,
            "regimen": regimen,
        },
        "quantities": quantities,
        "unitCosts": unit_costs,
        "costs": cost_values,
        "costEffectiveness": cost_effectiveness,
        "status": status,
    }
    econ["summaryRows"] = build_summary_rows(econ)
    econ["summaryTable"] = econ["summaryRows"]
    status["isComplete"] = not status["missingInputs"] and not status["notCalculated"]
    return econ


def run_health_economics_for_config(
    config: dict[str, Any],
    econ_config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    out = run_scenario_with_do_nothing(config)
    econ = run_health_economics(out, econ_config)
    econ["scenario"] = {
        "modelVersion": out["results"].get("modelVersion"),
        "backend": out["results"].get("backend"),
        "scenarioLabel": out["results"].get("interfaceConfig", {}).get("scenarioLabel"),
    }
    return econ


def _has_event_ledger(results: dict[str, Any]) -> bool:
    if isinstance(results.get("eventLedger"), dict):
        return True
    if isinstance(results.get("technical"), dict) and isinstance(results["technical"].get("eventLedger"), dict):
        return True
    if isinstance(results.get("bundle"), dict):
        return _has_event_ledger(results["bundle"])
    return False


def summary_metric(results: dict[str, Any], metric_name: str) -> float | None:
    rows = _summary_rows(results)
    for row in rows:
        if str(row.get("Metric", row.get("metric", ""))) == metric_name:
            return _coerce_number(row.get("Median", row.get("median")))
    return None


def get_config_text(results: dict[str, Any], field_name: str) -> str:
    for config in _candidate_interface_configs(results):
        value = config.get(field_name)
        if value not in (None, "", []):
            return str(value)
    return ""


def selected_test_cost(
    econ_config: dict[str, Any],
    test_type: str,
    status: dict[str, Any],
) -> float | None:
    if not test_type:
        status["missingInputs"].append("results.interfaceConfig.testType")
        return None
    costs = econ_config.get("costs", {}).get("test", {})
    value = _coerce_number(costs.get(str(test_type).upper()))
    if value is None:
        status["missingInputs"].append(f"costs.test.{test_type}")
    return value


def selected_regimen_cost(
    econ_config: dict[str, Any],
    regimen: str,
    status: dict[str, Any],
) -> float | None:
    if not regimen:
        status["missingInputs"].append("results.interfaceConfig.regimen")
        return None
    field = _regimen_cost_field(regimen)
    if field is None:
        status["missingInputs"].append(f"Unknown regimen label: {regimen}")
        return None
    costs = econ_config.get("costs", {}).get("regimen", {})
    value = _coerce_number(costs.get(field))
    if value is None:
        status["missingInputs"].append(f"costs.regimen.{field}")
    return value


def selected_converted_test_cost(
    cost_items: dict[str, dict[str, Any]],
    test_type: str,
    status: dict[str, Any],
) -> float | None:
    key = {"IGRA": "test_igra", "TST": "test_tst"}.get(str(test_type).upper())
    if key is None:
        status["missingInputs"].append(f"costItems test cost for {test_type}")
        return None
    value = cost_item_value(cost_items, key)
    if value is None:
        status["missingInputs"].append(f"costItems.{key}.convertedTargetYearCost")
    return value


def selected_converted_regimen_cost(
    cost_items: dict[str, dict[str, Any]],
    regimen: str,
    status: dict[str, Any],
) -> float | None:
    key = {
        "3HP": "regimen_3hp",
        "4R": "regimen_4r",
        "3HR": "regimen_3hr",
        "6H": "regimen_6h",
        "9H": "regimen_9h",
    }.get(str(regimen).upper())
    if key is None:
        status["missingInputs"].append(f"Unknown regimen label: {regimen}")
        return None
    value = cost_item_value(cost_items, key)
    if value is None:
        status["missingInputs"].append(f"costItems.{key}.convertedTargetYearCost")
    return value


def cost_item_value(cost_items: dict[str, dict[str, Any]], cost_item_id: str) -> float | None:
    return valid_converted_cost(cost_items.get(cost_item_id, {}))


def baseline_intervention_cases(results: dict[str, Any]) -> tuple[float | None, float | None]:
    derived = _do_nothing_derived(results)
    if derived is not None:
        baseline = _median_from_rows_or_frame(derived, "nActiveBy20y_DoNothing")
        intervention = _median_from_rows_or_frame(derived, "nActiveBy20y_AfterStrategy")
        if baseline is not None or intervention is not None:
            return baseline, intervention

    dynamic = _dynamic_comparison(results)
    if dynamic:
        baseline = _coerce_number(dynamic.get("cumulative_baseline_active_tb_cases"))
        intervention = _coerce_number(
            dynamic.get("cumulative_intervention_active_tb_cases")
        )
        if baseline is not None or intervention is not None:
            return baseline, intervention

    raw = _raw_frame(results)
    if raw is not None and {"nActiveBy20y", "nPreventedActiveTB"}.issubset(raw.columns):
        baseline_values = pd.to_numeric(raw["nActiveBy20y"], errors="coerce")
        intervention_values = baseline_values - pd.to_numeric(
            raw["nPreventedActiveTB"],
            errors="coerce",
        )
        return _median_series(baseline_values), _median_series(intervention_values)

    return None, None


def multiply_if_available(a: Any, b: Any) -> float | None:
    a_num = _coerce_number(a)
    b_num = _coerce_number(b)
    if a_num is None or b_num is None:
        return None
    return a_num * b_num


def subtract_if_available(a: Any, b: Any) -> float | None:
    a_num = _coerce_number(a)
    b_num = _coerce_number(b)
    if a_num is None or b_num is None:
        return None
    return a_num - b_num


def divide_if_available(
    a: Any,
    b: Any,
    field_name: str,
    status: dict[str, Any],
) -> float | None:
    a_num = _coerce_number(a)
    b_num = _coerce_number(b)
    if a_num is None or b_num is None or b_num == 0:
        status["notCalculated"].append(field_name)
        status["messages"].append(
            f"{field_name} not calculated because the denominator is missing, zero, or NaN."
        )
        return None
    return a_num / b_num


def sum_required_costs(values: list[Any]) -> float | None:
    total = 0.0
    for value in values:
        value_num = _coerce_number(value)
        if value_num is None:
            return None
        total += value_num
    return total


def zero_if_empty(value: Any) -> float:
    value_num = _coerce_number(value)
    return 0.0 if value_num is None else value_num


def build_summary_rows(econ: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    metadata = econ.get("metadata", {})
    currency = metadata.get("currencyCode", "")
    price_year = metadata.get("priceYear")
    specs = [
        ("Costs", "testingCost", econ["costs"].get("testingCost")),
        ("Costs", "treatmentCost", econ["costs"].get("treatmentCost")),
        (
            "Costs",
            "falsePositiveIncrementalCost",
            econ["costs"].get("falsePositiveIncrementalCost"),
        ),
        ("Costs", "programSetupCost", econ["costs"].get("programSetupCost")),
        ("Costs", "programRunningCost", econ["costs"].get("programRunningCost")),
        ("Costs", "baselineTBDiseaseCost", econ["costs"].get("baselineTBDiseaseCost")),
        (
            "Costs",
            "interventionTBDiseaseCost",
            econ["costs"].get("interventionTBDiseaseCost"),
        ),
        ("Costs", "tbDiseaseCostsAverted", econ["costs"].get("tbDiseaseCostsAverted")),
        ("Costs", "totalProgramCost", econ["costs"].get("totalProgramCost")),
        ("Costs", "netCostVsBaseline", econ["costs"].get("netCostVsBaseline")),
        (
            "Cost-effectiveness",
            "costPerInfectionCured",
            econ["costEffectiveness"].get("costPerInfectionCured"),
        ),
        (
            "Cost-effectiveness",
            "costPerTBCasePrevented",
            econ["costEffectiveness"].get("costPerTBCasePrevented"),
        ),
    ]
    for section, metric, value in specs:
        rows.append(
            {
                "Section": section,
                "Metric": metric,
                "Value": _coerce_number(value),
                "Calculated": _coerce_number(value) is not None,
                "CurrencyCode": currency,
                "PriceYear": price_year,
                "Notes": "",
            }
        )
    return rows


def _issue(field: str, message: str, severity: str = "error") -> dict[str, str]:
    return {
        "field": field,
        "severity": severity,
        "message": message,
    }


def _report(
    errors: list[dict[str, str]],
    warnings: list[dict[str, str]],
) -> dict[str, Any]:
    return {
        "isValid": not errors,
        "hasWarnings": bool(warnings),
        "errors": errors,
        "warnings": warnings,
        "fatalFieldNames": [issue["field"] for issue in errors],
        "warningFieldNames": [issue["field"] for issue in warnings],
    }


def _validate_text_field(
    errors: list[dict[str, str]],
    parent: dict[str, Any],
    field: str,
    full_name: str,
) -> None:
    value = parent.get(field)
    if value in (None, "", []):
        return
    if not isinstance(value, str):
        errors.append(_issue(full_name, f"{full_name} must be text."))


def _validate_optional_numeric_scalar(
    errors: list[dict[str, str]],
    parent: dict[str, Any],
    field: str,
    full_name: str,
) -> None:
    value = parent.get(field)
    if value in (None, "", []):
        return
    if _coerce_number(value) is None:
        errors.append(_issue(full_name, f"{full_name} must be a finite scalar."))


def _validate_optional_cost(
    errors: list[dict[str, str]],
    parent: dict[str, Any],
    field: str,
    full_name: str,
) -> None:
    value = parent.get(field)
    if value in (None, "", []):
        return
    number = _coerce_number(value)
    if number is None or number < 0:
        errors.append(
            _issue(full_name, f"{full_name} must be empty or a non-negative scalar.")
        )


def _empty_status(validation_report: dict[str, Any]) -> dict[str, Any]:
    return {
        "isComplete": False,
        "missingInputs": [],
        "partialCalculations": [],
        "notCalculated": [],
        "messages": [],
        "validationReport": validation_report,
    }


def _optional_cost(parent: dict[str, Any], field: str) -> float | None:
    return _coerce_number(parent.get(field))


def _note_missing_if_empty(
    status: dict[str, Any],
    value: Any,
    name: str,
    message: str,
) -> None:
    if _is_missing_number(value):
        status["notCalculated"].append(name)
        status["messages"].append(message)


def _is_missing_number(value: Any) -> bool:
    return _coerce_number(value) is None


def _coerce_number(value: Any) -> float | None:
    if value in (None, "", []):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    if number.is_integer():
        return int(number)
    return number


def _normalise_empty_to_none(value: Any) -> Any:
    if value == []:
        return None
    if isinstance(value, dict):
        return {key: _normalise_empty_to_none(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_normalise_empty_to_none(item) for item in value]
    return value


def _get_nested(config: dict[str, Any], path: tuple[str, ...]) -> Any:
    current: Any = config
    for key in path:
        if not isinstance(current, dict):
            return None
        current = current.get(key)
    return current


def _set_nested(config: dict[str, Any], path: tuple[str, ...], value: Any) -> None:
    current = config
    for key in path[:-1]:
        child = current.get(key)
        if not isinstance(child, dict):
            child = {}
            current[key] = child
        current = child
    current[path[-1]] = value


def _summary_rows(results: dict[str, Any]) -> list[dict[str, Any]]:
    if "results" in results and isinstance(results["results"], dict):
        return _summary_rows(results["results"])
    summary = results.get("summary")
    if isinstance(summary, pd.DataFrame):
        return summary.to_dict(orient="records")
    if isinstance(summary, list):
        return summary
    headline = results.get("headline")
    if isinstance(headline, dict):
        rows = headline.get("summaryRows") or headline.get("keyMetricsRows")
        if isinstance(rows, list):
            return rows
    return []


def _candidate_interface_configs(results: dict[str, Any]) -> list[dict[str, Any]]:
    configs = []
    if "results" in results and isinstance(results["results"], dict):
        configs.extend(_candidate_interface_configs(results["results"]))
    for candidate in (
        results.get("interfaceConfig"),
        results.get("technical", {}).get("interfaceConfig")
        if isinstance(results.get("technical"), dict)
        else None,
    ):
        if isinstance(candidate, dict):
            configs.append(candidate)
    return configs


def _do_nothing_derived(results: dict[str, Any]) -> Any:
    do_nothing = results.get("doNothing")
    if isinstance(do_nothing, dict):
        return do_nothing.get("derived")
    return None


def _dynamic_comparison(results: dict[str, Any]) -> dict[str, Any]:
    if "bundle" in results and isinstance(results["bundle"], dict):
        dynamic = _dynamic_comparison(results["bundle"])
        if dynamic:
            return dynamic
    technical = results.get("technical")
    if isinstance(technical, dict) and isinstance(technical.get("dynamicComparison"), dict):
        return technical["dynamicComparison"]
    return {}


def _raw_frame(results: dict[str, Any]) -> pd.DataFrame | None:
    if "results" in results and isinstance(results["results"], dict):
        return _raw_frame(results["results"])
    raw = results.get("raw")
    if isinstance(raw, pd.DataFrame):
        return raw
    if isinstance(raw, list):
        return pd.DataFrame(raw)
    return None


def _median_from_rows_or_frame(rows_or_frame: Any, column: str) -> float | None:
    if isinstance(rows_or_frame, pd.DataFrame):
        if column not in rows_or_frame.columns:
            return None
        return _median_series(pd.to_numeric(rows_or_frame[column], errors="coerce"))
    if isinstance(rows_or_frame, list):
        values = [
            _coerce_number(row.get(column))
            for row in rows_or_frame
            if isinstance(row, dict)
        ]
        values = [value for value in values if value is not None]
        if not values:
            return None
        return _median_series(pd.Series(values, dtype="float64"))
    return None


def _median_series(values: pd.Series) -> float | None:
    clean = values.dropna()
    if clean.empty:
        return None
    value = float(clean.median())
    if value.is_integer():
        return int(value)
    return value


def _regimen_cost_field(regimen: str) -> str | None:
    normalised = str(regimen).strip().upper()
    mapping = {
        "3HP": "x3HP",
        "4R": "x4R",
        "3HR": "x3HR",
        "6H": "x6H",
        "9H": "x9H",
    }
    return mapping.get(normalised)
