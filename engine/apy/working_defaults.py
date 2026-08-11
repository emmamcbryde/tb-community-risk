from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from typing import Any

from engine.apy.config import build_default_config
from engine.apy.economics import build_economics_preset_dale2019_aud
from engine.apy.scenario import (
    DEFAULT_POPULATION_PRESET_ID,
    build_scenario_contract,
    config_updates_from_scenario,
)


UNIFIED_WORKING_DEFAULT_PRESET_ID = "apy_sa_health_working_defaults"
UNIFIED_WORKING_DEFAULT_PRESET_VERSION = "apy_sa_health_working_defaults_v1"
UNIFIED_WORKING_DEFAULT_LABEL = "APY / SA Health working defaults"
PRIMARY_DALY_OUTCOME_PRESET = "Primary acute TB DALYs"
HIGHER_BURDEN_DALY_OUTCOME_PRESET = "Higher TB burden - include post-TB health loss"


def build_unified_working_default_preset() -> dict[str, Any]:
    """Return the standard clinician-interface working-default bundle."""
    scenario = build_scenario_contract(DEFAULT_POPULATION_PRESET_ID)
    config = build_default_config()
    config.update(config_updates_from_scenario(scenario))
    config.update(
        {
            "scenarioLabel": "APY / SA Health working-default analysis",
            "analysisMethod": "expected_value",
            "analysisMethodLabel": "Expected outcomes",
            "simulationMode": "standard",
            "simulationModeLabel": "Standard analysis: use the current validated repository standard",
            "workingDefaultPresetId": UNIFIED_WORKING_DEFAULT_PRESET_ID,
            "workingDefaultPresetVersion": UNIFIED_WORKING_DEFAULT_PRESET_VERSION,
            "workingDefaultPresetLabel": UNIFIED_WORKING_DEFAULT_LABEL,
        }
    )
    economics_config = build_economics_preset_dale2019_aud(config.get("regimen", "3HP"))
    _add_sa_health_pathway_assumptions(economics_config)
    _apply_outcome_preset(economics_config, PRIMARY_DALY_OUTCOME_PRESET)
    preset = {
        "contractVersion": "apy_unified_working_default_preset_v1",
        "presetId": UNIFIED_WORKING_DEFAULT_PRESET_ID,
        "presetVersion": UNIFIED_WORKING_DEFAULT_PRESET_VERSION,
        "label": UNIFIED_WORKING_DEFAULT_LABEL,
        "sourceComponentPresets": [
            "APY demonstration demography",
            "APY TB epidemiology defaults",
            "Dale 2019 AUD working defaults",
            "SA Health local pathway working assumptions",
            PRIMARY_DALY_OUTCOME_PRESET,
        ],
        "workingDefault": True,
        "referenceStatus": "working_default_not_final_apy_evidence",
        "config": config,
        "economicsConfig": economics_config,
        "provisionalAssumptions": _provisional_assumptions(config, economics_config),
        "unresolvedAssumptions": _unresolved_assumptions(config, economics_config),
    }
    preset["configurationHash"] = _hash_payload(
        {
            "config": preset["config"],
            "economicsConfig": preset["economicsConfig"],
            "presetVersion": preset["presetVersion"],
        }
    )
    config["workingDefaultPresetHash"] = preset["configurationHash"]
    economics_config.setdefault("metadata", {})["workingDefaultPresetHash"] = preset["configurationHash"]
    economics_config["metadata"]["workingDefaultPresetId"] = UNIFIED_WORKING_DEFAULT_PRESET_ID
    economics_config["metadata"]["workingDefaultPresetVersion"] = UNIFIED_WORKING_DEFAULT_PRESET_VERSION
    return preset


def apply_program_setup_preset(economics_config: dict[str, Any], preset_id: str) -> dict[str, Any]:
    out = deepcopy(economics_config)
    presets = ((out.get("localSAHealthPathwayCosts") or {}).get("programSetupPresets") or {})
    preset = presets.get(preset_id)
    if not preset:
        raise ValueError(f"Unknown programme setup preset: {preset_id}")
    value = float(preset["value"])
    out.setdefault("localSAHealthPathwayCosts", {})["selectedProgramSetupPreset"] = preset_id
    out.setdefault("metadata", {})["programSetupPreset"] = preset_id
    for item in out.get("costItems") or []:
        if item.get("costItemId") == "program_setup":
            item["originalCost"] = value
            item["sourceCitation"] = preset.get("source", "SA Health local working assumption")
            item["notes"] = preset.get("notes", "")
            item["conversionStatus"] = "not_converted"
            item["conversionApplied"] = False
            item["costRecordType"] = "source"
            item["warnings"] = []
    for row in out.get("assumptionEvidenceRegistry") or []:
        if row.get("assumptionId") == "cost.program_setup":
            row["currentValue"] = value
            row["sourceCitation"] = preset.get("source", "SA Health local working assumption")
            row["reviewStatus"] = "configured_reviewed"
            row["provisional"] = False
            row["inclusionStatus"] = "included"
            row["notes"] = preset.get("notes", "")
            row["unresolvedReason"] = ""
    return out


def apply_outcome_preset(economics_config: dict[str, Any], preset_name: str) -> dict[str, Any]:
    out = deepcopy(economics_config)
    _apply_outcome_preset(out, preset_name)
    return out


def _apply_outcome_preset(economics_config: dict[str, Any], preset_name: str) -> None:
    daly = economics_config.setdefault("dalyAssumptions", {})
    daly["outcomePreset"] = preset_name
    daly["standardOutcomeMetric"] = "DALYs"
    daly["qalyStandardInterfaceVisible"] = False
    if preset_name == HIGHER_BURDEN_DALY_OUTCOME_PRESET:
        daly["includePostTBSequelae"] = True
        daly["postTBSequelaeStatus"] = "configured_reviewed"
        daly["postTBSequelaeExclusionRationale"] = ""
        daly["postTBDALYsPerActiveTBCase"] = {
            "value": 5.8,
            "unit": "DALYs per active-TB case",
            "source": "Named higher-burden working scenario",
            "status": "configured_reviewed",
            "notes": "Post-TB DALY sensitivity option; DALYs and QALYs are not combined.",
            "provisional": False,
        }
    else:
        daly["includePostTBSequelae"] = False
        daly["postTBSequelaeStatus"] = "reviewed_exclusion"
        daly.setdefault(
            "postTBSequelaeExclusionRationale",
            "Primary analysis estimates acute active-TB disability and mortality.",
        )
        daly.pop("postTBDALYsPerActiveTBCase", None)
    economics_config.setdefault("metadata", {})["outcomePreset"] = preset_name


def _add_sa_health_pathway_assumptions(economics_config: dict[str, Any]) -> None:
    economics_config["localSAHealthPathwayCosts"] = {
        "currency": "AUD",
        "priceYear": "2019",
        "status": "local_working_assumption_not_final_apy_evidence",
        "returnForResults": {
            "value": 50.0,
            "basis": "per_person_screened",
            "engineMapping": "not_separately_mapped",
            "notes": "Do not double count if already included in the selected test cost.",
        },
        "clinicalReview": {
            "value": 50.0,
            "basis": "per_tpt_start",
            "engineMapping": "not_separately_mapped",
            "notes": "Do not double count if already included in the selected regimen cost.",
        },
        "activeTBExclusionWorkup": {
            "value": 150.0,
            "basis": "per_tpt_start",
            "engineMapping": "not_separately_mapped",
            "notes": "CXR/laboratory work-up assumption; not added unless a non-overlapping mapping is configured.",
        },
        "travelOutreachStaffSupport": {
            "value": 0.0,
            "basis": "per_person_screened",
            "status": "local_working_assumption_not_yet_locally_costed",
            "notes": "Not yet locally costed.",
        },
        "programRunning": {
            "value": 0.0,
            "basis": "annual_during_screening_window",
            "status": "local_working_assumption_not_yet_locally_costed",
            "notes": "Not yet locally costed.",
        },
        "programSetupPresets": {
            "existing_service_base": {
                "label": "Existing-service/base scenario",
                "value": 0.0,
                "source": "SA Health local working assumption",
                "notes": "Existing-service/base working scenario.",
            },
            "new_programme_implementation": {
                "label": "New-programme implementation scenario",
                "value": 500000.0,
                "source": "SA Health local working assumption",
                "notes": "Alternative setup-cost scenario for a new programme.",
            },
        },
        "selectedProgramSetupPreset": "existing_service_base",
    }
    economics_config.setdefault("metadata", {})["saHealthPathwayCosts"] = "working_defaults_v1"
    economics_config["metadata"]["costFramework"] = "2019 AUD"


def _provisional_assumptions(config: dict[str, Any], economics_config: dict[str, Any]) -> list[str]:
    out = []
    ltbi = config.get("ltbiStateAssumptions") or {}
    if ltbi.get("provisional"):
        out.append("epi.baseline_recent_ltbi_proportion")
    for row in economics_config.get("assumptionEvidenceRegistry") or []:
        if row.get("provisional"):
            out.append(str(row.get("assumptionId")))
    return sorted(set(out))


def _unresolved_assumptions(config: dict[str, Any], economics_config: dict[str, Any]) -> list[str]:
    out = []
    ltbi = config.get("ltbiStateAssumptions") or {}
    if str(ltbi.get("baselineRecentLTBIProportionStatus") or "").startswith("unresolved"):
        out.append("epi.baseline_recent_ltbi_proportion")
    for row in economics_config.get("assumptionEvidenceRegistry") or []:
        if row.get("reviewStatus") == "unresolved" or row.get("unresolvedReason"):
            out.append(str(row.get("assumptionId")))
    return sorted(set(out))


def _hash_payload(payload: dict[str, Any]) -> str:
    text = json.dumps(payload, sort_keys=True, default=str, separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()
