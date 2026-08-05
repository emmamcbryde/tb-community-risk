from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from typing import Any


SCENARIO_CONTRACT_VERSION = "ltbi_screening_scenario_v1"
DEFAULT_POPULATION_PRESET_ID = "apy_demonstration"
DEFAULT_COMPARATOR = "current practice / no additional systematic LTBI screening"
DEFAULT_INTERVENTION = "targeted LTBI screening and preventive treatment"
DIRECT_EFFECTS_SCOPE_STATEMENT = (
    "The current static and individual-based analyses estimate direct benefits, "
    "harms and costs for screened and treated individuals. Transmission-mediated "
    "benefits are not yet included."
)

REPO_ROOT = Path(__file__).resolve().parents[2]
POPULATION_PRESETS_PATH = REPO_ROOT / "data" / "population_presets.json"


def load_population_presets(path: str | Path | None = None) -> dict[str, dict[str, Any]]:
    source = Path(path) if path is not None else POPULATION_PRESETS_PATH
    payload = json.loads(source.read_text(encoding="utf-8"))
    presets = payload.get("presets")
    if not isinstance(presets, list):
        raise ValueError("Population preset file must contain a presets list.")
    out: dict[str, dict[str, Any]] = {}
    for preset in presets:
        if not isinstance(preset, dict):
            continue
        preset_id = str(preset.get("populationPresetId", "")).strip()
        if not preset_id:
            raise ValueError("Population preset is missing populationPresetId.")
        out[preset_id] = deepcopy(preset)
    return out


def load_population_preset(
    population_preset_id: str = DEFAULT_POPULATION_PRESET_ID,
    path: str | Path | None = None,
) -> dict[str, Any]:
    presets = load_population_presets(path)
    try:
        return deepcopy(presets[population_preset_id])
    except KeyError as exc:
        raise ValueError(f"Unknown population preset: {population_preset_id}") from exc


def build_scenario_contract(
    population_preset_id: str = DEFAULT_POPULATION_PRESET_ID,
    overrides: dict[str, Any] | None = None,
) -> dict[str, Any]:
    preset = load_population_preset(population_preset_id)
    scenario = {
        "contractVersion": SCENARIO_CONTRACT_VERSION,
        "scenarioName": preset.get("scenarioName", ""),
        "scenarioVersion": preset.get("scenarioVersion", "m1"),
        "populationPresetId": preset.get("populationPresetId", population_preset_id),
        "populationSize": preset.get("populationSize"),
        "ageProfileSource": deepcopy(preset.get("ageProfileSource", {})),
        "ltbiPrevalenceAssumptions": deepcopy(
            preset.get("ltbiPrevalenceAssumptions", {})
        ),
        "riskFactorAssumptions": deepcopy(preset.get("riskFactorAssumptions", {})),
        "targetingCriteria": deepcopy(preset.get("targetingCriteria", {})),
        "eligible": deepcopy(preset.get("eligible", {})),
        "screened": deepcopy(preset.get("screened", {})),
        "screeningWindowYears": preset.get("screeningWindowYears"),
        "earlyProgressionPeriodYears": preset.get("earlyProgressionPeriodYears", 2),
        "activeTBCalibrationHorizonYears": preset.get("activeTBCalibrationHorizonYears", 2),
        "followUpHorizonYears": preset.get("followUpHorizonYears"),
        "ltbiStateAssumptions": deepcopy(
            preset.get(
                "ltbiStateAssumptions",
                {
                    "baselineRecentLTBIProportion": None,
                    "recentToRemoteTransitionRatePerYear": 0.2,
                    "transitionModel": "continuous_markov_recent_remote",
                    "stateDefinition": (
                        "fast/recent latent state with mean residence time of "
                        "five years before transition to remote latent infection"
                    ),
                    "source": (
                        "Transition structure from older static/transmission-dynamic "
                        "architecture; APY-specific baseline recent fraction unresolved."
                    ),
                    "status": "unresolved",
                    "notes": "",
                    "developmentCompatibilityMode": False,
                    "provisional": True,
                },
            )
        ),
        "comparator": {
            "name": DEFAULT_COMPARATOR,
            "notes": (
                "Comparator includes ordinary background clinical care already "
                "represented in the model; it is not absence of all TB care."
            ),
        },
        "intervention": {
            "name": DEFAULT_INTERVENTION,
            "notes": "Systematic targeted LTBI screening with preventive treatment.",
        },
        "sourcesAndNotes": list(preset.get("sourcesAndNotes", [])),
        "scopeStatement": DIRECT_EFFECTS_SCOPE_STATEMENT,
    }
    if overrides:
        _deep_update(scenario, overrides)
    return scenario


def config_updates_from_scenario(scenario: dict[str, Any]) -> dict[str, Any]:
    screened = scenario.get("screened") or {}
    ltbi = scenario.get("ltbiPrevalenceAssumptions") or {}
    risk = scenario.get("riskFactorAssumptions") or {}
    targeting = scenario.get("targetingCriteria") or {}
    age_source = scenario.get("ageProfileSource") or {}
    updates: dict[str, Any] = {
        "scenarioLabel": scenario.get("scenarioName"),
        "populationPresetId": scenario.get("populationPresetId"),
        "N": scenario.get("populationSize"),
        "screenWindow": scenario.get("screeningWindowYears"),
        "screeningWindowYears": scenario.get("screeningWindowYears"),
        "earlyProgressionPeriodYears": scenario.get("earlyProgressionPeriodYears", 2),
        "activeTBCalibrationHorizonYears": scenario.get("activeTBCalibrationHorizonYears", 2),
        "followHorizon": scenario.get("followUpHorizonYears"),
        "followUpHorizonYears": scenario.get("followUpHorizonYears"),
        "screenCoverage": _coverage_from_screened(screened, scenario.get("populationSize")),
        "screeningStrategy": targeting.get("strategy"),
        "ltbiPrevalence": ltbi.get("value"),
        "ltbiStateAssumptions": deepcopy(scenario.get("ltbiStateAssumptions", {})),
        "scenario": deepcopy(scenario),
    }
    age_path = age_source.get("path")
    if age_path:
        updates["ageDistributionFile"] = age_path
    if isinstance(risk.get("riskPrev"), dict):
        updates["riskPrev"] = deepcopy(risk["riskPrev"])
    if isinstance(risk.get("diseaseOR"), dict):
        updates["diseaseOR"] = deepcopy(risk["diseaseOR"])
    return {key: value for key, value in updates.items() if value is not None}


def scenario_from_config(config: dict[str, Any]) -> dict[str, Any]:
    existing = config.get("scenario")
    if isinstance(existing, dict):
        scenario = deepcopy(existing)
    else:
        scenario = build_scenario_contract(
            str(config.get("populationPresetId") or DEFAULT_POPULATION_PRESET_ID)
        )
    scenario["scenarioName"] = config.get("scenarioLabel", scenario.get("scenarioName"))
    scenario["populationPresetId"] = config.get(
        "populationPresetId", scenario.get("populationPresetId")
    )
    scenario["populationSize"] = config.get("N", scenario.get("populationSize"))
    scenario["screeningWindowYears"] = config.get(
        "screeningWindowYears", config.get("screenWindow", scenario.get("screeningWindowYears"))
    )
    scenario["earlyProgressionPeriodYears"] = config.get(
        "earlyProgressionPeriodYears", scenario.get("earlyProgressionPeriodYears", 2)
    )
    scenario["activeTBCalibrationHorizonYears"] = config.get(
        "activeTBCalibrationHorizonYears",
        scenario.get("activeTBCalibrationHorizonYears", 2),
    )
    scenario["followUpHorizonYears"] = config.get(
        "followUpHorizonYears", config.get("followHorizon", scenario.get("followUpHorizonYears"))
    )
    scenario.setdefault("screened", {})
    scenario["screened"]["proportion"] = config.get(
        "screenCoverage", scenario["screened"].get("proportion")
    )
    scenario.setdefault("targetingCriteria", {})
    scenario["targetingCriteria"]["strategy"] = config.get(
        "screeningStrategy", scenario["targetingCriteria"].get("strategy")
    )
    scenario.setdefault("ltbiPrevalenceAssumptions", {})
    scenario["ltbiPrevalenceAssumptions"]["value"] = config.get(
        "ltbiPrevalence", scenario["ltbiPrevalenceAssumptions"].get("value")
    )
    scenario["ltbiStateAssumptions"] = deepcopy(
        config.get("ltbiStateAssumptions", scenario.get("ltbiStateAssumptions", {}))
    )
    scenario.setdefault("riskFactorAssumptions", {})
    scenario["riskFactorAssumptions"]["riskPrev"] = deepcopy(config.get("riskPrev", {}))
    scenario["riskFactorAssumptions"]["diseaseOR"] = deepcopy(config.get("diseaseOR", {}))
    scenario["scopeStatement"] = DIRECT_EFFECTS_SCOPE_STATEMENT
    return scenario


def _coverage_from_screened(screened: dict[str, Any], population_size: Any) -> Any:
    if screened.get("proportion") is not None:
        return screened["proportion"]
    if screened.get("number") is None or not population_size:
        return None
    return float(screened["number"]) / float(population_size)


def _deep_update(base: dict[str, Any], overrides: dict[str, Any]) -> None:
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            _deep_update(base[key], value)
        else:
            base[key] = deepcopy(value)
