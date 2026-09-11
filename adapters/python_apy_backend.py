from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

from adapters.backend import JsonDict
from adapters.serialization import to_json_like
from engine.apy.config import build_default_config
from engine.apy.economics import (
    build_default_economics_config,
    build_economics_preset_dale2019_aud,
    build_economics_preset_kwab150,
    run_health_economics,
    run_health_economics_for_config,
)
from engine.apy.expected_value import run_expected_value
from engine.apy.results_bundle import build_results_bundle
from engine.apy.summary import summarise_numeric_rows
from engine.apy.runner import run_scenario, run_scenario_with_do_nothing
from engine.apy.validation import collect_validation_issues
import pandas as pd


class PythonApyBackend:
    """Experimental pure-Python APY v9 backend adapter.

    This adapter intentionally does not import or call MATLAB. MATLAB remains
    the reference backend while Python parity validation is expanded.
    """

    def __init__(self, root: Path) -> None:
        self.root = Path(root)

    def status(self) -> JsonDict:
        return {
            "name": "python_apy",
            "started": True,
            "abm_path": "",
            "error": "",
            "experimental": True,
            "matlabRequired": False,
        }

    def default_config(self) -> JsonDict:
        return to_json_like(build_default_config())

    def validate_config(self, config: JsonDict) -> JsonDict:
        return to_json_like(collect_validation_issues(_matlab_empty_to_none(config)))

    def run_scenario(self, config: JsonDict) -> JsonDict:
        return to_json_like(run_scenario(_matlab_empty_to_none(config)))

    def build_results_bundle(
        self,
        results: JsonDict,
        validation_report: JsonDict | None = None,
        economics: JsonDict | None = None,
    ) -> JsonDict:
        bundle = build_results_bundle(results)
        if validation_report is not None:
            bundle["validation"] = {"report": validation_report}
        if economics is not None:
            bundle["economics"] = economics
        return to_json_like(bundle)

    def run_scenario_bundle(
        self,
        config: JsonDict,
        validation_report: JsonDict | None = None,
        progress_callback: Callable[[dict[str, Any]], None] | None = None,
    ) -> JsonDict:
        clean_config = _matlab_empty_to_none(config)
        if str(clean_config.get("analysisMethod") or "agent_based") == "expected_value":
            out = _run_expected_value_with_bundle(clean_config)
        else:
            out = run_scenario_with_do_nothing(
                clean_config,
                progress_callback=progress_callback,
            )
        bundle = out["bundle"]
        if validation_report is not None:
            bundle["validation"] = {"report": validation_report}
        return to_json_like(bundle)

    def save_scenario(
        self,
        config: JsonDict,
        path: str,
        economics_config: JsonDict | None = None,
    ) -> JsonDict:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "contractVersion": "ltbi_screening_scenario_v1",
            "backend": "python_apy",
            "scenarioLabel": config.get("scenarioLabel", ""),
            "scenario": to_json_like(config.get("scenario")),
            "config": to_json_like(config),
            "economics": to_json_like(economics_config),
        }
        target.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        return {
            "filename": str(target),
            "saved": True,
            "backend": "python_apy",
        }

    def load_scenario(self, path: str) -> tuple[JsonDict, JsonDict, JsonDict]:
        source = Path(path)
        payload = json.loads(source.read_text(encoding="utf-8"))
        config = payload.get("config")
        if not isinstance(config, dict):
            raise ValueError("Scenario JSON does not contain a config object.")
        report = self.validate_config(config)
        load_info = {
            "filename": str(source),
            "contractVersion": payload.get("contractVersion", ""),
            "scenarioLabel": payload.get("scenarioLabel", ""),
            "backend": payload.get("backend", ""),
        }
        economics = payload.get("economics")
        return to_json_like(config), report, {**load_info, "economics": economics}

    def default_economics_config(self) -> JsonDict:
        return to_json_like(build_default_economics_config())

    def economics_preset_kwab150(self) -> JsonDict:
        return to_json_like(build_economics_preset_kwab150())

    def economics_preset_dale2019_aud(self, regimen: str | None = None) -> JsonDict:
        return to_json_like(build_economics_preset_dale2019_aud(regimen))

    def run_economics(self, results: JsonDict, economics_config: JsonDict) -> JsonDict:
        return to_json_like(
            run_health_economics(
                _matlab_empty_to_none(results),
                _matlab_empty_to_none(economics_config),
            )
        )

    def run_economics_for_config(
        self,
        config: JsonDict,
        economics_config: JsonDict,
    ) -> JsonDict:
        return to_json_like(
            run_health_economics_for_config(
                _matlab_empty_to_none(config),
                _matlab_empty_to_none(economics_config),
            )
        )


def _matlab_empty_to_none(value):
    if value == []:
        return None
    if isinstance(value, dict):
        return {
            key: _matlab_empty_to_none(item)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [_matlab_empty_to_none(item) for item in value]
    return value


def _run_expected_value_with_bundle(config: JsonDict) -> dict[str, Any]:
    from engine.apy.natural_history import run_do_nothing_summary

    results = run_expected_value(config)
    raw = _expected_raw_from_event_ledger(results)
    results["raw"] = raw
    results["summary"] = summarise_numeric_rows(raw)
    do_nothing = run_do_nothing_summary(results)
    bundle = build_results_bundle(results, do_nothing=do_nothing)
    return {
        "results": results,
        "doNothing": do_nothing,
        "bundle": bundle,
    }


def _expected_raw_from_event_ledger(results: dict[str, Any]) -> pd.DataFrame:
    ledger = results.get("eventLedger") or {}
    totals = ledger.get("replicateTotals")
    annual = ledger.get("annualEvents")
    if not isinstance(totals, pd.DataFrame) or totals.empty:
        raise ValueError("Expected-value result did not return event-ledger totals.")

    def value(event: str, arm: str = "intervention") -> float:
        subset = totals[(totals["eventName"] == event) & (totals["arm"] == arm)]
        if subset.empty:
            return 0.0
        return float(subset["value"].iloc[0])

    def annual_value(event: str, max_year: int, arm: str = "comparator") -> float:
        if not isinstance(annual, pd.DataFrame) or annual.empty:
            return 0.0
        subset = annual[
            (annual["eventName"] == event)
            & (annual["arm"] == arm)
            & (annual["modelYear"].astype(float) < float(max_year))
        ]
        return float(subset["value"].sum()) if not subset.empty else 0.0

    screened = value("screened")
    started = value("tpt_started_total")
    completed = value("tpt_completed_total")
    effectively_treated = value("infection_effectively_treated_total")
    prevented = value("active_tb_cases_prevented")
    active_20 = value("active_tb_cases", arm="comparator")
    active_2 = annual_value("active_tb_cases", 2, arm="comparator")
    false_positive_treated = value("tpt_started_false_positive")
    raw = {
        "rep": 1,
        "seed": None,
        "nScreened": screened,
        "nInfected": value("infected_at_baseline", arm="comparator"),
        "nRecentLTBIAtBaseline": value("recent_ltbi_at_baseline", arm="comparator"),
        "nRemoteLTBIAtBaseline": value("remote_ltbi_at_baseline", arm="comparator"),
        "nLatentAtScreen": value("latent_infected_at_screen"),
        "nRecentLatentAtScreen": value("recent_latent_at_screen"),
        "nRemoteLatentAtScreen": value("remote_latent_at_screen"),
        "nActiveAtScreen": value("active_tb_at_screen"),
        "nTruePositiveLatent": value("true_positive_latent"),
        "nTruePositiveRecent": value("true_positive_recent"),
        "nTruePositiveRemote": value("true_positive_remote"),
        "nTestPositiveActive": value("test_positive_active"),
        "nFalseNegativeLatent": value("false_negative_latent"),
        "nTestNegativeActive": value("test_negative_active"),
        "nTrueNegative": value("true_negative"),
        "nTPTEligible": value("tpt_eligible"),
        "nTPTStartedTruePositive": value("tpt_started_true_positive"),
        "nTPTStartedRecent": value("tpt_started_recent"),
        "nTPTStartedRemote": value("tpt_started_remote"),
        "nTPTStartedFalsePositive": false_positive_treated,
        "nTPTCompletedTruePositive": value("tpt_completed_true_positive"),
        "nTPTCompletedRecent": value("tpt_completed_recent"),
        "nTPTCompletedRemote": value("tpt_completed_remote"),
        "nTPTCompletedFalsePositive": value("tpt_completed_false_positive"),
        "nTPTADRstopTruePositive": value("tpt_adr_stop_true_positive"),
        "nTPTADRstopFalsePositive": value("tpt_adr_stop_false_positive"),
        "nTPTStoppedOtherTruePositive": value("tpt_other_stop_true_positive"),
        "nTPTStoppedOtherFalsePositive": value("tpt_other_stop_false_positive"),
        "nTestPositive": value("test_positive_total"),
        "nTestPositiveNonActive": value("test_positive_total") - value("test_positive_active"),
        "nIGRApos": value("test_positive_total") - value("test_positive_active"),
        "nFalsePositiveTests": value("false_positive"),
        "nFalsePositiveTestsBCG": value("false_positive_bcg"),
        "nFalsePositiveTestsNoBCG": value("false_positive_no_bcg"),
        "nFalsePositiveTreated": false_positive_treated,
        "nFalsePositiveTreatedBCG": 0.0,
        "nFalsePositiveTreatedNoBCG": false_positive_treated,
        "nFalsePositiveCompleted": value("tpt_completed_false_positive"),
        "nFalsePositiveCompletedBCG": 0.0,
        "nFalsePositiveCompletedNoBCG": value("tpt_completed_false_positive"),
        "nExcessFalsePositiveTestsDueToBCG": value("false_positive_due_to_bcg"),
        "nExcessCoursesStartedDueToBCG": 0.0,
        "nExcessCoursesCompletedDueToBCG": 0.0,
        "nExcessCoursesDueToBCG": 0.0,
        "nScreenedBCG": 0.0,
        "nScreenedNoBCG": screened,
        "nUninfectedScreenedBCG": 0.0,
        "nUninfectedScreenedNoBCG": value("uninfected_screened"),
        "nStartTPT": started,
        "nCompleteTPT": completed,
        "nADRstop": value("tpt_adr_stop_total"),
        "nStoppedOther": value("tpt_other_stop_total"),
        "nPartialCourses": value("tpt_partial_course_total"),
        "nTotalCoursesStarted": started,
        "nTotalCoursesCompleted": completed,
        "nCuredInfection": effectively_treated,
        "nCuredInfectionRecent": value("infection_effectively_treated_recent"),
        "nCuredInfectionRemote": value("infection_effectively_treated_remote"),
        "nCuredInfectionFull": value("infection_effectively_treated_full"),
        "nCuredInfectionPartial": value("infection_effectively_treated_partial"),
        "nPreventedActiveTB": prevented,
        "nPreventedActiveTBRecent": value("active_tb_cases_prevented_recent"),
        "nPreventedActiveTBRemote": value("active_tb_cases_prevented_remote"),
        "nPreventedActiveTBFull": 0.0,
        "nPreventedActiveTBPartial": 0.0,
        "nActiveBy2y": active_2,
        "nActiveBy20y": active_20,
        "completionRateObserved": _safe_fraction(completed, started),
        "adrRateObserved": _safe_fraction(value("tpt_adr_stop_total"), started),
        "partialCourseRateObserved": _safe_fraction(value("tpt_partial_course_total"), started),
        "propCoursesFalsePositive": _safe_fraction(false_positive_treated, started),
        "falsePositiveRateObservedBCG": None,
        "falsePositiveRateObservedNoBCG": _safe_fraction(value("false_positive_no_bcg"), value("uninfected_screened")),
        "bcgAttributableFalsePositiveRateObserved": None,
        "propFalsePositiveTestsDueToBCGAmongBCG": None,
        "propCoursesDueToBCGAmongFalsePositiveCoursesBCG": None,
        "NNS_cureInfection": _safe_fraction(screened, effectively_treated),
        "NNS_preventActiveTB": _safe_fraction(screened, prevented),
        "NNS_falsePositiveTreated": _safe_fraction(screened, false_positive_treated),
        "NNT_started_cureInfection": _safe_fraction(started, effectively_treated),
        "NNT_started_preventActiveTB": _safe_fraction(started, prevented),
    }
    return pd.DataFrame([raw])


def _safe_fraction(numerator: float, denominator: float) -> float | None:
    if denominator == 0:
        return None
    return float(numerator) / float(denominator)
