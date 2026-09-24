"""Lazy provenance package for the general application (built only on request)."""

from __future__ import annotations

import csv
from datetime import date
import io
import json
from pathlib import Path
import platform
import re
import subprocess
import sys
from typing import Any
import zipfile

from app.general.profile_editing import risk_factors_csv
from engine.profiles.effect_measures import crosswalk_rows, effect_warnings
from engine.profiles.population_profile import PopulationProfile, Provenance, user_override_fields
from engine.who_incidence.trend import TrendResult


REPO_ROOT = Path(__file__).resolve().parents[2]
LIMITATIONS = """# Limitations of this analysis

- Epidemiological, risk-factor, intervention, cost and DALY inputs are demonstration working
  defaults unless marked otherwise; they are not evidence for any particular country.
- Country incidence data describe the TB disease burden and trend. They are not yet used to infer
  infection pressure or transmission, and do not change the model's epidemiology.
- Results estimate direct effects for people screened and treated; transmission benefits are excluded.
- Stochastic simulation intervals show variation between simulated populations only. Trend
  uncertainty intervals reflect published WHO bounds only. Neither includes parameter or structural
  uncertainty.
- Risk-factor effect estimates are applied as progression-hazard multipliers whatever their declared
  measure type; odds ratios and risk ratios are not converted.
- WHO has not reviewed or endorsed this application or its use of the data.
"""


def package_filename(profile: PopulationProfile, analysis_label: str, today: date | None = None) -> str:
    today = today or date.today()
    place = profile.location.name if profile.location.kind.value != "demonstration" else "demonstration"
    if profile.incidence.provenance is Provenance.WHO_SNAPSHOT:
        vintage = f"who-gtb{dict(profile.incidence.source_detail).get('reportYear', '')}"
    elif profile.incidence.provenance is Provenance.LOCAL_UPLOAD:
        vintage = "local-incidence"
    else:
        vintage = "no-incidence"
    return f"{_slug(place)}_{_slug(analysis_label)}_{vintage}_{today.isoformat()}.zip"


def build_provenance_package(
    *,
    profile: PopulationProfile,
    trend: TrendResult | None,
    snapshot_manifest: dict[str, Any] | None,
    engine_config: dict[str, Any] | None,
    results_bundle: dict[str, Any] | None,
    economics_config: dict[str, Any] | None,
    economics_results: dict[str, Any] | None,
) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("population_profile.json", profile.to_json())
        archive.writestr("incidence_used.csv", _incidence_csv(profile))
        if trend is not None:
            archive.writestr("trend_fit.csv", _trend_csv(trend))
            payload = trend.to_dict()
            archive.writestr("trend_settings.json", _json({"method": payload["method"], "settings": payload["settings"]}))
            archive.writestr("trend_diagnostics.json", _json({k: v for k, v in payload.items() if k not in {"points", "settings"}}))
        if snapshot_manifest is not None and profile.incidence.provenance is Provenance.WHO_SNAPSHOT:
            archive.writestr("incidence_snapshot_manifest.json", _json(snapshot_manifest))
        archive.writestr("risk_factors.csv", risk_factors_csv(profile))
        archive.writestr("effect_measure_crosswalk.csv", _rows_csv(crosswalk_rows(profile)))
        archive.writestr("override_audit.csv", _override_csv(profile))
        archive.writestr("environment.json", _json(environment_record()))
        if engine_config is not None:
            archive.writestr("model_configuration.json", _json(engine_config))
        if results_bundle:
            headline = results_bundle.get("headline") or {}
            archive.writestr("results_key_metrics.csv", _rows_csv(headline.get("keyMetricsRows") or []))
            archive.writestr("results_metadata.json", _json(results_bundle.get("metadata") or {}))
        if economics_config is not None:
            archive.writestr("economics_assumptions.json", _json({"costItems": economics_config.get("costItems"), "discounting": economics_config.get("discounting"), "dalyAssumptions": economics_config.get("dalyAssumptions")}))
        if economics_results:
            archive.writestr("economics_summary.csv", _rows_csv([row for row in economics_results.get("summaryRows") or [] if row.get("discountProfile") == "primary"]))
        warnings = [item.to_dict() for item in effect_warnings(profile)]
        archive.writestr("LIMITATIONS.md", LIMITATIONS + ("\n## Effect-measure warnings\n\n" + "\n".join(f"- {w['message']}" for w in warnings) + "\n" if warnings else ""))
    return buffer.getvalue()


def environment_record() -> dict[str, Any]:
    versions = {}
    for name in ("streamlit", "pandas", "numpy", "scipy", "altair"):
        try:
            from importlib import metadata

            versions[name] = metadata.version(name)
        except Exception:
            versions[name] = None
    return {
        "softwareCommit": _git(["rev-parse", "HEAD"]),
        "workingTreeModified": bool(_git(["status", "--porcelain"])) if _git(["rev-parse", "HEAD"]) else None,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
        "packages": versions,
    }


def _incidence_csv(profile: PopulationProfile) -> str:
    rows = [
        {"year": p.year, "estimate_per_100k": p.estimate, "lower": p.lower, "upper": p.upper, "provenance": profile.incidence.provenance.value, "source": profile.incidence.source, "snapshot_id": profile.incidence.snapshot_id or ""}
        for p in profile.incidence.series
    ]
    return _rows_csv(rows, ["year", "estimate_per_100k", "lower", "upper", "provenance", "source", "snapshot_id"])


def _trend_csv(trend: TrendResult) -> str:
    rows = [
        {"year": p.year, "observed": p.observed, "lower": p.observed_lower, "upper": p.observed_upper, "fitted": p.fitted, "used_in_fit": p.used_in_fit, "in_period": p.in_window, "covid_disruption_year": p.covid_disruption}
        for p in trend.points
    ]
    return _rows_csv(rows, ["year", "observed", "lower", "upper", "fitted", "used_in_fit", "in_period", "covid_disruption_year"])


def _override_csv(profile: PopulationProfile) -> str:
    rows = [{"input": item, "status": "user-defined or local file"} for item in user_override_fields(profile)]
    return _rows_csv(rows, ["input", "status"])


def _rows_csv(rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> str:
    buffer = io.StringIO()
    names = fieldnames or (list(rows[0].keys()) if rows else ["empty"])
    writer = csv.DictWriter(buffer, fieldnames=names, lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow({key: "" if row.get(key) is None else row.get(key) for key in names})
    return buffer.getvalue()


def _json(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, default=str) + "\n"


def _slug(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "-", str(text)).strip("-").lower() or "profile"


def _git(args: list[str]) -> str | None:
    try:
        return subprocess.run(["git", "-C", str(REPO_ROOT), *args], check=True, capture_output=True, text=True, timeout=5).stdout.strip()
    except Exception:
        return None
