"""Preview and deliberately apply country (WHO snapshot) or local incidence data.

Applying country data changes only fields supported by WHO evidence: location
identity, the incidence series and its bounds, national population metadata and
source provenance. LTBI prevalence, age distribution, test accuracy, treatment
cascade, costs, DALYs and risk factors are never changed.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

from engine.profiles.demonstration import DEMONSTRATION_PROFILE_ID, DEMONSTRATION_PROFILE_LABEL
from engine.profiles.population_profile import (
    IncidenceData,
    Location,
    LocationKind,
    PopulationProfile,
    Provenance,
)
from engine.who_incidence.snapshot import IncidenceSnapshot


INCIDENCE_LINK_NOTE = (
    "Country incidence data currently describe the TB disease burden and trend. "
    "They are not yet used to infer infection pressure or transmission."
)
UNCHANGED_BY_COUNTRY_DATA = (
    "Simulated population size",
    "Age distribution",
    "LTBI prevalence",
    "Risk-factor prevalence and effects",
    "Screening test accuracy and treatment cascade",
    "Costs and DALY assumptions",
)
KEEP_CURRENT = "keep_current"
USE_NEW = "use_new"


class ConflictRequiresChoice(ValueError):
    """Raised when applying data would replace user-supplied values without a choice."""

    def __init__(self, fields: list[str]) -> None:
        self.fields = fields
        super().__init__("A choice is required for user-supplied value(s): " + ", ".join(fields))


@dataclass(frozen=True)
class ProposedChange:
    key: str
    label: str
    current: str
    proposed: str
    source: str
    conflict: bool = False

    def to_row(self) -> dict[str, Any]:
        return {
            "Field": self.label,
            "Current": self.current,
            "After applying": self.proposed,
            "Source": self.source,
            "Needs your choice": "Yes" if self.conflict else "",
        }


def preview_country_application(profile: PopulationProfile, snapshot: IncidenceSnapshot, iso3: str) -> list[ProposedChange]:
    summary = snapshot.country_summary(iso3)
    series = snapshot.series_for(iso3)
    source = f"WHO {snapshot.report_year}-round snapshot"
    user_incidence = profile.incidence.provenance in {Provenance.LOCAL_UPLOAD, Provenance.USER_DEFINED}
    user_location = profile.location.kind is LocationKind.SUBNATIONAL
    latest = summary["latest"]
    changes = [
        ProposedChange(
            "location",
            "Country or area",
            _location_text(profile.location),
            f"{summary['name']} ({iso3}), WHO region {summary['region']}",
            source,
            conflict=user_location,
        ),
        ProposedChange(
            "incidence",
            "Estimated TB incidence series",
            _series_text(profile.incidence),
            f"{len(series)} years" + (f", {series[0].year}-{series[-1].year}" if series else " (no estimates published)"),
            source,
            conflict=user_incidence,
        ),
        ProposedChange(
            "national_population",
            "National population (context only)",
            _population_text(profile.location),
            "Not published" if latest is None else f"{latest['population']:,.0f} ({latest['year']})",
            source,
        ),
    ]
    return changes


def apply_snapshot_country(
    profile: PopulationProfile,
    snapshot: IncidenceSnapshot,
    iso3: str,
    *,
    resolutions: dict[str, str] | None = None,
) -> PopulationProfile:
    """Return ``profile`` with the country's WHO incidence applied.

    Raises ``ConflictRequiresChoice`` when user-supplied incidence or location would
    be replaced and no resolution was given for that field.
    """
    resolutions = resolutions or {}
    changes = preview_country_application(profile, snapshot, iso3)
    unresolved = [change.label for change in changes if change.conflict and change.key not in resolutions]
    if unresolved:
        raise ConflictRequiresChoice(unresolved)
    summary = snapshot.country_summary(iso3)
    rows = snapshot.rows_for(iso3)
    latest_row = rows[-1]
    manifest = snapshot.manifest
    provenance = snapshot.provenance_summary()
    base_id, base_name = _base_identity(profile)
    updated = profile
    if resolutions.get("location", USE_NEW) == USE_NEW:
        updated = replace(
            updated,
            profile_id=f"{iso3.lower()}-{base_id}",
            name=f"{summary['name']}: {base_name}",
            location=Location(
                name=summary["name"],
                kind=LocationKind.COUNTRY,
                iso3=iso3,
                who_region=summary["region"],
                national_population=latest_row.number("population"),
                national_population_year=latest_row.year,
                national_population_source=f"WHO {snapshot.report_year}-round snapshot (UN Population Division estimates)",
            ),
        )
    if resolutions.get("incidence", USE_NEW) == USE_NEW:
        updated = replace(
            updated,
            incidence=IncidenceData(
                source=f"WHO TB burden estimates, Global Tuberculosis Report {snapshot.report_year} round",
                snapshot_id=snapshot.snapshot_id,
                provenance=Provenance.WHO_SNAPSHOT,
                series=snapshot.series_for(iso3),
                notes=INCIDENCE_LINK_NOTE,
                source_detail=tuple(
                    sorted(
                        {
                            "snapshotId": snapshot.snapshot_id,
                            "reportYear": snapshot.report_year,
                            "accessDate": provenance["accessDate"],
                            "dataSha256": snapshot.data_sha256,
                            "sourceUrl": provenance["sourceUrl"],
                            "citation": manifest.get("citation"),
                            "iso3": iso3,
                        }.items()
                    )
                ),
            ),
            data_vintage=f"WHO Global Tuberculosis Report {snapshot.report_year} round; snapshot {snapshot.snapshot_id}",
        )
    return updated


def remove_country_data(profile: PopulationProfile) -> PopulationProfile:
    """Return to the demonstration location without touching other inputs."""
    from engine.profiles.demonstration import build_demonstration_profile

    demo = build_demonstration_profile()
    return replace(
        profile,
        profile_id=demo.profile_id,
        name=demo.name,
        location=demo.location,
        incidence=demo.incidence,
        data_vintage=demo.data_vintage,
    )


def _base_identity(profile: PopulationProfile) -> tuple[str, str]:
    base_id, base_name = profile.profile_id, profile.name
    previous = profile.location.iso3
    if previous and base_id.startswith(f"{previous.lower()}-"):
        base_id = base_id[len(previous) + 1 :]
    if DEMONSTRATION_PROFILE_ID in base_id and not base_id.startswith(DEMONSTRATION_PROFILE_ID):
        base_id = base_id[base_id.index(DEMONSTRATION_PROFILE_ID) :]
    if base_name.startswith(f"{profile.location.name}: "):
        base_name = base_name[len(profile.location.name) + 2 :]
    if profile.demonstration and not base_id:
        base_id, base_name = DEMONSTRATION_PROFILE_ID, DEMONSTRATION_PROFILE_LABEL
    return base_id, base_name


def _location_text(location: Location) -> str:
    if location.kind is LocationKind.DEMONSTRATION:
        return "Demonstration (no country)"
    suffix = f" ({location.iso3})" if location.iso3 else ""
    return f"{location.name}{suffix}"


def _series_text(incidence: IncidenceData) -> str:
    if not incidence.series:
        return "None linked"
    first, last = incidence.data_year_range
    labels = {
        Provenance.WHO_SNAPSHOT: "WHO snapshot",
        Provenance.LOCAL_UPLOAD: "local file",
        Provenance.USER_DEFINED: "local file",
        Provenance.BUNDLED: "bundled",
    }
    return f"{len(incidence.series)} years, {first}-{last} ({labels[incidence.provenance]})"


def _population_text(location: Location) -> str:
    if location.national_population is None:
        return "Not recorded"
    return f"{location.national_population:,.0f} ({location.national_population_year})"
