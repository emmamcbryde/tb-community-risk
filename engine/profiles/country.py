"""Attach country or local incidence series to a population profile."""

from __future__ import annotations

from dataclasses import replace

from engine.profiles.population_profile import (
    IncidenceData,
    IncidencePoint,
    Location,
    LocationKind,
    PopulationProfile,
    Provenance,
)
from engine.who_incidence.schema import IncidenceRecord
from engine.who_incidence.snapshot import IncidenceSnapshot


INCIDENCE_LINK_NOTE = (
    "Displayed for review. In this version the incidence series does not yet drive "
    "infection pressure, calibration or the dynamic model."
)


def apply_snapshot_country(profile: PopulationProfile, snapshot: IncidenceSnapshot, iso3: str) -> PopulationProfile:
    """Return ``profile`` with the country's bundled incidence series attached.

    Other inputs are left unchanged, so demonstration values remain labelled as
    demonstration values rather than country-specific evidence.
    """
    records = snapshot.records_for(iso3)
    if not records:
        raise KeyError(f"{iso3} is not available in snapshot {snapshot.snapshot_id}.")
    name = records[0].country or iso3
    manifest = snapshot.manifest
    base_id, base_name = profile.profile_id, profile.name
    previous = profile.location.iso3
    if previous and base_id.startswith(f"{previous.lower()}-"):
        base_id = base_id[len(previous) + 1 :]
    if previous and base_name.startswith(f"{profile.location.name}: "):
        base_name = base_name[len(profile.location.name) + 2 :]
    return replace(
        profile,
        profile_id=f"{iso3.lower()}-{base_id}",
        name=f"{name}: {base_name}",
        location=Location(name=name, kind=LocationKind.COUNTRY, iso3=iso3),
        incidence=IncidenceData(
            source=f"{manifest.get('sourceDataset')} (report year {manifest.get('sourceReportYear')})",
            snapshot_id=snapshot.snapshot_id,
            provenance=Provenance.BUNDLED,
            series=snapshot.series_for(iso3),
            notes=INCIDENCE_LINK_NOTE,
        ),
        data_vintage=f"{manifest.get('sourceDataset')} {manifest.get('sourceReportYear')}; snapshot {snapshot.snapshot_id}",
    )


def apply_user_incidence(profile: PopulationProfile, records: list[IncidenceRecord], *, source_label: str) -> PopulationProfile:
    """Attach a validated user-uploaded (local or subnational) incidence series."""
    areas = {(record.iso3, record.country) for record in records}
    if len(areas) != 1:
        raise ValueError("An uploaded incidence file must describe exactly one area.")
    iso3, area = next(iter(areas))
    series = tuple(
        IncidencePoint(
            year=record.year,
            estimate=record.incidence_per_100k,
            lower=record.incidence_per_100k_lo,
            upper=record.incidence_per_100k_hi,
        )
        for record in sorted(records, key=lambda item: item.year)
    )
    kind = LocationKind.COUNTRY if iso3 and not area else LocationKind.SUBNATIONAL
    return replace(
        profile,
        location=Location(name=area or iso3 or "User-defined area", kind=kind, iso3=iso3 or None),
        incidence=IncidenceData(
            source=f"User-defined file: {source_label}",
            snapshot_id=None,
            provenance=Provenance.USER_DEFINED,
            series=series,
            notes=INCIDENCE_LINK_NOTE,
        ),
        data_vintage=f"User-defined incidence file: {source_label}",
    )
