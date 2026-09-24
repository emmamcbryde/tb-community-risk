"""Snapshot contract v2 for country-level WHO TB incidence.

Values are stored exactly as published text (no re-rounding or float
reformatting). Field definitions are recorded in every manifest.
"""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Any


SNAPSHOT_CONTRACT_VERSION = "who_incidence_snapshot_v2"
IMPORTER_VERSION = "who_incidence_importer_v2"

ISO3_PATTERN = re.compile(r"^[A-Z]{3}$")
ISO2_PATTERN = re.compile(r"^[A-Z]{2}$")
WHO_REGIONS = {"AFR", "AMR", "EMR", "EUR", "SEA", "WPR"}

# (snapshot column, WHO public CSV column or None when derived, kind, definition)
FIELDS: tuple[tuple[str, str | None, str, str], ...] = (
    ("iso3", "iso3", "code", "ISO 3166-1 alpha-3 country or area code (WHO)."),
    ("iso2", "iso2", "code", "ISO 3166-1 alpha-2 code (WHO)."),
    ("iso_numeric", "iso_numeric", "code", "ISO 3166-1 numeric code as published (text, leading zeros kept)."),
    ("country", "country", "text", "Country or area name as published by WHO."),
    ("who_region", "g_whoregion", "code", "WHO region: AFR, AMR, EMR, EUR, SEA or WPR."),
    ("year", "year", "year", "Calendar year of the estimate."),
    ("population", "e_pop_num", "number", "Estimated total population number (WHO, from UN Population Division)."),
    ("incidence_per_100k", "e_inc_100k", "number", "Estimated incidence of TB disease (all forms) per 100 000 population."),
    ("incidence_per_100k_lo", "e_inc_100k_lo", "number", "Low bound of estimated incidence per 100 000 population."),
    ("incidence_per_100k_hi", "e_inc_100k_hi", "number", "High bound of estimated incidence per 100 000 population."),
    ("incident_cases", "e_inc_num", "number", "Estimated number of incident TB cases (all forms)."),
    ("incident_cases_lo", "e_inc_num_lo", "number", "Low bound of estimated number of incident cases."),
    ("incident_cases_hi", "e_inc_num_hi", "number", "High bound of estimated number of incident cases."),
    (
        "hiv_positive_incidence_per_100k",
        "e_inc_tbhiv_100k",
        "number",
        "Estimated incidence of TB among people living with HIV, per 100 000 total population.",
    ),
    ("hiv_positive_incidence_per_100k_lo", "e_inc_tbhiv_100k_lo", "number", "Low bound of HIV-positive TB incidence per 100 000."),
    ("hiv_positive_incidence_per_100k_hi", "e_inc_tbhiv_100k_hi", "number", "High bound of HIV-positive TB incidence per 100 000."),
    (
        "incidence_status",
        None,
        "code",
        "Derived: 'estimated' when all three incidence-rate values are published, 'not_estimated' when none are.",
    ),
)
COLUMNS: tuple[str, ...] = tuple(field[0] for field in FIELDS)
WHO_COLUMN_MAP = {who: column for column, who, _, _ in FIELDS if who is not None}
REQUIRED_WHO_COLUMNS: tuple[str, ...] = tuple(who for _, who, _, _ in FIELDS if who is not None)
NUMERIC_COLUMNS: tuple[str, ...] = tuple(column for column, _, kind, _ in FIELDS if kind == "number")
RATE_TRIPLE = ("incidence_per_100k", "incidence_per_100k_lo", "incidence_per_100k_hi")
COUNT_TRIPLE = ("incident_cases", "incident_cases_lo", "incident_cases_hi")
HIV_TRIPLE = ("hiv_positive_incidence_per_100k", "hiv_positive_incidence_per_100k_lo", "hiv_positive_incidence_per_100k_hi")

MANIFEST_REQUIRED_FIELDS: tuple[str, ...] = (
    "contractVersion",
    "snapshotId",
    "snapshotKind",
    "isCompleteDataset",
    "reportYear",
    "source",
    "sourceFiles",
    "crossCheck",
    "importerVersion",
    "importerCommit",
    "transformation",
    "dataFile",
    "coverage",
    "fields",
    "citation",
    "terms",
    "validationSummary",
    "knownExclusions",
    "volatile",
)


@dataclass(frozen=True)
class SnapshotRow:
    """One country-year row; ``raw`` keeps published text, parsed values are floats."""

    raw: tuple[tuple[str, str], ...]

    def text(self, column: str) -> str:
        return dict(self.raw)[column]

    def number(self, column: str) -> float | None:
        text = dict(self.raw)[column]
        return None if text == "" else float(text)

    @property
    def iso3(self) -> str:
        return self.text("iso3")

    @property
    def year(self) -> int:
        return int(self.text("year"))

    @property
    def is_estimated(self) -> bool:
        return self.text("incidence_status") == "estimated"


def field_definitions() -> list[dict[str, Any]]:
    return [
        {"name": column, "whoSourceColumn": who, "kind": kind, "definition": definition}
        for column, who, kind, definition in FIELDS
    ]


def rounding_half_unit(text: str) -> float:
    """Half of the last published digit's unit (e.g. '6.8' -> 0.05, '1300' -> 50)."""
    text = text.strip().lstrip("+-")
    if "e" in text.lower():
        mantissa, exponent = text.lower().split("e")
        return rounding_half_unit(mantissa) * 10 ** int(exponent)
    if "." in text:
        return 0.5 * 10 ** (-len(text.split(".", 1)[1]))
    stripped = text.rstrip("0")
    zeros = len(text) - len(stripped) if stripped else 0
    return 0.5 * 10**zeros
