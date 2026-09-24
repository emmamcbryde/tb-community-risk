"""Source adapters that transform local source files into application records.

Adapters never download data. The maintainer obtains source files explicitly and
passes local paths. Two audited sources are supported:

* ``WhoBurdenCsvAdapter`` - the public WHO TB burden-estimates CSV
  (``https://extranet.who.int/tme/generateCSV.asp?ds=estimates``), which is the
  current retrospective series for the latest report round and needs no
  credentials.
* ``GtbReportEstimatesAdapter`` - ``inc_mort/analysis/est.rda`` in a pinned local
  checkout of the Global TB Report repository. Reading ``.rda`` needs the
  optional ``pyreadr`` package, which the application itself does not require.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import subprocess
from typing import Any, Iterable, Protocol

from engine.who_incidence.schema import IncidenceRecord, ValidationReport, parse_rows
from engine.who_incidence.snapshot import sha256_file


WHO_BURDEN_CSV_URL = "https://extranet.who.int/tme/generateCSV.asp?ds=estimates"
WHO_DATA_DICTIONARY_URL = "https://extranet.who.int/tme/generateCSV.asp?ds=dictionary"
GTB_REPORT_UPSTREAM = "https://github.com/GTB-TME/gtbreport2025"
GTB_REPORT_ESTIMATES_FILE = "inc_mort/analysis/est.rda"
GTB_REPORT_COUNTRY_FILE = "data/gtb/other/cty.rda"

WHO_BURDEN_CSV_COLUMN_MAP = {
    "iso3": "iso3",
    "country": "country",
    "year": "year",
    "e_inc_100k": "incidence_per_100k",
    "e_inc_100k_lo": "incidence_per_100k_lo",
    "e_inc_100k_hi": "incidence_per_100k_hi",
    "e_inc_num": "incident_cases",
    "e_inc_num_lo": "incident_cases_lo",
    "e_inc_num_hi": "incident_cases_hi",
    "e_pop_num": "population",
}
GTB_REPORT_EST_COLUMN_MAP = {
    "iso3": "iso3",
    "year": "year",
    "inc": "incidence_per_100k",
    "inc.lo": "incidence_per_100k_lo",
    "inc.hi": "incidence_per_100k_hi",
    "inc.num": "incident_cases",
    "inc.lo.num": "incident_cases_lo",
    "inc.hi.num": "incident_cases_hi",
    "pop": "population",
    "source.inc": "estimation_method",
}


@dataclass(frozen=True)
class SourceFile:
    path: str
    sha256: str
    role: str


class IncidenceSourceAdapter(Protocol):
    source_id: str

    def read_rows(self) -> list[dict[str, Any]]:
        """Return rows in the application column vocabulary."""

    def source_files(self) -> list[SourceFile]:
        """Return the local files read, with checksums."""

    def upstream(self) -> dict[str, Any]:
        """Return upstream repository / commit / URL metadata."""


def transform(
    adapter: IncidenceSourceAdapter,
    *,
    iso3_filter: Iterable[str] | None = None,
    rate_decimals: int | None = None,
) -> tuple[list[IncidenceRecord], ValidationReport]:
    rows = adapter.read_rows()
    wanted = {code.upper() for code in iso3_filter} if iso3_filter else None
    if wanted is not None:
        rows = [row for row in rows if str(row.get("iso3") or "").upper() in wanted]
    if rate_decimals is not None:
        rows = [_round_row(row, rate_decimals) for row in rows]
    return parse_rows(rows)


class WhoBurdenCsvAdapter:
    source_id = "who_tb_burden_estimates_csv"

    def __init__(self, csv_path: Path | str, *, report_year: int | None = None) -> None:
        self.csv_path = Path(csv_path)
        self.report_year = report_year

    def read_rows(self) -> list[dict[str, Any]]:
        with self.csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            missing = [column for column in WHO_BURDEN_CSV_COLUMN_MAP if column not in (reader.fieldnames or [])]
            if missing:
                raise ValueError(f"WHO burden CSV schema changed; missing column(s): {', '.join(missing)}.")
            return [_rename(row, WHO_BURDEN_CSV_COLUMN_MAP) for row in reader]

    def source_files(self) -> list[SourceFile]:
        return [SourceFile(self.csv_path.name, sha256_file(self.csv_path), "incidence estimates")]

    def upstream(self) -> dict[str, Any]:
        return {"repository": None, "commit": None, "url": WHO_BURDEN_CSV_URL}


class GtbReportEstimatesAdapter:
    source_id = "gtbreport_inc_mort_est_rda"

    def __init__(self, repo_root: Path | str) -> None:
        self.repo_root = Path(repo_root)

    def read_rows(self) -> list[dict[str, Any]]:
        try:
            import pyreadr  # optional maintainer-only dependency
        except ImportError as exc:
            raise RuntimeError("Reading .rda files requires the optional 'pyreadr' package.") from exc
        est = pyreadr.read_r(str(self.repo_root / GTB_REPORT_ESTIMATES_FILE))["est"]
        cty = pyreadr.read_r(str(self.repo_root / GTB_REPORT_COUNTRY_FILE))["cty"]
        names = dict(zip(cty["iso3"], cty["country"]))
        missing = [column for column in GTB_REPORT_EST_COLUMN_MAP if column not in est.columns]
        if missing:
            raise ValueError(f"est.rda schema changed; missing column(s): {', '.join(missing)}.")
        rows = []
        for raw in est[list(GTB_REPORT_EST_COLUMN_MAP)].to_dict(orient="records"):
            row = _rename(raw, GTB_REPORT_EST_COLUMN_MAP)
            row["country"] = names.get(row["iso3"], "")
            rows.append(row)
        return rows

    def source_files(self) -> list[SourceFile]:
        return [
            SourceFile(GTB_REPORT_ESTIMATES_FILE, sha256_file(self.repo_root / GTB_REPORT_ESTIMATES_FILE), "incidence estimates"),
            SourceFile(GTB_REPORT_COUNTRY_FILE, sha256_file(self.repo_root / GTB_REPORT_COUNTRY_FILE), "country names"),
        ]

    def upstream(self) -> dict[str, Any]:
        return {"repository": GTB_REPORT_UPSTREAM, "commit": _git_head(self.repo_root), "url": None}


def parse_user_incidence_csv(text: str) -> tuple[list[IncidenceRecord], ValidationReport]:
    """Parse a user-supplied local or subnational incidence file (application schema).

    ``iso3`` may be blank for subnational areas; ``country`` then names the area.
    """
    reader = csv.DictReader(text.splitlines())
    return parse_rows(list(reader), require_iso3=False)


def _rename(row: dict[str, Any], mapping: dict[str, str]) -> dict[str, Any]:
    return {target: row.get(source) for source, target in mapping.items()}


def _round_row(row: dict[str, Any], decimals: int) -> dict[str, Any]:
    out = dict(row)
    for column in ("incidence_per_100k", "incidence_per_100k_lo", "incidence_per_100k_hi"):
        value = out.get(column)
        if isinstance(value, (int, float)) and value == value:
            out[column] = round(float(value), decimals)
    for column in ("incident_cases", "incident_cases_lo", "incident_cases_hi", "population"):
        value = out.get(column)
        if isinstance(value, (int, float)) and value == value:
            out[column] = float(round(float(value)))
    return out


def _git_head(path: Path) -> str | None:
    try:
        result = subprocess.run(
            ["git", "-C", str(path), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip() or None
