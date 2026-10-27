# Annual frozen data updates

## Explicit acquisition and publication

A maintainer acquires a report's source CSV and dictionary outside an analysis run and records report year, final observation year, acquisition date, source URL/repository and commit where applicable, original-file SHA-256, transformation/schema version, transformed checksum, citation and terms. A new report must have a separate content-addressed snapshot ID.

The pure validate_who_bundle checks acquired bytes using the same semantic rules as load_who_snapshot. freeze_who_bundle publishes validated bytes as a new immutable directory. It does not download anything or replace the current snapshot. The loader supports the existing flat WHO 2025 layout and explicit nested immutable folders.

Current WHO schema/importer support is who_incidence_snapshot_v2 / who_incidence_importer_v2. A changed schema or transformation requires a reviewed adapter/migration, not permissive parsing. Snapshot data must end in report_year−1, which is recorded independently of the report year. Baseline composition remains explicitly bound to GTR 2025 and year 2024 in this milestone.

## Compare before selection

compare_snapshots is a pure, non-mutating function. It reports country additions/removals, record-year additions/removals and per-country/year revisions to disease rates, counts, population and HIV-associated rates, including changed lower/upper uncertainty values and numerical differences. Missing/published transitions remain visible. It does not smooth, combine or migrate observations.

Before approving a new source, review historical revisions, territorial changes, uncertainty changes, availability and terms. Retain the old snapshot and comparison report. No anticipated future WHO report was downloaded in Milestone 2A.

## Preserve old analyses

An analysis records exact source IDs/hashes and the profile hash. Restoring it resolves those IDs, verifies their bytes and refuses an unsupported revision or missing source. Catalogue defaults do not affect restoration. Choosing a newer snapshot is an explicit new configuration and hash, with comparison available first.

Demographic updates follow the same immutable publication rule. Reviewed user counts can be imported with freeze_user_demography; official ABS/WPP import functions accept already acquired files and an explicit import date. Repeated identical imports cannot append records. Changed input or metadata cannot overwrite a prior snapshot; new user vintages have distinct provenance-based IDs.

Future demographic snapshots can add years or variants. Medium remains the explicit default projection scenario after 2024; alternatives remain named. No carry-forward or extrapolation fills a year absent from the installed frozen coverage.
