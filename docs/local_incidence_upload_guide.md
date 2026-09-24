# Local or subnational incidence upload

Status: **implemented and validated**. Code: `engine/profiles/local_incidence.py`.

Use this to add estimated TB disease incidence for a country, district or other
area that is not in the WHO snapshot, or to use a local estimate instead of WHO's.
The file stays in your browser session and is never sent to an external service.

## Template

Download the template from *Set up population > Use a local or subnational
incidence file*. Columns:

| Column | Required | Meaning |
| --- | --- | --- |
| `location` | yes | Area name. One location per file. |
| `iso3` | no | ISO3 code; leave blank for subnational areas. |
| `year` | yes | One row per year. |
| `measure` | yes | Must be `estimated_incidence`. Notifications are rejected. |
| `incidence_per_100k` | yes | Estimated incidence per 100,000 per year. |
| `lower`, `upper` | no | Uncertainty bounds; give both or neither. |
| `population` | no | Area population for that year. |
| `source` | yes | Citation or description of the estimate. |
| `notes` | no | Free text. |

## Validation

The file is rejected if it:
* is not UTF-8;
* lacks required columns;
* has more than one location, duplicate years or non-numeric values;
* has negative values, only one of the two bounds, or bounds on the wrong side of
  the estimate;
* is missing a source;
* declares a measure other than `estimated_incidence`.

Warnings are shown when:
* bounds are missing, so propagated trend uncertainty is unavailable;
* the source cites WHO. The file is still recorded as a local file and never as
  the checksummed WHO snapshot.

## After applying

* **Provenance:** the profile records `provenance = local_upload`, the file name,
  the file SHA-256, a content hash, and the cited sources. These go into exported
  profiles and provenance packages.
* **Replacing existing incidence:** if the profile already has WHO or local
  incidence data, you must choose to keep it or replace it.
* **Use:** local incidence, like WHO incidence, is descriptive only and does not
  change the model's epidemiology.
