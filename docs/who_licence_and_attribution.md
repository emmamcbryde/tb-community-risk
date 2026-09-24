# WHO data: licence, terms and attribution

Status: **reviewed; one question open** (see "Open question" below). Reviewed 2026-09-24.

## Sources consulted

| Source | Finding |
| --- | --- |
| WHO TB data page, `https://www.who.int/teams/global-programme-on-tuberculosis-and-lung-health/data` | Downloads are provided "in accordance with WHO's data policy and their use are subject to WHO's terms and conditions". No Creative Commons licence is named. |
| WHO data policy, `https://www.who.int/about/policies/publishing/data-policy` | Data are made available on terms that allow non-commercial, not-for-profit use for public health purposes. Detailed conditions are in the terms and conditions. |
| WHO dataset terms and conditions, `https://www.who.int/about/policies/publishing/data-policy/terms-and-conditions` | See the conditions below. |
| `GTB-TME/gtbreport2025` (fork `emmamcbryde/gtbreport2025`, commit `666088c`) | No LICENSE file, and no licence statement in the README. |

## Conditions that apply to this application (WHO dataset terms)

* **Permitted:** a royalty-free, worldwide, non-exclusive right to use, reproduce,
  extract, download, copy, distribute, display or include the datasets for public
  health purposes.
* **Attribution:** cite WHO, the dataset title, the year, the date of access, and
  acknowledge the countries that provided the underlying data.
* **Modification:** only minimal alteration of figures and tables, to match
  presentation style, is allowed. Any other alteration needs prior written WHO
  authorisation.
* **No commercial promotion:** the data must not be used with the promotion of a
  commercial enterprise, product or service.
* **No endorsement:** do not state or imply that WHO endorses the use, the product
  or any entity.
* **Name and emblem:** do not use the WHO name, any abbreviation of it, or the WHO
  emblem without prior written approval. This application uses the name only in
  source attribution and citation, and uses no WHO logo or emblem.
* **No warranty:** the data are provided without warranty, and the user is
  responsible for their use.

## How this application complies

* **Source:** the production data are the public WHO TB burden-estimates CSV
  (`https://extranet.who.int/tme/generateCSV.asp?ds=estimates`), downloaded
  2026-09-24.
* **Values unaltered:** published values are stored exactly as published text.
  The only changes are selecting columns, renaming them to the snapshot contract,
  sorting rows, and adding one derived status column. Nothing is imputed, smoothed
  or combined across report rounds. Trend estimates are separate analyses and are
  labelled as such.
* **Citation:** the manifest and the Evidence page carry the citation and the
  access date.
* **No endorsement:** the manifest and the Evidence page state that WHO has not
  reviewed or endorsed this application or its use of the data.
* **Non-commercial:** the application is a public, non-commercial planning tool.
* **No repository code reused:** no code from `gtbreport2025` is copied. The
  repository is used read-only, only to cross-check published values against the
  report-round analysis output (`inc_mort/analysis/est.rda`). The exported
  cross-check file is not committed.

## Required citation

> World Health Organization. WHO TB burden estimates (country level), Global Tuberculosis Report 2025 data.
> Geneva: WHO; 2025. Accessed 2026-09-24. Estimates are based on data reported by countries and areas to WHO.

## Absence of a repository licence

Without a licence, the code and any original content in `GTB-TME/gtbreport2025` are
"all rights reserved" by default. Being public on GitHub does not grant reuse
rights. This limits reuse of the repository's code and derived objects
independently of the WHO data terms. So:

* no repository code is used;
* the production snapshot is built from the public WHO CSV, not from `est.rda`;
* `est.rda` is only compared against locally for validation, and no values from it
  are redistributed.

## Open question requiring confirmation from WHO

Does storing a subset of published WHO columns under renamed headers (values
unchanged, one derived status column) count as "minimal alteration" under the
dataset terms? And is redistributing that subset inside this public non-commercial
application within the permitted right to "distribute ... or include the Datasets
... in other products for public health purposes"?

Our reading is yes to both, since the values are not altered. Written
confirmation would remove the ambiguity. Until then the manifest records
`terms.status = reuse_permitted_with_conditions` together with this open question.
