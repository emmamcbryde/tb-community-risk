# Frozen baseline worked examples

These are national data profiles at baseline 2024, using WHO report vintage 2025 and disease-incidence observations ending in 2024. They contain no modelled infection or intervention results.

| Country/area | Selected demographic source | Female | Male | Derived total | Year treatment |
| --- | --- | ---: | ---: | ---: | --- |
| Australia AUS | ABS Table 59, June 2024 release | 13,698,889 | 13,505,920 | 27,204,809 | Exact 2024 estimate, stock at 30 June |
| Viet Nam VNM | UN WPP 2024 Medium | 51,526,284 | 49,461,398 | 100,987,682 | Exact 2024 source value, classified as projection, stock at 1 July |
| PRK | UN WPP 2024 Medium | 13,390,111 | 13,108,724 | 26,498,835 | Exact 2024 source value, classified as projection, stock at 1 July |

Every profile retains 101 source age intervals per sex, including open 100+, and derives five display bands per sex. Neither interpolation nor year displacement was needed; numerical adjustment from the selected source age/sex counts is zero.

ABS female+male agrees exactly with published Persons in all age cells. WPP published sex columns have independent rounding: Viet Nam's sum of published both-sex ages is four people above the female+male total, while PRK's is one person below. Maximum cell differences are one person and median absolute cell differences zero. The active total always comes from female+male; the published total is not averaged in.

Australia's unselected WPP alternative is retained in the source catalogue: 26,713,198 derived people versus ABS 27,204,809. The difference is 491,611 people. This reflects distinct source vintages/definitions and is not a reconciliation error to average away. Australia has no supplied reviewed user alternative, so ABS wins the approved hierarchy.

AUS and VNM have complete 2000–2024 WHO incidence series with uncertainty. PRK has 25 country/year records but no publicly published incidence estimates; its preview shows missing central/lower/upper values and an explicit unavailable-incidence warning. No rate of zero is substituted.

ANT is a separate historical WHO entity with records only for 2000–2009; no 2024 demographic population is assigned to it. It is not merged into CUW or SXM. The country selector exposes the incomplete history, and requesting an unsupported demographic baseline fails.

For 2025–2030, the installed WPP Medium scenario can be selected explicitly and is the default for those years. An Australia 2025 request selects WPP because the installed ABS snapshot contains 2024 only; source-resolution diagnostics retain that change. A 2031 request fails until a new explicit frozen source supplies it.

Generate read-only CSV/JSON/Markdown previews with scripts/preview_country_tb_data.py. Exact output folders, hashes and test results are recorded in the milestone status and validation record.
