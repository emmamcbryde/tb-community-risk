# Individual TB-risk platform status

Updated2026-11-04: **Milestone2A implemented: independent audit and publication-only
Menzies2008 reference ready for review.** This is model reproduction and evidence
review, not clinical validation or external scientific sign-off.

Default behaviour remains `hazard_based_progression_v1` with clinical parameters
external-required. The separate `menzies_2008_compatibility_v1` must be selected
explicitly with `menzies_2008_published_reference_v1`. No historical-only package
is silently used by the modern model. `synthetic_demonstration_v1` contains invented
values; `contemporary_parameters_external_required_v1` contains no clinical defaults.

The independent M1 audit corrected five demonstrated gaps in unclassified CXR,
nested evidence status, conditional ranking intervals, typed-ID tie ordering and
complete parameter provenance. The modern numerical formulation and working
highest-RR default remain intact. Additive-excess and product remain sensitivities.

Historical arithmetic implements documented PPV/rate/linear lifetime equations and
retains raw invalid outputs separately. Exact BCG size mappings, original country/
NTM/immigration tables, IGRA, general concurrent-factor rules, conflicting recent
lifetime integration and treatment-benefit calculations remain unavailable.
All42 published disease outputs are accounted for; no full end-to-end reproduction
is claimed. The current public website is not assumed equivalent to the paper.

Both formulations concern infection present at baseline. Future acquisition,
reinfection, transmission, migration, recurrence, treatment strategies and changing
risk factors remain outside this milestone. No UI, country download, population
generation, Starsim, MATLAB, intervention or economic work is integrated.

Review: [audit](core_audit_m2a.md), [historical specification](historical_algorithm_m2a.md),
[evidence](parameter_evidence_m2a.csv), [reproduction](reproduction_m2a.md),
[formulation comparison](formulation_comparison_m2a.md),
[current public-information audit](current_tstin3d_audit_m2a.md),
[author questions](author_questions_m2a.md),
[future country-profile specification](country_profile_interface_m2a.md),
[validation](validation_m2a.md).

Next: scientific review of the unresolved historical mappings and publication
inconsistencies; obtain author clarification and licensed reference vectors;
then assess contemporary parameter evidence and independent cohort validation.
No author contact or contemporary parameter adoption is authorized automatically.
The model is for planning and sequencing, not denying care.
