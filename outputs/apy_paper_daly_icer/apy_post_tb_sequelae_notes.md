# APY Post-TB Sequelae Scenario Notes

Primary source: `docs/source_material/Post TB Mortality & Standardised Follow up_Byrne June 2026.pptx`.

Evidence summary from the Byrne June 2026 talk and linked literature:

- Global post-TB burden was presented as 58 million DALYs due to post-TB sequelae, around 47% of the combined incident-TB plus post-TB burden.
- The talk reports 12.1 DALYs per incident TB case in total, with 6.3 during active TB and 5.8 after TB treatment.
- Byrne highlighted increased long-term mortality after TB treatment, including a pooled SMR of 2.91 from Romanowski et al. and an Australian cause-mortality ratio of 1.16 with most deaths within 10 years.
- PTLD manifestations listed in the talk include fibrosis, bronchiectasis, emphysema/COPD, chronic respiratory failure, pulmonary hypertension, lung cancer, and chronic fungal disease.
- Disability domains highlighted in the talk include respiratory, mental health, renal, visual, hearing, and musculoskeletal impairment.

APY interpretation:

- The acute APY model remains unchanged. Post-TB morbidity and mortality are layered on as a scenario tail attached to prevented incident TB cases.
- Base post-TB assumptions remain uncertain for APY. PTLD prevalence is left unresolved by default and the workbook keeps that assumption editable.
- The default post-TB annual care cost is zero until APY-specific utilisation and costing inputs are agreed.

Default scenario settings: postTBDalysPerCase=5.8, postTBQalysLostPerCase=0.93, postTBDurationYears=10, postTBExcessMortalityMultiplier=1.16.

Illustrative prevention-targeted strategy under the default tail:

- Acute-only DALYs averted: 19.758.
- Acute plus post-TB DALY tail: 89.358.
- Acute-only QALYs gained: 15.84.
- Acute plus post-TB QALY tail: 27.

Caveats:

- These post-TB values are scenario extensions, not APY-specific measured outcomes.
- Excess post-TB mortality is recorded as an evidence input but is not used to alter the APY epidemiological simulation.
- Lifetime post-TB duration remains a sensitivity concept until a numeric APY-compatible duration assumption is agreed.
- The APY model should support planning and sequencing, not denying care.
