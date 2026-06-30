# Why the APY Strategies Are Cost-Saving Under Current Assumptions

This note is an audit trail for the Word report economics wording. It uses the generated offline APY DALY/QALY/ICER outputs and does not rerun the model.

## Direct Cost-Saving Check

The current KWAB150-informed active TB disease cost input is AUD $19,080 per active TB case. Current included programme costs are testing plus preventive treatment. Programme setup, programme running, outreach, travel, staff time, follow-up support, and false-positive incremental costs are currently set to zero unless supplied separately.

| Strategy | Active TB cases prevented | TB disease costs averted | Included programme cost | Net cost vs baseline | Break-even cases prevented | Actual minus break-even | Cost-saving? |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Random | 4.0 | AUD $76,318 | AUD $56,362 | AUD $-19,956 | 2.95 | 1.05 | yes |
| LTBI-targeted | 11.0 | AUD $209,876 | AUD $62,155 | AUD $-147,721 | 3.26 | 7.74 | yes |
| Cure-targeted | 11.0 | AUD $209,876 | AUD $62,155 | AUD $-147,721 | 3.26 | 7.74 | yes |
| Prevention-targeted | 12.0 | AUD $228,955 | AUD $61,327 | AUD $-167,628 | 3.21 | 8.79 | yes |

## Interpretation

The strategies are cost-saving because the modelled cost of screening and preventive treatment is lower than the active TB disease costs averted. In break-even terms, each strategy only needs to prevent roughly 3 active TB cases to cover included programme costs, while the generated outputs prevent 4 to 12 active TB cases per 1,500-person cohort depending on the sequencing strategy.

This differs from Dale et al. for several reasons:

- The APY analysis uses a different population with higher baseline TB burden and higher expected yield than the onshore Australian migrant setting.
- The APY model is a closed-cohort, person-level prevention model with no emigration process.
- Prevented cases are direct person-level cases prevented; downstream transmission benefits are not included.
- Programme setup/running, outreach, travel, staff time, follow-up support, and false-positive incremental costs are currently omitted or set to zero.
- Dale et al. used a lifetime Markov structure and explicitly explored treatment utility decrement sensitivity.
- The APY DALY results depend strongly on TB case fatality risk and years-of-life-lost assumptions, which still require review.

## What Would Make the Strategies No Longer Cost-Saving?

Cost-saving status would be reduced or removed if added implementation costs exceed the current margin between TB disease costs averted and included programme costs. Approximate additional cost headroom before losing cost-saving status equals the absolute value of net cost vs baseline.

| Strategy | Additional cost headroom before no longer cost-saving |
| --- | ---: |
| Random | AUD $19,956 |
| LTBI-targeted | AUD $147,721 |
| Cure-targeted | AUD $147,721 |
| Prevention-targeted | AUD $167,628 |

If APY-specific programme setup, travel, outreach, community engagement, staff time, follow-up support, or false-positive management costs exceed these headroom amounts, the relevant strategy would no longer be cost-saving, although it could still remain cost-effective depending on DALYs/QALYs gained and the willingness-to-pay threshold.
