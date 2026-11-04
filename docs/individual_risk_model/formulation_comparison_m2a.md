# Historical versus hazard arithmetic

These are controlled numerical comparisons, not adopted clinical parameters. The annual historical number is treated as a constant hazard only for this explicitly labelled mathematical sensitivity. Actual hazard-model examples keep the existing independent synthetic curves and values.

| Scenario | Horizon years | RR | Historical linear raw | Exponential mathematical sensitivity | Difference |
| --- | ---: | ---: | ---: | ---: | ---: |
| low_short | 1 | 1 | 0.0005 | 0.000499750083313 | 2.49916687496e-07 |
| low_long | 45 | 1 | 0.0225 | 0.0220012590835 | 0.00049874091655 |
| moderate | 20 | 1 | 0.1 | 0.090634623461 | 0.00936537653899 |
| high_rr_short | 1 | 74 | 0.074 | 0.0713283061587 | 0.00267169384129 |
| high_rr_long | 45 | 74 | 3.33 | 0.964206894932 | 2.36579310507 |
| near80 | 1 | 1 | 0.001 | 0.000999500166625 | 4.99833374992e-07 |
| above80_literal | -5 | 1 | -0.005 | - | - |
| multiple_highest_applicable_rr | 45 | 3.6 | 0.162 | 0.14955879546 | 0.0124412045402 |
| multiple_additive_excess_rr | 45 | 5.6 | 0.252 | 0.222755261931 | 0.0292447380689 |
| multiple_multiplicative_rr | 45 | 10.8 | 0.486 | 0.384918192686 | 0.101081807314 |

At low hazard and short horizons the difference is small; at larger hazard/horizon it is material. RR 74 over45 years gives raw3.33, outside probability bounds, while the same-number exponential calculation remains below1. No cap is applied. Negative age-to80 durations remain explicitly invalid rather than becoming zero.

For recent contact/conversion, the initial annual number, the two-year rate footnote and the printed lifetime results conflict. Neither a constant early multiplier over life nor a modern exponential conversion resolves that documentation gap. The raw five-year schedule at tau0 is .033, compared with .075 if the initial early rate were extended; the reference does not silently choose a long-horizon interpretation.

The historical general multiple-factor rule remains unverified. Modern highest, additive-excess and multiplicative values are separate sensitivities and are not labelled the historical default.

Both models apply only to infection already present at baseline. They do not account for new acquisition, future reinfection, transmission, migration, recurrence or changing risk factors. Disagreement is not evidence that either model has clinical validity.
