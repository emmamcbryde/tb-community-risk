# APY DALY/QALY Source Audit Notes

This audit separates morbidity and mortality sources for DALYs and QALYs.

DALYs use active TB YLD plus YLL from TB deaths prevented. YLD uses the active TB disability weight of 0.333 over 0.5 years. YLL uses the configurable TB case fatality risk and YLL per TB death.

Primary QALYs are Dale-compatible and mortality-adjusted. Morbidity uses the Dale/Bauer utility decrement from healthy/LTBI utility 0.8733 to active TB utility 0.8182 over the active TB duration. Mortality QALYs use TB deaths prevented multiplied by quality-adjusted life expectancy per TB death.

When qualityAdjustedLifeExpectancyPerTBDeath is not supplied, it is calculated as yllPerTBDeath multiplied by healthyUtility.

GBD-aligned QALYs are provided as a morbidity sensitivity. They use the active TB disability weight for morbidity and the same mortality QALY logic.

TB case fatality risk, YLL per TB death, and whether to prefer Dale/Bauer or GBD morbidity weights remain important sensitivity assumptions.
