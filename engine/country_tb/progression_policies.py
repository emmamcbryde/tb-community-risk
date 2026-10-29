"""Explicit progression effects, clock, mortality and exogenous-pressure policies."""
from copy import deepcopy
from math import log
import re

from engine.country_tb.progression_common import artifact,verify,require,metadata,number,num,canonical

EFFECT='country_tb.progression_effect_application.v1'
CLOCK='country_tb.reinfection_progression_clock.v1'
MORTALITY='country_tb.competing_mortality.v1'
FUTURE='country_tb.future_exogenous_pressure.v1'
RESET='most_recent_reinfection_resets_progression_clock_v1'
NO_RESET='no_reinfection_clock_reset_sensitivity_v1'
PARTIAL='partial_reinfection_clock_reset_sensitivity_v1'
PROTECTION='partial_protection_against_reinfection_sensitivity_v1'


def effect_application(*,source,provenance,review):
    metadata(source,provenance,review)
    return artifact(EFFECT,policy='provisional_rr_as_hazard_multiplier_v1',multimorbidity='highest_applicable_reviewed_rr_v1',
        formula='1-exp(-m*incremental_baseline_cumulative_hazard)',original_effect_measure='relative_risk_not_automatically_hazard_ratio',
        qualification='modelling_approximation_bounded_probability_may_not_reproduce_original_RR_requires_scientific_review',
        cap=None,source=source,provenance=provenance,review=review)


def reinfection_clock(*,source,provenance,review,policy=RESET,remaining_clock_fraction=None,reinfection_susceptibility=None):
    metadata(source,provenance,review);require(policy in {RESET,NO_RESET,PARTIAL,PROTECTION},'Unknown reinfection progression clock.')
    require((policy==PARTIAL)==(remaining_clock_fraction is not None),'Partial reset fraction must be explicit and used only by that sensitivity.')
    require((policy==PROTECTION)==(reinfection_susceptibility is not None),'Reinfection protection must be explicit, not inferred from TBI.')
    if remaining_clock_fraction is not None:number(remaining_clock_fraction,probability=True)
    if reinfection_susceptibility is not None:number(reinfection_susceptibility,probability=True)
    return artifact(CLOCK,policy=policy,calculation_only_sensitivity=policy!=RESET,
        remaining_clock_fraction=None if remaining_clock_fraction is None else num(remaining_clock_fraction),
        reinfection_susceptibility=None if reinfection_susceptibility is None else num(reinfection_susceptibility),
        partial_reset_scope='single_event_calculation_only' if policy==PARTIAL else None,
        identifiable_from_total_TBI=False,RR_multiplied_by_reinfection=False,source=source,provenance=provenance,review=review)


def reset_tau(clock,tau,*,already_infected):
    if not already_infected:return 0.
    if clock['policy']==NO_RESET:return tau
    if clock['policy']==PARTIAL:return tau*number(clock['remaining_clock_fraction'])
    return 0.


def competing_mortality(*,geography_id,source,provenance,review,mode,rows=None,licensing=None,snapshot_hash=None):
    metadata(source,provenance,review)
    require(mode in {'no_mortality_diagnostic_v1','external_age_sex_survival_v1','reviewed_national_life_table_v1'},'Explicit mortality mode required.')
    records=deepcopy(rows or [])
    if mode=='no_mortality_diagnostic_v1':require(not records,'No-mortality mode cannot consume survival data.')
    else:
        require(records,'Mortality survival records required; life expectancy is not an age-specific survival curve.')
        for r in records:
            require({'age_lower','age_upper_exclusive','sex','subgroup_id','survival_probability_one_year'}<=r.keys(),'Incomplete age/sex/subgroup survival record.')
            lo=number(r['age_lower']);hi=None if r['age_upper_exclusive'] is None else number(r['age_upper_exclusive'])
            require(lo>=0 and (hi is None or hi>lo) and r['sex'] in {'female','male','both'},'Invalid survival age/sex.')
            p=number(r['survival_probability_one_year'],probability=True);require(p>0,'Zero annual survival needs an explicit terminal-event model, not infinite hazard.')
            r.update(age_lower=num(lo),age_upper_exclusive=None if hi is None else num(hi),
                survival_probability_one_year=num(r['survival_probability_one_year']),hazard_per_year=num(-log(p)))
        for i,a in enumerate(records):
            for b in records[i+1:]:
                if a['subgroup_id']!=b['subgroup_id']:continue
                require(not((a['sex']==b['sex'] or 'both' in {a['sex'],b['sex']}) and
                    (a['age_upper_exclusive'] is None or number(b['age_lower'])<number(a['age_upper_exclusive'])) and
                    (b['age_upper_exclusive'] is None or number(a['age_lower'])<number(b['age_upper_exclusive']))),'Overlapping survival rows at the same subgroup scope.')
    if mode=='reviewed_national_life_table_v1':
        require(isinstance(licensing,dict) and licensing.get('status')=='reviewed' and licensing.get('url') and licensing.get('terms')
            and isinstance(snapshot_hash,str) and re.fullmatch('[a-f0-9]{64}',snapshot_hash),
            'National mortality needs licensing review and an immutable source snapshot hash.')
    return artifact(MORTALITY,geography_id=geography_id,mode=mode,rows=sorted(records,key=canonical),licensing=licensing,snapshot_hash=snapshot_hash,
        status='competing_mortality_not_modelled' if mode=='no_mortality_diagnostic_v1' else 'externally_supplied_competing_survival',
        bias_if_unmodelled='overestimates_long_horizon_first_episode_TB_risk_and_disease_free_survivors',
        assumptions='external_background_cause_specific_nonTB_mortality_independent_of_latent_progression_conditional_on_strata',
        source=source,provenance=provenance,review=review)


def mortality_rate(policy,age,sex,subgroup):
    if policy['mode']=='no_mortality_diagnostic_v1':return 0.
    candidates=[r for r in policy['rows'] if r['sex'] in {'both',sex} and r['subgroup_id'] in {None,subgroup}
        and number(r['age_lower'])<=age and (r['age_upper_exclusive'] is None or age<number(r['age_upper_exclusive']))]
    specific=[r for r in candidates if r['subgroup_id']==subgroup]
    chosen=specific or [r for r in candidates if r['subgroup_id'] is None]
    require(len(chosen)==1,'mortality_age_sex_survival_external_required')
    return number(chosen[0]['hazard_per_year'])


def future_pressure(*,geography_id,source,provenance,review,scenario,horizon_years,rows=None,baseline_rate=None,
    annual_change=None,rate_bounds=None):
    metadata(source,provenance,review);h=number(horizon_years);require(h>0,'Positive future-pressure coverage required.')
    require(scenario in {'zero','constant_at_2024','user_supplied','bounded_decline','bounded_increase'},'Unknown explicit future pressure scenario.')
    records=deepcopy(rows or [])
    if scenario=='user_supplied':require(records and baseline_rate is None and annual_change is None and rate_bounds is None,'User pressure requires explicit intervals only.')
    else:
        require(not records,'Generated future scenario cannot silently combine user intervals.')
        if scenario=='zero':require(baseline_rate is None and annual_change is None and rate_bounds is None,'Zero pressure has no nonzero scenario parameters.');rate=0.
        else:
            require(baseline_rate is not None,'Explicit reviewed 2024 infection-event rate required; WHO proxy is not extrapolated.')
            rate=number(baseline_rate);require(rate>=0,'Negative infection pressure.')
        if scenario in {'bounded_decline','bounded_increase'}:
            require(annual_change is not None and rate_bounds is not None and len(rate_bounds)==2,'Bounded trends need a signed change and explicit rate bounds.')
            change=number(annual_change);lower,upper=map(number,rate_bounds)
            require(0<=lower<=rate<=upper and (change<0 if scenario=='bounded_decline' else change>0),'Invalid signed trend/bounds.')
            # Explicit annual staircase, rather than an undisclosed continuous interpolation.
            t=0.
            while t<h:
                value=max(lower,min(upper,rate+change*int(t)))
                records.append({'start':num(t),'end':num(min(t+1,h)),'rate_per_year':num(value),'sex':'both','subgroup_id':None});t+=1
        else:
            require(annual_change is None and rate_bounds is None,'Unused pressure trend parameters.')
            records=[{'start':'0','end':num(h),'rate_per_year':num(rate),'sex':'both','subgroup_id':None}]
    for r in records:
        require(set(r)=={'start','end','rate_per_year','sex','subgroup_id'} and r['sex'] in {'female','male','both'},'Complete exogenous calendar-offset intervals required.')
        require(0<=number(r['start'])<number(r['end']) and number(r['rate_per_year'])>=0,'Invalid future pressure interval.')
        for k in ('start','end','rate_per_year'):r[k]=num(r[k])
    for i,a in enumerate(records):
        for b in records[i+1:]:
            if a['subgroup_id']==b['subgroup_id'] and (a['sex']==b['sex'] or 'both' in {a['sex'],b['sex']}):
                require(number(a['end'])<=number(b['start']) or number(b['end'])<=number(a['start']),'Ambiguous future pressure intervals.')
    return artifact(FUTURE,geography_id=geography_id,scenario=scenario,horizon_years=num(h),rows=sorted(records,key=canonical),
        baseline_rate=None if baseline_rate is None else num(baseline_rate),annual_change=None if annual_change is None else num(annual_change),
        rate_bounds=None if rate_bounds is None else [num(x) for x in rate_bounds],
        trend_interpolation='explicit_annual_piecewise_constant_staircase',exogenous=True,endogenous_transmission=False,
        automatic_WHO_extrapolation=False,source=source,provenance=provenance,review=review)


def pressure_rate(policy,time,sex,subgroup):
    rows=[r for r in policy['rows'] if r['sex'] in {sex,'both'} and r['subgroup_id'] in {subgroup,None}
        and number(r['start'])<=time<number(r['end'])]
    specific=[r for r in rows if r['subgroup_id']==subgroup];chosen=specific or [r for r in rows if r['subgroup_id'] is None]
    require(len(chosen)==1,'future_pressure_coverage_external_required')
    return number(chosen[0]['rate_per_year'])
