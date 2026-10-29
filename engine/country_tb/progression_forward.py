"""Absorbing first-episode TB with competing survival and exogenous resets.

Baseline-only kernels are analytic for constant hazard segments. Future pressure
uses an explicit Strang-split renewal grid; its refinement error is reported by
the validator, never presented as statistical uncertainty.
"""
from collections import defaultdict
from math import exp,expm1,fsum,floor,ceil

from engine.country_tb.progression_common import artifact,verify,require,metadata,number,num,quadrature,EXCLUSIONS,DataValidationError
from engine.country_tb.progression_natural_history import (
    NATURAL_HISTORY,select_curve,curve_hazard,curve_integral,interval_hazard,instantaneous_hazard,
)
from engine.country_tb.progression_policies import (
    EFFECT,CLOCK,MORTALITY,FUTURE,RESET,NO_RESET,PARTIAL,PROTECTION,mortality_rate,pressure_rate,reset_tau,
)
from engine.country_tb.progression_baseline import COHORT

OPERATOR='country_tb.active_tb_forward_operator.v1'
RESULT='country_tb.active_tb_projection_result.v1'
MODES={'baseline_cohort_only_v1','exogenous_future_infection_v1'}
TSI_BANDS=((0,1,'0-1'),(1,2,'1-2'),(2,5,'2-5'),(5,10,'5-10'),(10,25,'10-25'),(25,None,'25+'))


def forward_operator(*,natural_history,effect,clock,mortality,future,mode,step_years,horizons,
    source,provenance,review,diagnostic_only,recent_window_years='5'):
    metadata(source,provenance,review)
    require(type(diagnostic_only) is bool,'Explicit boolean diagnostic-only status required.')
    for value,family in ((natural_history,NATURAL_HISTORY),(effect,EFFECT),(clock,CLOCK),(mortality,MORTALITY),(future,FUTURE)):verify(value,family)
    require(mode in MODES and 0<number(step_years)<=1,'Explicit forward mode and bounded numerical step required.')
    require(number(recent_window_years)>0,'Positive recency window required.')
    require(effect['policy']=='provisional_rr_as_hazard_multiplier_v1' and effect['multimorbidity']=='highest_applicable_reviewed_rr_v1' and effect['cap'] is None,
        'Forward population policy must retain the uncapped highest-RR hazard approximation.')
    require(clock['policy']!=PARTIAL,'Partial reset is a single-event calculation-only sensitivity, not a population default.')
    require(isinstance(horizons,list) and horizons and all(number(h)>0 for h in horizons),'Positive projection horizons required.')
    times=sorted(set(number(h) for h in horizons));require(times[-1]<=100,'Focused forward contract supports at most 100 follow-up years.')
    require(times[-1]/number(step_years)<=10000,'Forward grid is bounded at 10,000 steps.')
    require(mortality['geography_id']==future['geography_id'],'Mortality and future pressure must match residence geography.')
    require(number(future['horizon_years'])>=times[-1],'Future pressure scenario does not cover the horizon.')
    if mode=='baseline_cohort_only_v1':require(future['scenario']=='zero','Baseline cohort progression forbids nonzero future pressure.')
    if not diagnostic_only:
        require(mortality['mode']!='no_mortality_diagnostic_v1','No-mortality calculations are mathematical diagnostics only.')
        require(not clock['calculation_only_sensitivity'],'Reinfection sensitivities cannot be central reviewed projections.')
        require(all(c is None or c['status']=='reviewed_external' for curves in [natural_history['group_curves'],*natural_history['sex_curves'].values()] for c in curves.values()),
            'Synthetic/surrogate curves require explicit diagnostic-only operation.')
    return artifact(OPERATOR,natural_history=natural_history,effect_application=effect,reinfection_clock=clock,
        competing_mortality=mortality,future_pressure=future,mode=mode,step_years=num(step_years),horizons=[num(x) for x in times],
        recent_window_years=num(recent_window_years),
        diagnostic_only=diagnostic_only,disease_event='first_active_TB_episode_after_baseline_absorbing',
        inverse_incidence_calibration='unavailable_pending_progression_and_risk_forward_operator',inverse_approval='independent_validation_required_not_unlocked_by_forward_code',
        numerical_method='analytic_piecewise_competing_hazard_kernel_positive_GL8_for_declining_tail_Strang_renewal_grid_for_future_infection',
        numerical_precision='binary64_full_precision_math_fsum_no_internal_reporting_rounding',
        exclusions=list(EXCLUSIONS),source=source,provenance=provenance,review=review)


def _breaks(op,row,tau,start,end):
    age=number(row['age']);points={start,end}
    points.update(x-age for x in range(int(floor(age+start))+1,int(floor(age+end))+1) if start<x-age<end)
    for r in op['competing_mortality']['rows']:
        for a in ('age_lower','age_upper_exclusive'):
            if r[a] is not None and start<number(r[a])-age<end:points.add(number(r[a])-age)
    if tau is not None:
        for curves in [op['natural_history']['group_curves'],*op['natural_history']['sex_curves'].values()]:
            for c in curves.values():
                if c:
                    points.update(number(k['time_years'])-tau for k in c['knots'] if start<number(k['time_years'])-tau<end)
        points.update(a-tau for a,_,_ in TSI_BANDS if start<a-tau<end)
    return sorted(points)


def _kernel(op,row,anchor,actual_anchor,start,end,cache):
    """Return survival, TB, deaths, person-years and event-age/time partitions."""
    age=number(row['age']);rr=number(row['effective_rr']);sex=row['sex'];group=row['subgroup_id']
    key=(age,rr,sex,group,anchor,actual_anchor,start,end)
    if key in cache:return cache[key]
    survive=1.;cases=deaths=py=0.;parts=[]
    points=set(_breaks(op,row,anchor,start,end))
    if actual_anchor is not None:
        points.update(_breaks(op,row,actual_anchor,start,end))
        cut=number(op['recent_window_years'])-actual_anchor
        if start<cut<end:points.add(cut)
    points=sorted(points)
    for a,b in zip(points,points[1:]):
        duration=b-a;mu=mortality_rate(op['competing_mortality'],age+a,sex,group)
        curve=None if anchor is None else select_curve(op['natural_history'],age+a,sex)
        t=None if anchor is None else anchor+a
        hazard=0. if curve is None else rr*curve_integral(curve,t,t+duration)
        varying=curve is not None and curve['tail']=='declining_to_remote_plateau_v1' and t>=number(curve['knots'][-1]['time_years'])
        if not varying:
            rate=hazard/duration+mu;event=-expm1(-rate*duration)
            y=exp(-rate*duration);years=event/rate if rate else duration
            tb=(hazard/duration)*years;death=mu*years
        else:
            def survival(x):return exp(-rr*curve_integral(curve,t,t+x)-mu*x)
            count=max(1,ceil(hazard+mu*duration));require(count<=10000,'Competing-risk quadrature budget exceeded; revise reviewed interval grid.')
            y=survival(duration);years=fsum(w*survival(x) for j in range(count)
                for x,w in quadrature(j*duration/count,(j+1)*duration/count))
            death=mu*years;tb=-expm1(-hazard-mu*duration)-death
            require(tb>=-1e-13 and death>=0,'Competing-risk quadrature failed positivity; refine interval.')
            require(tb>=0,'Negative TB mass cannot be clipped.')
        tau_mid=None if actual_anchor is None else actual_anchor+(a+b)/2
        clock_mid=None if anchor is None else anchor+(a+b)/2
        band='never_infected' if tau_mid is None else next(label for lo,hi,label in TSI_BANDS if lo<=tau_mid and (hi is None or tau_mid<hi))
        clock_band='never_infected' if clock_mid is None else next(label for lo,hi,label in TSI_BANDS if lo<=clock_mid and (hi is None or clock_mid<hi))
        state='never_infected' if tau_mid is None else 'recent' if tau_mid<number(op['recent_window_years']) else 'remote_only'
        parts.append((int(floor(age+(a+b)/2)),band,clock_band,state,survive*tb))
        cases+=survive*tb;deaths+=survive*death;py+=survive*years;survive*=y
    require(abs(survive+cases+deaths-1)<1e-10,'First-event/death/survival mass does not reconcile.')
    value=(survive,cases,deaths,py,parts)
    if len(cache)<200000:cache[key]=value
    return value


def _future_points(op,row,horizon):
    step=number(op['step_years']);points={0.,horizon,*[min(horizon,i*step) for i in range(1,int(horizon/step)+2)]}
    for r in op['future_pressure']['rows']:
        points.update(number(r[k]) for k in ('start','end') if 0<number(r[k])<horizon)
    return sorted(points)


def _preflight(op,row,horizon):
    reasons=list(row['unresolved_reasons'])
    known_age=number(row['age']) if row['age'] is not None else row.get('source_age_lower')
    if known_age is not None and known_age<15:
        try:select_curve(op['natural_history'],known_age,row['sex'])
        except DataValidationError as exc:reasons.append(str(exc))
    if reasons:return reasons
    age=number(row['age']);sex=row['sex'];group=row['subgroup_id'];tau=None if row['tau'] is None else number(row['tau'])
    scope=row['effect_age_applicability']
    if scope['age_upper_exclusive'] is not None and age+horizon>number(scope['age_upper_exclusive']):
        raise DataValidationError('effect_age_applicability_coverage_external_required')
    points=_future_points(op,row,horizon)
    for a,b in zip(points,points[1:]):pressure_rate(op['future_pressure'],(a+b)/2,sex,group)
    for a,b in zip(_breaks(op,row,None,0,horizon),_breaks(op,row,None,0,horizon)[1:]):mortality_rate(op['competing_mortality'],age+(a+b)/2,sex,group)
    positive_future=op['mode']=='exogenous_future_infection_v1' and any(pressure_rate(op['future_pressure'],(a+b)/2,sex,group)>0 for a,b in zip(points,points[1:]))
    if tau is not None or positive_future:
        interval_hazard(op['natural_history'],tau=tau or 0.,age_current=age,age_at_infection=age-(tau or 0.),sex=sex,duration=horizon)
    return []


def _trajectory(op,row,horizon,cache):
    tau=None if row['tau'] is None else number(row['tau'])
    initial_origin='uninfected' if tau is None else 'baseline_infection'
    states={(tau,initial_origin,tau):1.};out={'cases':0.,'deaths':0.,'person_years':0.,'first_infection_events':0.,'reinfection_events':0.}
    breakdown=defaultdict(float);age_parts=defaultdict(float);time_parts=defaultdict(float);clock_parts=defaultdict(float);state_parts=defaultdict(float)
    nofuture=op['mode']=='baseline_cohort_only_v1' or op['future_pressure']['scenario']=='zero'
    points=[0.,horizon] if nofuture else _future_points(op,row,horizon)
    def advance(a,b):
        updated={}
        for (anchor,origin,actual_anchor),mass in sorted(states.items(),key=lambda x:str(x[0])):
            y,tb,death,py,parts=_kernel(op,row,anchor,actual_anchor,a,b,cache)
            updated[anchor,origin,actual_anchor]=mass*y;out['cases']+=mass*tb;out['deaths']+=mass*death;out['person_years']+=mass*py
            breakdown[origin]+=mass*tb
            for age,band,clock_band,state,count in parts:
                age_parts[age]+=mass*count;time_parts[band]+=mass*count;clock_parts[clock_band]+=mass*count;state_parts[state]+=mass*count
        return updated
    for a,b in zip(points,points[1:]):
        if nofuture:states=advance(a,b);continue
        mid=(a+b)/2;states=advance(a,mid);updated=defaultdict(float)
        rate=pressure_rate(op['future_pressure'],mid,row['sex'],row['subgroup_id']);dt=b-a
        clock=op['reinfection_clock'];r1=rate*(number(clock['reinfection_susceptibility']) if clock['policy']==PROTECTION else 1.)
        for (anchor,origin,actual_anchor),mass in sorted(states.items(),key=lambda x:str(x[0])):
            infected=anchor is not None;r=r1 if infected else rate;none=exp(-r*dt)
            updated[anchor,origin,actual_anchor]+=mass*none
            if not infected:
                first=-expm1(-rate*dt)
                if r1==0:exactly_one=first
                elif rate==r1:exactly_one=rate*dt*exp(-rate*dt)
                elif abs((r1-rate)*dt)<100:exactly_one=rate*exp(-rate*dt)*(-expm1(-(r1-rate)*dt))/(r1-rate)
                else:exactly_one=rate/(r1-rate)*(exp(-rate*dt)-exp(-r1*dt))
                multiple=first-exactly_one
                require(multiple>=-1e-14,'Infection-event probability cancellation; refine interval.')
                if rate==0:events=0.
                else:events=first+r1*(dt-first/rate)
                out['first_infection_events']+=mass*first;out['reinfection_events']+=mass*(events-first)
                updated[-mid,'following_future_first_infection',-mid]+=mass*exactly_one
                if multiple>0:updated[-mid,'following_reinfection',-mid]+=mass*multiple
            else:
                new_tau=reset_tau(clock,anchor+mid,already_infected=True)
                updated[new_tau-mid,'following_reinfection',-mid]+=mass*(-expm1(-r*dt))
                out['reinfection_events']+=mass*r*dt
        states={k:v for k,v in updated.items() if v>0};states=advance(mid,b)
    out['surviving_disease_free']=fsum(states.values());out['by_infection_origin']=dict(breakdown)
    out['by_attained_age']=dict(age_parts);out['by_time_since_infection_band']=dict(time_parts)
    out['by_progression_clock_band']=dict(clock_parts);out['by_event_tbi_state']=dict(state_parts)
    out['population_accounting_residual']=1-out['cases']-out['deaths']-out['surviving_disease_free']
    require(abs(out['population_accounting_residual'])<1e-9,'Absorbing risk-set depletion failed.')
    return out


def _concentration(rows,population):
    available=[r for r in rows if r['status']=='resolved'];covered=fsum(number(r['risk_set_population']) for r in available)
    total=fsum(number(r['expected_cases']) for r in available);ranked=sorted(available,key=lambda r:(-number(r['cumulative_risk']),r['cohort_id']))
    out=[]
    for fraction in (.01,.05,.10):
        remaining=covered*fraction;cases=0.
        for r in ranked:
            take=min(remaining,number(r['risk_set_population']));cases+=take*number(r['cumulative_risk']);remaining-=take
            if remaining<=0:break
        out.append({'population_fraction':num(fraction),'population':num(covered*fraction),'expected_cases':num(cases),
            'fraction_of_expected_cases':num(cases/total) if total else None})
    return {'scope':'resolved_baseline_disease_free_population_ranked_by_expected_individual_risk',
        'covered_population':num(covered),'rows':out,'warning':'small_group_dominates_expected_cases' if total and number(out[0]['fraction_of_expected_cases'])>.5 else None,
        'warning_rule':'top_1_percent_of_covered_population_produces_more_than_half_of_expected_cases_no_RR_cap'}


def project_active_tb(cohort,operator):
    verify(cohort,COHORT);verify(operator,OPERATOR)
    require(all(isinstance(r.get('effect_age_applicability'),dict) for r in cohort['rows']),
        'Cohort lacks explicit effect-followup age applicability; rebuild from the validated source contract.')
    require(cohort['geography_id']==operator['competing_mortality']['geography_id'],'Forward geography mismatch; no Australian mortality fallback.')
    require(cohort['recent_window_years']==operator['recent_window_years'],'Baseline and follow-up recent-window policies must match.')
    require(not cohort['diagnostic_only'] or operator['diagnostic_only'],'Diagnostic baseline cannot become a reviewed country estimate.')
    op=operator;cache={};results=[];baseline_rates=[]
    for row in cohort['rows']:
        n=number(row['population']);risk=n-number(row['baseline_active_tb_population']);reasons=list(row['unresolved_reasons'])
        rate=None
        try:
            if not reasons:
                rate=0. if row['tau'] is None else number(row['effective_rr'])*instantaneous_hazard(op['natural_history'],tau=number(row['tau']),
                    age_current=number(row['age']),age_at_infection=number(row['age_at_most_recent_infection']),sex=row['sex'])
        except DataValidationError as exc:reasons.append(str(exc))
        baseline_rates.append({'cohort_id':row['cohort_id'],'risk_set_population':num(risk),'instantaneous_rate_per_person_year':None if rate is None else num(rate),
            'expected_cases_per_year_at_baseline':None if rate is None else num(risk*rate),'unresolved_reasons':reasons})
    for horizon in map(number,op['horizons']):
        rows=[];decomp={key:defaultdict(float) for key in ('subgroup','sex','baseline_tbi_state','event_tbi_state','effective_rr','highest_rr_responsible_factor',
            'attained_age','time_since_infection_band','progression_clock_band','infection_origin')}
        for source_row in cohort['rows']:
            risk=number(source_row['population'])-number(source_row['baseline_active_tb_population'])
            reasons=[]
            try:reasons=_preflight(op,source_row,horizon)
            except DataValidationError as exc:reasons.append(str(exc))
            base={k:source_row[k] for k in ('cohort_id','subgroup_id','sex','baseline_tbi_state','effective_rr')}
            base.update(risk_set_population=num(risk),baseline_active_TB_excluded=source_row['baseline_active_tb_population'])
            if reasons:
                rows.append({**base,'status':'unresolved','unresolved_reasons':sorted(set(reasons)),
                    'expected_cases':None,'cumulative_risk':None,'surviving_disease_free':None,'competing_deaths':None,'person_years':None})
                continue
            if risk==0:trajectory={'cases':0.,'deaths':0.,'person_years':0.,'surviving_disease_free':1.,'by_infection_origin':{},'by_attained_age':{},'by_time_since_infection_band':{},
                'by_progression_clock_band':{},'by_event_tbi_state':{},'first_infection_events':0.,'reinfection_events':0.,'population_accounting_residual':0.}
            else:trajectory=_trajectory(op,source_row,horizon,cache)
            cases=risk*trajectory['cases'];years=risk*trajectory['person_years']
            row={**base,'status':'resolved','unresolved_reasons':[],'expected_cases':num(cases),'cumulative_risk':num(trajectory['cases']),
                'surviving_disease_free':num(risk*trajectory['surviving_disease_free']),'competing_deaths':num(risk*trajectory['deaths']),
                'person_years':num(years),'incidence_per_100000_person_years':num(cases/years*100000) if years else None,
                'first_infection_events':num(risk*trajectory['first_infection_events']),'reinfection_events':num(risk*trajectory['reinfection_events']),
                'population_accounting_residual':num(risk*trajectory['population_accounting_residual']),
                'prior_remote_plus_recent_baseline_subset_population':num(number(source_row['population'])*number(source_row['prior_remote_plus_recent_fraction']))}
            rows.append(row)
            for dimension,field in [('subgroup','subgroup_id'),('sex','sex'),('baseline_tbi_state','baseline_tbi_state'),('effective_rr','effective_rr')]:decomp[dimension][source_row[field]]+=cases
            for factor,weight in source_row['responsibility_weights'].items():decomp['highest_rr_responsible_factor'][factor]+=cases*number(weight)
            if not source_row['responsibility_weights']:decomp['highest_rr_responsible_factor']['no_factor_or_not_attributable']+=cases
            for dimension,key in [('attained_age','by_attained_age'),('time_since_infection_band','by_time_since_infection_band'),
                ('progression_clock_band','by_progression_clock_band'),('event_tbi_state','by_event_tbi_state'),('infection_origin','by_infection_origin')]:
                for k,v in trajectory[key].items():decomp[dimension][str(k)]+=risk*v
        available=[r for r in rows if r['status']=='resolved'];unresolved=fsum(number(r['risk_set_population']) for r in rows if r['status']=='unresolved')
        cases=fsum(number(r['expected_cases']) for r in available);years=fsum(number(r['person_years']) for r in available)
        covered=fsum(number(r['risk_set_population']) for r in available);deaths=fsum(number(r['competing_deaths']) for r in available);survive=fsum(number(r['surviving_disease_free']) for r in available)
        results.append({'horizon_years':num(horizon),'rows':rows,'covered_risk_set_population':num(covered),'unresolved_population':num(unresolved),
            'complete_result_available':unresolved==0,'expected_cases_covered':num(cases),'expected_cases_full':num(cases) if unresolved==0 else None,
            'cumulative_risk_covered':num(cases/covered) if covered else None,'cumulative_risk_full':num(cases/covered) if covered and unresolved==0 else None,
            'surviving_disease_free_covered':num(survive),'competing_deaths_covered':num(deaths),'disease_free_person_years_covered':num(years),
            'incidence_per_100000_person_years_covered':num(cases/years*100000) if years else None,
            'incidence_per_100000_person_years_full':num(cases/years*100000) if years and unresolved==0 else None,
            'incidence_denominator':'integrated_alive_first_episode_disease_free_person_time_not_initial_population_times_horizon',
            'expected_case_decomposition':{k:[{'category':label,'expected_cases':num(v)} for label,v in sorted(values.items())] for k,values in decomp.items()},
            'factor_attribution':'nonexclusive_all_tied_highest_contributors_counts_must_not_be_added_to_unique_total',
            'extreme_risk_concentration':_concentration(rows,number(cohort['population']))})
    common=set.intersection(*[{r['cohort_id'] for r in p['rows'] if r['status']=='resolved'} for p in results])
    for p in results:
        shared=[r for r in p['rows'] if r['cohort_id'] in common]
        p['common_horizon_covered_population']=num(fsum(number(r['risk_set_population']) for r in shared))
        p['common_horizon_expected_cases']=num(fsum(number(r['expected_cases']) for r in shared))
    covered=fsum(number(r['risk_set_population']) for r in baseline_rates if r['instantaneous_rate_per_person_year'] is not None)
    cases=fsum(number(r['expected_cases_per_year_at_baseline']) for r in baseline_rates if r['expected_cases_per_year_at_baseline'] is not None)
    unresolved=fsum(number(r['risk_set_population']) for r in baseline_rates if r['instantaneous_rate_per_person_year'] is None)
    return artifact(RESULT,baseline_cohort_hash=cohort['profile_hash'],operator=op,operator_hash=op['profile_hash'],geography_id=cohort['geography_id'],
        population=cohort['population'],baseline_active_TB_excluded=cohort['baseline_active_TB_excluded'],diagnostic_only=op['diagnostic_only'],
        instantaneous_baseline={'rows':baseline_rates,'covered_risk_set_population':num(covered),'unresolved_population':num(unresolved),
            'expected_cases_per_year_covered':num(cases),'rate_per_100000_person_years_covered':num(cases/covered*100000) if covered else None,
            'rate_per_100000_person_years_full':num(cases/covered*100000) if covered and unresolved==0 else None},
        projections=results,complete_country_estimate=False,
        coverage_interpretation='covered_sets_can_shrink_with_horizon_use_common_horizon_counts_for_longitudinal_comparisons_no_cumulative_cases_added',
        infection_time_interpretation='time_since_infection_is_actual_most_recent_event_progression_clock_reported_separately_for_sensitivities',
        country_estimate_status='independent_validation_historical_TB_depletion_mortality_and_curve_transport_review_required',
        mortality_status=op['competing_mortality']['status'],exclusions=list(EXCLUSIONS),
        historical_conditioning_limitation=cohort['prior_TB_history'],
        horizon_interpretation='cumulative_first_episode_cases_from_same_absorbing_baseline_cohort_never_sum_across_horizons',
        inverse_incidence_calibration='unavailable_pending_progression_and_risk_forward_operator')


def single_reinfection_sensitivity(history,clock,*,tau,age,sex,multiplier,horizon):
    verify(clock,CLOCK);verify(history,NATURAL_HISTORY)
    new=reset_tau(clock,number(tau),already_infected=True)
    hazard=interval_hazard(history,tau=new,age_current=number(age),age_at_infection=number(age)-number(tau),sex=sex,duration=number(horizon))
    return {'policy':clock['policy'],'calculation_only_sensitivity':clock['calculation_only_sensitivity'],
        'effective_clock_after_event':num(new),'conditional_on_one_reinfection_risk':num(-expm1(-number(multiplier)*hazard)),
        'event_susceptibility':clock['reinfection_susceptibility'],'not_a_population_projection':True}
