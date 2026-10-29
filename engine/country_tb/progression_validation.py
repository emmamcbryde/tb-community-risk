"""Explicit-RNG event-time validation and incidence comparison without fitting."""
from bisect import bisect_right
from collections import Counter
from math import exp,fsum,log1p,sqrt

from engine.country_tb.progression_common import artifact,verify,require,number,num,canonical
from engine.country_tb.data_common import digest
from engine.country_tb.progression_baseline import COHORT
from engine.country_tb.progression_forward import OPERATOR,RESULT,_preflight,_breaks
from engine.country_tb.progression_natural_history import interval_hazard,instantaneous_hazard
from engine.country_tb.progression_policies import mortality_rate,pressure_rate,reset_tau,PROTECTION

COMPARISON='country_tb.incidence_validation_comparison.v1'
MC='country_tb.active_tb_monte_carlo_validation.v1'


def _exponential(rng):
    value=rng.random();require(0<=value<1,'RNG must return values in [0,1).')
    return -log1p(-value)


def _one_person(row,op,horizon,rng,*,event_budget,recent_window):
    age=number(row['age']);m=number(row['effective_rr']);sex=row['sex'];group=row['subgroup_id'];clock=op['reinfection_clock']
    anchor=None if row['tau'] is None else number(row['tau']);actual_last=None if anchor is None else -anchor
    baseline_infected=anchor is not None;first=reinfections=0;prior=rng.random()<number(row['prior_remote_plus_recent_fraction'])
    if row['baseline_tbi_state']=='remote_only':prior=True
    t=0.;event_log=[];outcome='surviving_disease_free';outcome_time=horizon;case_tau=None
    for _ in range(event_budget):
        if t>=horizon:break
        candidates=_breaks(op,row,anchor,t,horizon)[1:]
        candidates += [number(r[k]) for r in op['future_pressure']['rows'] for k in ('start','end') if t<number(r[k])<=horizon]
        end=min(candidates);duration=end-t
        mu=mortality_rate(op['competing_mortality'],age+t,sex,group)
        rate=pressure_rate(op['future_pressure'],(t+end)/2,sex,group)
        if anchor is not None and clock['policy']==PROTECTION:rate*=number(clock['reinfection_susceptibility'])
        def integrated(dt):
            h=0. if anchor is None else m*interval_hazard(op['natural_history'],tau=anchor+t,age_current=age+t,
                age_at_infection=age+t-(anchor+t),sex=sex,duration=dt)
            return h+(mu+rate)*dt
        threshold=_exponential(rng);total=integrated(duration)
        if threshold>=total:t=end;continue
        lo,hi=0.,duration
        for _ in range(55):
            mid=(lo+hi)/2
            if integrated(mid)<threshold:lo=mid
            else:hi=mid
        t+=(lo+hi)/2
        tb=0. if anchor is None else m*instantaneous_hazard(op['natural_history'],tau=anchor+t,age_current=age+t,
            age_at_infection=age-anchor,sex=sex)
        draw=rng.random()*(tb+mu+rate)
        if draw<tb:
            outcome='first_active_TB_episode';outcome_time=t;case_tau=t-actual_last;break
        if draw<tb+mu:outcome='competing_death';outcome_time=t;break
        infected=anchor is not None
        if infected:
            prior=prior or t-actual_last>=recent_window;reinfections+=1
        else:first+=1
        new_tau=reset_tau(clock,0. if anchor is None else anchor+t,already_infected=infected)
        anchor=new_tau-t;actual_last=t
        event_log.append({'time_after_baseline':num(t),'event':'reinfection' if infected else 'first_infection',
            'progression_clock_after_event':num(new_tau)})
    else:raise ValueError('Monte Carlo event budget exhausted; no events silently discarded.')
    origin='following_reinfection' if reinfections else 'following_future_first_infection' if first else 'baseline_infection' if baseline_infected else 'uninfected'
    return {'cohort_id':row['cohort_id'],'sex':sex,'subgroup_id':group,'outcome':outcome,'outcome_time':num(outcome_time),
        'most_recent_infection_time_after_baseline':None if actual_last is None else num(actual_last),
        'time_since_most_recent_infection_at_case':None if case_tau is None else num(case_tau),
        'progression_clock_at_outcome':None if anchor is None else num(anchor+outcome_time),
        'event_tbi_state_at_case':None if case_tau is None else 'recent' if case_tau<recent_window else 'remote_only',
        'prior_remote_history':prior,'infection_events':event_log,'case_origin':origin if outcome=='first_active_TB_episode' else None,
        'effective_rr':row['effective_rr'],'highest_rr_responsibility_weights':row['responsibility_weights'],
        'first_TB_episode_count':int(outcome=='first_active_TB_episode')}


def monte_carlo_validate(cohort,operator,projection,*,rng,sample_size,horizon,standard_errors=5.,numerical_absolute_tolerance=.001,event_budget=10000):
    verify(cohort,COHORT);verify(operator,OPERATOR);verify(projection,RESULT)
    require(projection['baseline_cohort_hash']==cohort['profile_hash'] and projection['operator_hash']==operator['profile_hash'],'Monte Carlo and deterministic inputs must be identical.')
    require(rng is not None and callable(getattr(rng,'random',None)) and callable(getattr(rng,'getstate',None)),'An explicit stateful RNG is required; global RNG is forbidden.')
    require(type(sample_size) is int and 0<sample_size<=100000 and standard_errors>=3 and numerical_absolute_tolerance>=0,'Predeclared bounded Monte Carlo size/tolerance required.')
    horizon=number(horizon);p=next((p for p in projection['projections'] if number(p['horizon_years'])==horizon),None)
    require(p is not None,'Requested validation horizon is not projected.')
    eligible={r['cohort_id'] for r in p['rows'] if r['status']=='resolved' and number(r['risk_set_population'])>0}
    candidates=[r for r in cohort['rows'] if r['cohort_id'] in eligible]
    weights=[number(r['population'])-number(r['baseline_active_tb_population']) for r in candidates];total=fsum(weights)
    require(total>0,'No resolved disease-free population for validation.')
    cumulative=[];accumulated=0.
    for w in weights:accumulated+=w;cumulative.append(accumulated)
    initial=digest(canonical(rng.getstate()));audit=[]
    for _ in range(sample_size):
        index=min(len(candidates)-1,bisect_right(cumulative,rng.random()*total))
        audit.append(_one_person(candidates[index],operator,horizon,rng,event_budget=event_budget,recent_window=number(cohort['recent_window_years'])))
    checks=[]
    outcomes=Counter(r['outcome'] for r in audit)
    for name,field in [('first_active_TB_episode','expected_cases_covered'),('competing_death','competing_deaths_covered'),('surviving_disease_free','surviving_disease_free_covered')]:
        expected=number(p[field])/total;observed=outcomes[name]/sample_size
        se=sqrt(expected*(1-expected)/sample_size)
        # A one-count floor makes degenerate probabilities testable without inventing uncertainty.
        tolerance=standard_errors*max(se,1/sample_size)+numerical_absolute_tolerance
        checks.append({'outcome':name,'deterministic_probability':num(expected),'sampled_probability':num(observed),'MC_standard_error':num(se),
            'predeclared_absolute_tolerance':num(tolerance),'passed':abs(expected-observed)<=tolerance})
    years=[number(r['outcome_time']) for r in audit];mean=fsum(years)/sample_size
    se=sqrt(fsum((y-mean)**2 for y in years)/(sample_size*max(1,sample_size-1)))
    expected_years=number(p['disease_free_person_years_covered'])/total
    tolerance=standard_errors*max(se,horizon/sample_size)+numerical_absolute_tolerance*horizon
    checks.append({'outcome':'disease_free_person_time','deterministic_mean_years':num(expected_years),'sampled_mean_years':num(mean),
        'MC_standard_error':num(se),'predeclared_absolute_tolerance':num(tolerance),'passed':abs(expected_years-mean)<=tolerance})
    return artifact(MC,baseline_cohort_hash=cohort['profile_hash'],operator_hash=operator['profile_hash'],projection_hash=projection['profile_hash'],
        horizon_years=num(horizon),sample_size=sample_size,covered_population=num(total),checks=checks,passed=all(r['passed'] for r in checks),
        RNG_type=type(rng).__name__,RNG_initial_state_sha256=initial,RNG_final_state_sha256=digest(canonical(rng.getstate())),
        global_RNG_used=False,event_time_method='continuous_integrated_competing_hazard_inversion_not_renewal_grid_event_rounding',
        tolerance_policy={'MC_standard_errors':num(standard_errors),'numerical_absolute_tolerance':num(numerical_absolute_tolerance),
            'one_count_probability_floor':num(1/sample_size),'numerical_tolerance_is_not_a_confidence_interval':True},rows=audit)


def incidence_comparison(projection,target_contract,*,observation_id,prediction_measure='instantaneous_baseline',horizon=None):
    from engine.country_tb.active_tb_targets import validate_active_tb_targets
    verify(projection,RESULT);validate_active_tb_targets(target_contract)
    targets=[r for r in target_contract['rows'] if r['observation_id']==observation_id];require(len(targets)==1,'Exactly one supplied incidence observation must be selected.')
    target=targets[0];require(target['geography_id']==projection['geography_id'],'Incidence comparison geography mismatch.')
    warnings=['first_episode_operator_excludes_recurrence_prior_TB_depletion_imported_disease_migration_and_transmission',
        'WHO_incidence_and_notifications_are_distinct_no_observation_model_fitted','diagnostic_population_may_not_match_national_age_sex_risk_structure',
        'validation_only_not_calibration_inverse_TBI_inference_remains_blocked']
    if target['active_tb_target_type']!='incident_disease':warnings.append('target_is_not_incident_disease_incompatible_observation_process')
    if target['rate'] is not None:
        require(target['unit'] in {'cases_per_100000_person_year','cases_per_person_year'},'Prevalent disease cannot be compared as an incidence rate.')
        value=number(target['rate'])*(100000 if target['unit']=='cases_per_person_year' else 1)
    else:
        duration=number(target['observation_period']['duration_years']);require(duration>0,'Counts need an observation duration, not prevalent disease.')
        value=number(target['count'])/number(target['denominator'])/duration*100000
    if prediction_measure=='instantaneous_baseline':
        predicted=projection['instantaneous_baseline']['rate_per_100000_person_years_full'];unresolved=projection['instantaneous_baseline']['unresolved_population']
        partial=projection['instantaneous_baseline']['rate_per_100000_person_years_covered'];warnings.append('instantaneous_baseline_rate_is_not_calendar_year_average')
    else:
        require(prediction_measure=='followup_person_time_rate' and horizon is not None,'Explicit prediction units/horizon required.')
        p=next((p for p in projection['projections'] if number(p['horizon_years'])==number(horizon)),None);require(p is not None,'Unknown comparison horizon.')
        predicted=p['incidence_per_100000_person_years_full'];partial=p['incidence_per_100000_person_years_covered'];unresolved=p['unresolved_population']
        warnings.append('followup_cohort_person_time_rate_is_not_calendar_year_population_incidence')
    pred=None if predicted is None else number(predicted)
    return artifact(COMPARISON,status='validation_only_not_calibration',target=target,target_contract_hash=target_contract['profile_hash'],
        projection_hash=projection['profile_hash'],target_incidence_per_100000_person_year=num(value),predicted_incidence_per_100000_person_year=predicted,
        predicted_covered_population_rate=partial,difference=None if pred is None else num(pred-value),
        ratio=None if pred is None or value==0 else num(pred/value),ratio_unavailable_reason='zero_target_or_unresolved_prediction' if pred is None or value==0 else None,
        uncertainty=target['uncertainty'],scenario_hash=projection['operator_hash'],prediction_measure=prediction_measure,
        unresolved_population=unresolved,compatibility_warnings=warnings,model_exclusions=projection['exclusions'],
        fitted_parameters=None,changes_to_TBI_or_pressure=None,inverse_incidence_calibration='unavailable_pending_progression_and_risk_forward_operator')


def refinement_comparison(coarse,fine,*,absolute_risk_tolerance):
    verify(coarse,RESULT);verify(fine,RESULT)
    require(coarse['baseline_cohort_hash']==fine['baseline_cohort_hash'],'Refinement needs identical baseline support.')
    checks=[]
    for a,b in zip(coarse['projections'],fine['projections']):
        require(a['horizon_years']==b['horizon_years'],'Refinement horizons differ.')
        x=a['cumulative_risk_full'];y=b['cumulative_risk_full']
        error=None if x is None or y is None else abs(number(x)-number(y))
        checks.append({'horizon_years':a['horizon_years'],'absolute_risk_difference':None if error is None else num(error),
            'passed':None if error is None else error<=number(absolute_risk_tolerance)})
    assessed=[x for x in checks if x['passed'] is not None]
    return {'method':'explicit_numerical_step_refinement_not_statistical_uncertainty','coarse_hash':coarse['profile_hash'],'fine_hash':fine['profile_hash'],
        'predeclared_absolute_risk_tolerance':num(absolute_risk_tolerance),'checks':checks,'passed':all(x['passed'] for x in assessed) if assessed else None}
