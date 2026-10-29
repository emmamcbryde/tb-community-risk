"""Offline, invented M3B forward diagnostics. NOT national TB estimates."""
import argparse
from copy import deepcopy
from math import exp,expm1,fsum
from pathlib import Path
import random
import time

from engine.country_tb.progression_common import artifact,num,number,canonical
from engine.country_tb.data_common import REPO_ROOT,publish_immutable,digest,csv_bytes,DataValidationError
from engine.country_tb.progression_natural_history import audited_adult_candidate,curve_specification,natural_history,cumulative_hazard
from engine.country_tb.progression_policies import effect_application,reinfection_clock,competing_mortality,future_pressure,NO_RESET,PARTIAL,PROTECTION
from engine.country_tb.progression_baseline import baseline_cohort,cohort_from_allocation
from engine.country_tb.progression_forward import forward_operator,project_active_tb,single_reinfection_sensitivity
from engine.country_tb.progression_validation import monte_carlo_validate,incidence_comparison,refinement_comparison

SOURCE={'source_id':'synthetic:M3B','citation':'Invented forward-operator diagnostics; not country estimates'}
PROVENANCE={'source_id':'synthetic:M3B','diagnostic_only':True,'real_TBI_or_risk_prevalence_adopted':False}
REVIEW={'status':'reviewed','source_id':'synthetic:M3B','rationale':'Software diagnostic assumptions reviewed as invented inputs, not clinical evidence.'}
META={'source':SOURCE,'provenance':PROVENANCE,'review':REVIEW}


def diagnostic_row(identifier,*,population='1000',age='40',tau='.5',state='recent',rr='1',sex='female',subgroup='g'):
    return {'cohort_id':identifier,'subgroup_id':subgroup,'sex':sex,'population':population,'age':age,
        'baseline_tbi_state':state,'tau':tau,'effective_rr':rr,'responsibility_weights':{} if rr=='1' else {'synthetic_effect_'+rr:'1'},
        'prior_remote_plus_recent_fraction':'0','baseline_active_tb_population':'0','effect_status':'diagnostic_provisional',
        'pediatric_effect_approved':False,'unresolved_reasons':[]}


def make_cohort(rows,geography='who:AUS'):
    return baseline_cohort(rows,population=num(fsum(number(r['population']) for r in rows)),geography_id=geography,
        baseline_time='2024.5',diagnostic_only=True,active_TB_exclusion_policy='assumed_zero_baseline_active_TB_diagnostic_v1',**META)


def make_operator(history,*,geography='who:AUS',rate=None,clock=None,mortality=None,step='.25',horizons=None):
    future=future_pressure(geography_id=geography,scenario='zero' if rate is None else 'constant_at_2024',
        baseline_rate=rate,horizon_years='20',**META)
    return forward_operator(natural_history=history,effect=effect_application(**META),clock=clock or reinfection_clock(**META),
        mortality=mortality or competing_mortality(geography_id=geography,mode='no_mortality_diagnostic_v1',**META),
        future=future,mode='baseline_cohort_only_v1' if rate is None else 'exogenous_future_infection_v1',
        step_years=step,horizons=horizons or ['1','5','10','20'],diagnostic_only=True,**META)


def invented_conditional_support(allocation):
    # These are explicitly supplied diagnostic point ages/times, not inferred
    # midpoints or a reconstruction from a mean infection-time summary.
    support={}
    for group in allocation.get('subgroups',[allocation]):
        for r in group['rows']:
            if number(r['population']['value'])==0:continue
            for state,tau in [('never_infected',None),('recent','.5'),('remote_only','10' if r['age_lower']<45 else '30')]:
                key=group['subgroup_id']+':'+r['sex']+':'+str(r['age_lower'])+':'+state
                support[key]={**META,'nodes':[{'age':num(r['age_lower']+.25),'tau':tau,'weight':'1'}],
                    'assumption':'explicit_synthetic_conditional_age_and_time_support_selected_for_diagnostics'}
    return support


def _vietnam_allocation():
    from engine.country_tb.composition_common import population_demography
    from engine.country_tb.baseline_tbi_source import resolve_baseline_tbi
    from engine.country_tb.preview_risk_allocation import allocate_example,prevalence_record
    from engine.country_tb.risk_factor_catalog import progression_risk_factor_catalog
    from scripts.preview_country_tb_sources import direct_selection
    demo=population_demography([{'age_lower':age,'age_upper_exclusive':age+1,'sex':s,'population':'250'}
        for age in (30,60) for s in ('female','male')],geography_id='who:VNM',geography_level='national',population='1000',**META)
    states={'record_id':'invented_states','age_lower':0,'age_upper_exclusive':None,'sex':'both',
        'age_allocation_policy':'piecewise_constant_within_band_v1','any_tbi':'.2','recent':'.05','remote_only':'.15',
        'never_infected':'.8','prior_remote_plus_recent':'.01'}
    base=resolve_baseline_tbi(direct_selection('synthetic_vietnam',demo,[states]),demo);catalog=progression_risk_factor_catalog()
    records=[{**prevalence_record(catalog,'heavy_smoking',p,sex=s),**META} for s,p in [('female','.01'),('male','.20')]]
    return allocate_example(base,{'heavy_smoking':None},catalog=catalog,records=records)


def _who_target():
    from engine.country_tb.incidence import load_who_snapshot,WHO_2025
    from engine.country_tb.active_tb_targets import active_tb_targets
    snapshot=load_who_snapshot(WHO_2025);r=next(r for r in snapshot['rows'] if r['iso3']=='VNM' and r['year']=='2024')
    source={'source_id':WHO_2025,'citation':'WHO Global Tuberculosis Report 2025 frozen country incidence estimates, 2024',
        'url':'https://www.who.int/teams/global-programme-on-tuberculosis-and-lung-health/tb-reports/global-tuberculosis-report-2025'}
    provenance={'source_snapshot_hash':snapshot['source_data_hash'],'retained_raw_row':r}
    observation={'observation_id':'WHO:VNM:2024','active_tb_target_type':'incident_disease','rate':r['incidence_per_100k'],'count':None,
        'unit':'cases_per_100000_person_year','denominator':r['population'],'denominator_unit':'persons',
        'observation_period':{'start_time':'2024','end_time':'2025'},'calendar_year':2024,'age_lower':0,
        'age_upper_exclusive':None,'sex':'both','geography_id':'who:VNM','geography_level':'national',
        'uncertainty':{'kind':'interval','lower':r['incidence_per_100k_lo'],'upper':r['incidence_per_100k_hi'],'method':'WHO_reported_estimate_bounds_not_TBI_CI'},
        'source':source,'source_class':'WHO_estimate','ascertainment':{'status':'reviewed','description':'WHO estimate used only for diagnostic comparison; not notifications'},
        'review':REVIEW,'provenance':provenance,'observation_model':None}
    return active_tb_targets([observation],source=source,provenance=provenance,review=REVIEW)


def diagnostics(*,mc_sample_size=4000):
    from engine.country_tb.preview_risk_allocation import diagnostics as risk_diagnostics
    comparison_target=_who_target()
    constant_curve=curve_specification([{'time_years':'0','cumulative_risk':'0'},
        {'time_years':'10','cumulative_risk':num(-expm1(-.2))}],tail='fixed_remote_hazard_v1',
        tail_parameters={'hazard_per_year':'.02'},status='diagnostic_synthetic',**META)
    constant_nh=natural_history({'0-4':None,'5-14':None,'15+':constant_curve},**META)
    analytic_cohort=make_cohort([diagnostic_row('constant_hazard')]);analytic=project_active_tb(analytic_cohort,make_operator(constant_nh))
    central=audited_adult_candidate();tail=audited_adult_candidate(tail='fixed_remote_hazard_v1',tail_parameters={'hazard_per_year':'.00038'})
    simple=make_cohort([diagnostic_row('recent',population='500',tau='.1'),diagnostic_row('remote',population='500',tau='10',state='remote_only')])
    recent_remote=project_active_tb(simple,make_operator(tail))
    risk=make_cohort([diagnostic_row('RR_'+rr,population='250',tau='0',rr=rr) for rr in ('1','3.6','25.3','74')])
    risks=project_active_tb(risk,make_operator(tail))
    rare=make_cohort([diagnostic_row('rare_transplant',population='5',tau='0',rr='74'),
        diagnostic_row('never',population='995',tau=None,state='never_infected')])
    extremes=project_active_tb(rare,make_operator(tail))
    remote=make_cohort([diagnostic_row('remote',tau='10',state='remote_only')])
    reset=project_active_tb(remote,make_operator(tail,rate='.05',step='.125'))
    noreset=project_active_tb(remote,make_operator(tail,rate='.05',clock=reinfection_clock(policy=NO_RESET,**META),step='.125'))
    protection=project_active_tb(remote,make_operator(tail,rate='.05',clock=reinfection_clock(policy=PROTECTION,reinfection_susceptibility='.2',**META),step='.125'))
    finer=project_active_tb(remote,make_operator(tail,rate='.05',step='.0625'))
    refine=refinement_comparison(reset,finer,absolute_risk_tolerance='.0001')
    mortality=competing_mortality(geography_id='who:AUS',mode='external_age_sex_survival_v1',rows=[
        {'age_lower':0,'age_upper_exclusive':65,'sex':'both','subgroup_id':None,'survival_probability_one_year':num(exp(-.005))},
        {'age_lower':65,'age_upper_exclusive':None,'sex':'both','subgroup_id':None,'survival_probability_one_year':num(exp(-.08))}],**META)
    mc_projection=project_active_tb(remote,make_operator(tail,rate='.05',mortality=mortality,step='.125',horizons=['5']))
    mc=monte_carlo_validate(remote,mc_projection['operator'],mc_projection,rng=random.Random(31829),sample_size=mc_sample_size,horizon='5',
        standard_errors=5,numerical_absolute_tolerance=.0001)
    older=make_cohort([diagnostic_row('older',age='70',tau='30',state='remote_only')])
    old_without=project_active_tb(older,make_operator(tail));old_with=project_active_tb(older,make_operator(tail,mortality=mortality))
    allocation=risk_diagnostics()['examples']['Australian_mixture']
    australian=cohort_from_allocation(allocation,conditional_support=invented_conditional_support(allocation),
        active_TB_exclusion_policy='assumed_zero_baseline_active_TB_diagnostic_v1',**META)
    au=project_active_tb(australian,make_operator(tail,mortality=mortality))
    au_unextrapolated=project_active_tb(australian,make_operator(central))
    vnm_allocation=_vietnam_allocation();vnm=cohort_from_allocation(vnm_allocation,conditional_support=invented_conditional_support(vnm_allocation),
        active_TB_exclusion_policy='assumed_zero_baseline_active_TB_diagnostic_v1',**META)
    viet=project_active_tb(vnm,make_operator(tail,geography='who:VNM'));comparison=incidence_comparison(viet,comparison_target,observation_id='WHO:VNM:2024')
    children=make_cohort([diagnostic_row('young',population='250',age='3',tau='.5'),diagnostic_row('school',population='250',age='10',tau='.5'),
        diagnostic_row('adult',population='500',age='30',tau='.5')])
    child_result=project_active_tb(children,make_operator(central))
    # Hypothetical historical first-TB removal, NOT a fitted historical epidemic.
    original=make_cohort([diagnostic_row('ordinary',population='9500',age='30',tau='.5'),
        diagnostic_row('high',population='500',age='40',tau='10',state='remote_only',rr='74')])
    depleted=[]
    for r in original['rows']:
        H=cumulative_hazard(tail,tau=number(r['tau']),age_at_infection=number(r['age_at_most_recent_infection']),age_current=number(r['age']),sex=r['sex'])
        depleted.append({**r,'population':num(number(r['population'])*exp(-number(r['effective_rr'])*H))})
    conditional=make_cohort(depleted)
    history_sensitivity={'assumption':'hypothetical_remove_every_prebaseline_TB_event_no_treatment_survival_migration_or_frailty_model_not_a_country_correction',
        'original_population':original['population'],'hypothetical_never_TB_surviving_population':conditional['population'],
        'hypothetical_prior_TB_removed':num(number(original['population'])-number(conditional['population'])),
        'unconditioned':project_active_tb(original,make_operator(tail)), 'historical_depletion_sensitivity':project_active_tb(conditional,make_operator(tail))}
    tail_cases={}
    for policy,pars in [('constant_final_segment_hazard_v1',None),('declining_to_remote_plateau_v1',{'plateau_hazard_per_year':'.0001','decay_rate_per_year':'.1'}),
        ('fixed_remote_hazard_v1',{'hazard_per_year':'.00038'}),('no_extrapolation_v1',None)]:
        tail_cases[policy]=project_active_tb(remote,make_operator(audited_adult_candidate(tail=policy,tail_parameters=pars)))
    return artifact('country_tb.active_tb_forward_diagnostic_set.v1',diagnostic_only=True,
        warning='INVENTED TBI/RISK/TIME SUPPORT; NOT AUSTRALIAN OR VIET NAM INCIDENCE ESTIMATES; NO INVERSE CALIBRATION',
        natural_history_audit=central,projections={'constant_hazard_analytic':analytic,'recent_remote':recent_remote,'risk_multipliers':risks,'rare_high_RR':extremes,
            'reset':reset,'no_reset_sensitivity':noreset,'reinfection_protection_sensitivity':protection,
            'older_no_mortality':old_without,'older_synthetic_mortality':old_with,'Australian_synthetic_mixture':au,
            'Australian_no_extrapolation':au_unextrapolated,'Viet_Nam_diagnostic':viet,'children_external_required':child_result,
            'stochastic_competing_mortality':mc_projection},
        baseline_cohorts={'Australia':australian,'Viet_Nam':vnm,'reinfection':remote,'children':children,'analytic':analytic_cohort,
            'recent_remote':simple,'risk_multipliers':risk,'rare_high_RR':rare,'older':older,'historical_original':original,'historical_depleted':conditional},
        constant_hazard_reference={num(t):num(-expm1(-.02*t)) for t in (1,5,10,20)},
        incidence_comparison=comparison,tail_scenarios=tail_cases,numerical_refinement=refine,monte_carlo=mc,
        partial_reset_single_event=single_reinfection_sensitivity(tail,reinfection_clock(policy=PARTIAL,remaining_clock_fraction='.5',**META),
            tau='10',age='40',sex='female',multiplier='1',horizon='1'),historical_conditioning_sensitivity=history_sensitivity,
        source=SOURCE,provenance=PROVENANCE,review=REVIEW)


def summary_rows(value):
    return [{'scenario':name,'horizon_years':p['horizon_years'],'population':projection['population'],
        'expected_cases_covered':p['expected_cases_covered'],'expected_cases_full':p['expected_cases_full'],
        'unresolved_population':p['unresolved_population'],'competing_deaths_covered':p['competing_deaths_covered'],
        'incidence_per_100000_disease_free_person_years':p['incidence_per_100000_person_years_full'],'diagnostic_only':True}
        for name,projection in sorted(value['projections'].items()) for p in projection['projections']]


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=REPO_ROOT/'abm/output/country_tb_m3b_preview')
    parser.add_argument('--mc-samples',type=int,default=4000);args=parser.parse_args(argv)
    output=args.output.resolve()
    if not output.is_relative_to((REPO_ROOT/'abm/output').resolve()):raise DataValidationError('Developer diagnostics must stay under abm/output.')
    start=time.perf_counter();value=diagnostics(mc_sample_size=args.mc_samples);rows=summary_rows(value)
    if not value['numerical_refinement']['passed'] or not value['monte_carlo']['passed']:raise DataValidationError('Diagnostic numerical or Monte Carlo validation failed.')
    files={'diagnostics.json':canonical(value)+b'\n','summary.csv':csv_bytes(rows,tuple(rows[0])),
        'README.md':b'# Synthetic forward diagnostics\n\nNOT COUNTRY ESTIMATES. No incidence fitting. Child curves external-required. Adult four-knot curve is a diagnostic surrogate.\n'}
    files['manifest.json']=canonical(artifact('country_tb.active_tb_forward_export_receipt.v1',profile_hash_reference=value['profile_hash'],
        files={k:digest(v) for k,v in sorted(files.items())}))+b'\n'
    print(publish_immutable(output,'progression-diagnostics-'+value['profile_hash'][:12],files))
    for r in rows:print(r)
    print('WHO validation only:',{k:value['incidence_comparison'][k] for k in ('target_incidence_per_100000_person_year','predicted_incidence_per_100000_person_year','difference','ratio')})
    print('Monte Carlo passed:',value['monte_carlo']['passed'],'refinement:',value['numerical_refinement']['passed'])
    print(f'Runtime seconds: {time.perf_counter()-start:.3f}')


if __name__=='__main__':main()
