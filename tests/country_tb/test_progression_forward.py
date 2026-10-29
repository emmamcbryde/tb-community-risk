"""Bounded, offline M3B mathematics and first-event process tests."""
from copy import deepcopy
from decimal import Decimal
from math import exp,expm1,fsum,log1p
from pathlib import Path
import ast
import json
import random
import socket
import subprocess

import pytest

from engine.country_tb.data_common import DataValidationError,profile_hash,canonical,REPO_ROOT,digest
from engine.country_tb.progression_common import number
from engine.country_tb.progression_common import IMPLEMENTATION,artifact
from engine.country_tb.progression_natural_history import (
    curve_specification,natural_history,audited_adult_candidate,risk_to_hazard,curve_hazard,curve_integral,
    cumulative_hazard,interval_hazard,instantaneous_hazard,
)
from engine.country_tb.progression_policies import (
    effect_application,reinfection_clock,competing_mortality,future_pressure,mortality_rate,pressure_rate,
    RESET,NO_RESET,PARTIAL,PROTECTION,reset_tau,
)
from engine.country_tb.progression_baseline import baseline_cohort,cohort_from_allocation
from engine.country_tb.progression_forward import forward_operator,project_active_tb,single_reinfection_sensitivity
from engine.country_tb.progression_validation import monte_carlo_validate,incidence_comparison,refinement_comparison

SOURCE={'source_id':'synthetic:M3B','citation':'Invented software validation only'}
PROVENANCE={'source_id':'synthetic:M3B','diagnostic_only':True}
REVIEW={'status':'reviewed','source_id':'synthetic:M3B','rationale':'Software fixture, no empirical effect approval.'}
META={'source':SOURCE,'provenance':PROVENANCE,'review':REVIEW}
GEO='who:AUS'


def constant_history(h=.02,*,children=False):
    def c(rate):return curve_specification([{'time_years':'0','cumulative_risk':'0'},
        {'time_years':'10','cumulative_risk':str(-expm1(-rate*10))}],tail='fixed_remote_hazard_v1',
        tail_parameters={'hazard_per_year':str(rate)},status='diagnostic_synthetic',**META)
    return natural_history({'0-4':c(.08) if children else None,'5-14':c(.03) if children else None,'15+':c(h)},**META)


def row(identifier='adult',*,n='1000',age='30',tau='.5',rr='1',state='recent',active='0',**kwargs):
    return {'cohort_id':identifier,'subgroup_id':'g','sex':'female','population':n,'age':age,'baseline_tbi_state':state,
        'tau':tau,'effective_rr':rr,'responsibility_weights':{} if rr=='1' else {'synthetic_highest_factor':'1'},
        'prior_remote_plus_recent_fraction':'0','baseline_active_tb_population':active,'effect_status':'diagnostic_provisional',
        'pediatric_effect_approved':False,'unresolved_reasons':[],**kwargs}


def cohort(rows=None):
    rows=rows or [row()]
    return baseline_cohort(rows,population=str(fsum(number(r['population']) for r in rows)),geography_id=GEO,
        baseline_time='2024.5',diagnostic_only=True,active_TB_exclusion_policy='assumed_zero_baseline_active_TB_diagnostic_v1',**META)


def operator(history=None,*,horizons=None,rate=None,mu=None,clock=None,step='.25'):
    mortality=competing_mortality(geography_id=GEO,mode='no_mortality_diagnostic_v1' if mu is None else 'external_age_sex_survival_v1',
        rows=None if mu is None else [{'age_lower':0,'age_upper_exclusive':None,'sex':'both','subgroup_id':None,
            'survival_probability_one_year':str(exp(-mu))}],**META)
    future=future_pressure(geography_id=GEO,scenario='zero' if rate is None else 'constant_at_2024',
        baseline_rate=None if rate is None else str(rate),horizon_years='20',**META)
    return forward_operator(natural_history=history or constant_history(),effect=effect_application(**META),
        clock=clock or reinfection_clock(**META),mortality=mortality,future=future,
        mode='baseline_cohort_only_v1' if rate is None else 'exogenous_future_infection_v1',
        step_years=step,horizons=horizons or ['1','5','10','20'],diagnostic_only=True,**META)


@pytest.mark.parametrize('risk',[0,.038,.05,.066,.072,.079,.999999])
def test_cumulative_risk_conversion(risk):
    assert number(risk_to_hazard(str(risk)))==pytest.approx(-float((Decimal(1)-Decimal(str(risk))).ln()),abs=1e-14)


@pytest.mark.parametrize('bad',[1,'1.1','-.01',None,True,'NaN'])
def test_invalid_cumulative_risks_not_clipped(bad):
    with pytest.raises(DataValidationError):risk_to_hazard(bad)


def test_adult_candidate_four_knots_and_unverified_two_year_rejected():
    nh=audited_adult_candidate();c=nh['group_curves']['15+']
    assert [k['time_years'] for k in c['knots']]==['0','1','5','10','25']
    assert nh['provenance']['two_year_005_knot']=='unverified_not_adopted'
    assert nh['group_curves']['0-4'] is None and nh['group_curves']['5-14'] is None
    previous=0.
    for knot in c['knots']:
        h=curve_integral(c,0,number(knot['time_years']))
        assert h>=previous and -expm1(-h)==pytest.approx(number(knot['cumulative_risk']),abs=1e-14);previous=h
    for segment in c['segments']:
        delta=number(segment['end'])-number(segment['start']);rate=number(segment['hazard_per_year'])
        assert rate>=0 and curve_integral(c,number(segment['start']),number(segment['end']))==pytest.approx(rate*delta)
    assert number(c['segments'][0]['hazard_per_year'])==pytest.approx(.03874082831643056)


@pytest.mark.parametrize('tail,pars',[('constant_final_segment_hazard_v1',None),
    ('fixed_remote_hazard_v1',{'hazard_per_year':'.00038'}),
    ('declining_to_remote_plateau_v1',{'plateau_hazard_per_year':'.0001','decay_rate_per_year':'.1'})])
def test_long_tail_explicit_positive_monotone(tail,pars):
    c=audited_adult_candidate(tail=tail,tail_parameters=pars)['group_curves']['15+']
    assert curve_integral(c,0,45)>=curve_integral(c,0,25)
    assert curve_hazard(c,45)>=0
    if tail=='declining_to_remote_plateau_v1':assert .0001<curve_hazard(c,45)<curve_hazard(c,25)


def test_no_extrapolation_and_malformed_curves_fail():
    c=audited_adult_candidate()['group_curves']['15+']
    with pytest.raises(DataValidationError,match='tail'):curve_integral(c,20,25.1)
    with pytest.raises(DataValidationError):curve_specification([{'time_years':0,'cumulative_risk':0},{'time_years':1,'cumulative_risk':.3},
        {'time_years':2,'cumulative_risk':.2}],tail='no_extrapolation_v1',status='diagnostic_synthetic',**META)
    with pytest.raises(DataValidationError):curve_specification(c['knots'],tail='declining_to_remote_plateau_v1',status='diagnostic_synthetic',**META)


def test_age_at_infection_and_attained_age_trajectory_and_common_sexes():
    nh=constant_history(children=True)
    for sex in ('female','male'):
        assert cumulative_hazard(nh,tau=20,age_at_infection=0,age_current=20,sex=sex)==pytest.approx(.08*5+.03*10+.02*5)
        assert interval_hazard(nh,tau=14,age_current=14,age_at_infection=0,sex=sex,duration=2)==pytest.approx(.03+.02)
    with pytest.raises(DataValidationError):cumulative_hazard(nh,tau=20,age_at_infection=1,age_current=20,sex='female')


@pytest.mark.parametrize('rr',['1','3.6','25.3','74'])
@pytest.mark.parametrize('horizon',['1','5','10','20'])
def test_analytic_constant_hazard_highest_multiplier_and_depletion(rr,horizon):
    result=project_active_tb(cohort([row(rr=rr)]),operator(horizons=[horizon]))
    p=result['projections'][0];risk=-expm1(-.02*float(rr)*float(horizon))
    assert number(p['cumulative_risk_full'])==pytest.approx(risk,abs=2e-13)
    assert 0<=number(p['cumulative_risk_full'])<=1
    assert number(p['expected_cases_full'])+number(p['surviving_disease_free_covered'])==pytest.approx(1000)
    assert number(p['incidence_per_100000_person_years_full'])==pytest.approx(.02*float(rr)*100000,rel=1e-11)


def test_recent_vs_remote_and_conditioning_is_incremental_hazard():
    nh=audited_adult_candidate(tail='fixed_remote_hazard_v1',tail_parameters={'hazard_per_year':'.00038'})
    result=project_active_tb(cohort([row('recent',tau='.1'),row('remote',tau='10',state='remote_only')]),operator(nh))
    p=result['projections'][0];a={r['cohort_id']:r for r in p['rows']}
    assert number(a['recent']['cumulative_risk'])>number(a['remote']['cumulative_risk'])
    expected=-expm1(-interval_hazard(nh,tau=10,age_current=30,age_at_infection=20,sex='female',duration=1))
    assert number(a['remote']['cumulative_risk'])==pytest.approx(expected)
    totals=[number(p['expected_cases_full']) for p in result['projections']]
    assert totals==sorted(totals) and totals[-1]<=2000


def test_competing_mortality_closed_form_and_reconciliation():
    result=project_active_tb(cohort([row(rr='3.6')]),operator(mu=.1))
    for p in result['projections']:
        t=number(p['horizon_years']);rate=.072+.1;event=-expm1(-rate*t)
        assert number(p['expected_cases_full'])==pytest.approx(1000*.072/rate*event)
        assert number(p['competing_deaths_covered'])==pytest.approx(1000*.1/rate*event)
        assert number(p['surviving_disease_free_covered'])==pytest.approx(1000*exp(-rate*t))
    without=project_active_tb(cohort([row(rr='3.6')]),operator())
    assert number(result['projections'][-1]['expected_cases_full'])<number(without['projections'][-1]['expected_cases_full'])
    assert without['mortality_status']=='competing_mortality_not_modelled'


def test_mortality_scope_licence_gaps_and_subgroup_override():
    records=[{'age_lower':0,'age_upper_exclusive':None,'sex':'both','subgroup_id':None,'survival_probability_one_year':'.99'},
        {'age_lower':0,'age_upper_exclusive':None,'sex':'female','subgroup_id':'g','survival_probability_one_year':'.9'}]
    policy=competing_mortality(geography_id=GEO,mode='external_age_sex_survival_v1',rows=records,**META)
    assert mortality_rate(policy,30,'female','g')>mortality_rate(policy,30,'male','g')
    with pytest.raises(DataValidationError):competing_mortality(geography_id=GEO,mode='reviewed_national_life_table_v1',rows=records,**META)
    with pytest.raises(DataValidationError):competing_mortality(geography_id=GEO,mode='external_age_sex_survival_v1',rows=[{**records[0],'survival_probability_one_year':'0'}],**META)


def test_children_external_required_block_complete_incidence_no_adult_fallback():
    pop=cohort([row('adult',n='600'),row('child',n='200',age='3'),row('school',n='200',age='10')])
    result=project_active_tb(pop,operator())
    for p in result['projections']:
        assert p['expected_cases_full'] is None and number(p['unresolved_population'])==400
        assert any('external_required_natural_history_0-4' in r['unresolved_reasons'] for r in p['rows'])
    with pytest.raises(DataValidationError,match='Adult'):cohort([row(age='3',rr='3.6')])
    synthetic=project_active_tb(pop,operator(constant_history(children=True)))
    assert synthetic['projections'][0]['expected_cases_full'] is not None and synthetic['diagnostic_only']


def test_baseline_active_TB_excluded_and_no_future_recurrence():
    result=project_active_tb(cohort([row(active='100')]),operator())
    assert result['baseline_active_TB_excluded']=='100'
    p=result['projections'][-1]
    assert number(p['expected_cases_full'])+number(p['competing_deaths_covered'])+number(p['surviving_disease_free_covered'])==pytest.approx(900)
    assert result['operator']['disease_event']=='first_active_TB_episode_after_baseline_absorbing'


@pytest.mark.parametrize('state,tau',[('never_infected',None),('recent','.5'),('remote_only','10')])
def test_all_uninfected_recent_remote_and_zero_future(state,tau):
    result=project_active_tb(cohort([row(state=state,tau=tau)]),operator())
    for p in result['projections']:
        if state=='never_infected':assert number(p['expected_cases_full'])==0
        assert number(p['surviving_disease_free_covered'])+number(p['expected_cases_full'])==pytest.approx(1000)
    zero=project_active_tb(cohort([row(state=state,tau=tau)]),operator(rate=0))
    assert [float(p['expected_cases_full']) for p in zero['projections']]==pytest.approx([float(p['expected_cases_full']) for p in result['projections']],abs=1e-9)


def test_future_first_infection_closed_form_and_origin_decomposition():
    population=cohort([row(tau=None,state='never_infected')]);op=operator(rate=.1,horizons=['1','5'],step='.0625')
    result=project_active_tb(population,op)
    for p in result['projections']:
        t=number(p['horizon_years']);a=.1;b=.02
        expected=1-(b*exp(-a*t)-a*exp(-b*t))/(b-a)
        assert number(p['cumulative_risk_full'])==pytest.approx(expected,abs=2e-6)
        assert fsum(number(r['expected_cases']) for r in p['expected_case_decomposition']['infection_origin'])==pytest.approx(number(p['expected_cases_full']))
        assert number(p['rows'][0]['first_infection_events'])<=1000 and number(p['rows'][0]['reinfection_events'])>0


def test_reinfection_reset_no_reset_partial_and_protection_are_explicit():
    history=audited_adult_candidate(tail='fixed_remote_hazard_v1',tail_parameters={'hazard_per_year':'.00038'})
    pop=cohort([row(tau='10',state='remote_only')])
    central=project_active_tb(pop,operator(history,rate=.1,horizons=['1','5']))
    noreset=project_active_tb(pop,operator(history,rate=.1,clock=reinfection_clock(policy=NO_RESET,**META),horizons=['1','5']))
    protected=project_active_tb(pop,operator(history,rate=.1,clock=reinfection_clock(policy=PROTECTION,reinfection_susceptibility='.2',**META),horizons=['1','5']))
    assert number(central['projections'][-1]['expected_cases_full'])>number(noreset['projections'][-1]['expected_cases_full'])
    assert number(central['projections'][-1]['expected_cases_full'])>number(protected['projections'][-1]['expected_cases_full'])
    partial=reinfection_clock(policy=PARTIAL,remaining_clock_fraction='.5',**META)
    r=single_reinfection_sensitivity(history,partial,tau='10',age='30',sex='female',multiplier='3.6',horizon='1')
    assert number(r['effective_clock_after_event'])==5 and r['not_a_population_projection']
    with pytest.raises(DataValidationError,match='single-event'):operator(history,rate=.1,clock=partial)
    assert not central['operator']['reinfection_clock']['RR_multiplied_by_reinfection']


def test_numerical_refinement_and_horizon_prefix_consistency():
    history=audited_adult_candidate(tail='fixed_remote_hazard_v1',tail_parameters={'hazard_per_year':'.00038'})
    pop=cohort([row(tau='10',state='remote_only')])
    coarse=project_active_tb(pop,operator(history,rate=.05,horizons=['1','5'],step='.25'))
    fine=project_active_tb(pop,operator(history,rate=.05,horizons=['1','5'],step='.125'))
    assert refinement_comparison(coarse,fine,absolute_risk_tolerance='.0001')['passed']
    single=project_active_tb(pop,operator(history,rate=.05,horizons=['5'],step='.25'))
    assert coarse['projections'][-1]['expected_cases_full']==single['projections'][0]['expected_cases_full']


def test_sex_subgroup_age_and_highest_factor_decompositions_and_ties():
    pop=cohort([row('a',n='800',sex='female'),row('b',n='200',sex='male',subgroup_id='other',rr='3.6',
        responsibility_weights={'tie_a':'1','tie_b':'1'})])
    result=project_active_tb(pop,operator())
    for p in result['projections']:
        for dimension in ('sex','subgroup','baseline_tbi_state','effective_rr','attained_age','time_since_infection_band','infection_origin'):
            assert fsum(number(r['expected_cases']) for r in p['expected_case_decomposition'][dimension])==pytest.approx(number(p['expected_cases_full']))
        tied={r['category']:r for r in p['expected_case_decomposition']['highest_rr_responsible_factor']}
        assert tied['tie_a']['expected_cases']==tied['tie_b']['expected_cases']


def test_extreme_risk_concentration_without_cap():
    result=project_active_tb(cohort([row('rare_high',n='5',tau='0',rr='74'),row('others',n='995',tau=None,state='never_infected')]),operator())
    p=result['projections'][0];top=p['extreme_risk_concentration']
    assert top['warning']=='small_group_dominates_expected_cases'
    assert number(top['rows'][0]['fraction_of_expected_cases'])==pytest.approx(1)
    assert result['operator']['effect_application']['cap'] is None


def test_future_pressure_explicit_bounds_coverage_and_no_WHO_extrapolation():
    for scenario,change in [('bounded_decline','-.02'),('bounded_increase','.02')]:
        p=future_pressure(geography_id=GEO,scenario=scenario,baseline_rate='.1',annual_change=change,rate_bounds=['.01','.2'],horizon_years=20,**META)
        assert all(.01<=number(r['rate_per_year'])<=.2 for r in p['rows']) and not p['automatic_WHO_extrapolation']
        assert pressure_rate(p,19.5,'female','g')==pytest.approx(.01 if change[0]=='-' else .2)
    with pytest.raises(DataValidationError):future_pressure(geography_id=GEO,scenario='constant_at_2024',horizon_years=20,**META)
    with pytest.raises(DataValidationError):future_pressure(geography_id=GEO,scenario='user_supplied',horizon_years=20,rows=[],**META)


def test_exact_event_time_MC_constant_and_future_reinfection():
    pop=cohort([row('remote',n='600',tau='10',state='remote_only'),row('uninfected',n='400',tau=None,state='never_infected')])
    op=operator(rate=.05,mu=.03,horizons=['5'],step='.125');prediction=project_active_tb(pop,op)
    mc=monte_carlo_validate(pop,op,prediction,rng=random.Random(731),sample_size=2500,horizon='5',standard_errors=5,numerical_absolute_tolerance=.001)
    assert mc['passed'] and not mc['global_RNG_used']
    assert any(len(r['infection_events'])>=2 for r in mc['rows'])
    assert all(r['first_TB_episode_count'] in (0,1) for r in mc['rows'])
    assert any(r['outcome']=='competing_death' for r in mc['rows'])
    with pytest.raises(DataValidationError):monte_carlo_validate(pop,op,prediction,rng=None,sample_size=100,horizon=5)


def test_deterministic_hash_row_order_and_offline_isolation(monkeypatch):
    def no_network(*args,**kwargs):raise AssertionError('No network in forward operator')
    monkeypatch.setattr(socket,'socket',no_network);monkeypatch.setattr(socket,'create_connection',no_network)
    records=[row('a',n='800'),row('b',n='200',sex='male')]
    a=project_active_tb(cohort(records),operator());b=project_active_tb(cohort(list(reversed(records))),operator())
    assert a['profile_hash']==b['profile_hash']==profile_hash(a)
    assert json.loads(canonical(a))==a and a['inverse_incidence_calibration']=='unavailable_pending_progression_and_risk_forward_operator'
    for folder in ('app','ui','general_ui','explicit_tbi_ui','engine/apy'):
        for p in (REPO_ROOT/folder).rglob('*.py'):
            for node in ast.walk(ast.parse(p.read_text(encoding='utf-8-sig'))):
                assert not(isinstance(node,ast.ImportFrom) and (node.module or '').startswith('engine.country_tb'))


def test_incidence_comparison_does_not_unlock_inverse_or_fit():
    from scripts.preview_country_tb_sources import incidence_selection
    from engine.country_tb.preview_risk_allocation import toy_demography
    from engine.country_tb.baseline_tbi_source import resolve_baseline_tbi
    demo=toy_demography();selection=incidence_selection('g',demo)
    result=project_active_tb(cohort(),operator())
    comparison=incidence_comparison(result,selection['configuration']['targets'],observation_id='synthetic:active-disease')
    assert comparison['status']=='validation_only_not_calibration' and comparison['fitted_parameters'] is None
    assert number(comparison['predicted_incidence_per_100000_person_year'])==pytest.approx(2000)
    assert number(comparison['difference'])==pytest.approx(1990)
    blocked=resolve_baseline_tbi(selection,demo)
    assert blocked['status']=='unavailable_pending_progression_and_risk_forward_operator' and blocked['tbi_estimate'] is None


def test_adapter_direct_times_required_and_marginals_preserved():
    from scripts.preview_country_tb_sources import direct_selection
    from engine.country_tb.preview_risk_allocation import toy_demography,allocate_example
    from engine.country_tb.baseline_tbi_source import resolve_baseline_tbi
    demo=toy_demography();record={'record_id':'states','age_lower':0,'age_upper_exclusive':None,'sex':'both',
        'age_allocation_policy':'piecewise_constant_within_band_v1','any_tbi':'.3','recent':'.1','remote_only':'.2',
        'never_infected':'.7','prior_remote_plus_recent':'.04'}
    base=resolve_baseline_tbi(direct_selection('g',demo,[record]),demo);allocation=allocate_example(base,{'diabetes':'.2'})
    before=deepcopy(allocation)
    unresolved=cohort_from_allocation(allocation,active_TB_exclusion_policy='assumed_zero_baseline_active_TB_diagnostic_v1',**META)
    assert all(r['unresolved_reasons'] for r in unresolved['rows'])
    support={}
    for cell in allocation['rows']:
        if not number(cell['population']['value']):continue
        for state,tau in [('never_infected',None),('recent','.5'),('remote_only','10')]:
            support['g:'+cell['sex']+':'+str(cell['age_lower'])+':'+state]={**META,'nodes':[
                {'age':str(cell['age_lower']+.25),'tau':tau,'weight':'1'}]}
    compiled=cohort_from_allocation(allocation,conditional_support=support,
        active_TB_exclusion_policy='assumed_zero_baseline_active_TB_diagnostic_v1',**META)
    assert fsum(number(r['population']) for r in compiled['rows'])==pytest.approx(1000)
    for state,p in [('never_infected',.7),('recent',.1),('remote_only',.2)]:
        assert fsum(number(r['population']) for r in compiled['rows'] if r['baseline_tbi_state']==state)==pytest.approx(1000*p)
    assert fsum(number(r['population'])*number(r['prior_remote_plus_recent_fraction']) for r in compiled['rows'])==pytest.approx(40)
    assert allocation==before and compiled['source_allocation']['profile_hash']==allocation['profile_hash']


def test_conditional_last_event_support_not_mean_time_or_display_midpoint():
    from scripts.preview_country_tb_sources import overall_selection
    from engine.country_tb.preview_risk_allocation import toy_demography,allocate_example
    from engine.country_tb.baseline_tbi_source import resolve_baseline_tbi
    demo=toy_demography();base=resolve_baseline_tbi(overall_selection('g',demo,'.2'),demo)
    allocation=allocate_example(base,{'heavy_smoking':'0'})
    compiled=cohort_from_allocation(allocation,active_TB_exclusion_policy='assumed_zero_baseline_active_TB_diagnostic_v1',**META)
    assert len(compiled['rows'])>100 and all(r['tau'] is not None for r in compiled['rows'] if r['baseline_tbi_state']!='never_infected')
    assert number(compiled['population'])==1000
    assert all(number(r['age'])%1!=.5 for r in compiled['rows'])


def test_open_tail_and_missing_curve_stay_visible():
    unresolved=row('tail',n='20',age=None,tau=None,state='unresolved',rr=None,
        unresolved_reasons=['unresolved_100plus_age_distribution'],source_age_lower=100,source_age_upper_exclusive=None)
    result=project_active_tb(cohort([row('adult',n='980'),unresolved]),operator())
    assert number(result['population'])==1000 and number(result['projections'][0]['unresolved_population'])==20
    assert result['projections'][0]['expected_cases_full'] is None


def test_protected_state_and_all_retained_2a_3a_files():
    p=REPO_ROOT/'abm/output/country_tb_m3b_audit_20261029/starting_state.json'
    if not p.is_file():pytest.skip('Live protected-state audit unavailable outside milestone workspace.')
    start=json.loads(p.read_text())
    for path,h in {**start['frozen_files'],**start['retained_2a_3a_files']}.items():assert digest((REPO_ROOT/path).read_bytes())==h,path
    for name,old in start['protected'].items():
        folder=REPO_ROOT.parent/name
        def git(*args):return subprocess.check_output(['git','-C',str(folder),*args],stderr=subprocess.DEVNULL)
        assert git('rev-parse','HEAD').decode().strip()==old['head'] and git('status','--short').decode().splitlines()==old['status']
        assert digest(git('diff','--binary','HEAD'))==old['tracked_diff_sha256'] and digest(git('diff','--cached','--binary'))==old['index_diff_sha256']
        for file,h in old['changed_files'].items():assert digest((folder/file).read_bytes())==h
    refs=dict(x.split(' ',1) for x in subprocess.check_output(['git','for-each-ref','--format=%(refname) %(objectname)']).decode().splitlines())
    allowed={'refs/heads/feature/country-subnational-tb-model-v1','refs/remotes/origin/feature/country-subnational-tb-model-v1'}
    assert {k:v for k,v in refs.items() if k not in allowed}=={k:v for k,v in start['refs'].items() if k not in allowed}


def test_observed_cumulative_incidence_not_silently_a_cause_specific_hazard():
    knots=[{'time_years':0,'cumulative_risk':0},{'time_years':1,'cumulative_risk':.05}]
    with pytest.raises(DataValidationError,match='interpretation'):curve_specification(knots,tail='no_extrapolation_v1',**META)
    with pytest.raises(DataValidationError):curve_specification(knots,tail='no_extrapolation_v1',
        risk_interpretation='diagnostic_cumulative_risk_surrogate_pending_cause_specific_review_v1',**META)
    accepted=curve_specification(knots,tail='no_extrapolation_v1',risk_interpretation='net_risk_without_competing_mortality_v1',**META)
    assert accepted['status']=='reviewed_external'


def test_no_reset_still_updates_actual_recency_and_most_recent_time():
    history=audited_adult_candidate(tail='fixed_remote_hazard_v1',tail_parameters={'hazard_per_year':'.00038'})
    result=project_active_tb(cohort([row(tau='10',state='remote_only')]),
        operator(history,rate=.1,horizons=['5'],clock=reinfection_clock(policy=NO_RESET,**META)))
    d=result['projections'][0]['expected_case_decomposition']
    state={x['category']:number(x['expected_cases']) for x in d['event_tbi_state']}
    clock={x['category']:number(x['expected_cases']) for x in d['progression_clock_band']}
    actual={x['category']:number(x['expected_cases']) for x in d['time_since_infection_band']}
    assert state['recent']>0 and clock['10-25']==pytest.approx(number(result['projections'][0]['expected_cases_full']))
    assert any(actual.get(k,0)>0 for k in ['0-1','1-2','2-5'])


def test_partial_coverage_has_common_horizon_monotonic_comparator():
    pop=cohort([row('recent'),row('remote',tau='10',state='remote_only')])
    result=project_active_tb(pop,operator(audited_adult_candidate()))
    counts=[number(p['common_horizon_expected_cases']) for p in result['projections']]
    assert counts==sorted(counts) and len({p['common_horizon_covered_population'] for p in result['projections']})==1
    assert number(result['projections'][-1]['unresolved_population'])==1000


def test_declining_tail_and_strong_mortality_remain_bounded():
    nh=audited_adult_candidate(tail='declining_to_remote_plateau_v1',tail_parameters={'plateau_hazard_per_year':'.0001','decay_rate_per_year':'.1'})
    result=project_active_tb(cohort([row(age='70',tau='30',state='remote_only',rr='74')]),operator(nh,mu=2.))
    for p in result['projections']:
        assert 0<=number(p['expected_cases_full'])<=1000
        assert number(p['expected_cases_full'])+number(p['competing_deaths_covered'])+number(p['surviving_disease_free_covered'])==pytest.approx(1000,abs=1e-9)


def test_hashes_independent_of_caller_decimal_context_and_explicit_RNG():
    from decimal import localcontext
    with localcontext() as c:
        c.prec=12;a=project_active_tb(cohort(),operator())
    with localcontext() as c:
        c.prec=60;b=project_active_tb(cohort(),operator())
    assert a['profile_hash']==b['profile_hash']
    pop=cohort();op=operator(horizons=['1']);prediction=project_active_tb(pop,op)
    left=monte_carlo_validate(pop,op,prediction,rng=random.Random(443),sample_size=100,horizon=1)
    right=monte_carlo_validate(pop,op,prediction,rng=random.Random(443),sample_size=100,horizon=1)
    assert left['profile_hash']==right['profile_hash'] and left['rows']==right['rows']


def test_finite_reviewed_RR_age_scope_not_silently_extended():
    scoped=row(age='60',rr='3.6',effect_status='reviewed',effect_age_applicability={
        'age_lower':'15','age_upper_exclusive':'65','policy':'synthetic_reviewed_age_range'})
    result=project_active_tb(cohort([scoped]),operator())
    assert result['projections'][0]['expected_cases_full'] is not None
    assert result['projections'][1]['expected_cases_full'] is not None
    assert result['projections'][2]['expected_cases_full'] is None
    assert 'effect_age_applicability_coverage_external_required' in result['projections'][2]['rows'][0]['unresolved_reasons']


def test_pediatric_factor_scope_needs_review_and_stops_at_its_boundary():
    r=row(age='3',rr='2',pediatric_effect_approved=True,effect_status='reviewed_pediatric')
    with pytest.raises(DataValidationError,match='applicability'):cohort([r])
    r['effect_age_applicability']={'age_lower':'0','age_upper_exclusive':'5','policy':'synthetic_reviewed_pediatric_effect'}
    result=project_active_tb(cohort([r]),operator(constant_history(children=True)))
    assert result['projections'][0]['expected_cases_full'] is not None
    assert result['projections'][1]['expected_cases_full'] is None


@pytest.fixture(scope='module')
def diagnostic_bundle():
    from engine.country_tb.preview_progression import diagnostics
    return diagnostics(mc_sample_size=250)


def test_complete_diagnostic_suite_is_synthetic_and_WHO_comparison_only(diagnostic_bundle):
    d=diagnostic_bundle
    assert d['monte_carlo']['passed'] and d['numerical_refinement']['passed']
    assert d['incidence_comparison']['target_incidence_per_100000_person_year']=='182'
    assert d['incidence_comparison']['status']=='validation_only_not_calibration'
    assert d['incidence_comparison']['fitted_parameters'] is None
    au=d['projections']['Australian_synthetic_mixture']
    assert number(au['population'])==10000 and not au['complete_country_estimate']
    child=d['projections']['children_external_required']
    assert number(child['projections'][0]['unresolved_population'])==500
    for name,p in d['projections'].items():
        assert p['diagnostic_only'] and not p['complete_country_estimate']
        for horizon in p['projections']:
            for dimension,rows in horizon['expected_case_decomposition'].items():
                if dimension!='highest_rr_responsible_factor':
                    assert fsum(number(r['expected_cases']) for r in rows)==pytest.approx(number(horizon['expected_cases_covered']))


def test_twelve_versioned_contract_families_offline_schema(diagnostic_bundle,monkeypatch):
    from jsonschema import Draft202012Validator
    schema=json.loads((REPO_ROOT/'schemas/country_tb/v1/progression_forward.schema.json').read_text())
    Draft202012Validator.check_schema(schema);validator=Draft202012Validator(schema)
    def forbidden(*args,**kwargs):raise AssertionError('Schema cannot access network')
    monkeypatch.setattr(socket,'socket',forbidden);monkeypatch.setattr(socket,'create_connection',forbidden)
    samples={}
    def walk(x):
        if isinstance(x,dict):
            if x.get('implementation_version')==IMPLEMENTATION:samples[x['profile_hash']]=x
            for v in x.values():walk(v)
        elif isinstance(x,list):
            for v in x:walk(v)
    walk(diagnostic_bundle)
    walk(artifact('country_tb.active_tb_forward_export_receipt.v1',profile_hash_reference=diagnostic_bundle['profile_hash'],files={'diagnostics.json':'a'*64}))
    assert len({x['schema_id'] for x in samples.values()})==12
    for x in samples.values():validator.validate(x);assert x['profile_hash']==profile_hash(x)
    bad=deepcopy(diagnostic_bundle['natural_history_audit']);bad['group_curves']['15+']['knots'][1]['cumulative_risk']='3.8'
    assert not validator.is_valid(bad)
    bad=deepcopy(diagnostic_bundle['projections']['reset']['operator']);bad['mode']='inverse_fit'
    assert not validator.is_valid(bad)
    bad=deepcopy(diagnostic_bundle['incidence_comparison']);bad['status']='calibrated'
    assert not validator.is_valid(bad)
