"""Milestone 2D: invented source/subgroup fixtures, no empirical country inputs."""
from copy import deepcopy
from decimal import Decimal,localcontext
import ast
import json
import socket
import subprocess

import pytest

from engine.country_tb.composition_common import population_demography,canonical,profile_summary
from engine.country_tb.data_common import DataValidationError,seal,profile_hash,REPO_ROOT,digest
from engine.country_tb.baseline_tbi_source import (
    baseline_tbi_source,resolve_baseline_tbi,interface_specification,INCIDENCE,OBSERVED,OVERALL,
    AGE_ANCHORS,DIRECT_STATES,
)
from engine.country_tb.active_tb_targets import active_tb_targets,BLOCKED,TYPES
from engine.country_tb.direct_tbi_profile import direct_tbi_profile,CONDITIONAL_RECENCY,UNAVAILABLE_RECENCY
from engine.country_tb.population_subgroups import population_subgroup,population_mixture,validate_resolution
from engine.country_tb.catalytic_calibration import calibration_policy
from engine.country_tb.catalytic_core import evaluate_age
from engine.country_tb.tbi_observations import observations_profile,assay_operator
from scripts.preview_country_tb_sources import (
    REVIEW,SOURCE,PROVENANCE,GEO,overall_selection,direct_selection,incidence_selection,
    synthetic_history,synthetic_demography,comparison_model,make_subgroup,make_mixture,membership,
    synthetic_scenarios,direct_operator,
)


def small_demo(population=1000,tail=0):
    with localcontext() as c:
        c.prec=80
        rows=[{'age_lower':a,'age_upper_exclusive':a+1 if a<100 else None,'sex':s,
            'population':str((Decimal(population)-Decimal(tail))/4 if a<100 else Decimal(tail)/2)}
            for s in ('female','male') for a in (20,60,100)]
    return population_demography(rows,geography_id=GEO,geography_level='national',population=population,
        source=SOURCE,provenance=PROVENANCE,review=REVIEW)


def total_record(value='.2',**changes):
    row={'record_id':'all','age_lower':0,'age_upper_exclusive':None,'sex':'both',
         'age_allocation_policy':'piecewise_constant_within_band_v1','any_tbi':value}
    row.update(changes);return row


def selection_kwargs(selection):
    return {k:deepcopy(selection[k]) for k in ('mode','geography_id','geography_level','subgroup_id',
        'source','provenance','review','profile_kind','configuration','baseline_year')}


@pytest.fixture(scope='module')
def demo():return small_demo()


@pytest.fixture(scope='module')
def model(demo):return comparison_model(demo)


@pytest.fixture(scope='module')
def scenarios():return synthetic_scenarios()


@pytest.mark.parametrize('mode',[None,'unknown',[],[INCIDENCE,OBSERVED],{'mode':OBSERVED},''])
def test_exactly_one_known_primary_mode_required(demo,mode):
    kwargs=selection_kwargs(overall_selection('g',demo,'.2'));kwargs['mode']=mode
    with pytest.raises(DataValidationError):baseline_tbi_source(**kwargs)


def test_no_fallback_or_operator_bypass(demo):
    kwargs=selection_kwargs(incidence_selection('g',demo))
    kwargs['configuration']['anchor']={'prevalence':'.2'}
    with pytest.raises(DataValidationError,match='no prevalence fallback'):baseline_tbi_source(**kwargs)
    kwargs=selection_kwargs(incidence_selection('g',demo));kwargs['configuration']['forward_operator']={'approved':True}
    with pytest.raises(DataValidationError,match='No forward operator'):baseline_tbi_source(**kwargs)
    kwargs=selection_kwargs(overall_selection('g',demo,'.2'));kwargs['configuration']['direct_profile']={}
    with pytest.raises(DataValidationError,match='mixed'):baseline_tbi_source(**kwargs)


@pytest.mark.parametrize('kind',TYPES)
def test_incident_notification_prevalent_imported_observations_stay_distinct_and_blocked(demo,kind):
    kwargs=selection_kwargs(incidence_selection('g',demo))
    raw=deepcopy(kwargs['configuration']['targets']['rows'][0]['receipt']['source_record'])
    raw['active_tb_target_type']=kind
    if 'prevalent' in kind:raw['unit']='cases_per_100000_persons'
    targets=active_tb_targets([raw],source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    kwargs['configuration']['targets']=targets
    result=resolve_baseline_tbi(baseline_tbi_source(**kwargs),demo)
    assert result['status']==BLOCKED and result['tbi_estimate'] is None
    assert not result['usable_for_downstream_simulation'] and result['calibration'] is None
    assert result['summary']['unresolved_population']=='1000'
    assert all(r['any_tbi'] is None for r in result['rows'])
    assert targets['rows'][0]['active_tb_target_type']==kind
    assert targets['rows'][0]['eligible_for_future_pressure_calibration']==(kind=='incident_disease')
    assert len(result['missing_approved_components'])==8


@pytest.mark.parametrize('change',[{'denominator':None},{'denominator':'0'},{'observation_period':None},
    {'observation_period':{'start_time':'2024','end_time':'2024'}},{'rate':'10','count':'1'},
    {'rate':None,'count':None},{'unit':'prevalence_percent'},{'active_tb_target_type':'unclear_typo'}])
def test_bad_active_targets_fail_instead_of_becoming_rates(demo,change):
    raw=deepcopy(incidence_selection('g',demo)['configuration']['targets']['rows'][0]['receipt']['source_record'])
    raw.update(change)
    with pytest.raises(DataValidationError):active_tb_targets([raw],source=SOURCE,provenance=PROVENANCE,review=REVIEW)


def test_count_and_who_ascertainment_provenance(demo):
    raw=deepcopy(incidence_selection('g',demo)['configuration']['targets']['rows'][0]['receipt']['source_record'])
    raw.update(rate=None,count=12,unit='cases',source_class='WHO_estimate')
    result=active_tb_targets([raw],source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    assert result['rows'][0]['receipt']['source_record']==raw
    assert result['rows'][0]['count']=='12' and result['rows'][0]['rate'] is None
    raw['ascertainment']['status']='not_reviewed'
    assert not active_tb_targets([raw],source=SOURCE,provenance=PROVENANCE,review=REVIEW)['rows'][0]['eligible_for_future_pressure_calibration']


def test_overall_anchor_covered_population_recovery_and_tail():
    demo=small_demo(tail=20);selection=overall_selection('single',demo,'.25')
    original=deepcopy(selection);result=resolve_baseline_tbi(selection,demo)
    assert selection==original and result['calibration']['fit_status']=='accepted'
    assert float(result['summary']['components']['any_tbi']['prevalence_covered'])==pytest.approx(.25,abs=1e-12)
    assert result['summary']['estimated_population']=='980' and result['summary']['unresolved_population']=='20'
    assert result['summary']['components']['any_tbi']['prevalence_full'] is None
    assert not result['usable_for_downstream_simulation']
    assert all(r['any_tbi'] is None for r in result['rows'] if r['age_upper_exclusive'] is None)
    assert result['calibration']['observations']['rows'][0]['numerator'] is None
    assert result['recent_remote_composition']=='mixed_or_partial'


def test_overall_known_scale_recovery(demo):
    h=synthetic_history();true='.013'
    with localcontext() as ctx:
        ctx.prec=80
        target=sum(Decimal(r['population'])*Decimal(evaluate_age(h,alpha=true,age=r['age_lower'],age_upper_exclusive=r['age_upper_exclusive'],
            age_policy='uniform_within_single_year',sex=r['sex'],time='2024.5',include_event_times=False)['any_tbi'])
            for r in demo['rows'] if r['age_upper_exclusive'] is not None)/1000
    result=resolve_baseline_tbi(overall_selection('known',demo,str(target)),demo)
    assert abs(Decimal(result['calibration']['alpha'])-Decimal(true))<Decimal('1e-17')
    assert all(r['most_recent_infection_time'] is not None for r in result['rows'] if Decimal(r['population'])>0)


@pytest.mark.parametrize('p',['-0.1','1.01','20%','NaN','1'])
def test_invalid_or_infeasible_prevalences_are_not_clipped(demo,p):
    with pytest.raises(DataValidationError):resolve_baseline_tbi(overall_selection('g',demo,p),demo)


def test_numerator_denominator_consistency_and_no_invented_counts(demo):
    kwargs=selection_kwargs(overall_selection('g',demo,'.2'));a=kwargs['configuration']['anchor']
    a.update(numerator=20,denominator=100)
    result=resolve_baseline_tbi(baseline_tbi_source(**kwargs),demo)
    assert result['calibration']['observations']['rows'][0]['numerator']==20
    a['numerator']=None
    with pytest.raises(DataValidationError):baseline_tbi_source(**kwargs)
    a.update(numerator=30,denominator=100)
    with pytest.raises(DataValidationError,match='Counts/prevalence'):baseline_tbi_source(**kwargs)


@pytest.mark.parametrize('sexes',[('both',),('female',),('male',),('female','male')])
def test_age_specific_one_scale_anchors_and_assays(demo,sexes):
    h=synthetic_history();op=assay_operator(assay_class='IGRA',name='Invented assay',positivity_threshold='synthetic',
        sensitivity='.85',specificity='.95',review=REVIEW,provenance=PROVENANCE)
    records=[]
    for age in (20,60):
        for sex in sexes:
            p=Decimal(evaluate_age(h,alpha='.011',age=age,sex=sex,include_event_times=False)['any_tbi'])
            records.append({'observation_id':f'{sex}:{age}','subgroup_id':'age_anchors','geography_id':GEO,'geography_level':'national',
                'observation_policy':'reviewed_age_specific_anchor','age_lower':str(age),'age_upper_exclusive':str(age+1),'sex':sex,
                'survey_year':2024,'survey_time':'2024','prevalence':str(Decimal('.85')*p+Decimal('.05')*(1-p)),
                'numerator':None,'denominator':None,'review':REVIEW,'source':SOURCE,'sampling_frame':'invented sample',
                'survey_design':{'kind':'synthetic_only'},'uncertainty':{'kind':'exact'},'assay_operator':op,
                'weight_source':PROVENANCE,'strata':[{'age_lower':str(age),'age_upper_exclusive':None,
                    'age_policy':'exact_age','sex':sex,'weight':'1'}]})
    observations=observations_profile(records,review=REVIEW)
    policy=calibration_policy(alpha_upper='2',observation_weights={r['observation_id']:'1' for r in records},
        residual_tolerance='.00001',review=REVIEW,method='weighted_least_squares')
    selection=baseline_tbi_source(mode=OBSERVED,geography_id=GEO,geography_level='national',subgroup_id='age_anchors',
        source=SOURCE,provenance=PROVENANCE,review=REVIEW,profile_kind='diagnostic',configuration={'submode':AGE_ANCHORS,
            'baseline_time':'2024.5','recent_window_years':'5','pressure_history':h,'calibration_policy':policy,'observations':observations})
    result=resolve_baseline_tbi(selection,demo)
    assert float(result['calibration']['alpha'])==pytest.approx(.011,abs=1e-14)
    assert result['calibration']['policy']['free_parameters']==['alpha']
    assert len(result['calibration']['rows'])==len(records)


def test_direct_profile_original_values_and_model_comparison_preserved(demo,model):
    records=[total_record('0.2000')];before=deepcopy(records);m=deepcopy(model)
    selection=direct_selection('direct',demo,records,comparison=model)
    result=resolve_baseline_tbi(selection,demo)
    assert records==before and model==m and result['model_derived_profile']==model
    assert result['selection']['configuration']['direct_profile']['receipt']['source_records']==before
    assert result['summary']['components']['any_tbi']['prevalence_full']=='0.2'
    assert result['recent_remote_composition']=='unavailable' and not result['usable_for_downstream_simulation']
    assert all(r['recent'] is None and r['remote_only'] is None for r in result['rows'])


def test_direct_recent_remote_retained_and_subset_not_added(demo):
    record=total_record('.3',recent='.1',remote_only='.2',never_infected='.7',prior_remote_plus_recent='.04')
    kwargs=selection_kwargs(direct_selection('states',demo,[record],profile_kind='user_defined'))
    kwargs['configuration']['submode']=DIRECT_STATES
    result=resolve_baseline_tbi(baseline_tbi_source(**kwargs),demo)
    assert result['usable_for_downstream_simulation']
    assert result['summary']['components']['any_tbi']['complete_count']=='300'
    assert result['summary']['components']['prior_remote_plus_recent']['complete_count']=='40'
    assert all(r['recent']=='0.1' and r['remote_only']=='0.2' for r in result['rows'])
    validate_resolution(result)


@pytest.mark.parametrize('changes',[{'recent':'.1'},{'recent':'.1','remote_only':'.3','never_infected':'.7'},
    {'recent':'.1','remote_only':'.1','never_infected':'.8','prior_remote_plus_recent':'.2'},
    {'any_tbi':'1.1'},{'age_allocation_policy':'interpolate'},{'age_allocation_policy':'exact_single_year_values_v1'}])
def test_invalid_direct_states_and_ages_fail(changes):
    with pytest.raises(DataValidationError):direct_tbi_profile([total_record(**changes)],source=SOURCE,provenance=PROVENANCE,review=REVIEW)


def test_overlap_and_gaps_are_explicit(demo):
    rows=[total_record(record_id='both'),total_record(record_id='female',sex='female')]
    with pytest.raises(DataValidationError,match='Overlapping'):direct_tbi_profile(rows,source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    records=[total_record(age_lower=20,age_upper_exclusive=21,age_allocation_policy='exact_single_year_values_v1')]
    result=resolve_baseline_tbi(direct_selection('gaps',demo,records),demo)
    assert result['status']=='partial_direct_profile_coverage'
    assert result['summary']['unresolved_population']=='500'
    assert result['summary']['components']['any_tbi']['prevalence_full'] is None
    assert result['direct_profile_gaps']


def test_open_tail_direct_value_and_reviewed_tail_weights():
    demo=small_demo(tail=20)
    direct=resolve_baseline_tbi(direct_selection('tail',demo,[total_record('.3')]),demo)
    assert direct['summary']['unresolved_population']=='0'
    assert direct['summary']['components']['any_tbi']['complete_count']=='300'
    kwargs=selection_kwargs(overall_selection('tail_model',demo,'.2'))
    kwargs['configuration']['tail_distributions']={s+':100':{'review':REVIEW,'provenance':PROVENANCE,
        'ages':[{'age':'101','weight':'1'}]} for s in ('female','male')}
    modeled=resolve_baseline_tbi(baseline_tbi_source(**kwargs),demo)
    assert modeled['summary']['unresolved_population']=='0'
    assert float(modeled['summary']['components']['any_tbi']['prevalence_full'])==pytest.approx(.2,abs=1e-12)


def test_explicit_conditional_recency_matches_formula(demo,model):
    default=resolve_baseline_tbi(direct_selection('default',demo,[total_record('.3')],comparison=model),demo)
    hybrid=resolve_baseline_tbi(direct_selection('hybrid',demo,[total_record('.3')],comparison=model,recency_policy=CONDITIONAL_RECENCY),demo)
    assert default['profile_hash']!=hybrid['profile_hash']
    for out,m in zip(sorted(hybrid['rows'],key=lambda r:(r['sex'],r['age_lower'])),sorted(model['rows'],key=lambda r:(r['sex'],r['age_lower']))):
        assert float(out['recent'])==pytest.approx(.3*float(m['recent'])/float(m['any_tbi']),rel=1e-14)
        assert float(out['recent'])+float(out['remote_only'])==pytest.approx(.3)
        assert Decimal(out['prior_remote_plus_recent'])<=Decimal(out['recent'])
    assert hybrid['recent_remote_composition']=='model_derived'
    assert all(r['most_recent_infection_time'] is None for r in hybrid['rows'])


def test_conditional_recency_needs_model_and_never_overwrites_direct_states(demo,model):
    with pytest.raises(DataValidationError,match='requires a comparison'):
        resolve_baseline_tbi(direct_selection('bad',demo,[total_record()],recency_policy=CONDITIONAL_RECENCY),demo)
    record=total_record('.3',recent='.02',remote_only='.28',never_infected='.7')
    result=resolve_baseline_tbi(direct_selection('states',demo,[record],comparison=model,recency_policy=CONDITIONAL_RECENCY),demo)
    assert all(r['recent']=='0.02' for r in result['rows'])


def test_supplied_within_band_values_and_no_interpolation(demo):
    children=[{'age_lower':a,'age_upper_exclusive':a+1 if a<100 else None,'sex':'both','any_tbi':p}
              for a,p in ((20,'.1'),(60,'.3'),(100,'.4'))]
    record=total_record('.2',age_allocation_policy='supplied_within_band_distribution_v1',within_band_values=children)
    result=resolve_baseline_tbi(direct_selection('distributed',demo,[record]),demo)
    assert result['summary']['components']['any_tbi']['prevalence_full']=='0.2'
    assert {r['any_tbi'] for r in result['rows']}=={'0.1','0.3','0.4'}
    record['any_tbi']='.9'
    with pytest.raises(DataValidationError,match='parent-band'):resolve_baseline_tbi(direct_selection('bad',demo,[record]),demo)


def test_source_row_order_hashes_and_canonical_serialization(demo):
    records=[total_record('.1',record_id='young',age_upper_exclusive=40),total_record('.3',record_id='older',age_lower=40)]
    a=direct_selection('g',demo,records)
    b=direct_selection('g',demo,list(reversed(records)))
    assert a['profile_hash']==b['profile_hash']
    assert resolve_baseline_tbi(a,demo)['profile_hash']==resolve_baseline_tbi(b,demo)['profile_hash']
    assert json.loads(canonical(a))==a and a['profile_hash']==profile_hash(a)
    before=deepcopy(demo);shuffled=deepcopy(demo['rows']);shuffled.reverse()
    d=population_demography(shuffled,geography_id=GEO,geography_level='national',population=demo['population'],source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    assert d['profile_hash']==demo['profile_hash'] and demo==before


def test_subgroup_metadata_and_weighted_mixture(demo):
    groups=[]
    for name,pop,p in [('a',1000,'.1'),('b',3000,'.3')]:
        d=small_demo(pop)
        states=total_record(p,recent=str(Decimal(p)/2),remote_only=str(Decimal(p)/2),never_infected=str(1-Decimal(p)))
        selection=direct_selection(name,d,[states],profile_kind='user_defined')
        groups.append(make_subgroup(name,d,selection,country_or_region_of_birth={'geography_id':'who:VNM','geography_level':'national'}))
    result=make_mixture(groups,subgroup_weights={'a':'.25','b':'.75'})
    assert result['summary']['population']=='4000'
    assert result['summary']['components']['any_tbi']['prevalence_full']=='0.25'
    assert result['summary']['components']['recent']['prevalence_full']=='0.125'
    assert result['usable_for_downstream_simulation']
    assert result['summary']['components']['prior_remote_plus_recent']['prevalence_full'] is None
    assert result['subgroups'][0]['residence_metadata']['exposure_history_geography'] is None
    assert result['subgroups'][0]['current_residence_geography']!=result['subgroups'][0]['residence_metadata']['country_or_region_of_birth']['geography_id']
    assert result['profile_hash']==make_mixture(list(reversed(groups)),subgroup_weights={'b':'.75','a':'.25'})['profile_hash']
    for field in ('population','estimated_population','unresolved_population'):assert result['summary'][field]==profile_summary([r for g in groups for r in g['tbi_profile']['rows']])[field]


def test_incomplete_subgroup_blocks_full_mixture_but_keeps_denominator(demo):
    a=make_subgroup('a',demo,direct_selection('a',demo,[total_record('.2')]))
    b=make_subgroup('b',demo,incidence_selection('b',demo))
    result=make_mixture([a,b])
    assert result['total_population']=='2000' and result['summary']['unresolved_population']=='1000'
    assert result['summary']['components']['any_tbi']['estimated_count']=='200'
    assert result['summary']['components']['any_tbi']['prevalence_covered']=='0.2'
    assert result['summary']['components']['any_tbi']['prevalence_full'] is None
    assert result['complete_mixture_result'] is None and not result['usable_for_downstream_simulation']
    assert any(r['subgroup_id']=='b' and r['population']=='1000' for r in result['unresolved_subgroups'])


@pytest.mark.parametrize('case',['ID','category','person','denominator','weights','age_grid','reference'])
def test_invalid_mixtures_fail(demo,case):
    a=make_subgroup('a',demo,direct_selection('a',demo,[total_record()]))
    b=make_subgroup('b',demo,direct_selection('b',demo,[total_record()]))
    kwargs={}
    if case=='ID':b=a
    elif case=='category':b['membership']['category_id']='a';b=seal(b)
    elif case=='person':
        for g in (a,b):g['membership']['person_ids']=[str(i) for i in range(1000)]
        a,b=seal(a),seal(b)
    elif case=='denominator':kwargs['total_population']=3000
    elif case=='weights':kwargs['subgroup_weights']={'a':'.1','b':'.9'}
    elif case=='age_grid':
        d=small_demo();r=deepcopy(d['rows'])
        for cell in r:
            if cell['age_lower']==20:cell['age_lower']=21;cell['age_upper_exclusive']=22
        d=population_demography(r,geography_id=GEO,geography_level='national',population=1000,source=SOURCE,provenance=PROVENANCE,review=REVIEW)
        b=make_subgroup('b',d,direct_selection('b',d,[total_record()]))
    else:kwargs['reference_age_sex_distribution']=[{**r,'population':'1'} for r in demo['rows']]
    with pytest.raises(DataValidationError):make_mixture([a,b],**kwargs)


def test_subgroup_population_and_risk_placeholder_validation(demo):
    with pytest.raises(DataValidationError,match='denominator'):
        population_subgroup(subgroup_id='g',label='g',description='g',current_residence_geography=GEO,geography_level='national',
            population=999,demography=demo,tbi_source=direct_selection('g',demo,[total_record()]),membership=membership('g'),source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    with pytest.raises(DataValidationError,match='deferred'):
        make_subgroup('g',demo,direct_selection('g',demo,[total_record()]),risk_factor_profile={'smoking':'.2'})


def test_exposure_proxy_requires_explicit_review_and_preserves_original(demo):
    h=synthetic_history('who:VNM');selection=overall_selection('migrant',demo,'.2',history=h)
    with pytest.raises(DataValidationError,match='Exposure-history geography'):resolve_baseline_tbi(selection,demo)
    kwargs=selection_kwargs(selection);kwargs['configuration']['history_application_policy']={
        'mode':'reviewed_subgroup_exposure_proxy_v1','source':SOURCE,'provenance':PROVENANCE,'review':REVIEW,
        'rationale':'Invented proxy transport; not a migration reconstruction.'}
    result=resolve_baseline_tbi(baseline_tbi_source(**kwargs),demo)
    assert result['calibration']['pressure_history']['shape']['provenance']['source_pressure_history']==h
    assert result['geography_id']==GEO and h['geography_id']=='who:VNM'


def test_unreviewed_source_does_not_fall_back(demo):
    kwargs=selection_kwargs(direct_selection('unreviewed',demo,[total_record()]))
    kwargs['review']['status']='not_reviewed'
    result=resolve_baseline_tbi(baseline_tbi_source(**kwargs),demo)
    assert result['status']=='unavailable_pending_source_review' and result['summary']['unresolved_population']=='1000'


def test_six_invented_scenarios(scenarios):
    assert scenarios['diagnostic_only'] and 'NOT AUSTRALIAN' in scenarios['warning']
    s=scenarios['scenarios'];mix=s['2_Australian_mixture']
    assert float(s['1_overall_anchor']['summary']['components']['any_tbi']['prevalence_covered'])==pytest.approx(.2,abs=1e-12)
    assert mix['total_population']=='10000' and float(mix['summary']['components']['any_tbi']['prevalence_full'])==pytest.approx(.08,abs=1e-12)
    assert s['5_incidence_mode_blocked']['status']==BLOCKED
    assert s['6_incomplete_mixture']['summary']['unresolved_population']=='8000'
    assert s['4_total_only_default']['summary']['components']['recent']['prevalence_full'] is None
    assert s['4_explicit_conditional_recency']['summary']['recency_complete']
    assert s['3_complete_direct_override']['model_derived_profile'] is not None
    assert s['3_complete_direct_override']['summary']['unresolved_population']=='0'
    assert not mix['usable_for_downstream_simulation']


def test_executable_schema_all_contracts_offline_and_invalid_mode(monkeypatch,scenarios,demo):
    from jsonschema import Draft202012Validator,ValidationError
    from engine.country_tb.composition_common import IMPLEMENTATION
    def forbidden(*args,**kwargs):raise AssertionError('Schema validation must be offline')
    monkeypatch.setattr(socket,'socket',forbidden);monkeypatch.setattr(socket,'create_connection',forbidden)
    schema=json.loads((REPO_ROOT/'schemas/country_tb/v1/baseline_source_composition.schema.json').read_text())
    Draft202012Validator.check_schema(schema);validator=Draft202012Validator(schema)
    checked={}
    def visit(value):
        if isinstance(value,dict):
            if value.get('implementation_version')==IMPLEMENTATION and 'schema_id' in value:
                if value['profile_hash'] not in checked:
                    validator.validate(value)
                    assert value['profile_hash']==profile_hash(value)
                    checked[value['profile_hash']]=value['schema_id']
            for k,v in value.items():
                if k!='receipt':visit(v)
        elif isinstance(value,list):
            for v in value:visit(v)
    visit(scenarios)
    assert len(set(checked.values()))>=9
    bad=deepcopy(overall_selection('bad',demo,'.2'));bad['mode']='automatic_fallback'
    with pytest.raises(ValidationError):validator.validate(bad)
    bad=deepcopy(scenarios['scenarios']['5_incidence_mode_blocked']);bad['tbi_estimate']='0.2'
    with pytest.raises(ValidationError):validator.validate(bad)


def test_age_specific_binomial_anchors_use_counts(demo):
    records=[]
    for sex in ('female','male'):
        records.append({'observation_id':sex,'subgroup_id':'counts','geography_id':GEO,'geography_level':'national',
            'observation_policy':'synthetic_independent_binomial','age_lower':'20','age_upper_exclusive':'21','sex':sex,
            'survey_year':2024,'survey_time':'2024','prevalence':'.2','numerator':200,'denominator':1000,
            'review':REVIEW,'source':SOURCE,'sampling_frame':'invented sample','survey_design':{'kind':'synthetic_only'},
            'uncertainty':{'kind':'unspecified'},'assay_operator':direct_operator(),'weight_source':PROVENANCE,
            'strata':[{'age_lower':'20','age_upper_exclusive':None,'age_policy':'exact_age','sex':sex,'weight':'1'}]})
    obs=observations_profile(records,review=REVIEW)
    policy=calibration_policy(alpha_upper='2',observation_weights={r['observation_id']:'1' for r in records},
        residual_tolerance='.001',review=REVIEW,method='binomial')
    selection=baseline_tbi_source(mode=OBSERVED,geography_id=GEO,geography_level='national',subgroup_id='counts',source=SOURCE,
        provenance=PROVENANCE,review=REVIEW,profile_kind='diagnostic',configuration={'submode':AGE_ANCHORS,'baseline_time':'2024.5',
            'recent_window_years':'5','pressure_history':synthetic_history(),'calibration_policy':policy,'observations':obs})
    result=resolve_baseline_tbi(selection,demo)
    assert float(result['calibration']['alpha'])==pytest.approx(float(-Decimal('.8').ln()/20))
    assert result['calibration']['diagnostics']['likelihood_claim']
    assert all(r['numerator']==200 for r in result['calibration']['rows'])


def test_receipts_and_usability_cannot_be_tampered_with(demo):
    source=direct_selection('g',demo,[total_record()])
    source['configuration']['direct_profile']['receipt']['source_records'][0]['any_tbi']='.9'
    with pytest.raises(DataValidationError,match='receipt'):resolve_baseline_tbi(source,demo)
    bad=deepcopy(demo);bad['receipt']['source_rows'][0]['population']='9000'
    with pytest.raises(DataValidationError,match='receipt'):resolve_baseline_tbi(direct_selection('g',demo,[total_record()]),bad)
    result=resolve_baseline_tbi(direct_selection('g',demo,[total_record()]),demo)
    result['usable_for_downstream_simulation']=True
    with pytest.raises(DataValidationError,match='usability'):validate_resolution(seal(result))


def test_incidence_mode_never_calls_catalytic_or_progression_functions(monkeypatch,demo):
    import engine.country_tb.catalytic_calibration as calibration
    import engine.country_tb.tbi_baseline as baseline
    def forbidden(*args,**kwargs):raise AssertionError('A blocked incidence selection cannot calculate TBI')
    selection=incidence_selection('no_fallback',demo)
    kwargs=selection_kwargs(selection);kwargs['configuration']['pressure_history']=synthetic_history()
    selection=baseline_tbi_source(**kwargs)
    monkeypatch.setattr(calibration,'calibrate_scale',forbidden)
    monkeypatch.setattr(baseline,'baseline_tbi_profile',forbidden)
    assert resolve_baseline_tbi(selection,demo)['status']==BLOCKED


def test_direct_separate_sex_values_and_age_sex_reconciliation(demo):
    rows=[total_record('.1',record_id='female',sex='female'),total_record('.3',record_id='male',sex='male')]
    r=resolve_baseline_tbi(direct_selection('sexes',demo,rows),demo)
    assert r['summary']['components']['any_tbi']['prevalence_full']=='0.2'
    assert {x['any_tbi'] for x in r['rows'] if x['sex']=='female'}=={'0.1'}
    assert {x['any_tbi'] for x in r['rows'] if x['sex']=='male'}=={'0.3'}


def test_incompatible_age_anchors_return_diagnostics_without_fallback(demo):
    records=[]
    for age,p in ((20,'.5'),(60,'.01')):
        records.append({'observation_id':str(age),'subgroup_id':'incompatible','geography_id':GEO,'geography_level':'national',
            'observation_policy':'synthetic_exact_age','age_lower':str(age),'age_upper_exclusive':str(age+1),'sex':'both',
            'survey_year':2024,'survey_time':'2024','prevalence':p,'numerator':None,'denominator':None,'review':REVIEW,
            'source':SOURCE,'sampling_frame':'invented sample','survey_design':{'kind':'synthetic_only'},'uncertainty':{'kind':'exact'},
            'assay_operator':direct_operator(),'weight_source':PROVENANCE,'strata':[{'age_lower':str(age),'age_upper_exclusive':None,
            'age_policy':'exact_age','sex':'both','weight':'1'}]})
    observations=observations_profile(records,review=REVIEW)
    policy=calibration_policy(alpha_upper='2',observation_weights={r['observation_id']:'1' for r in records},
        residual_tolerance='.001',review=REVIEW,method='weighted_least_squares')
    selection=baseline_tbi_source(mode=OBSERVED,geography_id=GEO,geography_level='national',subgroup_id='incompatible',
        source=SOURCE,provenance=PROVENANCE,review=REVIEW,profile_kind='diagnostic',configuration={'submode':AGE_ANCHORS,
            'baseline_time':'2024.5','recent_window_years':'5','pressure_history':synthetic_history(),'calibration_policy':policy,'observations':observations})
    result=resolve_baseline_tbi(selection,demo)
    assert result['status']=='incompatible_anchors' and result['model_derived_profile'] is None
    assert result['summary']['unresolved_population']=='1000'
    assert any(not r['compatible'] for r in result['calibration']['rows'])


def test_explicit_subnational_containment_and_migration_metadata(demo):
    subdemo=population_demography(demo['rows'],geography_id='custom:AUS:NSW',geography_level='subnational',
        population=1000,source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    a=make_subgroup('national',demo,direct_selection('national',demo,[total_record()]))
    b=make_subgroup('local',subdemo,direct_selection('local',subdemo,[total_record()]),
        previous_residence=[{'geography_id':'who:IND','geography_level':'national'}],arrival_year=2010,
        arrival_date='2010-01-01',duration_in_current_setting_years='14.5')
    with pytest.raises(DataValidationError,match='containment'):make_mixture([a,b])
    policy={'mode':'reviewed_containment_v1','parent_geography_id':GEO,'contained_geography_ids':['custom:AUS:NSW'],
            'source':SOURCE,'provenance':PROVENANCE,'review':REVIEW}
    mixture=make_mixture([a,b],geography_policy=policy)
    assert mixture['total_population']=='2000'
    assert b['residence_metadata']['previous_residence'][0]['geography_id']=='who:IND'


def test_population_rows_not_relabelled_and_sex_grids_not_invented(demo):
    rows=deepcopy(demo['rows']);rows[0]['year']=2023
    with pytest.raises(DataValidationError,match='relabelled'):
        population_demography(rows,geography_id=GEO,geography_level='national',population=1000,source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    rows=deepcopy(demo['rows']);rows[0]['age_lower']=21;rows[0]['age_upper_exclusive']=22
    with pytest.raises(DataValidationError,match='age grids'):
        population_demography(rows,geography_id=GEO,geography_level='national',population=1000,source=SOURCE,provenance=PROVENANCE,review=REVIEW)


def test_direct_profile_must_not_split_coarse_population_cell():
    demo=population_demography([{'age_lower':20,'age_upper_exclusive':25,'sex':'both','population':'100'}],
        geography_id=GEO,geography_level='national',population=100,source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    with pytest.raises(DataValidationError,match='cuts a demographic cell'):
        resolve_baseline_tbi(direct_selection('boundary',demo,[total_record(age_lower=20,age_upper_exclusive=21)]),demo)


def test_zero_anchor_and_caller_decimal_context_are_stable(demo):
    selection=overall_selection('zero',demo,'0')
    selection_kwargs_value=selection_kwargs(selection)
    selection_kwargs_value['configuration']['anchor'].update(numerator=0,denominator=100)
    selection=baseline_tbi_source(**selection_kwargs_value)
    normal=resolve_baseline_tbi(selection,demo)
    assert normal['calibration']['alpha']=='0'
    assert normal['summary']['components']['any_tbi']['complete_count']=='0'
    with localcontext() as context:
        context.prec=6
        other=resolve_baseline_tbi(selection,demo)
    assert normal['profile_hash']==other['profile_hash']


def test_direct_empirical_recency_classification_and_no_reallocation(demo):
    record=total_record('.3',recent='.02',remote_only='.28',never_infected='.7')
    selection=direct_selection('empirical',demo,[record],profile_kind='empirical',recency_policy=CONDITIONAL_RECENCY)
    assert selection['recent_remote_composition']=='observed'
    result=resolve_baseline_tbi(selection,demo)
    assert result['recent_remote_composition']=='observed'
    assert all(r['recent']=='0.02' for r in result['rows'])


def test_offline_isolated_interface_and_imports(monkeypatch,demo):
    def forbidden(*args,**kwargs):raise AssertionError('Network is forbidden')
    monkeypatch.setattr(socket,'socket',forbidden);monkeypatch.setattr(socket,'create_connection',forbidden)
    result=resolve_baseline_tbi(direct_selection('offline',demo,[total_record()]),demo)
    assert result['profile_hash']==profile_hash(result)
    spec=interface_specification()
    assert spec['question']=='How should baseline TBI prevalence be determined?'
    assert not next(r for r in spec['options'] if r['mode']==INCIDENCE)['operational']
    for directory in ('app','ui','general_ui','explicit_tbi_ui','engine/apy'):
        for file in (REPO_ROOT/directory).rglob('*.py'):
            assert all(x not in file.read_text(encoding='utf-8-sig') for x in ('baseline_tbi_source','population_subgroups','preview_country_tb_sources'))
    for name in ('baseline_tbi_source','composition_common','active_tb_targets','direct_tbi_profile','population_subgroups'):
        tree=ast.parse((REPO_ROOT/'engine/country_tb'/f'{name}.py').read_text())
        imports=[n.module or '' for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)]
        assert not any(x.startswith(('engine.apy','streamlit','engine.starsim','engine.static','engine.dynamic')) for x in imports)


def test_protected_state_integrity():
    receipt=REPO_ROOT/'abm/output/country_tb_m2d_audit_20261028/starting_state.json'
    if not receipt.is_file():pytest.skip('Session protection receipt unavailable on this checkout.')
    state=json.loads(receipt.read_text());assert state['starting_head']=='ef3269064bbbac32e250bea3953a6a55bce957a8'
    for f,h in state['frozen_files'].items():assert digest((REPO_ROOT/f).read_bytes())==h,f
    for name,old in state['protected'].items():
        folder=REPO_ROOT.parent/name
        def git(*args):return subprocess.check_output(['git','-C',str(folder),*args],stderr=subprocess.DEVNULL)
        assert git('rev-parse','HEAD').decode().strip()==old['head']
        assert git('status','--short').decode().splitlines()==old['status']
        assert digest(git('diff','--binary','HEAD'))==old['tracked_diff_sha256']
        assert digest(git('diff','--cached','--binary'))==old['index_diff_sha256']
        for f,h in old['changed_files'].items():assert digest((folder/f).read_bytes())==h
    refs=dict(x.split(' ',1) for x in subprocess.check_output(['git','for-each-ref','--format=%(refname) %(objectname)']).decode().splitlines())
    allowed={'refs/heads/feature/country-subnational-tb-model-v1','refs/remotes/origin/feature/country-subnational-tb-model-v1'}
    assert {k:v for k,v in refs.items() if k not in allowed}=={k:v for k,v in state['refs'].items() if k not in allowed}
