"""Bounded offline M3D source, comparison and scientific-boundary tests."""
from copy import deepcopy
from decimal import Decimal
from math import exp,expm1,log
from pathlib import Path
import ast
import hashlib
import json
import socket
from uuid import uuid4

import pytest

from engine.country_tb.data_common import DataValidationError,canonical,profile_hash,REPO_ROOT,seal
from engine.country_tb.validation_data import *
from engine.country_tb.empirical_validation import *
from engine.country_tb.review_parameters_3d import progression_sets,review_parameter_set,assay_scenario,evidence_registry
from engine.country_tb.validation_numerics import one_year_kernel,cohort_reference
from engine.country_tb.preview_incidence_inversion import (
    make_mapping,make_operator,constant_history,observation,synthetic_gate,DMETA,stratum,
)
from engine.country_tb.progression_natural_history import audited_adult_candidate,curve_integral
from engine.country_tb.incidence_inversion import PUBLIC_STATUS

SNAPSHOT='tbi-validation-fe05d6d1ab8c'
MORTALITY='abs-life2022-2024-6785449a5200'


@pytest.fixture(scope='module')
def inputs():return load_validation_snapshot(SNAPSHOT)


def selected(inputs,**fields):
    return [r for r in inputs['rows'] if all(r[k]==v for k,v in fields.items())]


def raw(inputs):
    return deepcopy(next(r for r in inputs['source_records'] if r['source_id']=='PMC10446657' and r['assay']=='QFT-Plus' and r['age_lower']==30))


def test_frozen_sources_and_exact_counts(inputs):
    assert len(inputs['rows'])==96 and len(inputs['sources'])==3
    assert not inputs['production_available'] and not inputs['externally_validated']
    r=next(r for r in inputs['rows'] if r['source_id']=='PMC10446657' and r['assay']=='QFT-Plus' and r['age_lower']==30 and r['survey_year']==2020)
    assert (r['numerator'],r['denominator'])==(85,217)
    assert float(r['observed_positivity'])==pytest.approx(85/217,abs=1e-16)
    assert r['source_preserved']


def test_weighted_survey_does_not_invent_counts(inputs):
    rows=selected(inputs,source_id='PMC4633161')
    assert rows and all(r['numerator'] is None and r['denominator'] is None for r in rows)
    assert all(r['lower'] is not None and r['upper'] is not None for r in rows)


def test_migration_assay_and_source_cells_preserved(inputs):
    for original in inputs['source_records']:
        transformed=next(r for r in inputs['rows'] if r['record_id']==original['record_id'])
        assert transformed['migration']==original['migration']
        assert transformed['assay']==original['assay'] and transformed['threshold']==original['threshold']
        assert transformed['source_record_sha256']==digest(canonical(original))
        assert original['source_cells'] and original['citation'] and original['url']


@pytest.mark.parametrize('field',['citation','url','migration','threshold','source_cells','review_status'])
def test_missing_source_fields_rejected(inputs,field):
    r=raw(inputs);r.pop(field)
    with pytest.raises(DataValidationError):normalize_observation(r)


@pytest.mark.parametrize('value',[-1,101,'NaN',True])
def test_invalid_percent_rejected(inputs,value):
    r=raw(inputs);r.update(numerator=None,denominator=None,reported_prevalence_percent=value)
    with pytest.raises(DataValidationError):normalize_observation(r)


def test_zero_distinct_from_missing_and_counts_checked(inputs):
    r=raw(inputs);r.update(numerator=0,denominator=10,reported_prevalence_percent='0')
    assert normalize_observation(r)['observed_positivity']=='0'
    for change in [dict(numerator=None,denominator=None,reported_prevalence_percent=None),
        dict(numerator=11),dict(numerator=None),dict(denominator=0),dict(reported_prevalence_percent=90)]:
        with pytest.raises(DataValidationError):normalize_observation(dict(r,**change))


def test_canonical_row_order_and_hash(inputs):
    a=validation_inputs(inputs['source_records'],inputs['sources'])
    b=validation_inputs(list(reversed(inputs['source_records'])),list(reversed(inputs['sources'])))
    assert canonical(a)==canonical(b) and a['profile_hash']==b['profile_hash']
    c=deepcopy(a);c['rows'][0]['observed_positivity']='.999'
    with pytest.raises(DataValidationError):verify(c,OBSERVATIONS)


def test_immutable_snapshot_reimport_and_checksum(inputs):
    # pytest's Windows private temp ACL is inaccessible to restricted tokens.
    # Ordinary inherited workspace ACL; retain every fixture, no cleanup/delete.
    tmp_path=REPO_ROOT/'abm/output/country_tb_m3d_test_fixtures'/uuid4().hex
    tmp_path.mkdir(parents=True,exist_ok=False)
    path=freeze_validation_snapshot(inputs['source_records'],inputs['sources'],directory=tmp_path)
    assert freeze_validation_snapshot(inputs['source_records'],inputs['sources'],directory=tmp_path)==path
    assert load_validation_snapshot(path.name,directory=tmp_path)==inputs
    # Tamper only a new test copy, never the repository snapshot.
    with (path/'source_observations.json').open('ab') as stream:stream.write(b' ')
    with pytest.raises(DataValidationError,match='checksum'):load_validation_snapshot(path.name,directory=tmp_path)
    with pytest.raises(DataValidationError,match='Immutable'):freeze_validation_snapshot(inputs['source_records'],inputs['sources'],directory=tmp_path)


def test_unlicensed_source_not_published(inputs):
    sources=deepcopy(inputs['sources']);sources[0]['licence']['redistribution']='unknown'
    with pytest.raises(DataValidationError):validation_inputs(inputs['source_records'],sources)


def test_split_rejects_reuse_paired_assay_parent_and_nativity_overlap(inputs):
    v=selected(inputs,source_id='PMC10446657',survey_year=2020,age_lower=30,view='age')
    a=next(r for r in v if r['assay']=='QFT-Plus');b=next(r for r in v if r['assay']=='TST')
    for second in [a,b,next(r for r in inputs['rows'] if r['source_id']=='PMC10446657' and r['survey_year']==2020 and r['view']=='overall' and r['assay']=='QFT-Plus')]:
        with pytest.raises(DataValidationError):validation_split(inputs,[a['record_id']],[second['record_id']],kind='held_out_age_records',rationale='fixture')
    us=selected(inputs,source_id='PMC4633161',age_lower=15,view='age')
    allr=next(r for r in us if r['subgroup_id']=='all');born=next(r for r in us if r['subgroup_id']=='us_born')
    with pytest.raises(DataValidationError):validation_split(inputs,[allr['record_id']],[born['record_id']],kind='held_out_age_records',rationale='fixture')


def test_split_sorted_and_cluster_caveat(inputs):
    rows=selected(inputs,source_id='PMC4633161',subgroup_id='us_born',view='age')
    ids=[r['record_id'] for r in rows];a=validation_split(inputs,ids[:1],ids[1:],kind='held_out_age_records',rationale='fixture')
    b=validation_split(inputs,ids[:1],list(reversed(ids[1:])),kind='held_out_age_records',rationale='fixture')
    assert a==b and not a['independence_claim'] and a['cluster_caveat']
    with pytest.raises(DataValidationError):validation_split(inputs,ids[:1],ids[1:],kind='independent_survey',rationale='false_independence')


@pytest.mark.parametrize('identifier',[None,'central','unknown'])
def test_explicit_parameter_selection_no_fallback(identifier):
    with pytest.raises(DataValidationError):review_parameter_set(identifier)


@pytest.mark.parametrize('identifier',['conservative','central_candidate','higher_progression'])
def test_explicit_parameter_selection_still_pending(identifier):
    assert review_parameter_set(identifier)['scientific_approval']=='pending'


@pytest.mark.parametrize('product,se,sp',[('QFT-Plus',.89,.98),('QFT-GIT',.81,.99)])
def test_assay_transformation_fixed_and_product_specific(product,se,sp):
    a=assay_scenario(product)
    for p in [0,.2,1]:assert float(assay_probability(p,a['operator']))==pytest.approx(se*p+(1-sp)*(1-p))
    assert a['pediatric_performance'] is None and not a['production_available']
    assert a['provenance']['scenario'] and a['operator']['name']==product
    assert assay_scenario(product,variant='lower_specificity')['profile_hash']!=a['profile_hash']
    with pytest.raises(DataValidationError):assay_scenario('TST')


def test_metric_units_errors_missing_and_no_threshold_approval(inputs):
    r=selected(inputs,source_id='PMC4633161',subgroup_id='us_born',age_lower=15,view='age')[0]
    r=deepcopy(r);r['observed_positivity']='.2';r['lower']='.1';r['upper']='.3'
    w={r['record_id']:1};p={r['record_id']:{'predicted_positivity':'.25'}}
    v=comparison_metrics([r],p,weights=w,weight_policy='explicit_fixture',source_hash='a'*64,scenario_hash='b'*64)
    a=v['rows'][0]
    assert a['absolute_error']==pytest.approx(.05) and a['relative_error']==pytest.approx(.25)
    assert a['calibration_ratio']==pytest.approx(1.25) and v['weighted_RMSE']==pytest.approx(.05)
    assert v['interval_inclusion_fraction']==1 and v['unresolved_tested_population'] is None
    assert not v['thresholds_approved']
    v=comparison_metrics([r],{r['record_id']:{'predicted_positivity':None,'unavailable_reason':'external_required'}},
        weights=w,weight_policy='explicit_fixture',source_hash='a'*64,scenario_hash='b'*64)
    assert v['weighted_RMSE'] is None and v['covered_metric_weight']==0
    with pytest.raises(DataValidationError):comparison_metrics([r],{},weights=w,weight_policy='fixture',source_hash='a'*64,scenario_hash='b'*64)


def test_metric_zero_and_sex_residuals_and_row_order(inputs):
    a=deepcopy(inputs['rows'][0]);b=deepcopy(a)
    a.update(record_id='female',sex='female',observed_positivity='0',denominator=10)
    b.update(record_id='male',sex='male',observed_positivity='.2',denominator=20)
    p={'female':{'predicted_positivity':'.1'},'male':{'predicted_positivity':'.3'}}
    kwargs=dict(weights={'female':10,'male':20},weight_policy='synthetic_tested_count',source_hash='0'*64,scenario_hash='1'*64)
    v=comparison_metrics([a,b],p,**kwargs);w=comparison_metrics([b,a],p,**kwargs)
    assert v==w and v['rows'][0]['calibration_ratio'] is None
    assert v['sex_specific']['male']['weighted_mean_residual']==pytest.approx(.1)
    assert v['sex_specific']['female']['weighted_mean_residual']==pytest.approx(.1)


def test_tbi_anchor_real_input_has_heldout_and_unresolved_ages(inputs):
    from scripts.validate_country_tb_empirical_3d import _history
    rows=selected(inputs,source_id='PMC10446657',survey_year=2020,view='age',assay='QFT-Plus')
    anchor=next(r for r in rows if r['age_lower']==30);valid=[r for r in rows if r!=anchor]
    split=validation_split(inputs,[anchor['record_id']],[r['record_id'] for r in valid],kind='held_out_age_records',rationale='prespecified_test')
    before=canonical(inputs);history,_=_history('VNM')
    v=tbi_anchored_validation(inputs,split,history,assay_scenario('QFT-Plus'),
        metric_weights={r['record_id']:r['denominator'] for r in valid},weight_policy='tested_sample_descriptive',scenario_id='fixture')
    assert float(v['fitted_scale'])==pytest.approx(.011211954245899264,abs=1e-10)
    assert v['held_out_metrics']['weighted_RMSE']==pytest.approx(.08254881625,abs=1e-9)
    assert any(r['unavailable_reason'] for r in v['held_out_metrics']['rows'])
    assert v['sources_unchanged'] and canonical(inputs)==before
    assert v['incidence_comparison_status']=='unavailable_no_population_compatible_first_episode_target'


def test_unresolved_anchor_and_cross_country_refused(inputs):
    from scripts.validate_country_tb_empirical_3d import _history
    rows=selected(inputs,source_id='PMC10446657',survey_year=2020,view='age',assay='QFT-Plus')
    child=next(r for r in rows if r['age_lower']==0);adult=next(r for r in rows if r['age_lower']==30)
    split=validation_split(inputs,[child['record_id']],[adult['record_id']],kind='held_out_age_records',rationale='fixture')
    history,_=_history('VNM');kwargs=dict(metric_weights={adult['record_id']:1},weight_policy='fixture',scenario_id='fixture')
    v=tbi_anchored_validation(inputs,split,history,assay_scenario('QFT-Plus'),**kwargs)
    assert v['status']=='blocked_anchor_inputs' and v['fitted_scale'] is None
    split=validation_split(inputs,[adult['record_id']],[child['record_id']],kind='held_out_age_records',rationale='fixture')
    history,_=_history('USA')
    with pytest.raises(DataValidationError):tbi_anchored_validation(inputs,split,history,assay_scenario('QFT-Plus'),
        metric_weights={child['record_id']:1},weight_policy='fixture',scenario_id='fixture')


def synthetic_input(inputs):
    r=raw(inputs);r.update(record_id='synthetic_independent_TBI',source_id='synthetic',citation='Invented mathematical fixture',
        url='https://example.org/synthetic-not-empirical',geography_id='synthetic:closed',geography_level='custom',
        subgroup_id='g',participant_pool='fictional_independent_survey',membership_group='all',survey_year=2024,
        age_lower=40,age_upper_exclusive=41,numerator=50,denominator=100,reported_prevalence_percent='50',
        source_cells=['invented:50/100'],source_locator='synthetic_fixture',uncertainty={'kind':'not_reported'})
    s={'source_id':'synthetic','citation':r['citation'],'url':r['url'],'access_date':r['access_date'],
        'source_bytes_sha256':'0'*64,'licence':{'redistribution':'permitted_with_attribution','url':r['url']},'extraction':'invented_fixture'}
    return validation_inputs([r],[s])


def test_incidence_direction_synthetic_recovery_and_independent_TBI_not_fit(inputs):
    v=synthetic_input(inputs);m=make_mapping();obs=observation(m.incidence(.017));g=synthetic_gate(m,obs);ids=[v['rows'][0]['record_id']]
    result=incidence_anchored_validation(m,obs,g,v,ids,assay_scenario('QFT-Plus'),lower=0,upper=.2,
        metric_weights={ids[0]:1},weight_policy='fictional',independence_provenance={'TBI_observations_used_to_construct_incidence_target':False},**DMETA)
    assert result['fitted_scale']==pytest.approx(.017,abs=1e-10)
    assert not result['TBI_observations_used_in_incidence_fit'] and result['held_out_metrics']['weighted_RMSE']>0
    assert result['diagnostic_status']=='diagnostic_only_not_approved_for_production'
    assert v['rows'][0]['numerator']==50 and not result['production_available']
    with pytest.raises(DataValidationError):incidence_anchored_validation(m,obs,g,v,ids,assay_scenario('QFT-Plus'),lower=0,upper=.2,
        metric_weights={ids[0]:1},weight_policy='fictional',independence_provenance={'TBI_observations_used_to_construct_incidence_target':True},**DMETA)


@pytest.mark.parametrize('target_type',['notifications','baseline_prevalent_disease','imported_active_disease'])
def test_incidence_direction_incompatible_observation_stays_blocked(inputs,target_type):
    v=synthetic_input(inputs);m=make_mapping();obs=observation(100,target_type=target_type);ids=[v['rows'][0]['record_id']]
    result=incidence_anchored_validation(m,obs,synthetic_gate(m,obs),v,ids,assay_scenario('QFT-Plus'),lower=0,upper=.2,
        metric_weights={ids[0]:1},weight_policy='fictional',independence_provenance={'TBI_observations_used_to_construct_incidence_target':False},**DMETA)
    assert result['fitted_scale'] is None and result['held_out_metrics']['weighted_RMSE'] is None
    assert result['incidence_calibration']['public_status']==PUBLIC_STATUS


def test_mortality_snapshot_and_no_country_fallback():
    v=load_mortality_validation_snapshot(MORTALITY,geography_id='who:AUS')
    assert len(v['source_rows'])==202 and len(v['rows'])==200
    r=next(r for r in v['rows'] if r['age_lower']==0 and r['sex']=='male')
    assert Decimal(r['survival_probability_one_year'])==Decimal('1')-Decimal('.00353')
    assert max(r['age_upper_exclusive'] for r in v['rows'])==100
    assert v['metadata']['licence']['name']=='CC BY4.0'
    with pytest.raises(DataValidationError):load_mortality_validation_snapshot(MORTALITY,geography_id='who:VNM')


def test_parameter_sets_are_separate_not_approved_or_defaults():
    sets=progression_sets();assert [s['parameter_set_id'] for s in sets]==['conservative','central_candidate','higher_progression']
    assert len({s['profile_hash'] for s in sets})==3
    for s in sets:
        assert not s['selected_by_default'] and s['central_adult_curve'] is None and s['central_remote_tail'] is None
        assert s['pediatric_curves']=={'0-4':None,'5-14':None} and s['scientific_approval']=='pending'
        assert all(float(r['hazard_per_year'])>=0 for r in s['segment_hazards'])
        knots=s['cumulative_risk_knots'];assert [r['time_years'] for r in knots]==['0','1','5','10','25']
        for a,b in zip(knots,knots[1:]):assert float(b['cumulative_hazard'])>=float(a['cumulative_hazard'])


def test_factor_evidence_qualification_and_extreme_warning():
    e=evidence_registry();assert len(e['factor_rows'])==11
    assert all(not r['applicable_as_progression_hazard_multiplier'] for r in e['factor_rows'])
    assert next(r for r in e['factor_rows'] if r['factor_id']=='cxr_granuloma')['status']=='incompatible_definition'
    assert e['RR74_warning']['one_year_probability']==pytest.approx(1-.962**74)
    assert not e['externally_validated'] and e['no_independent_cohort_double_counting']


@pytest.mark.parametrize('tau',[None,0,4.5,30])
@pytest.mark.parametrize('rr',[1,74])
def test_independent_competing_event_kernel(tau,rr):
    curve=audited_adult_candidate(tail='fixed_remote_hazard_v1',tail_parameters={'hazard_per_year':'.00038'})['group_curves']['15+']
    v=one_year_kernel(curve,tau,rr=rr,mortality_rate=.08,future_pressure=.05)
    assert v['cases']+v['deaths']+v['surviving_disease_free']==pytest.approx(1,abs=1e-11)
    assert min(v.values())>=0 and 0<v['person_years']<=1
    if tau is not None:
        base=one_year_kernel(curve,tau,rr=rr)
        expected=-expm1(-rr*curve_integral(curve,tau,tau+1))
        assert base['cases']==pytest.approx(expected,abs=2e-14)


def test_independent_last_event_quadrature_analytic():
    m=make_mapping();curve=m.operator['natural_history']['group_curves']['15+']
    v=cohort_reference(m.history,[stratum()],curve,.017,instantaneous=True)
    assert v['incidence']==pytest.approx(2000*(-expm1(-.017*40)),abs=1e-10)


def test_discrepancy_and_reviewer_response():
    names=['WHO_estimation_uncertainty','observation_incompatibility','migration','imported_disease','recurrence','prior_TB',
        'unobserved_risk_factors','assay_error','historical_incidence_revision','pressure_shape_misspecification','endogenous_transmission','age_sex_aggregation','mortality']
    rows=[dict(component=n,status='unresolved',direction='unknown',magnitude=None,evidence='missing',required_data='review') for n in names]
    a=model_discrepancy(source_hash='0'*64,scenario_hash='1'*64,entries=rows)
    b=model_discrepancy(source_hash='0'*64,scenario_hash='1'*64,entries=list(reversed(rows)))
    assert a==b and not a['discrepancy_fitted_into_alpha']
    with pytest.raises(DataValidationError):model_discrepancy(source_hash='0'*64,scenario_hash='1'*64,entries=rows[:-1])
    r=reviewer_response(['remote_tail','adult']);assert not r['automatically_approved']
    assert all(a['decision'] is None and a['reviewer'] is None for a in r['rows'])


def test_all_m3d_schemas_validate_offline(inputs):
    from jsonschema import Draft202012Validator
    s=json.loads((REPO_ROOT/'schemas/country_tb/v1/empirical_validation.schema.json').read_text(encoding='utf-8'))
    Draft202012Validator.check_schema(s);validator=Draft202012Validator(s)
    values=[inputs,*progression_sets(),assay_scenario('QFT-Plus'),evidence_registry(),reviewer_response(['adult']),
        load_mortality_validation_snapshot(MORTALITY,geography_id='who:AUS')]
    for v in values:
        validator.validate(v);assert v['profile_hash']==profile_hash(v)
    invalid=deepcopy(inputs);invalid['production_available']=True
    assert list(validator.iter_errors(invalid))


def test_no_tail_not_zero_and_children_explicit_in_candidates():
    from engine.country_tb.progression_natural_history import curve_hazard
    c=audited_adult_candidate()['group_curves']['15+']
    with pytest.raises(DataValidationError,match='external_required'):curve_hazard(c,26)
    assert all(s['pediatric_curves']=={'0-4':None,'5-14':None} for s in progression_sets())


def test_offline_loading_and_no_production_imports(monkeypatch):
    monkeypatch.setattr(socket,'create_connection',lambda *a,**k:pytest.fail('Network access forbidden'))
    assert load_validation_snapshot(SNAPSHOT)['rows']
    assert load_mortality_validation_snapshot(MORTALITY,geography_id='who:AUS')['rows']
    for root in ['app','apps','engine/apy','explicit_tbi_ui']:
        for p in (REPO_ROOT/root).rglob('*.py'):
            if p.is_file():assert 'country_tb.empirical_validation' not in p.read_text(encoding='utf-8')
    for name in ['validation_data','empirical_validation','validation_numerics','review_parameters_3d']:
        tree=ast.parse((REPO_ROOT/'engine/country_tb'/f'{name}.py').read_text(encoding='utf-8'))
        for node in ast.walk(tree):
            if isinstance(node,ast.Import):assert not any(a.name.split('.')[0] in {'requests','urllib','socket','streamlit'} for a in node.names)


def test_retained_public_status_and_protection_manifest():
    assert PUBLIC_STATUS=='unavailable_pending_progression_and_risk_forward_operator'
    start=json.loads((REPO_ROOT/'docs/country_subnational_model/protected_fingerprints_3c.json').read_text(encoding='utf-8'))
    assert start['passed'] and (REPO_ROOT/'engine/country_tb/baseline_tbi_source.py').is_file()
    for path,h in start['retained_2a_3b_files'].items():
        assert hashlib.sha256((REPO_ROOT/path).read_bytes()).hexdigest()==h
