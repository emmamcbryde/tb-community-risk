"""Scientific boundaries and pure adapters for the isolated3E interface."""
from copy import deepcopy
from decimal import Decimal
from pathlib import Path
import hashlib
import inspect
import json
import socket
import time
import pytest

from engine.country_tb.data_common import canonical,DataValidationError,profile_hash,REPO_ROOT
from engine.country_tb.risk_common import ratio
from country_model_review_ui.state import *
from country_model_review_ui.model import *
from country_model_review_ui.forward import selected_history,_conditional_hazard,RR74_WARNING
from country_model_review_ui.gates import gate_panel,prepare_diagnostic,audit_diagnostic,execute_diagnostic,PUBLIC_STATUS


@pytest.fixture(scope='module')
def example():return working_example()


@pytest.fixture(scope='module')
def run_example(example):return execute_user_run(example)[0]


def entry(factor='diabetes',p='.1',state='value',sex='both'):
    return dict(record_id=factor+sex,subgroup_id='all',factor_id=factor,age_lower=15,age_upper_exclusive=None,sex=sex,
        state=state,prevalence=p,allocation_policy='piecewise_constant_within_supplied_band_v1',
        supplied_exposure_definition=catalogue_records(catalog())[factor]['exposure_definition'],source='Synthetic compatible exposure',review_status='diagnostic',reason='Explicit test input')


def profile_csv(config,records):
    base=dict(subgroup_id='main',recent='',remote_only='',never_infected='',prior_remote_plus_recent='',source='Synthetic TBI table',review_status='diagnostic',
        geography_id=geography(config),baseline_year=2024,allocation_policy='piecewise_constant_within_band_v1')
    return csv_text([{**base,**r} for r in records],TBI_COLUMNS)


def test_user_overall_catalytic(example,run_example):
    s=run_example['baseline']['summary'];assert float(s['components']['any_tbi']['prevalence_full'])==pytest.approx(.2,abs=1e-12)
    assert run_example['baseline']['subgroups'][0]['tbi_profile']['calibration']['alpha']
    assert not run_example['production_available'] and not run_example['externally_validated']


def test_total_only_does_not_invent_recency(example):
    c=deepcopy(example);c['tbi']['derive']=False;r,_=execute_user_run(c)
    assert all(x['recent'] is None and x['remote_only'] is None and x['prior_remote_plus_recent'] is None for x in r['baseline']['rows'])
    assert all(x['any_tbi']=='0.2' for x in r['baseline']['rows'])
    assert r['projection']['instantaneous']['full_population_rate_per100000PY'] is None


@pytest.mark.parametrize('value',['-1','101','nan','Infinity','20%'])
def test_bad_percent_no_clipping(value):
    with pytest.raises(DataValidationError):percent_input(value)


def test_zero_is_not_missing():
    assert percent_input('0')=='0' and percent_input('',required=False) is None
    with pytest.raises(DataValidationError):percent_input('')
    assert percent_input('12.3456789012345678901234567890123456')=='0.123456789012345678901234567890123456'


@pytest.mark.parametrize('bysex',[False,True])
def test_direct_age_sex_preservation(example,bysex):
    c=deepcopy(example);c['tbi'].update(input_kind='Age-by-sex total profile' if bysex else 'Age-specific total profile',derive=False)
    records=[dict(record_id='f',age_lower=0,age_upper_exclusive='',sex='female',any_tbi='0.1000'),
        dict(record_id='m',age_lower=0,age_upper_exclusive='',sex='male',any_tbi='0.3000')] if bysex else [
        dict(record_id='young',age_lower=0,age_upper_exclusive=45,sex='both',any_tbi='0.1000'),
        dict(record_id='older',age_lower=45,age_upper_exclusive='',sex='both',any_tbi='0.3000')]
    c['tbi']['csv']=profile_csv(c,records);before=c['tbi']['csv'];r,_=execute_user_run(c)
    assert r['configuration']['tbi']['csv']==before
    assert float(r['baseline']['summary']['components']['any_tbi']['prevalence_full'])==pytest.approx(.2)
    assert all(x['recent'] is None for x in r['baseline']['rows'])
    assert {'0.1','0.3'}=={x['any_tbi'] for x in r['baseline']['rows']}


def test_direct_recent_remote_no_double_count(example):
    c=deepcopy(example);c['tbi'].update(input_kind='Direct recent and remote',recent_percent='5',remote_percent='15',prior_percent='2',derive=False)
    r,_=execute_user_run(c);s=r['baseline']['summary']['components']
    assert Decimal(s['any_tbi']['complete_count'])==200 and Decimal(s['recent']['complete_count'])==50
    assert Decimal(s['remote_only']['complete_count'])==150 and Decimal(s['prior_remote_plus_recent']['complete_count'])==20
    assert r['projection']['instantaneous']['full_population_rate_per100000PY'] is None
    c['tbi']['conditional_clock']=True;r,_=execute_user_run(c)
    assert r['projection']['instantaneous']['full_population_rate_per100000PY'] is not None
    assert r['baseline']['summary']['components']['recent']['complete_count']=='50'


def test_invalid_direct_states(example):
    c=deepcopy(example);c['tbi'].update(input_kind='Direct recent and remote',recent_percent='60',remote_percent='50')
    with pytest.raises(DataValidationError):build_baseline(c,select_population(c))
    c['tbi'].update(recent_percent='5',remote_percent='15',prior_percent='6')
    with pytest.raises(DataValidationError):build_baseline(c,select_population(c))


def test_explicit_hybrid_retains_total(example):
    c=deepcopy(example);c['tbi'].update(derive=False,hybrid_recency=True)
    r,_=execute_user_run(c)
    assert all(x['any_tbi']=='0.2' and x['recent'] is not None for x in r['baseline']['rows'])
    assert all(x['recency_composition']=='model_derived' for x in r['baseline']['rows'])


def test_mixture_population_weights_and_missing_no_fallback():
    c=working_example(mixture=True);r,_=execute_user_run(c)
    assert float(r['baseline']['summary']['components']['any_tbi']['prevalence_full'])==pytest.approx(.1,abs=1e-12)
    c['population']['subgroups'][1]['tbi_percent']='';r,_=execute_user_run(c)
    assert r['baseline']['summary']['components']['any_tbi']['prevalence_full'] is None
    assert r['baseline']['summary']['unresolved_population']=='200'
    assert r['projection']['instantaneous']['full_population_rate_per100000PY'] is None
    assert r['baseline']['unresolved_subgroups'][0]['population']=='200'


@pytest.mark.parametrize('source_mode',['Unresolved','Incidence inference — blocked'])
def test_incomplete_subgroup_blocks_complete_mixture(source_mode):
    c=working_example(mixture=True);c['population']['subgroups'][0]['source_mode']=source_mode;r,_=execute_user_run(c)
    assert r['baseline']['declared_population']=='1000' and r['baseline']['summary']['unresolved_population']=='700'
    assert not r['baseline']['complete_mixture_available']


def test_no_subgroup_overlap_or_wrong_population():
    c=working_example(mixture=True);c['population']['subgroups'][1]['subgroup_id']=c['population']['subgroups'][0]['subgroup_id']
    with pytest.raises(DataValidationError):build_baseline(c,select_population(c))
    c=working_example(mixture=True);c['population']['partition_confirmed']=False
    with pytest.raises(DataValidationError):build_baseline(c,select_population(c))


def test_csv_schema_provenance_missing_and_assay_distinction(example):
    c=deepcopy(example);text=profile_csv(c,[dict(record_id='all',age_lower=0,age_upper_exclusive='',sex='both',any_tbi='0.12345678901234567890123456789')])
    x=parse_csv(text.encode('utf-8'),'tbi',geography_id=geography(c),subgroup_ids=['main'])
    assert x['source_rows'][0]['any_tbi']=='0.12345678901234567890123456789' and x['rows'][0]['recent'] is None
    assert x['source_sha256']==digest(text.encode('utf-8'))
    for bad in [text.replace('Synthetic TBI table',''),text.replace('0.12345678901234567890123456789','20'),text.replace('custom:synthetic_review','who:AUS'),text.replace('record_id,','wrong_column,')]:
        with pytest.raises(DataValidationError):parse_csv(bad,'tbi',geography_id=geography(c))
    c['tbi'].update(input_kind='Age-specific total profile',csv=text,assay='QFT-GIT')
    with pytest.raises(DataValidationError,match='latent TBI'):build_baseline(c,select_population(c))


def test_direct_profile_gaps_and_overlaps(example):
    c=deepcopy(example);c['tbi'].update(input_kind='Age-specific total profile',derive=False)
    one=dict(record_id='young',age_lower=0,age_upper_exclusive=45,sex='both',any_tbi='.1')
    c['tbi']['csv']=profile_csv(c,[one]);b=build_baseline(c,select_population(c))
    assert b['summary']['unresolved_population']=='500'
    c['tbi']['csv']=profile_csv(c,[one,dict(one,record_id='overlap')])
    with pytest.raises(DataValidationError):build_baseline(c,select_population(c))


def test_sex_specific_risk_highest_rr_and_definition_mismatch(example):
    c=deepcopy(example);c['risk'].update(no_factors=False,legacy_effects=True,
        records=[entry('diabetes','.2'),entry('heavy_smoking','.3',sex='male'),entry('heavy_smoking','.01',sex='female'),entry('haemodialysis','.05')])
    r,_=execute_user_run(c);assert not r['risk_allocation'].get('input_error')
    for row in r['risk_allocation']['rows']:
        assert all(float(ratio(x['effective_rr']))<=25.3 for x in row['effective_rr_distribution'])
        assert row['effective_rr_distribution'][-1]['effective_rr']['value']=='25.3'
    c['risk']['records']=[dict(entry('heavy_smoking','.2'),supplied_exposure_definition='Current smoking')]
    r,_=execute_user_run(c)
    assert r['risk_allocation'].get('input_error') or all(x['any_risk_prevalence'] is None for x in r['risk_allocation']['rows'])
    assert '94.31%' in RR74_WARNING


@pytest.mark.parametrize('control',['-1','0','.4','1'])
def test_association_preserves_all_marginals(example,control):
    c=deepcopy(example);c['risk'].update(no_factors=False,legacy_effects=True,records=[entry('diabetes','.2'),entry('heavy_smoking','.3')])
    base,_=execute_user_run(c);c['association']['global_control']=control;r,_=execute_user_run(c)
    for a,b in zip(base['risk_allocation']['rows'],r['risk_allocation']['rows']):
        assert a['factor_marginals']==b['factor_marginals'] and a['any_risk_prevalence']==b['any_risk_prevalence']
        assert a['tbi_prevalence']==b['tbi_prevalence'] and a['tbi_state_prevalences']==b['tbi_state_prevalences']
        assert Decimal(b['association']['normalized_association_input'])==Decimal(control)
        assert b['association']['achieved_phi'] is not None
        for field in ['recent','remote_only','any_tbi','prior_remote_plus_recent']:
            assert sum(ratio(x['tbi_state_counts'][field]) for x in b['compressed_cells'])==ratio(b['population'])*ratio(b['tbi_state_prevalences'][field])
    c['association']['global_control']='1.1'
    with pytest.raises(DataValidationError):allocate_risk(c,r['baseline'])


@pytest.mark.parametrize('state,p,known',[('value','0',True),('missing',None,False),('not_reviewed',None,False),('excluded',None,True),('not_applicable',None,True)])
def test_risk_input_states_distinct(example,state,p,known):
    c=deepcopy(example);c['risk'].update(no_factors=False,records=[entry(p=p,state=state)])
    b=build_baseline(c,select_population(c));a=allocate_risk(c,b)
    assert all((x['any_risk_prevalence'] is not None)==known for x in a['rows'])


def test_no_tail_unresolved_not_zero(example):
    c=deepcopy(example);c['natural']['tail']='No extrapolation';r,_=execute_user_run(c);p=r['projection']['instantaneous']
    assert p['unresolved_population']>0 and p['covered_instantaneous_case_flow_per_year']>0
    assert p['full_population_rate_per100000PY'] is None
    assert any('remote_tail_external_required_not_zero' in x['unresolved_reasons'] for x in p['rows'])


def test_children_and_open_age_retained(example):
    c=deepcopy(example)
    rows=[dict(age_lower=a,age_upper_exclusive=a+1 if a<100 else '',sex=s,population=n,source='Synthetic',review_status='diagnostic',geography_id=geography(c),baseline_year=2024)
        for a,n in [(3,'50'),(20,'400'),(100,'50')] for s in ('female','male')]
    c['population']['csv']=csv_text(rows,POP_COLUMNS);r,_=execute_user_run(c)
    assert r['context']['open_age_population']=='100' and r['baseline']['summary']['unresolved_population']=='100'
    child=[x for x in r['projection']['instantaneous']['rows'] if x['age_lower']==3]
    assert all(x['covered_instantaneous_case_flow_per_year']==0 and x['unresolved_population']>0 for x in child)
    assert r['projection']['instantaneous']['full_population_rate_per100000PY'] is None


def test_exact_cdf_instantaneous_matches_native_3b_support(run_example):
    from engine.country_tb.progression_baseline import _support
    from engine.country_tb.progression_natural_history import instantaneous_hazard
    group=run_example['risk_allocation']['subgroups'][0];model=group['resolved_baseline']['model_derived_profile'];cell=group['rows'][-1]
    history=selected_history(run_example['configuration']);curve=history['group_curves']['15+']
    for state in ['recent','remote_only']:
        h,gap=_conditional_hazard(model,cell,state,curve)
        nodes=_support(model,cell,state,max_width=1)
        expected=sum(w*instantaneous_hazard(history,tau=t,age_current=a,age_at_infection=a-t,sex=cell['sex']) for a,t,w in nodes)
        assert h==pytest.approx(expected,abs=1e-13) and gap==0


def test_stale_hash_and_exports_preserve_precision(example,run_example):
    from country_model_review_ui.presentation import prepare_exports,percent
    s={'review_config':deepcopy(example),'review_result':run_example};assert not stale(s)
    before=canonical(run_example);display=percent(run_example['baseline']['rows'][0]['any_tbi']);assert display.endswith('%')
    export=prepare_exports(example,run_example);v=json.loads(export['json']);assert canonical(v['run'])==before
    assert 'exact' in export['json'].decode('utf8') and 'run_scientific_hash' in export['csv'].decode('utf8')
    s['review_config']['display']['sexes']='Combined display';assert not stale(s)
    s['review_config']['tbi']['overall_percent']='21';assert stale(s)
    with pytest.raises(DataValidationError):prepare_exports(s['review_config'],run_example)
    assert canonical(run_example)==before


def make_review_bundle(config):
    from engine.country_tb.preview_incidence_inversion import make_mapping,make_operator,constant_history,stratum,observation,synthetic_gate
    from engine.country_tb.incidence_inversion import incidence_observation
    from engine.country_tb.incidence_mapping import IncidenceMapping
    geo=geography(config)
    rows=[dict(stratum(s+str(a),age=a+.25,sex=s,population=250),subgroup_id='main',
        age_support_policy='reviewed_within_cell_age_distribution',age_distribution_review={'status':'reviewed','source_id':'explicit_synthetic_age_support'})
        for s in ['female','male'] for a in [20,60]]
    m0=make_mapping(rows=rows,history=constant_history(geo),operator=make_operator(geography=geo),geography=geo)
    args=deepcopy(m0._constructor);args['subgroup_id']='main';args['kernel_width']=.5;m=IncidenceMapping(**args)
    obs0=observation(m.incidence(.017),geography=geo);obsargs={k:v for k,v in obs0.items() if k in inspect.signature(incidence_observation).parameters}
    obsargs.update(subgroup_id='main',provenance_record=obs0['provenance']);obs=incidence_observation(**obsargs);gate=synthetic_gate(m,obs)
    return dict(review_binding_hash=review_binding_hash(config),component_reviews={i:dict(status='Resolved and reviewed',source_id='synthetic:M3E_mathematics_only',reviewer='Fictional test protocol',rationale='Invented fixture,not country scientific approval') for i in GATE_ITEMS},
        mapping_inputs=args,observation=obs,unresolved_threshold_approval=gate['unresolved_threshold_approval'],
        scientific_review=gate['scientific_review'],concentration_policy=gate['concentration_policy'],
        numerical_protocol=dict(status='reviewed',source_id='synthetic_test_only',policy='retained_3C_bounded_grid_all_roots_and_refinement_v1',absolute_refinement_tolerance=.001),
        assay_decision=dict(status='reviewed',source_id='synthetic_direct_latent_fixture'))


@pytest.fixture(scope='module')
def diagnostic_fixture(example):
    c=deepcopy(example);c['tbi']['mode']=SOURCE_MODES[1];c['diagnostic']['enabled']=True
    context=select_population(c);bundle=make_review_bundle(c);c['diagnostic']['bundle']=json.dumps(bundle)
    return c,context,bundle


def test_incidence_disabled_no_fallback(example,monkeypatch):
    c=deepcopy(example);c['tbi']['mode']=SOURCE_MODES[1]
    with pytest.raises(DataValidationError):execute_user_run(c)
    with pytest.raises(DataValidationError,match='disabled'):prepare_diagnostic(c,select_population(c))
    c['diagnostic']['enabled']=True
    with pytest.raises(DataValidationError,match='Upload'):prepare_diagnostic(c,select_population(c))
    assert gate_panel(c)['public_status']==PUBLIC_STATUS


@pytest.mark.parametrize('item',GATE_ITEMS)
def test_each_unresolved_gate_prevents_mapping(diagnostic_fixture,item,monkeypatch):
    c,context,bundle=deepcopy(diagnostic_fixture);bundle['component_reviews'][item]['status']='Unresolved';c['diagnostic']['bundle']=json.dumps(bundle)
    import country_model_review_ui.gates as gate_module
    monkeypatch.setattr(gate_module,'IncidenceMapping',lambda **kwargs:pytest.fail('Mapping constructed before science gate resolution'))
    with pytest.raises(DataValidationError,match='Unresolved scientific gate'):prepare_diagnostic(c,context)


def test_guarded_uploaded_synthetic_diagnostic_and_no_direct_bypass(diagnostic_fixture):
    c,context,bundle=deepcopy(diagnostic_fixture)
    with pytest.raises(DataValidationError,match='unresolved'):execute_diagnostic(c,context)
    report=audit_diagnostic(c,context);assert report['passed'] and not report['fitting_accepted']
    result=execute_diagnostic(c,context,numerical_clearance=report)
    assert result['diagnostic']['fitted_scale']==pytest.approx(.017,abs=1e-10)
    assert result['public_status']==PUBLIC_STATUS and not result['production_available']
    c['association']['global_control']='.3'
    with pytest.raises(DataValidationError,match='unresolved'):execute_diagnostic(c,context,numerical_clearance=report)


def test_prk_no_invented_incidence(example):
    c=deepcopy(example);c['population']['country']='PRK';context=select_population(c)
    assert context['observed_incidence']['central'] is None
    c['diagnostic']['enabled']=True
    with pytest.raises(DataValidationError,match='unavailable'):prepare_diagnostic(c,context)


def test_offline_sources_and_retained_fingerprints(example,monkeypatch):
    monkeypatch.setattr(socket,'create_connection',lambda *a,**k:pytest.fail('Unexpected network call'))
    assert execute_user_run(example)[0]['baseline']['summary']['population']=='1000'
    manifest=json.loads((REPO_ROOT/'docs/country_subnational_model/protected_fingerprints_3d.json').read_text(encoding='utf8'))
    for field in ['frozen_files','retained_2a_3c_files']:
        for p,h in manifest[field].items():assert hashlib.sha256((REPO_ROOT/p).read_bytes()).hexdigest()==h


def test_ui_schemas_and_no_production_activation(example,run_example):
    from jsonschema import Draft202012Validator
    schema=json.loads((REPO_ROOT/'schemas/country_tb/v1/review_ui.schema.json').read_text(encoding='utf8'))
    Draft202012Validator.check_schema(schema);v=Draft202012Validator(schema)
    for payload in [example,run_example,run_example['baseline'],run_example['risk_allocation'],run_example['projection'],gate_panel(example)]:v.validate(payload)
    fake=deepcopy(run_example);fake['production_available']=True;assert list(v.iter_errors(fake))


def test_source_population_cache_cannot_be_reused_after_edit(example):
    c=deepcopy(example);context=select_population(c);c['population']['label']='Changed population definition'
    with pytest.raises(DataValidationError,match='stale'):build_baseline(c,context)


def test_uploaded_diagnostic_cannot_replace_selected_risk_association(diagnostic_fixture):
    c,context,bundle=deepcopy(diagnostic_fixture)
    bundle['mapping_inputs']['strata'][0]['association_control']=.5;c['diagnostic']['bundle']=json.dumps(bundle)
    with pytest.raises(DataValidationError,match='association differs'):prepare_diagnostic(c,context)
    c,context,bundle=deepcopy(diagnostic_fixture)
    bundle['mapping_inputs']['strata'][0]['factors']=[{'factor_id':'diabetes'}];c['diagnostic']['bundle']=json.dumps(bundle)
    with pytest.raises(DataValidationError,match='RR1'):prepare_diagnostic(c,context)


def test_example_files_validate_as_actual_uploads():
    folder=REPO_ROOT/'country_model_review_ui/examples'
    for name,kind in [('population','population'),('tbi_age_sex','tbi'),('tbi_states','tbi'),('risk','risk'),('association','association')]:
        value=parse_csv((folder/(name+'.csv')).read_bytes(),kind,subgroup_ids=['main'])
        assert value['source_rows'] and value['source_sha256']


def test_explicit_scenario_comparison_includes_tail_and_mortality(example):
    c=deepcopy(example);c['natural']['compare_scenarios']=True;r,_=execute_user_run(c)
    rows=r['projection']['scenarios'];assert len(rows)==6
    no_tail=next(s for s in rows if s['parameter_set_id']=='Remote tail: No extrapolation')
    assert no_tail['full_population_rate_per100000PY'] is None and no_tail['unresolved_population']>0
    assert any(s['parameter_set_id'].startswith('Mortality:') for s in rows)


def test_tiny_positive_unresolved_population_is_never_called_complete(example):
    c=deepcopy(example);c['tbi'].update(derive=False,overall_percent='0.0000000000000000000000000001')
    r,_=execute_user_run(c);v=r['projection']['instantaneous']
    assert 0<v['unresolved_population']<1e-9 and v['full_population_rate_per100000PY'] is None


def test_mortality_interval_scoping_matches_full_source_reference(example):
    c=deepcopy(example);c['natural'].update(horizon='One-year prospective (bounded)',
        mortality='ABS2022–2024 validation snapshot — sensitivity',mortality_proxy=True)
    r,_=execute_user_run(c);p=r['projection']['prospective'];assert p is not None
    assert float(p['projections'][0]['expected_cases_full'])==pytest.approx(.37144663879619993,abs=1e-12)
    mortality=p['operator']['competing_mortality']
    assert len(mortality['rows'])==8 and mortality['snapshot_hash']
    assert mortality['provenance']['source_interval_count']==200
