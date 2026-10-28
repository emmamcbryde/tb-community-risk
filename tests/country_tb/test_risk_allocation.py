"""Exact risk allocation using invented data and explicitly provisional diagnostics."""
from copy import deepcopy
from decimal import localcontext
from fractions import Fraction
import ast
import json
import socket
import subprocess

import pytest

from engine.country_tb.data_common import DataValidationError,canonical,profile_hash,REPO_ROOT,digest
from engine.country_tb.risk_common import ratio,amount,text_number,clinical_group
from engine.country_tb.risk_common import IMPLEMENTATION,artifact
from engine.country_tb.composition_common import population_demography
from engine.country_tb.risk_factor_catalog import progression_risk_factor_catalog,catalogue_records,applicable_effect,FACTOR_IDS
from engine.country_tb.risk_prevalence import (
    risk_factor_prevalence_profile,definition_mapping,resolve_risk_prevalence,DIRECT,CONSTANT,UNIFORM,PATTERN,
)
from engine.country_tb.risk_policies import (
    overlap_policy,association_policy,control_for_cell,feasible_association,phi,multimorbidity_policy,
    combine_risk_factors,HIGHEST,MULTIPLICATIVE,ADDITIVE,JOINT,
)
from engine.country_tb.risk_allocation import allocate_baseline_risk,combine_risk_allocations
from engine.country_tb.preview_risk_allocation import toy_demography,prevalence_record,allocate_example,diagnostics
from engine.country_tb.baseline_tbi_source import resolve_baseline_tbi
from scripts.preview_country_tb_sources import REVIEW,SOURCE,PROVENANCE,GEO,direct_selection,incidence_selection,make_subgroup,make_mixture


@pytest.fixture(scope='module')
def catalog():return progression_risk_factor_catalog()


@pytest.fixture(scope='module')
def resolved():
    demo=toy_demography()
    record={'record_id':'invented','age_lower':0,'age_upper_exclusive':None,'sex':'both',
        'age_allocation_policy':'piecewise_constant_within_band_v1','any_tbi':'.3','recent':'.1','remote_only':'.2',
        'never_infected':'.7','prior_remote_plus_recent':'.04'}
    return resolve_baseline_tbi(direct_selection('risk_test',demo,[record]),demo)


def profile(catalog,records,subgroup='risk_test'):
    return risk_factor_prevalence_profile(records,catalog=catalog,factor_ids=sorted({r['factor_id'] for r in records}),
        subgroup_id=subgroup,geography_id=GEO,source=SOURCE,provenance=PROVENANCE,review=REVIEW)


def policy(method=HIGHEST,**kwargs):
    return multimorbidity_policy(source=SOURCE,provenance=PROVENANCE,review=REVIEW,method=method,**kwargs)


def allocate(resolved,catalog,records,*,c='0',overlap=None,progression=None,diagnostic=True):
    return allocate_baseline_risk(resolved,catalog=catalog,prevalence=profile(catalog,records,resolved['subgroup_id']),
        overlap=overlap or overlap_policy(source=SOURCE,provenance=PROVENANCE,review=REVIEW),
        association=association_policy(global_control=c,source=SOURCE,provenance=PROVENANCE,review=REVIEW),
        progression_policy=progression or policy(effect_use='diagnostic_provisional_adult_effects_v1'),
        source=SOURCE,provenance=PROVENANCE,review=REVIEW,diagnostic_only=diagnostic)


def test_catalogue_values_definitions_provenance_and_separate_review(catalog):
    expected={'cxr_granuloma':'2','heavy_smoking':'3','underweight':'3','diabetes':'3.6','tnf_inhibitor':'4',
        'cxr_fibronodular':'6','head_neck_carcinoma':'16','haemodialysis':'25.3','silicosis':'30','hiv':'50','transplant_immunosuppression':'74'}
    assert set(FACTOR_IDS)==set(expected)
    for factor,row in catalogue_records(catalog).items():
        assert row['candidate_effect']==expected[factor] and row['review']['status']=='not_reviewed'
        assert row['transcription_review']['status']=='reviewed'
        assert row['source']['url']=='https://www.tstin3d.com/thinking3D_2008.pdf'
        assert row['age_applicability']['age_lower']==15 and row['definition_id'] and row['exposure_definition']
        assert row['effect_type']=='candidate_relative_risk_not_causal_HR' and row['profile_hash']==profile_hash(row)


@pytest.mark.parametrize('factor',['recent_infection','close_contact','aids','current_smoking'])
def test_recency_contact_and_unknown_factors_rejected(catalog,factor):
    with pytest.raises(DataValidationError):risk_factor_prevalence_profile([],catalog=catalog,factor_ids=[factor],
        subgroup_id='risk_test',geography_id=GEO,source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    with pytest.raises(DataValidationError):combine_risk_factors([factor],{factor:{'rr':'15','state':'reviewed'}},policy())


@pytest.mark.parametrize('state,p',[('value','0'),('missing',None),('not_applicable',None),('excluded',None),('not_reviewed','.2')])
def test_value_states_distinguish_zero_missing_excluded_and_unreviewed(catalog,resolved,state,p):
    r=prevalence_record(catalog,'diabetes',p,state=state,reason='explicit fixture reason')
    result=allocate(resolved,catalog,[r])
    first=result['rows'][0]['factor_marginals'][0]
    assert first['state']==state
    if state=='value':assert ratio(first['prevalence'])==0 and ratio(result['summary']['any_risk_prevalence']['prevalence_full'])==0
    else:assert first['prevalence'] is None
    if state in ('missing','not_reviewed'):assert ratio(result['summary']['unresolved_population'])==1000


@pytest.mark.parametrize('value',['-0.1','1.01','20%','NaN',True,None])
def test_invalid_prevalences_fail_without_clipping(catalog,value):
    with pytest.raises(DataValidationError):profile(catalog,[prevalence_record(catalog,'diabetes',value)])


def test_definition_mismatch_and_reviewed_mapping(catalog,resolved):
    raw=prevalence_record(catalog,'heavy_smoking','.2',definition='Current smoking')
    original=deepcopy(raw);missing=allocate(resolved,catalog,[raw])
    assert missing['complete_result'] is None and raw==original
    assert missing['rows'][0]['factor_marginals'][0]['state']=='unavailable_due_to_definition_mismatch'
    mapping=definition_mapping(catalog,factor_id='heavy_smoking',supplied_definition='Current smoking',
        source=SOURCE,provenance=PROVENANCE,review=REVIEW,rationale='Synthetic mapping demonstration, not real evidence of equivalence.')
    raw['definition_mapping']=mapping;result=allocate(resolved,catalog,[raw])
    assert ratio(result['summary']['any_risk_prevalence']['prevalence_full'])==Fraction(1,5)
    assert result['rows'][0]['factor_marginals'][0]['warnings']


@pytest.mark.parametrize('factor,definition',[('underweight','BMI below 18.5'),('haemodialysis','All chronic kidney disease'),
    ('head_neck_carcinoma','All cancer'),('transplant_immunosuppression','All immunosuppression')])
def test_incompatible_exposure_definitions_remain_unavailable(catalog,resolved,factor,definition):
    out=allocate(resolved,catalog,[prevalence_record(catalog,factor,'.1',definition=definition)])
    assert out['complete_result'] is None


def test_sex_age_specific_band_and_direct_values(catalog,resolved):
    records=[prevalence_record(catalog,'heavy_smoking',p,sex=s) for s,p in [('female','.01'),('male','.2')]]
    result=allocate(resolved,catalog,records)
    for row in result['rows']:assert ratio(row['any_risk_prevalence'])==Fraction('0.01' if row['sex']=='female' else '0.2')
    assert ratio(result['summary']['any_risk_prevalence']['prevalence_full'])==Fraction('0.105')
    records=[prevalence_record(catalog,'diabetes',p,identifier=f'age:{lo}',age_lower=lo,age_upper_exclusive=hi,policy=CONSTANT)
             for lo,hi,p in [(15,45,'.1'),(45,None,'.3')]]
    result=allocate(resolved,catalog,records)
    assert ratio(result['summary']['any_risk_prevalence']['prevalence_full'])==Fraction('.2')
    exact=prevalence_record(catalog,'diabetes','.1',age_lower=20,age_upper_exclusive=21,policy=DIRECT)
    assert profile(catalog,[exact])['rows'][0]['allocation_policy']==DIRECT


def test_reviewed_pattern_scaled_to_exact_requested_marginal(catalog,resolved):
    pattern={'source':SOURCE,'provenance':PROVENANCE,'review':REVIEW,'rows':[
        {'age_lower':15,'age_upper_exclusive':45,'sex':'both','relative_weight':'1'},
        {'age_lower':45,'age_upper_exclusive':None,'sex':'both','relative_weight':'2'}]}
    record=prevalence_record(catalog,'diabetes','.2',policy=PATTERN,pattern=pattern)
    result=allocate(resolved,catalog,[record]);check=result['pattern_reconciliation'][0]
    assert ratio(check['achieved_prevalence'])==ratio(check['requested_prevalence'])==Fraction(1,5)
    assert ratio(check['scale'])==Fraction(2,15)
    assert {ratio(r['any_risk_prevalence']) for r in result['rows']}=={Fraction(2,15),Fraction(4,15)}
    record['prevalence']='.9'
    with pytest.raises(DataValidationError,match='infeasible'):allocate(resolved,catalog,[record])


def test_missing_pattern_is_not_zero_and_policy_not_implicit(catalog,resolved):
    pattern={'source':SOURCE,'provenance':PROVENANCE,'review':REVIEW,'rows':[
        {'age_lower':15,'age_upper_exclusive':45,'sex':'both','relative_weight':'1'}]}
    with pytest.raises(DataValidationError,match='cover every'):allocate(resolved,catalog,[prevalence_record(catalog,'diabetes','.1',policy=PATTERN,pattern=pattern)])
    record=prevalence_record(catalog,'diabetes','.1');record.pop('allocation_policy')
    with pytest.raises(DataValidationError):profile(catalog,[record])


def test_full_precision_and_decimal_context_independence(catalog,resolved):
    value='0.'+'1234567890'*12
    r=prevalence_record(catalog,'diabetes',value);p=profile(catalog,[r]);assert p['rows'][0]['prevalence']==value.rstrip('0')
    a=allocate(resolved,catalog,[r])
    with localcontext() as context:
        context.prec=6;b=allocate(resolved,catalog,[r])
    assert a['profile_hash']==b['profile_hash']
    assert ratio(a['summary']['any_risk_prevalence']['prevalence_full'])==Fraction(value)


def test_independent_union_factor_counts_and_multimorbidity(catalog,resolved):
    result=allocate(resolved,catalog,[prevalence_record(catalog,k,p) for k,p in [('diabetes','.2'),('heavy_smoking','.1'),('haemodialysis','.05')]])
    s=result['summary']
    assert ratio(s['any_risk_prevalence']['prevalence_full'])==Fraction('.316')
    assert ratio(s['multimorbidity_prevalence']['prevalence_full'])==Fraction('.033')
    assert ratio(s['expected_factor_count']['prevalence_full'])==Fraction('.35')
    dist={ratio(r['effective_rr']):ratio(r['expected_population'])/1000 for r in s['effective_rr_distribution']}
    assert dist=={Fraction(1):Fraction('.684'),Fraction(3):Fraction('.076'),Fraction('3.6'):Fraction('.19'),Fraction('25.3'):Fraction('.05')}
    assert sum(ratio(x['expected_population']) for x in s['factor_count_distribution']['rows'])==1000


def test_highest_rr_ties_and_calculation_only_sensitivities(catalog):
    cell={'age_lower':20,'age_upper_exclusive':21,'sex':'female'}
    effects={k:applicable_effect(v,cell,effect_use='diagnostic_provisional_adult_effects_v1',diagnostic_only=True) for k,v in catalogue_records(catalog).items()}
    present=['diabetes','heavy_smoking','haemodialysis']
    h=combine_risk_factors(present,effects,policy(effect_use='diagnostic_provisional_adult_effects_v1'))
    assert ratio(h['effective_rr'])==Fraction('25.3') and h['responsible_factors']==['haemodialysis']
    assert h['all_factors_present']==sorted(present)
    tied=combine_risk_factors(['heavy_smoking','underweight'],effects,policy(effect_use='diagnostic_provisional_adult_effects_v1'))
    assert tied['tied_highest_contributors']==['heavy_smoking','underweight']
    for method,expected in [(MULTIPLICATIVE,'273.24'),(ADDITIVE,'29.9')]:
        out=combine_risk_factors(present,effects,policy(method,effect_use='diagnostic_provisional_adult_effects_v1'))
        assert ratio(out['effective_rr'])==Fraction(expected) and out['calculation_only_sensitivity'] and out['cap_policy']=='none'
    assert combine_risk_factors(present,effects,policy())['effective_rr'] is None
    with pytest.raises(DataValidationError):policy(cap_policy='none',cap='10')


def test_sensitivities_cannot_be_population_default(catalog,resolved):
    with pytest.raises(DataValidationError,match='calculation-only'):
        allocate(resolved,catalog,[prevalence_record(catalog,'diabetes','.1')],progression=policy(MULTIPLICATIVE))
    with pytest.raises(DataValidationError,match='diagnostic'):
        allocate(resolved,catalog,[prevalence_record(catalog,'diabetes','.1')],diagnostic=False)
    pending=allocate(resolved,catalog,[prevalence_record(catalog,'diabetes','.1')],progression=policy(),diagnostic=False)
    assert pending['complete_result'] is None and not pending['ready_as_input_for_future_operator']


@pytest.mark.parametrize('p,q',[('.3','.2'),('.8','.7'),('0','.2'),('1','.2'),('.3','0'),('.3','1')])
@pytest.mark.parametrize('c',['-1','0','.5','1'])
def test_feasible_joint_endpoints_and_phi(p,q,c):
    out=feasible_association(p,q,c);p,q,c=Fraction(p),Fraction(q),Fraction(c)
    j=ratio(out['joint_prevalence']);lo=max(0,p+q-1);hi=min(p,q)
    assert lo<=j<=hi and sum(ratio(v) for v in out['cells'].values())==1
    if c==-1:assert j==lo
    if c==0:assert j==p*q and (out['achieved_phi']=='0' or out['achieved_phi'] is None)
    if c==1:assert j==hi
    if p in (0,1) or q in (0,1):assert out['achieved_phi'] is None and out['phi_unavailable_reason']=='undefined_degenerate_marginal'
    assert not out['control_is_achieved_correlation']


@pytest.mark.parametrize('c',['-1','0','.37','1'])
def test_TBI_factor_and_recent_remote_marginals_preserved(catalog,resolved,c):
    out=allocate(resolved,catalog,[prevalence_record(catalog,k,p) for k,p in [('diabetes','.2'),('heavy_smoking','.1')]],c=c)
    for row in out['rows']:
        n=ratio(row['population'])
        for state,p in [('any_tbi','.3'),('recent','.1'),('remote_only','.2'),('never_infected','.7'),('prior_remote_plus_recent','.04')]:
            assert sum(ratio(x['tbi_state_counts'][state]) for x in row['compressed_cells'])==n*Fraction(p)
        for f in row['individual_factor_state_counts']:
            p=Fraction('.2' if f['factor_id']=='diabetes' else '.1')
            assert ratio(f['present_population'])==n*p
            for state,total in [('any_tbi','.3'),('recent','.1'),('remote_only','.2')]:
                assert ratio(f['present_tbi_state_counts'][state])+ratio(f['absent_tbi_state_counts'][state])==n*Fraction(total)
    assert ratio(out['summary']['tbi_prevalence']['prevalence_full'])==Fraction('.3')
    assert ratio(out['summary']['any_risk_prevalence']['prevalence_full'])==Fraction('.28')
    assert resolved['profile_hash']==out['resolved_baseline_hash']


def test_association_scope_precedence_and_bad_controls():
    def record(name,scope,c,sex='both',**ages):return {'control_id':name,'scope':scope,'control':c,'sex':sex,
        'source':SOURCE,'provenance':PROVENANCE,'review':REVIEW,**ages}
    controls=[record('sex','sex','.1','male'),record('band','age_band','.2',age_lower=15,age_upper_exclusive=65),
        record('band_sex','age_band','.3','male',age_lower=15,age_upper_exclusive=65),
        record('exact','exact_age_sex','.4',age_lower=20,age_upper_exclusive=21),
        record('exact_sex','exact_age_sex','.5','male',age_lower=20,age_upper_exclusive=21)]
    p=association_policy(global_control='0',overrides=controls,source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    for sex,age,expected in [('male',20,'.5'),('female',20,'.4'),('male',30,'.3'),('female',30,'.2'),('male',80,'.1'),('female',80,'0')]:
        c,_=control_for_cell(p,{'age_lower':age,'age_upper_exclusive':age+1,'sex':sex});assert c==Fraction(expected)
    for c in ('1.1','NaN',True):
        with pytest.raises(DataValidationError):association_policy(global_control=c,source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    with pytest.raises(DataValidationError):association_policy(global_control=0,overrides=[controls[0],{**controls[0],'control_id':'duplicate'}],source=SOURCE,provenance=PROVENANCE,review=REVIEW)


def test_reviewed_joint_distribution_and_incompatible_marginals(catalog,resolved):
    records=[prevalence_record(catalog,'diabetes','.2'),prevalence_record(catalog,'heavy_smoking','.1')]
    joint=[{'age_lower':15,'age_upper_exclusive':None,'sex':'both','patterns':[
        {'factor_ids':[],'probability':'.75'},{'factor_ids':['diabetes'],'probability':'.15'},
        {'factor_ids':['heavy_smoking'],'probability':'.05'},{'factor_ids':['heavy_smoking','diabetes'],'probability':'.05'}]}]
    p=overlap_policy(method=JOINT,joint_rows=joint,source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    out=allocate(resolved,catalog,records,overlap=p)
    assert ratio(out['summary']['any_risk_prevalence']['prevalence_full'])==Fraction('.25')
    assert ratio(out['summary']['multimorbidity_prevalence']['prevalence_full'])==Fraction('.05')
    records[0]['prevalence']='.3'
    with pytest.raises(DataValidationError,match='marginals'):allocate(resolved,catalog,records,overlap=p)


def test_children_no_adult_effects_and_future_reviewed_pediatric_effect():
    demo=toy_demography(1200,children=True)
    rows=[{'record_id':'all','age_lower':0,'age_upper_exclusive':None,'sex':'both',
        'age_allocation_policy':'piecewise_constant_within_band_v1','any_tbi':'.2','recent':'.2','remote_only':'0','never_infected':'.8','prior_remote_plus_recent':'0'}]
    resolved=resolve_baseline_tbi(direct_selection('children_test',demo,rows),demo);catalog=progression_risk_factor_catalog()
    result=allocate(resolved,catalog,[prevalence_record(catalog,'diabetes','.1')])
    for row in result['rows']:
        if row['age_lower']<15:
            assert row['factor_marginals'][0]['state']=='not_applicable_pending_pediatric_review'
            assert row['factor_marginals'][0]['effect']['rr'] is None and row['any_risk_prevalence'] is None
    assert ratio(result['summary']['unresolved_population'])==800
    reviewed=progression_risk_factor_catalog(reviewed_effects=[{'factor_id':'diabetes','age_lower':0,'age_upper_exclusive':15,
        'sex':'both','rr':'2','source':SOURCE,'evidence_vintage':'synthetic_future_pediatric_fixture','review':REVIEW,'provenance':PROVENANCE}])
    effects=applicable_effect(catalogue_records(reviewed)['diabetes'],{'age_lower':3,'age_upper_exclusive':4,'sex':'female'},effect_use='reviewed_effects_only_v1',diagnostic_only=False)
    assert effects['rr']=='2' and effects['state']=='reviewed'


def test_compressed_eleven_factor_bound_and_row_order_hash(catalog,resolved):
    records=[prevalence_record(catalog,k,'.1') for k in FACTOR_IDS]
    a=allocate(resolved,catalog,records);b=allocate(resolved,catalog,list(reversed(records)))
    assert a['profile_hash']==b['profile_hash'] and json.loads(canonical(a))==a
    assert all(len(r['compressed_cells'])<=144 for r in a['rows'])
    assert a['complexity']['no_independent_pattern_enumeration']
    assert not a['ready_as_input_for_future_operator'] and a['stochastic_allocator'].startswith('deferred')


def test_subgroup_weighting_incomplete_and_diagnostics():
    result=diagnostics();examples=result['examples']
    au=examples['Australian_mixture'];s=au['summary']
    assert ratio(s['population'])==10000 and abs(float(ratio(s['tbi_prevalence']['prevalence_full']))-.08)<1e-12
    assert ratio(s['any_risk_prevalence']['prevalence_full'])==Fraction('.187007')
    incomplete=examples['incomplete_incidence_selected_mixture']
    assert ratio(incomplete['summary']['unresolved_population'])==8000 and incomplete['complete_result'] is None
    assert incomplete['baseline_mixture']['subgroups'][0]['tbi_profile']['status']=='unavailable_pending_progression_and_risk_forward_operator'
    assert ratio(examples['current_smoking_definition_mismatch']['summary']['unresolved_population'])==1000
    assert ratio(examples['children_adult_effects_unavailable']['summary']['unresolved_population'])==800
    assert result['diagnostic_only'] and 'NOT AUSTRALIAN' in result['warning']
    for allocation in examples.values():assert allocation['forward_operator'] is None


def test_structural_association_is_not_average_control(catalog):
    demo=toy_demography();records=[{'record_id':str(lo),'age_lower':lo,'age_upper_exclusive':hi,'sex':'both',
        'age_allocation_policy':'piecewise_constant_within_band_v1','any_tbi':p,'recent':p,'remote_only':'0',
        'never_infected':str(1-Fraction(p)),'prior_remote_plus_recent':'0'} for lo,hi,p in [(0,45,'.1'),(45,None,'.8')]]
    # Source accepts decimal probabilities, not Fraction-formatted strings.
    for r in records:r['never_infected']=text_number(Fraction(r['never_infected']))
    base=resolve_baseline_tbi(direct_selection('structural',demo,records),demo)
    risks=[prevalence_record(catalog,'diabetes',p,identifier=str(lo),age_lower=lo,age_upper_exclusive=hi,policy=CONSTANT)
           for lo,hi,p in [(15,45,'.1'),(45,None,'.5')]]
    result=allocate(base,catalog,risks)
    agg=result['summary']['aggregate_association']
    assert ratio(agg['within_stratum_covariance'])==0 and ratio(agg['structural_age_sex_subgroup_covariance'])>0
    assert float(agg['achieved_phi'])>0 and agg['requested_controls_are_not_averaged']


def test_offline_and_no_production_or_disease_outputs(monkeypatch,catalog,resolved):
    def forbidden(*args,**kwargs):raise AssertionError('Network prohibited')
    monkeypatch.setattr(socket,'socket',forbidden);monkeypatch.setattr(socket,'create_connection',forbidden)
    out=allocate(resolved,catalog,[prevalence_record(catalog,'diabetes','.1')])
    assert out['forward_operator'] is None
    forbidden_names={'active_tb_incidence','active_tb_prevalence','cases_averted','DALYs','QALYs','costs'}
    assert not forbidden_names&out.keys()
    assert resolved['profile_hash']==out['resolved_baseline_hash']
    for folder in ('app','ui','general_ui','explicit_tbi_ui','engine/apy'):
        for file in (REPO_ROOT/folder).rglob('*.py'):
            tree=ast.parse(file.read_text(encoding='utf-8-sig'))
            imports=[a.name for n in ast.walk(tree) if isinstance(n,ast.Import) for a in n.names]
            imports += [n.module or '' for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)]
            assert not any(x.startswith('engine.country_tb') for x in imports),file


def test_protected_state_integrity_and_retained_contracts():
    receipt=REPO_ROOT/'abm/output/country_tb_m3a_audit_20261028/starting_state.json'
    if not receipt.is_file():pytest.skip('Live-session protection receipt unavailable.')
    state=json.loads(receipt.read_text());assert state['starting_head']=='5f56715b4f338338abcdbe04faabbf5a46a7471a'
    for f,h in {**state['frozen_files'],**state['retained_2a_2d_files']}.items():assert digest((REPO_ROOT/f).read_bytes())==h,f
    for name,old in state['protected'].items():
        folder=REPO_ROOT.parent/name
        def git(*args):return subprocess.check_output(['git','-C',str(folder),*args],stderr=subprocess.DEVNULL)
        assert git('rev-parse','HEAD').decode().strip()==old['head'] and git('status','--short').decode().splitlines()==old['status']
        assert digest(git('diff','--binary','HEAD'))==old['tracked_diff_sha256'] and digest(git('diff','--cached','--binary'))==old['index_diff_sha256']
        for f,h in old['changed_files'].items():assert digest((folder/f).read_bytes())==h
    refs=dict(x.split(' ',1) for x in subprocess.check_output(['git','for-each-ref','--format=%(refname) %(objectname)']).decode().splitlines())
    allowed={'refs/heads/feature/country-subnational-tb-model-v1','refs/remotes/origin/feature/country-subnational-tb-model-v1'}
    assert {k:v for k,v in refs.items() if k not in allowed}=={k:v for k,v in state['refs'].items() if k not in allowed}


def test_open_tail_retained_unresolved_or_explicit_direct_value(catalog):
    demo=population_demography([
        {'age_lower':20,'age_upper_exclusive':21,'sex':sex,'population':'450'}
        for sex in ('female','male')]+[
        {'age_lower':100,'age_upper_exclusive':None,'sex':sex,'population':'50'}
        for sex in ('female','male')],geography_id=GEO,geography_level='national',population=1000,
        source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    from scripts.preview_country_tb_sources import overall_selection
    base=resolve_baseline_tbi(overall_selection('tail',demo,'.2'),demo)
    records=[prevalence_record(catalog,'diabetes','.1')]
    out=allocate(base,catalog,records)
    assert ratio(out['summary']['population'])==1000
    assert ratio(out['summary']['tbi_prevalence']['unresolved_population'])==100
    assert out['summary']['tbi_prevalence']['prevalence_full'] is None and out['complete_result'] is None
    supplied=[{'record_id':'explicit_including_100plus','age_lower':0,'age_upper_exclusive':None,'sex':'both',
        'age_allocation_policy':'piecewise_constant_within_band_v1','any_tbi':'.3','recent':'.1',
        'remote_only':'.2','never_infected':'.7','prior_remote_plus_recent':'.04'}]
    direct=resolve_baseline_tbi(direct_selection('tail',demo,supplied),demo)
    result=allocate(direct,catalog,records)
    assert ratio(result['summary']['population'])==1000 and ratio(result['summary']['unresolved_population'])==0
    assert all(r['age_upper_exclusive'] is None for r in result['rows'] if r['age_lower']==100)


def test_default_effect_review_and_sensitivity_cannot_leak_into_allocation(catalog,resolved):
    records=[prevalence_record(catalog,'diabetes','.1')]
    default=allocate(resolved,catalog,records,progression=policy())
    assert default['complete_result'] is None
    assert any(r['effective_rr'] is None for r in default['summary']['effective_rr_distribution'])
    for method in (MULTIPLICATIVE,ADDITIVE):
        with pytest.raises(DataValidationError,match='calculation-only'):
            allocate(resolved,catalog,records,progression=policy(method))
    effects={'diabetes':{'rr':'3.6','state':'reviewed'},'heavy_smoking':{'rr':'3','state':'reviewed'}}
    capped=combine_risk_factors(list(effects),effects,policy(MULTIPLICATIVE,cap_policy='explicit_reviewed_cap_v1',cap='5'))
    assert ratio(capped['effective_rr'])==5 and capped['cap_policy']=='explicit_reviewed_cap_v1'


def test_all_contract_families_validate_offline_and_reject_invalid_structure(catalog,monkeypatch):
    from jsonschema import Draft202012Validator
    schema=json.loads((REPO_ROOT/'schemas/country_tb/v1/risk_allocation.schema.json').read_text())
    Draft202012Validator.check_schema(schema);validator=Draft202012Validator(schema)
    def forbidden(*args,**kwargs):raise AssertionError('Schema validation must remain offline')
    monkeypatch.setattr(socket,'socket',forbidden)
    monkeypatch.setattr(socket,'create_connection',forbidden)
    result=diagnostics();samples=[]
    def walk(value):
        if isinstance(value,dict):
            if value.get('implementation_version')==IMPLEMENTATION:samples.append(value)
            for child in value.values():walk(child)
        elif isinstance(value,list):
            for child in value:walk(child)
    walk(result)
    samples += [definition_mapping(catalog,factor_id='heavy_smoking',supplied_definition='Synthetic mapped definition',
        source=SOURCE,provenance=PROVENANCE,review=REVIEW,rationale='Invented transport fixture only'),
        artifact('country_tb.risk_allocation_export_receipt.v1',profile_hash_reference=result['profile_hash'],files={'diagnostics.json':'a'*64})]
    assert len({x['schema_id'] for x in samples})==11
    for sample in samples:
        validator.validate(sample)
        assert sample['profile_hash']==profile_hash(sample)
    bad=deepcopy(catalog);bad['rows'][0]['candidate_effect']='not a number'
    assert not validator.is_valid(bad)
    bad=deepcopy(result['examples']['multimorbidity']);bad['rows'][0]['sex']='both'
    assert not validator.is_valid(bad)
    bad=deepcopy(samples[-1]);bad['schema_id']='country_tb.unknown.v1'
    assert not validator.is_valid(bad)
