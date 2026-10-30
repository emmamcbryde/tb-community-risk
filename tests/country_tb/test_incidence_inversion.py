"""Bounded independent mathematical and guardrail checks for M3C."""
from copy import deepcopy
from decimal import Decimal,localcontext
from math import exp,expm1,log
from pathlib import Path
import ast
import hashlib
import json
import socket

import pytest

from engine.country_tb.data_common import DataValidationError,canonical,profile_hash,REPO_ROOT,seal
from engine.country_tb.incidence_inversion import *
from engine.country_tb.incidence_mapping import IncidenceMapping,_mean_offset,_association_cells
from engine.country_tb.progression_evidence_review import evidence_review
from engine.country_tb.preview_incidence_inversion import (
    make_mapping,make_operator,synthetic_natural_history,constant_history,stratum,factor,
    observation,synthetic_gate,fit_synthetic,META,DMETA,SYNTHETIC_THRESHOLD,CONCENTRATION,
)
from engine.country_tb.progression_natural_history import (
    risk_to_hazard,audited_adult_candidate,curve_hazard,curve_integral,curve_specification,natural_history,
)


def audit(f,target,lo=0,hi=1,n=64):
    return assess_identifiability(f,lower=lo,upper=hi,target=target,grid_intervals=n,scenario_hash='0'*64,**DMETA)


@pytest.fixture(scope='module')
def mapping():return make_mapping()


@pytest.mark.parametrize('risk',[0,.038,.066,.072,.079,.999999999])
def test_independent_hazard_transformation(risk):
    with localcontext() as ctx:
        ctx.prec=70
        expected=float(-(1-Decimal(str(risk))).ln())
    assert float(risk_to_hazard(risk))==pytest.approx(expected,abs=2e-15)


def test_independent_3b_segments_and_extreme_warning():
    curve=audited_adult_candidate()['group_curves']['15+']
    assert [float(s['hazard_per_year']) for s in curve['segments']]==pytest.approx(
        [.0387408283,.0073845031,.0012889411,.0005047798],abs=5e-11)
    assert -expm1(-74*float(risk_to_hazard('.038')))==pytest.approx(.9431205575309726)
    review=evidence_review()
    assert review['adult_curve_status']=='sensitivity_only'
    assert all(r['proposed_central_use']=='excluded_pending_factor_specific_review' for r in review['rows'])
    assert len(review['rows'])==11


@pytest.mark.parametrize('tau',[25,26,80])
def test_no_tail_is_not_zero(tau):
    c=audited_adult_candidate()['group_curves']['15+']
    with pytest.raises(DataValidationError,match='external_required'):curve_hazard(c,tau)
    with pytest.raises(DataValidationError,match='external_required'):curve_integral(c,24,tau+1)


def test_no_tail_population_remains_visible():
    m=make_mapping(operator=make_operator(audited_adult_candidate()))
    v=m.evaluate(.017)
    assert v['incidence'] is None and 0<v['unresolved_population']<1000
    assert m.tail_gaps and v['covered_numerator']>0
    assert not synthetic_gate(m,observation(1))['prefit_allowed']


@pytest.mark.parametrize('age',[0,3,10,14])
def test_children_unresolved(age):
    m=make_mapping(rows=[stratum(age=age)])
    v=m.evaluate(.01)
    assert v['incidence'] is None and v['unresolved_population']==1000
    g=synthetic_gate(m,observation(10))
    assert not g['checks']['all_material_age_curves'] and not g['prefit_allowed']


def test_100plus_unresolved_retained():
    m=make_mapping(rows=[stratum('adult',population=999),stratum('100plus',age=None,population=1)])
    v=m.evaluate(.01)
    assert v['population']==1000 and v['unresolved_population']==1 and v['incidence'] is None


def test_approved_small_exclusion_does_not_impute_unknown_forward_kernel():
    m=make_mapping(rows=[stratum('adult',population=999),stratum('100plus',age=None,population=1)])
    obs=observation(10);g=synthetic_gate(m,obs)
    approval=dict(SYNTHETIC_THRESHOLD,maximum_share=.001)
    g=inversion_gate(m,obs,unresolved_threshold_approval=approval,scientific_review=g['scientific_review'],concentration_policy=CONCENTRATION,**DMETA)
    assert g['checks']['unresolved_population'] and not g['checks']['forward_prediction_defined_for_target']
    result=diagnostic_calibration(m,obs,g,lower=0,upper=.2,**DMETA)
    assert result['status']=='blocked_scientific_gate' and result['fitted_scale'] is None


@pytest.mark.parametrize('kind,expected',[('notifications','validation_only'),('baseline_prevalent_disease','incompatible'),
    ('screen_detected_prevalent_disease','incompatible'),('imported_active_disease','incompatible'),('mixed_or_unclear','insufficiently_specified')])
def test_observation_types_not_silent_incidence(kind,expected):
    assert observation(10,target_type=kind)['compatibility']==expected


@pytest.mark.parametrize('episode,acquisition',[('new_and_recurrent','unknown'),('recurrent','represented_exposure_history'),('first','imported'),('first','mixed')])
def test_recurrent_imported_targets_not_compatible(episode,acquisition):
    assert observation(10,episode=episode,acquisition=acquisition)['compatibility']!='compatible'


def test_rate_count_consistency_and_required_person_time():
    obs=observation(200);args={k:v for k,v in obs.items() if k in {
        'geography_id','geography_level','subgroup_id','year','target_type','rate_per_100000_person_years',
        'count','denominator_person_years','observation_period_years','measure','age_lower','age_upper_exclusive',
        'sex','uncertainty','incident_or_prevalent','episode','acquisition','resident_compatibility','ascertainment','source','review'}}
    args['provenance_record']=obs['provenance'];args['count']=2
    assert incidence_observation(**args)['calibration_target']==200
    with pytest.raises(DataValidationError):incidence_observation(**dict(args,count=3))
    with pytest.raises(DataValidationError):incidence_observation(**dict(args,denominator_person_years=None))
    adj=dict(observation_model='reviewed notification ascertainment and episode/importation model',adjusted_rate=150,
        source=DMETA['source'],review=DMETA['review'],provenance=DMETA['provenance_record'],
        resolves=['incident_first_episode','represented_exposure_history','resident_population','ascertainment'])
    assert incidence_observation(**dict(args,target_type='notifications',adjustments=adj))['compatibility']=='compatible_after_reviewed_adjustment'
    with pytest.raises(DataValidationError):incidence_observation(**dict(args,target_type='baseline_prevalent_disease',adjustments=adj))


def test_threshold_requires_explicit_approval_and_country_cannot_use_test_review(mapping):
    target=observation(10);g=synthetic_gate(mapping,target)
    g2=inversion_gate(mapping,target,unresolved_threshold_approval=None,scientific_review=g['scientific_review'],concentration_policy=CONCENTRATION,**DMETA)
    assert not g2['prefit_allowed'] and 'threshold_approved' in g2['failed_criteria']
    country=make_mapping(kind='diagnostic_country_population')
    g3=inversion_gate(country,target,unresolved_threshold_approval=SYNTHETIC_THRESHOLD,
        scientific_review=g['scientific_review'],concentration_policy=CONCENTRATION,**DMETA)
    assert not g3['prefit_allowed'] and not g3['checks']['review_scope']


def test_no_fit_is_attempted_after_failed_gate(mapping,monkeypatch):
    target=observation(10,episode='new_and_recurrent');gate=synthetic_gate(mapping,target)
    monkeypatch.setattr(mapping,'incidence',lambda a:pytest.fail('Forward called after scientific gate failed'))
    result=diagnostic_calibration(mapping,target,gate,lower=0,upper=.2,**DMETA)
    assert result['fitted_scale'] is None and result['identifiability'] is None
    assert result['public_status']==PUBLIC_STATUS


def test_mortality_and_migration_missing_block(mapping):
    m=make_mapping(population_history={});g=synthetic_gate(m,observation(10))
    assert not g['checks']['migration_metadata_addressed']
    g=synthetic_gate(mapping,observation(10));receipt=deepcopy(g['scientific_review']);receipt['mortality']='pending'
    gate=inversion_gate(mapping,observation(10),unresolved_threshold_approval=SYNTHETIC_THRESHOLD,
        scientific_review=receipt,concentration_policy=CONCENTRATION,**DMETA)
    assert not gate['checks']['mortality_reviewed']


def test_legacy_rrs_block_even_when_finite():
    m=make_mapping(rows=[stratum(factors=[factor('transplant_immunosuppression',.001,74,status='legacy_sensitivity_only')])])
    assert not synthetic_gate(m,observation(10))['checks']['risk_effects_reviewed']
    assert m.evaluate(.01)['incidence'] is not None


@pytest.mark.parametrize('alpha',[0,1e-12,.001,.017,.2,1000])
def test_last_event_density_matches_poisson_analytic(mapping,alpha):
    v=mapping.evaluate(alpha);p=-expm1(-alpha*40)
    assert v['any_tbi']==pytest.approx(p,abs=1e-14)
    assert v['incidence']==pytest.approx(2000*p,abs=2e-10)
    r=v['rows'][0]
    assert r['recent']+r['remote_only']==pytest.approx(p,abs=1e-14)
    assert r['prior_remote_plus_recent']<=r['recent']+1e-15
    assert sum(x['population'] for x in v['risk_parts'])==pytest.approx(1000,abs=1e-10)


@pytest.mark.parametrize('control',[-1,0,.4,1])
def test_factor_and_tbi_marginals_and_recency_unchanged(control):
    factors=[factor('diabetes',.2,3.6),factor('heavy_smoking',.3,3)]
    m=make_mapping(rows=[stratum(factors=factors,control=control)]);v=m.evaluate(.02);r=v['rows'][0]
    assert r['any_risk']==pytest.approx(.44)
    assert r['any_tbi']==pytest.approx(-expm1(-.8))
    assert r['recent']==pytest.approx(-expm1(-.1))
    for f in factors:
        population=sum(x['population'] for x in v['risk_parts'] if f['factor_id'] in x['factor_ids'])
        assert population/1000==pytest.approx(f['prevalence'])
    assert max(x['effective_rr'] for x in v['risk_parts'])==3.6
    assert r['prior_remote_plus_recent']<=r['recent']


@pytest.mark.parametrize('p,q',[(1e-14,1e-14),(1-1e-7,1-1e-8),(.2,.3),(0,1),(1,0),(1,1)])
@pytest.mark.parametrize('control',[-1,0,.6,1])
def test_nearly_degenerate_association_no_negative_mass_or_clipping(p,q,control):
    cells=_association_cells(p,q,control)
    assert min(cells)>=0 and sum(cells)==pytest.approx(1,abs=1e-15)
    assert cells[0]+cells[1]==pytest.approx(p,abs=2e-16)
    assert cells[0]+cells[2]==pytest.approx(q,abs=2e-16)


def test_tiny_positive_risk_prevalence_not_lost_as_zero():
    m=make_mapping(rows=[stratum(factors=[factor('diabetes',1e-30,3.6)])])
    assert m.evaluate(.017)['rows'][0]['any_risk']==pytest.approx(1e-30,rel=1e-14)


def test_exact_synthetic_scale_recovery(mapping):
    result=fit_synthetic(mapping,mapping.incidence(.017))
    assert result['status']=='diagnostic_only_not_approved_for_production'
    assert result['fitted_scale']==pytest.approx(.017,abs=1e-10)
    assert result['production_available'] is False
    assert result['numerical_checks']['bracket_invariant']


@pytest.mark.parametrize('target,status',[(0,None),(1000,None),(3000,'target_outside_attainable_range')])
def test_boundaries_and_attainable_range(target,status):
    r=audit(lambda a:2000*(-expm1(-40*a)),target)
    assert r['failure_status']==status
    if target==0:assert r['all_roots'][0]['boundary'] and r['all_roots'][0]['alpha']==0


def test_multiple_root_actual_natural_history_operator():
    curve=curve_specification([dict(time_years=t,cumulative_risk=p) for t,p in
        [(0,0),(5,0),(10,-expm1(-.25)),(25,-expm1(-.25))]],tail='fixed_remote_hazard_v1',tail_parameters={'hazard_per_year':0},status='diagnostic_synthetic',**META)
    nh=natural_history({'0-4':None,'5-14':None,'15+':curve},**META)
    m=make_mapping(rows=[stratum(age=20)],operator=make_operator(nh))
    result=fit_synthetic(m,800,upper=1)
    r=result['identifiability']
    assert r['monotonicity']=='non_monotonic' and result['fitted_scale'] is None
    assert [p['alpha'] for p in r['all_roots']]==pytest.approx([-log(.8)/5,-log(.2)/5],abs=1e-9)


@pytest.mark.parametrize('function,target',[(lambda a:0,0),(lambda a:100,100),(lambda a:2000*(-expm1(-100*a)),2000)])
def test_flat_saturation_rejected(function,target):
    r=audit(function,target,lo=.5,hi=1)
    assert not r['identifiable_over_search_interval']


def test_stationary_tangent_root_is_nonidentifiable():
    r=audit(lambda a:1000*(a-.357)**2,0)
    assert not r['identifiable_over_search_interval']
    assert any(abs(x['alpha']-.357)<1e-6 for x in r['all_roots'])


def test_discontinuous_callback_fails_convergence():
    r=audit(lambda a:0 if a<.333 else 100,50)
    assert not r['numerical_convergence'] and not r['identifiable_over_search_interval']


def test_derivative_and_grid_bracket_invariance():
    f=lambda a:2000*(-expm1(-40*a));t=f(.017)
    a=audit(f,t,hi=.2,n=64);b=audit(f,t,hi=.1,n=128)
    r=a['all_roots'][0]
    assert r['derivative_incidence_per_alpha']==pytest.approx(80000*exp(-40*.017),rel=1e-7)
    assert r['alpha']==pytest.approx(b['all_roots'][0]['alpha'],abs=1e-9)
    assert a['grid_refinement']['coarse_crossings']==a['grid_refinement']['refined_crossings']==1


@pytest.mark.parametrize('rate',[1e-12,.1,100,10000])
def test_conditional_time_small_large_stable(rate):
    assert 0<=_mean_offset(rate,1)<=1


def test_prospective_depletion_mortality_and_future_first_infection():
    rate=.02;alpha=.017;p=-expm1(-40*alpha);mu=.08
    m=make_mapping(operator=make_operator(mortality_rate=mu),measure='prospective_per_100000_person_years')
    v=m.evaluate(alpha)
    cases=1000*p*rate/(rate+mu)*(-expm1(-(rate+mu)))
    py=1000*(p*(-expm1(-(rate+mu)))/(rate+mu)+(1-p)*(-expm1(-mu))/mu)
    assert v['covered_numerator']==pytest.approx(cases,abs=1e-10)
    assert v['incidence']==pytest.approx(cases/py*100000,abs=1e-8)
    f=make_mapping(rows=[stratum(age=20)],operator=make_operator(future_rate=.05),measure='prospective_per_100000_person_years',width=1)
    assert f.incidence(0)>0
    assert m.incidence(0)==0


def test_future_numerical_refinement_halves_time_and_infection_age_grids():
    m=make_mapping(rows=[stratum(age=20)],operator=make_operator(future_rate=.05),measure='prospective_per_100000_person_years',width=1)
    fine=m.refined()
    assert fine.width==m.width/2
    assert float(fine.operator['step_years'])==float(m.operator['step_years'])/2
    assert fine.specification['profile_hash']!=m.specification['profile_hash']


def test_user_tbi_preserved_validation_only(mapping):
    p=mapping.evaluate(.012);before=canonical(p)
    r=reconcile_incidence(p,observation(1),**DMETA)
    assert r['status']=='validation_only_not_calibration' and r['difference']>0
    assert canonical(p)==before and canonical(r['supplied_profile'])==before


def test_row_order_and_scientific_hashes():
    rows=[stratum('f',sex='female',factors=[factor('diabetes',.1,3.6),factor('heavy_smoking',.2,3)]),stratum('m',sex='male')]
    a=make_mapping(rows=rows)
    rows.reverse();rows[1]['factors'].reverse();b=make_mapping(rows=rows)
    assert a.specification['profile_hash']==b.specification['profile_hash']
    assert canonical(a.evaluate(.01))==canonical(b.evaluate(.01))
    assert evidence_review()['profile_hash']==evidence_review()['profile_hash']


@pytest.mark.parametrize('bad',[-.01,1.01,True,float('nan')])
def test_invalid_prevalence_no_clipping(bad):
    with pytest.raises(DataValidationError):make_mapping(rows=[stratum(factors=[factor('diabetes',bad,3.6)])])


def test_definition_mismatch_rejected():
    f=factor('heavy_smoking',.2,3);f['definition_compatibility']='current_smoking_not_heavy'
    with pytest.raises(DataValidationError):make_mapping(rows=[stratum(factors=[f])])


def test_concentration_rare_high_group_and_threshold():
    nh=audited_adult_candidate(tail='fixed_remote_hazard_v1',tail_parameters={'hazard_per_year':.00038})
    m=make_mapping(rows=[stratum(factors=[factor('transplant_immunosuppression',.001,74)])],operator=make_operator(nh))
    assert 0<m.concentration(.01,.01)['case_share']<=1
    assert CONCENTRATION['top_population_share']==.01 and CONCENTRATION['maximum_case_share']==.5


def test_compressed_scalar_preserves_exact_rr_union_mapping():
    fs=[factor('diabetes',.2,3.6),factor('heavy_smoking',.3,3),factor('haemodialysis',.01,25.3)]
    m=make_mapping(rows=[stratum(factors=fs,control=.6)])
    for alpha in (0,.001,.017,.5):
        assert m.incidence(alpha)==pytest.approx(m.evaluate(alpha)['incidence'],rel=1e-13)


def test_unreviewed_rare_group_can_drive_all_cases_and_is_blocked():
    m=make_mapping(rows=[stratum(control=1,factors=[factor('transplant_immunosuppression',.001,74,status='legacy_sensitivity_only')])])
    # TBI prevalence<0.001 and maximal feasible union association puts all TBI
    # in this0.1% group. The normalization control is not an achieved phi.
    c=m.concentration(.00001,.01)
    assert c['case_share']==pytest.approx(1)
    assert c['case_share']>CONCENTRATION['maximum_case_share']
    assert not synthetic_gate(m,observation(10))['checks']['risk_effects_reviewed']


def test_viet_nam_who_target_blocked_and_fictional_shape_fit_separate():
    from engine.country_tb.incidence import load_who_snapshot,WHO_2025,national_profile
    from engine.country_tb.pressure_shape import shape_from_incidence,apply_backcast
    national=national_profile(load_who_snapshot(WHO_2025),'VNM')
    history=apply_backcast(shape_from_incidence(national,review=DMETA['review']),start_time='1920')
    nh=audited_adult_candidate(tail='fixed_remote_hazard_v1',tail_parameters={'hazard_per_year':.00038})
    rows=[stratum(f'{s}-{a}',age=a,population=250,sex=s,factors=[factor('heavy_smoking',.01 if s=='female' else .2,3)])
        for s in ('female','male') for a in (30,60)]
    m=make_mapping(rows=rows,history=history,operator=make_operator(nh,geography='who:VNM'),geography='who:VNM')
    fictional=fit_synthetic(m,182)
    assert fictional['fitted_scale']==pytest.approx(.01819928768,abs=2e-9)
    assert fictional['profile']['any_tbi']==pytest.approx(.653314219,abs=1e-8)
    wr=next(r for r in national['rows'] if r['year']==2024)
    obs=observation(182,geography='who:VNM',empirical_who=True,denominator=wr['denominator'],episode='new_and_recurrent',
        acquisition='unknown',measure='prospective_per_100000_person_years')
    assert obs['source']['source_id']==WHO_2025 and obs['compatibility']=='validation_only'
    result=diagnostic_calibration(m,obs,synthetic_gate(m,obs),lower=0,upper=.2,**DMETA)
    assert result['status']=='blocked_scientific_gate' and result['fitted_scale'] is None


def test_australia_and_empirical_mortality_snapshot_gate():
    m=make_mapping(geography='who:AUS',kind='diagnostic_country_population',population_history={})
    obs=observation(10,geography='who:AUS',acquisition='unknown')
    receipt={'status':'pending','scenario_hash':m.specification['profile_hash'],'source_id':'independent_review',
        'provenance':{'missing_migration_and_importation':True},'protocol_scope':'reviewed_empirical_configuration'}
    g=inversion_gate(m,obs,unresolved_threshold_approval=None,scientific_review=receipt,concentration_policy=CONCENTRATION,**DMETA)
    assert not g['checks']['mortality_snapshot_reviewed'] and not g['checks']['migration_metadata_addressed']
    assert not g['prefit_allowed'] and g['production_available'] is False


def test_review_label_alone_does_not_approve_factor():
    with pytest.raises(DataValidationError,match='conditional-TBI'):
        make_mapping(rows=[stratum(factors=[factor('hiv',.01,50,status='reviewed')])])


def test_gate_checks_cannot_be_resealed_to_force_a_fit(mapping):
    obs=observation(10,target_type='notifications');gate=synthetic_gate(mapping,obs)
    gate['checks']['observation_compatible']=True;gate['failed_criteria']=[];gate['prefit_allowed']=True;gate=seal(gate)
    with pytest.raises(DataValidationError,match='must reproduce'):
        diagnostic_calibration(mapping,obs,gate,lower=0,upper=.2,**DMETA)


def test_six_contract_schema_exports(mapping):
    import jsonschema
    schema=json.loads((REPO_ROOT/'schemas/country_tb/v1/incidence_inversion.schema.json').read_text())
    result=fit_synthetic(mapping,mapping.incidence(.017))
    values=[result['observation'],result['gate'],result,result['identifiability'],evidence_review(),
        reconcile_incidence(mapping.evaluate(.017),observation(182),**DMETA)]
    assert len({v['schema_id'] for v in values})==6
    for value in values:
        jsonschema.Draft202012Validator(schema).validate(value)
        assert value['profile_hash']==profile_hash(value)


def test_offline_and_public_dispatch_unmodified(monkeypatch):
    monkeypatch.setattr(socket,'socket',lambda *a,**k:pytest.fail('Network access in pure diagnostics'))
    assert make_mapping().incidence(.01)>0
    text=(REPO_ROOT/'engine/country_tb/baseline_tbi_source.py').read_text()
    assert PUBLIC_STATUS in (REPO_ROOT/'engine/country_tb/active_tb_targets.py').read_text() and 'incidence_inversion' not in text
    for filename in ('incidence_inversion.py','incidence_mapping.py','progression_evidence_review.py','preview_incidence_inversion.py'):
        tree=ast.parse((REPO_ROOT/'engine/country_tb'/filename).read_text())
        imports=[n.module or '' for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)]
        imports += [v.name for n in ast.walk(tree) if isinstance(n,ast.Import) for v in n.names]
        assert not any(any(k in i.lower() for k in ('streamlit','starsim','requests','urllib','matlab')) for i in imports)


def test_frozen_2a_to_3b_source_bytes_unchanged():
    audit_file=REPO_ROOT/'abm/output/country_tb_m3c_audit_20261030/starting_state.json'
    if not audit_file.exists():pytest.skip('Session recovery receipt is local; tracked protection manifest is checked separately.')
    receipt=json.loads(audit_file.read_text())
    for name,expected in receipt['retained_2a_3b_files'].items():
        assert hashlib.sha256((REPO_ROOT/name).read_bytes()).hexdigest()==expected
