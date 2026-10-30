"""Offline, non-production M3C diagnostics. All fitted examples are fictional.

Run with python -m engine.country_tb.preview_incidence_inversion --output ...
Empirical country targets are validation-only and never acquire synthetic review.
"""
import argparse
from copy import deepcopy
from math import exp,expm1,fsum,log
from pathlib import Path
import json

from engine.country_tb.data_common import canonical,publish_immutable,seal
from engine.country_tb.incidence import WHO_2025,load_who_snapshot,national_profile
from engine.country_tb.pressure_shape import supplied_shape,apply_backcast,shape_from_incidence
from engine.country_tb.progression_natural_history import curve_specification,natural_history,audited_adult_candidate,curve_hazard,curve_integral
from engine.country_tb.progression_policies import effect_application,reinfection_clock,competing_mortality,future_pressure
from engine.country_tb.progression_forward import forward_operator
from engine.country_tb.incidence_mapping import IncidenceMapping
from engine.country_tb.incidence_inversion import (incidence_observation,inversion_gate,diagnostic_calibration,
    assess_identifiability,reconcile_incidence)
from engine.country_tb.progression_evidence_review import evidence_review

SOURCE={'source_id':'synthetic:M3C','citation':'Fictional closed-population mathematical validation; no empirical approval'}
REVIEW={'status':'reviewed','source_id':'synthetic:M3C','rationale':'Mathematical fixture only, not scientific approval of country assumptions.'}
PROVENANCE={'source_id':'synthetic:M3C','scope':'synthetic_test_only'}
META=dict(source=SOURCE,review=REVIEW,provenance=PROVENANCE)
DMETA=dict(source=SOURCE,review=REVIEW,provenance_record=PROVENANCE)
HISTORY_METADATA=dict(exposure_geography='explicit fictional exposure history',arrival_year_policy='closed_no_arrivals',
    imported_active_contribution=0,prior_tb_prevalence=0,recurrence_contribution=0,historical_depletion_policy='known_fictional_baseline_risk_set')
CONCENTRATION={'top_population_share':.01,'maximum_case_share':.5,'source_id':'M3B_top1_half_cases_warning',
    'review':'declared diagnostic warning threshold, not a validated clinical cutoff'}
SYNTHETIC_THRESHOLD={'maximum_share':0,'approved_by':'fictional_test_protocol','source_id':'synthetic:M3C',
    'provenance':{'scope':'unit_test_and_fictional_demonstration_only_not_user_approval'},'scope':'synthetic_test_only'}


def constant_history(geography='synthetic:closed'):
    shape=supplied_shape([{'start':'1900','end':'2025','shape':'1'}],geography_id=geography,
        review=REVIEW,provenance=PROVENANCE,normalization_time='2024')
    return apply_backcast(shape,start_time='1900')


def synthetic_natural_history(rate=.02,*,children=False):
    c=curve_specification([{'time_years':'0','cumulative_risk':'0'},
        {'time_years':'25','cumulative_risk':str(-expm1(-25*rate))}],tail='fixed_remote_hazard_v1',
        tail_parameters={'hazard_per_year':rate},status='diagnostic_synthetic',**META)
    return natural_history({'0-4':c if children else None,'5-14':c if children else None,'15+':c},**META)


def make_operator(history=None,*,geography='synthetic:closed',mortality_rate=0.,future_rate=0.,horizon=1):
    mortality=competing_mortality(geography_id=geography,mode='external_age_sex_survival_v1',
        rows=[{'age_lower':0,'age_upper_exclusive':None,'sex':'both','subgroup_id':None,
            'survival_probability_one_year':exp(-mortality_rate)}],**META)
    future=future_pressure(geography_id=geography,scenario='zero' if future_rate==0 else 'constant_at_2024',
        horizon_years=horizon,**({} if future_rate==0 else {'baseline_rate':future_rate}),**META)
    return forward_operator(natural_history=history or synthetic_natural_history(),effect=effect_application(**META),
        clock=reinfection_clock(**META),mortality=mortality,future=future,
        mode='baseline_cohort_only_v1' if future_rate==0 else 'exogenous_future_infection_v1',
        step_years=.1,horizons=[horizon],diagnostic_only=True,**META)


def stratum(identifier='adult',*,age=40,population=1000,sex='female',factors=None,control=0):
    return dict(stratum_id=identifier,age=age,sex=sex,subgroup_id='g',population=population,
        factors=factors or [],association_control=control,age_support_policy='exact_age_supplied' if age is not None else 'unresolved_open_age')


def factor(identifier,prevalence,rr,*,status='synthetic_mathematical'):
    return dict(factor_id=identifier,prevalence=prevalence,rr=rr,effect_status=status,definition_compatibility='exact',source=SOURCE)


def make_mapping(*,rows=None,history=None,operator=None,geography='synthetic:closed',kind='synthetic_closed_population',
    measure='instantaneous_baseline_per_100000_person_years',horizon=1,width=.5,population_history=None):
    return IncidenceMapping(history=history or constant_history(geography),strata=rows or [stratum()],
        operator=operator or make_operator(geography=geography),measure=measure,horizon_years=horizon,
        population_kind=kind,geography_id=geography,subgroup_id='g',sex_scope='both',age_lower=0,age_upper_exclusive=None,
        population_history=HISTORY_METADATA if population_history is None else population_history,kernel_width=width,**DMETA)


def observation(rate,*,geography='synthetic:closed',target_type='incident_disease',episode='first',
    acquisition='represented_exposure_history',measure='instantaneous_baseline_per_100000_person_years',uncertainty=None,
    empirical_who=False,denominator=1000):
    meta=DMETA if not empirical_who else dict(source={'source_id':WHO_2025,'citation':'WHO Global Tuberculosis Report2025, frozen2024 incidence estimate',
        'url':'https://www.who.int/teams/global-tuberculosis-programme/data'},review={'status':'reviewed','source_id':'M3C_source_transcription_not_observation_approval'},
        provenance_record={'immutable_snapshot':WHO_2025,'no_first_episode_adjustment_approved':True})
    return incidence_observation(geography_id=geography,geography_level='national' if empirical_who else 'custom',subgroup_id='g',year=2024,
        target_type=target_type,rate_per_100000_person_years=rate,count=None,denominator_person_years=denominator,
        observation_period_years=1,measure=measure,age_lower=0,age_upper_exclusive=None,sex='both',uncertainty=uncertainty,
        incident_or_prevalent='incident',episode=episode,acquisition=acquisition,resident_compatibility='compatible',
        ascertainment={'status':'reviewed','model':'WHO estimated incidence, not first-episode ascertainment' if empirical_who else 'complete fictional first-episode ascertainment'},**meta)


def synthetic_gate(mapping,target):
    require_synthetic=mapping.specification['population_kind']=='synthetic_closed_population'
    if not require_synthetic:raise ValueError('Synthetic protocol cannot review country populations.')
    receipt=dict(scenario_hash=mapping.specification['profile_hash'],status='reviewed',source_id='fictional_math_protocol',
        protocol_scope='synthetic_test_only',provenance=PROVENANCE,
        **{k:'reviewed' for k in ('progression_curves','remote_tail','risk_effect_application','mortality',
            'migration_importation_recurrence','prior_tb_and_depletion','future_pressure')})
    return inversion_gate(mapping,target,unresolved_threshold_approval=SYNTHETIC_THRESHOLD,scientific_review=receipt,
        concentration_policy=CONCENTRATION,**DMETA)


def fit_synthetic(mapping,rate,*,upper=.2,grid=64):
    target=observation(rate,geography=mapping.specification['geography_id'],measure=mapping.measure)
    return diagnostic_calibration(mapping,target,synthetic_gate(mapping,target),lower=0,upper=upper,grid_intervals=grid,**DMETA)


def _compact_fit(value):
    p=value.get('profile');rep=value.get('identifiability')
    return dict(status=value['status'],fitted_scale=value['fitted_scale'],any_tbi=None if p is None else p['any_tbi'],
        recent=None if p is None else p['recent'],remote_only=None if p is None else p['remote_only'],
        age_profiles=None if p is None else [{k:r[k] for k in ('stratum_id','age','sex','population','any_tbi','recent','remote_only')} for r in p['rows']],
        predicted_incidence=value['achieved_incidence'],condition_measures=None if rep is None else rep['all_roots'],
        numerical_checks=value.get('numerical_checks'),relative_condition_total_tbi=value.get('relative_condition_total_tbi'),
        failures=value['failures'],scientific_hash=value['profile_hash'])


def frozen_population_coverage():
    from engine.country_tb.demography import load_demography_snapshot
    out={}
    for geo,snapshot in [('AUS','abs-erp-jun2024-07ca11bcb092'),('VNM','un-wpp2024-medium-8028a7a3e163')]:
        value=load_demography_snapshot(snapshot)
        rows=[r for r in value['rows'] if r['geography_id']=='who:'+geo and str(r['year'])=='2024' and r['sex'] in {'female','male'}]
        # ABS geography is iso3166-1:AU in some frozen contexts; preserve actual source key.
        if not rows and geo=='AUS':rows=[r for r in value['rows'] if str(r['year'])=='2024' and r['sex'] in {'female','male'}]
        total=fsum(float(r['population']) for r in rows);child=fsum(float(r['population']) for r in rows if int(r['age_lower'])<15)
        terminal=fsum(float(r['population']) for r in rows if r['age_upper_exclusive'] is None)
        out[geo]={'snapshot_id':snapshot,'population':total,'pediatric_population':child,'pediatric_share':None if not total else child/total,
            'open_age_population':terminal,'open_age_share':None if not total else terminal/total,
            'interpretation':'demographic coverage audit only; no child or open-age incidence imputed'}
    return out


def diagnostics():
    closed=make_mapping();known=.017;closed_target=closed.incidence(known)
    recovery=fit_synthetic(closed,closed_target)
    # A non-monotonic map from the actual operator, not an arbitrary callback.
    knots=[(0,0),(5,0),(10,-expm1(-.25)),(25,-expm1(-.25))]
    curve=curve_specification([dict(time_years=t,cumulative_risk=p) for t,p in knots],
        tail='fixed_remote_hazard_v1',tail_parameters={'hazard_per_year':0},status='diagnostic_synthetic',**META)
    nonhistory=natural_history({'0-4':None,'5-14':None,'15+':curve},**META)
    nonmap=make_mapping(rows=[stratum(age=20)],operator=make_operator(nonhistory))
    nonident=fit_synthetic(nonmap,800,upper=1)
    # WHO temporal proxy is retained. All demographic/risk/curve assumptions below are fictional.
    national=national_profile(load_who_snapshot(WHO_2025),'VNM')
    shape=shape_from_incidence(national,review=REVIEW)
    pressure=apply_backcast(shape,start_time='1920')
    rows=[stratum(f'{sex}-{age}',age=age,sex=sex,population=250,
        factors=[factor('heavy_smoking',.01 if sex=='female' else .2,3,status='legacy_sensitivity_only')])
        for sex in ('female','male') for age in (30,60)]
    adult=audited_adult_candidate(tail='fixed_remote_hazard_v1',tail_parameters={'hazard_per_year':.00038})
    country=make_mapping(rows=rows,history=pressure,operator=make_operator(adult,geography='who:VNM'),geography='who:VNM',kind='diagnostic_country_population',population_history={})
    whorow=next(r for r in national['rows'] if r['year']==2024)
    target=observation(182,geography='who:VNM',episode='new_and_recurrent',acquisition='unknown',
        measure='prospective_per_100000_person_years',empirical_who=True,denominator=whorow['denominator'],
        uncertainty={'lower':whorow['lower'],'upper':whorow['upper'],'kind':'WHO_source_bounds_not_TBI_confidence_interval'})
    pending_review={'status':'pending','scenario_hash':country.specification['profile_hash'],'source_id':'M3C_country_review','provenance':{'no_approval':True},'protocol_scope':'reviewed_empirical_configuration'}
    gate=inversion_gate(country,target,unresolved_threshold_approval=None,scientific_review=pending_review,concentration_policy=CONCENTRATION,**DMETA)
    blocked=diagnostic_calibration(country,target,gate,lower=0,upper=.2,grid_intervals=64,**DMETA)
    grid=[i*.002 for i in range(101)]
    values=[{'alpha':a,'incidence':country.incidence(a)} for a in grid]
    brackets=[[a['alpha'],b['alpha']] for a,b in zip(values,values[1:]) if (a['incidence']-182)*(b['incidence']-182)<0]
    # Counterfactual closed mathematical population: same function, no claim of empirical validity.
    toyrows=deepcopy(rows)
    for r in toyrows:
        for f in r['factors']:f['effect_status']='synthetic_mathematical'
    toy=make_mapping(rows=toyrows,history=pressure,operator=make_operator(adult,geography='who:VNM'),geography='who:VNM')
    counterfactual=fit_synthetic(toy,182)
    scenarios=[]
    for name,tail,pars in [('constant_final','constant_final_segment_hazard_v1',None),
        ('declining_positive_plateau','declining_to_remote_plateau_v1',{'plateau_hazard_per_year':.0001,'decay_rate_per_year':.1}),
        ('declining_zero_plateau','declining_to_remote_plateau_v1',{'plateau_hazard_per_year':0,'decay_rate_per_year':.1}),
        ('fixed_supplied','fixed_remote_hazard_v1',{'hazard_per_year':.00038}),('no_estimate','no_extrapolation_v1',None)]:
        nh=audited_adult_candidate(tail=tail,tail_parameters=pars)
        m=make_mapping(rows=toyrows,history=pressure,operator=make_operator(nh,geography='who:VNM'),geography='who:VNM')
        c=nh['group_curves']['15+'];rfinal=float(c['segments'][-1]['hazard_per_year'])
        limit=None if tail=='no_extrapolation_v1' else -expm1(-float(c['knots'][-1]['cumulative_hazard'])-rfinal/.1) if name=='declining_zero_plateau' else 1.
        scenarios.append(dict(name=name,current_hazard_just_after25=None if tail=='no_extrapolation_v1' else curve_hazard(c,25),
            hazard_at50=None if tail=='no_extrapolation_v1' else curve_hazard(c,50),
            net_risk_in_infinite_time_without_death=limit,lifetime_risk_with_mortality=None,
            lifetime_warning='infinite-time net limit is not remaining-lifetime risk; no approved life table',
            incidence_at_known_scale=m.incidence(known),older_age_rows=m.evaluate(known)['rows'],
            fit=_compact_fit(fit_synthetic(m,182)) if not m.tail_gaps else {'status':'blocked_remote_tail','fitted_scale':None}))
    for name,multiplier,control,mortality in [('progression_half',.5,0,0),('progression_double',2,0,0),
        ('association_minus1',1,-1,0),('association_plus1',1,1,0),('no_risk_factors',1,0,0),('synthetic_mortality_008',1,0,.08)]:
        rs=deepcopy(toyrows)
        for r in rs:
            r['association_control']=control
            if name=='no_risk_factors':r['factors']=[]
        nh=deepcopy(adult)
        if multiplier!=1:
            knots=[{'time_years':k['time_years'],'cumulative_risk':-expm1(-multiplier*float(k['cumulative_hazard']))} for k in adult['group_curves']['15+']['knots']]
            curve=curve_specification(knots,tail='fixed_remote_hazard_v1',tail_parameters={'hazard_per_year':.00038*multiplier},status='diagnostic_synthetic',**META)
            nh=natural_history({'0-4':None,'5-14':None,'15+':curve},**META)
        m=make_mapping(rows=rs,history=pressure,operator=make_operator(nh,geography='who:VNM',mortality_rate=mortality),geography='who:VNM')
        scenarios.append({'name':name,'fit':_compact_fit(fit_synthetic(m,182)),
            'interpretation':'synthetic counterfactual only; instantaneous mortality sensitivity is exactly zero'})
    for name,rate in [('source_lower',next(r['lower'] for r in national['rows'] if r['year']==2024)),
        ('source_upper',next(r['upper'] for r in national['rows'] if r['year']==2024)),('imported_recurrent_adjustment_20percent',145.6)]:
        scenarios.append({'name':name,'target':rate,'fit':_compact_fit(fit_synthetic(toy,float(rate))),
            'interpretation':'numerical hypothetical target; no empirical first-episode adjustment approved'})
    for name,hist in [('bounded_pre2000',apply_backcast(shape,start_time='1920',method='bounded_backward_trend',
        parameters={'slope_per_year':.005,'minimum_shape':.001,'maximum_shape':100,'review':REVIEW})),
        ('WHO_lower_shape',apply_backcast(shape_from_incidence(national,review=REVIEW,source_bound='lower'),start_time='1920')),
        ('WHO_upper_shape',apply_backcast(shape_from_incidence(national,review=REVIEW,source_bound='upper'),start_time='1920'))]:
        m=make_mapping(rows=toyrows,history=hist,operator=make_operator(adult,geography='who:VNM'),geography='who:VNM')
        scenarios.append({'name':name,'fit':_compact_fit(fit_synthetic(m,182)),
            'interpretation':'source/history scenario, not a TBI confidence interval'})
    prospective=make_mapping(rows=toyrows,history=pressure,operator=make_operator(adult,geography='who:VNM'),
        geography='who:VNM',measure='prospective_per_100000_person_years',width=.25)
    prospective_target=prospective.incidence(known)
    for name,mu,pressure_rate in [('prospective_known_scale',0,0),('prospective_mortality_008',.08,0),('prospective_future_pressure_005',0,.05)]:
        m=prospective if name=='prospective_known_scale' else make_mapping(rows=toyrows,history=pressure,
            operator=make_operator(adult,geography='who:VNM',mortality_rate=mu,future_rate=pressure_rate),
            geography='who:VNM',measure='prospective_per_100000_person_years',width=.25)
        scenarios.append({'name':name,'target':prospective_target,'fit':_compact_fit(fit_synthetic(m,prospective_target)),
            'interpretation':'fictional prospective target; mortality/future scenarios fixed independently of historical alpha'})
    extreme_rows=[stratum(factors=[factor('diabetes',.1,3.6,status='legacy_sensitivity_only'),
        factor('heavy_smoking',.2,3,status='legacy_sensitivity_only'),factor('haemodialysis',.01,25.3,status='legacy_sensitivity_only'),
        factor('head_neck_carcinoma',.005,16,status='legacy_sensitivity_only'),factor('silicosis',.005,30,status='legacy_sensitivity_only'),
        factor('hiv',.005,50,status='legacy_sensitivity_only'),factor('transplant_immunosuppression',.001,74,status='legacy_sensitivity_only')])]
    # Concentration is a forward diagnostic, not a failed-gate inversion.
    concentration=[]
    for omit in [None,'transplant_immunosuppression','hiv','silicosis','haemodialysis','head_neck_carcinoma']:
        rs=deepcopy(extreme_rows)
        if omit:rs[0]['factors']=[f for f in rs[0]['factors'] if f['factor_id']!=omit]
        m=make_mapping(rows=rs,operator=make_operator(adult));v=m.evaluate(known);rrs={}
        for part in v['risk_parts']:
            key=str(part['effective_rr']);a=rrs.setdefault(key,{'population':0.,'instantaneous_expected_cases_per_year_at_baseline':0.})
            a['population']+=part['population'];a['instantaneous_expected_cases_per_year_at_baseline']+=part['numerator']
        fit_rows=deepcopy(rs)
        for r in fit_rows:
            for f in r['factors']:f['effect_status']='synthetic_mathematical'
        fit_map=make_mapping(rows=fit_rows,operator=make_operator(adult))
        concentration.append({'omitted_factor':omit,'incidence':v['incidence'],'rr_distribution':rrs,
            'top_shares':[m.concentration(known,s) for s in (.01,.05,.1)],
            'fictional_target_328_fit':_compact_fit(fit_synthetic(fit_map,328.05131036202033)),
            'warning':'legacy effects unreviewed; factor omission is an analysis sensitivity, not evidence of absence'})
    return seal(dict(schema_id='country_tb.m3c_diagnostics.internal.v1',production_available=False,
        evidence=evidence_review(),synthetic_recovery={'known_alpha':known,'result':recovery},
        nonidentifiable=nonident,viet_nam={'label':'DIAGNOSTIC ASSUMPTIONS; NOT A NATIONAL ESTIMATE',
            'actual_target_gate':blocked,'mapping':values,'target_182_in_scanned_range':min(v['incidence'] for v in values)<=182<=max(v['incidence'] for v in values),
            'unadopted_crossing_brackets':brackets,'fictional_closed_counterfactual':counterfactual,
            'validation_only':reconcile_incidence(country.evaluate(.005),target,**DMETA)},
        practical_scenarios=scenarios,extreme_factor_exclusions=concentration,
        population_coverage=frozen_population_coverage(),australia={'status':'blocked_migration_importation_prior_TB_and_natural_history',
            'preferred_pathway':'user_TBI_and_subgroup_forward_validation','empirical_subgroup_values_invented':False},
        unresolved_threshold={'actual_configuration_approval':'pending_user_approval','candidates':[0,.001,.01],
            'synthetic_fixtures':'zero unresolved, fictional protocol only'},
        dynamic_boundary='external historical pressure is not endogenous transmission; diagnostic initialization does not replace a dynamic model'))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True);args=parser.parse_args()
    result=diagnostics();folder=publish_immutable(Path(args.output),'incidence-inversion-'+result['profile_hash'][:12],
        {'diagnostics.json':canonical(result)+b'\n'})
    print(json.dumps({'directory':str(folder),'scientific_hash':result['profile_hash'],
        'known_alpha_recovered':result['synthetic_recovery']['result']['fitted_scale'],
        'country_gate_status':result['viet_nam']['actual_target_gate']['status']}))


if __name__=='__main__':main()
