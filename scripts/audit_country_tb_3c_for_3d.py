"""Bounded independent3C audit and prospective numerical-resolution diagnosis."""
import argparse
from copy import deepcopy
import json
from math import exp,expm1
from pathlib import Path
import sys
from time import perf_counter

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from engine.country_tb.data_common import seal,canonical,publish_immutable,DataValidationError
from engine.country_tb.preview_incidence_inversion import (
    make_mapping,make_operator,stratum,factor,fit_synthetic,observation,synthetic_gate,
    META,DMETA,CONCENTRATION,frozen_population_coverage,
)
from engine.country_tb.progression_natural_history import audited_adult_candidate,curve_specification,natural_history,curve_hazard
from engine.country_tb.progression_forward import _trajectory
from engine.country_tb.incidence_inversion import inversion_gate,diagnostic_calibration,assess_identifiability,PUBLIC_STATUS
from engine.country_tb.incidence import load_who_snapshot,WHO_2025,national_profile
from engine.country_tb.pressure_shape import shape_from_incidence,apply_backcast
from engine.country_tb.validation_numerics import cohort_reference,one_year_kernel
from engine.country_tb.validation_data import artifact,COMPARISON


def run(*,future_steps=(.1,.05,.025),widths=(.25,.125,.0625,.03125)):
    start=perf_counter();closed=make_mapping();recovery=fit_synthetic(closed,closed.incidence(.017));assert abs(recovery['fitted_scale']-.017)<1e-10
    c=curve_specification([dict(time_years=t,cumulative_risk=p) for t,p in
        [(0,0),(5,0),(10,-expm1(-.25)),(25,-expm1(-.25))]],tail='fixed_remote_hazard_v1',tail_parameters={'hazard_per_year':0},status='diagnostic_synthetic',**META)
    non=make_mapping(rows=[stratum(age=20)],operator=make_operator(natural_history({'0-4':None,'5-14':None,'15+':c},**META)))
    multiple=fit_synthetic(non,800,upper=1);assert multiple['fitted_scale'] is None and len(multiple['identifiability']['all_roots'])==2
    infeasible=assess_identifiability(lambda a:2000*(-expm1(-40*a)),lower=0,upper=.2,target=3000,grid_intervals=64,
        scenario_hash=closed.specification['profile_hash'],**DMETA);assert not infeasible['target_attainable']
    incidence=national_profile(load_who_snapshot(WHO_2025),'VNM');shape=shape_from_incidence(incidence,review=META['review']);history=apply_backcast(shape,start_time='1920')
    rows=[stratum(f'{sex}-{age}',age=age,sex=sex,population=250,factors=[factor('heavy_smoking',.01 if sex=='female' else .2,3)])
        for sex in ('female','male') for age in (30,60)]
    nh=audited_adult_candidate(tail='fixed_remote_hazard_v1',tail_parameters={'hazard_per_year':.00038});curve=nh['group_curves']['15+']
    instant=make_mapping(rows=rows,history=history,operator=make_operator(nh,geography='who:VNM'),geography='who:VNM')
    fictional=fit_synthetic(instant,182);assert abs(fictional['fitted_scale']-.01819928768)<2e-9
    wr=next(r for r in incidence['rows'] if r['year']==2024)
    target=observation(182,geography='who:VNM',empirical_who=True,denominator=wr['denominator'],episode='new_and_recurrent',
        acquisition='unknown',measure='prospective_per_100000_person_years')
    actual=diagnostic_calibration(instant,target,synthetic_gate(instant,target),lower=0,upper=.2,grid_intervals=64,**DMETA)
    assert actual['status']=='blocked_scientific_gate'
    aus=make_mapping(geography='who:AUS',kind='diagnostic_country_population',population_history={});obs=observation(5,geography='who:AUS',acquisition='unknown')
    receipt={'status':'pending','scenario_hash':aus.specification['profile_hash'],'source_id':'M3D_audit','provenance':{'migration_missing':True},'protocol_scope':'reviewed_empirical_configuration'}
    ag=inversion_gate(aus,obs,unresolved_threshold_approval=None,scientific_review=receipt,concentration_policy=CONCENTRATION,**DMETA)
    assert not ag['prefit_allowed']
    try:curve_hazard(audited_adult_candidate()['group_curves']['15+'],30)
    except DataValidationError:tail_blocked=True
    else:tail_blocked=False
    assert tail_blocked
    no_tail=make_mapping(rows=rows,history=history,operator=make_operator(audited_adult_candidate(),geography='who:VNM'),geography='who:VNM').evaluate(.017)
    assert no_tail['incidence'] is None and no_tail['unresolved_population']>0
    rare=make_mapping(rows=[stratum(control=1,factors=[factor('transplant_immunosuppression',.001,74,status='legacy_sensitivity_only')])])
    concentration=rare.concentration(.00001,.01);assert abs(concentration['case_share']-1)<1e-12
    resolution=[]
    for mu in (0.,.08):
        ref=cohort_reference(history,rows,curve,.017,mortality_rate=mu,width=.5)
        independent=cohort_reference(history,rows,curve,.017,mortality_rate=mu,width=.25)
        assert abs(ref['incidence']-independent['incidence'])<1e-8
        for width in widths:
            m=make_mapping(rows=rows,history=history,operator=make_operator(nh,geography='who:VNM',mortality_rate=mu),geography='who:VNM',
                measure='prospective_per_100000_person_years',width=width)
            value=m.incidence(.017)
            resolution.append({'mortality_rate':mu,'kernel_width':width,'known_alpha':.017,'incidence':value,
                'independent_GL8_reference':ref['incidence'],'absolute_error':abs(value-ref['incidence']),
                'scientific_gate':'still_fictional_not_country_approved'})
    future=[]
    for mu in (0.,.08):
        for tau in (.5,4.5,30.):
            for rr in (1.,3.):
                ref=one_year_kernel(curve,tau,rr=rr,mortality_rate=mu,future_pressure=.05)
                for dt in future_steps:
                    op=make_operator(nh,mortality_rate=mu,future_rate=.05);op['step_years']=str(dt);op=seal(op)
                    y=_trajectory(op,{'age':40.,'sex':'female','subgroup_id':'g','effective_rr':rr,'tau':tau},1.,{})
                    future.append({'tau':tau,'rr':rr,'mortality_rate':mu,'step_years':dt,'cases':y['cases'],'reference_cases':ref['cases'],
                        'absolute_case_probability_error':abs(y['cases']-ref['cases']),'person_years_error':abs(y['person_years']-ref['person_years'])})
    return artifact(COMPARISON,direction='independent_3C_numerical_and_scientific_gate_audit',status='audited',
        synthetic_recovery=recovery,multiple_root_rejection=multiple,infeasible_target=infeasible,
        fictional_viet_nam=fictional,fictional_label='Viet_Nam_shaped_fictional_closed_population_not_a_national_estimate',
        actual_WHO_gate=actual,australia_gate=ag,demographic_coverage=frozen_population_coverage(),
        remote_tail_blocked=tail_blocked,no_tail_population_result=no_tail,rare_group_concentration=concentration,
        prospective_response_resolution=resolution,future_time_resolution=future,public_status=PUBLIC_STATUS,
        diagnosis='Piecewise_linear_response_interpolation_and_Strang_time_resolution_are_tested_separately; see_measured_convergence_before_classifying')


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();start=perf_counter();v=run()
    target=publish_immutable(a.output,'audit3c-'+v['profile_hash'][:12],{'audit.json':canonical(v)+b'\n'})
    print(json.dumps({'directory':str(target),'elapsed_seconds':perf_counter()-start,'hash':v['profile_hash'],
        'resolution':v['prospective_response_resolution'],'future_max_error':max(r['absolute_case_probability_error'] for r in v['future_time_resolution'])}))


if __name__=='__main__':main()
