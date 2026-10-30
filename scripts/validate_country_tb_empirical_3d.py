"""Offline empirical TBI comparisons; no production inference activation."""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import sys
from time import perf_counter

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from engine.country_tb.data_common import canonical,publish_immutable,REPO_ROOT,seal
from engine.country_tb.validation_data import (
    load_validation_snapshot,load_mortality_validation_snapshot,validation_split,artifact,COMPARISON,require,
)
from engine.country_tb.empirical_validation import (
    tbi_anchored_validation,incidence_anchored_validation,model_discrepancy,reviewer_response,REVIEW,
)
from engine.country_tb.review_parameters_3d import progression_sets,assay_scenario,evidence_registry
from engine.country_tb.incidence import load_who_snapshot,national_profile,WHO_2025
from engine.country_tb.pressure_shape import shape_from_incidence,apply_backcast
from engine.country_tb.validation_numerics import band_instantaneous_reference
from engine.country_tb.preview_incidence_inversion import make_mapping,observation,DMETA,CONCENTRATION
from engine.country_tb.incidence_inversion import inversion_gate

SNAPSHOT='tbi-validation-fe05d6d1ab8c'
MORTALITY='abs-life2022-2024-6785449a5200'


def _history(country,*,alternative=False):
    incidence=national_profile(load_who_snapshot(WHO_2025),country)
    shape=shape_from_incidence(incidence,review=REVIEW)
    pars={'slope_per_year':'-.01','minimum_shape':'.01','maximum_shape':'100','review':REVIEW}
    return apply_backcast(shape,start_time='1900',method='bounded_backward_trend' if alternative else 'constant_earliest',
        parameters=pars if alternative else None),incidence


def run(snapshot=SNAPSHOT,mortality_snapshot=MORTALITY,*,include_forward=True):
    inputs=load_validation_snapshot(snapshot);mortality=load_mortality_validation_snapshot(mortality_snapshot,geography_id='who:AUS')
    before=canonical(inputs);sets=progression_sets();registry=evidence_registry();byid={r['record_id']:r for r in inputs['rows']}
    definitions=[('vnm2020','PMC10446657',2020,'household_contacts_2020','QFT-Plus',30),
        ('vnm2021','PMC10446657',2021,'household_contacts_2021','QFT-Plus',30),
        ('USA_USborn2012','PMC4633161',2012,'us_born','QFT-GIT',15)]
    cases=[];case_sources={};histories={};targets={}
    for case,sid,year,subgroup,product,anchor_age in definitions:
        rows=[r for r in inputs['rows'] if r['source_id']==sid and r['survey_year']==year and r['subgroup_id']==subgroup and r['assay']==product and r['view']=='age']
        fit=[r['record_id'] for r in rows if r['age_lower']==anchor_age];valid=[r['record_id'] for r in rows if r['record_id'] not in fit]
        split=validation_split(inputs,fit,valid,kind='held_out_age_records',rationale='Prespecified anchor age; disjoint age records,household/design correlation retained.')
        country='VNM' if sid=='PMC10446657' else 'USA';history,national=_history(country);histories[case]=history
        target=next(r for r in national['rows'] if r['year']==year);targets[case]=target
        weights={r['record_id']:r['denominator'] or 1 for r in rows if r['record_id'] in valid}
        policy='tested_sample_count_weighted_descriptive_not_national_population' if sid=='PMC10446657' else 'equal_age_record_weights_not_survey_population_weights'
        base=tbi_anchored_validation(inputs,split,history,assay_scenario(product),metric_weights=weights,weight_policy=policy,
            scenario_id=case+'_central_candidate',incidence_target=target)
        cases.append(base);case_sources[case]=(split,weights,policy,product,rows)
        for variant in ('lower_specificity','higher_specificity'):
            cases.append(tbi_anchored_validation(inputs,split,history,assay_scenario(product,variant=variant),metric_weights=weights,
                weight_policy=policy,scenario_id=case+'_'+variant,incidence_target=target))
        if case in {'vnm2020','USA_USborn2012'}:
            alternative,_=_history(country,alternative=True)
            cases.append(tbi_anchored_validation(inputs,split,alternative,assay_scenario(product),metric_weights=weights,
                weight_policy=policy,scenario_id=case+'_bounded_backcast',incidence_target=target))
    split2020=case_sources['vnm2020'][0]
    ids2021=[r['record_id'] for r in case_sources['vnm2021'][4]]
    cross=validation_split(inputs,split2020['fit_ids'],ids2021,kind='cross_cohort_transport_comparison',
        rationale='Different campaign districts and convenience subsets; not a national longitudinal validation.')
    cases.append(tbi_anchored_validation(inputs,cross,histories['vnm2020'],assay_scenario('QFT-Plus'),
        metric_weights={i:byid[i]['denominator'] for i in ids2021},weight_policy='tested_sample_count_weighted_transport_stress_test',
        scenario_id='vnm2020_to2021_transport',incidence_target=targets['vnm2021']))
    forward=[]
    if include_forward:
        for case,source,year,subgroup,product,anchor in definitions:
            result=next(r for r in cases if r['scenario_id']==case+'_central_candidate')
            if result['fitted_scale'] is None:continue
            split,weights,policy,product,rows=case_sources[case]
            eligible=[r for r in rows if r['age_lower']>=15 and r['age_upper_exclusive'] is not None]
            for pars in sets:
                curve=pars['natural_history']['group_curves']['15+'];predictions=[]
                for r in eligible:
                    value=band_instantaneous_reference(histories[case],curve,float(result['fitted_scale']),
                        age_lower=r['age_lower'],age_upper_exclusive=r['age_upper_exclusive'],time=float(year))
                    predictions.append({'record_id':r['record_id'],'instantaneous_rate_per100000_disease_free_PY':value,
                        'weight':r['denominator'] or 1})
                mean=sum(r['weight']*r['instantaneous_rate_per100000_disease_free_PY'] for r in predictions)/sum(r['weight'] for r in predictions)
                target=float(targets[case]['central'])
                forward.append({'case_id':case,'parameter_set_id':pars['parameter_set_id'],'fitted_TBI_scale_unchanged':result['fitted_scale'],
                    'rows':predictions,'covered_adult_descriptive_rate':mean,'WHO_same_year_context_rate':target,
                    'difference':mean-target,'ratio':mean/target if target else None,
                    'status':'validation_only_population_and_episode_incompatible',
                    'mortality':'not_required_for_instantaneous_rate; no_long_horizon_country_projection',
                    'risk_policy':'RR1_reference_only; subgroup_risk_prevalences_not_supplied; not_a_complete_population_risk_profile',
                    'unresolved_record_ids':[r['record_id'] for r in rows if r not in eligible],
                    'population_warning':'Raw contact sample or equal survey-age weights,not national incidence; both-sex common curve,no sex counts invented.'})
    b=[]
    for country,case in [('VNM','vnm2020'),('AUS',None)]:
        geo='who:'+country;m=make_mapping(geography=geo,kind='diagnostic_country_population',population_history={})
        wr=next(r for r in national_profile(load_who_snapshot(WHO_2025),country)['rows'] if r['year']==2024)
        obs=observation(float(wr['central']),geography=geo,empirical_who=True,denominator=wr['denominator'],episode='new_and_recurrent',
            acquisition='unknown',measure='prospective_per_100000_person_years')
        receipt={'status':'pending','scenario_hash':m.specification['profile_hash'],'source_id':'M3D_no_country_configuration_approved',
            'provenance':{'gate_only_placeholder_not_population_estimate':True},'protocol_scope':'reviewed_empirical_configuration'}
        gate=inversion_gate(m,obs,unresolved_threshold_approval=None,scientific_review=receipt,concentration_policy=CONCENTRATION,**DMETA)
        rows=[r for r in inputs['rows'] if r['geography_id']==geo and r['view']=='age' and (r['assay']=='QFT-Plus' if country=='VNM' else True)]
        rows=rows[:1] # Gate refusal demonstration; no empirical fit or pooled estimate.
        b.append(incidence_anchored_validation(m,obs,gate,inputs,[r['record_id'] for r in rows],assay_scenario('QFT-Plus'),
            lower=0,upper=.2,metric_weights={r['record_id']:1 for r in rows},weight_policy='gate_refusal_only_no_population_estimate',
            independence_provenance={'TBI_observations_used_to_construct_incidence_target':False,'incidence_source':WHO_2025,'TBI_source':snapshot},**DMETA))
    components=['WHO_estimation_uncertainty','observation_incompatibility','migration','imported_disease','recurrence','prior_TB',
        'unobserved_risk_factors','assay_error','historical_incidence_revision','pressure_shape_misspecification','endogenous_transmission','age_sex_aggregation','mortality']
    discrepancy=model_discrepancy(source_hash=inputs['profile_hash'],scenario_hash=sets[1]['profile_hash'],entries=[
        {'component':c,'status':'not_separately_identified','direction':'source_or_population_dependent','magnitude':None,
            'evidence':'see_evidence_registry_and_source_compatibility','required_data':'matched_outcomes_and_component_specific_data_not_an_alpha_offset'} for c in components])
    require(canonical(inputs)==before,'Validation altered source observations.')
    return artifact(COMPARISON,direction='M3D_internal_empirical_evidence_package',status='internal_review_package_not_external_validation',
        validation_input=inputs,snapshot_id=snapshot,evidence_registry=registry,parameter_sets=sets,TBI_anchored_cases=cases,
        forward_contextual_comparisons=forward,incidence_anchored_cases=b,mortality_validation_snapshot=mortality,
        australia={'observed_rows':[r for r in inputs['rows'] if r['geography_id']=='who:AUS'],
            'status':'descriptive_premigration_application_data_not_resident_TBI_validation',
            'locally_born':'data_required','long_term_migrant':'data_required','clinical_high_risk':'compatible_TBI_and_followup_data_required'},
        discrepancy=discrepancy,reviewer_response=reviewer_response(['adult','0-4','5-14','remote_tail','mortality','risk_factors','recurrence','imported_disease','assay','migration']),
        complete_country_estimate_available=False,central_parameter_set_approved=False,
        dynamic_handoff={'role':['initialization','diagnostic','comparison_target','bounded_prior'],
            'required':['source_and_parameter_hashes','age_sex_subgroup_population_and_unresolved_counts','latent_vs_assay_profiles',
                'most_recent_infection_distribution','RR_definitions_and_memberships','migration_and_exposure_metadata','prior_TB_and_active_baseline_exclusions'],
            'not_a_replacement_for_endogenous_transmission':True,'Summer_or_Starsim_modified':False})


def main():
    p=argparse.ArgumentParser();p.add_argument('--snapshot',default=SNAPSHOT);p.add_argument('--output',type=Path,required=True);p.add_argument('--skip-forward',action='store_true');a=p.parse_args()
    output=a.output.resolve();require(output.is_relative_to((REPO_ROOT/'abm/output').resolve()),'Generated outputs must remain under abm/output.')
    start=perf_counter();v=run(a.snapshot,include_forward=not a.skip_forward)
    target=publish_immutable(output,'empirical3d-'+v['profile_hash'][:12],{'validation.json':canonical(v)+b'\n'})
    print(json.dumps({'directory':str(target),'elapsed_seconds':perf_counter()-start,'scientific_hash':v['profile_hash'],
        'records':len(v['validation_input']['rows']),'cases':[{'id':r['scenario_id'],'status':r['status'],'alpha':r.get('fitted_scale'),
            'RMSE':r.get('held_out_metrics',{}).get('weighted_RMSE')} for r in v['TBI_anchored_cases']]}))


if __name__=='__main__':main()
