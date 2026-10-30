"""Separated empirical fit/validation directions and descriptive metrics.

Uses the retained2C calibration/assay contracts. No incidence adjustment,
clinical parameter approval or production source dispatch is performed here.
"""
from copy import deepcopy
from math import fsum,sqrt

from engine.country_tb.data_common import canonical,decimal,numeric_text,DataValidationError
from engine.country_tb.validation_data import (
    artifact,verify,require,probability,OBSERVATIONS,SPLIT,COMPARISON,DISCREPANCY,RESPONSE,
)
from engine.country_tb.catalytic_core import evaluate_age
from engine.country_tb.catalytic_calibration import calibration_policy,calibrate_scale
from engine.country_tb.tbi_observations import observations_profile,assay_probability
from engine.country_tb.incidence_inversion import diagnostic_calibration

REVIEW={'status':'reviewed','source_id':'M3D_internal_extraction_and_diagnostic_protocol',
    'rationale':'Internal transcription and mathematical diagnostic use only; no external scientific sign-off or country parameter approval.'}


def comparison_metrics(records,predictions,*,weights,weight_policy,source_hash,scenario_hash):
    """Weights must be explicit; complex-survey percentages never become counts."""
    require(weight_policy and set(weights)=={r['record_id'] for r in records},'Metric weights must cover exactly the declared comparison rows.')
    require(all(decimal(w)>0 for w in weights.values()),'Positive metric weights required.')
    rows=[];terms=[];total_weight=0.;covered_weight=0.;known_n=unresolved_n=0;counts_complete=True
    for r in sorted(records,key=lambda r:r['record_id']):
        identifier=r['record_id'];require(identifier in predictions,'Missing prediction status, not an implicit zero.')
        pred=predictions[identifier];value=pred.get('predicted_positivity')
        obs=float(probability(r['observed_positivity']));w=float(decimal(weights[identifier]));total_weight+=w
        if r['denominator'] is None:counts_complete=False
        else:known_n+=r['denominator']
        if value is None:
            require(pred.get('unavailable_reason'),'Unavailable prediction must identify missing input.')
            residual=absolute=relative=ratio=inside=None
            if r['denominator'] is not None:unresolved_n+=r['denominator']
        else:
            value=float(probability(value));residual=value-obs;absolute=abs(residual)
            relative=None if obs==0 else residual/obs;ratio=None if obs==0 else value/obs
            inside=None if r['lower'] is None else float(r['lower'])<=value<=float(r['upper'])
            terms.append(w*residual*residual);covered_weight+=w
        rows.append(dict(record_id=identifier,age_lower=r['age_lower'],age_upper_exclusive=r['age_upper_exclusive'],
            sex=r['sex'],subgroup_id=r['subgroup_id'],survey_year=r['survey_year'],assay=r['assay'],
            observed_positivity=obs,predicted_positivity=value,residual=residual,absolute_error=absolute,
            relative_error=relative,calibration_ratio=ratio,inside_reported_interval=inside,
            direction_of_bias=None if residual is None else 'overprediction' if residual>0 else 'underprediction' if residual<0 else 'equal',
            metric_weight=weights[identifier],unavailable_reason=pred.get('unavailable_reason'),source_compatibility=r['compatibility']))
    def group_summary(sex):
        selected=[r for r in rows if r['sex']==sex and r['residual'] is not None]
        den=fsum(float(r['metric_weight']) for r in selected)
        return {'records':len(selected),'weighted_mean_residual':None if not den else fsum(float(r['metric_weight'])*r['residual'] for r in selected)/den}
    intervals=[r['inside_reported_interval'] for r in rows if r['inside_reported_interval'] is not None]
    age_rows=sorted((r for r in rows if r['sex']=='both' and r['age_upper_exclusive'] is not None),key=lambda r:r['age_lower'])
    observed_decreases=[];predicted_decreases=[]
    for a,b in zip(age_rows,age_rows[1:]):
        if a['age_upper_exclusive']>b['age_lower']:continue
        if b['observed_positivity']<a['observed_positivity']:observed_decreases.append([a['record_id'],b['record_id']])
        if a['predicted_positivity'] is not None and b['predicted_positivity'] is not None and b['predicted_positivity']<a['predicted_positivity']:
            predicted_decreases.append([a['record_id'],b['record_id']])
    return artifact(COMPARISON,source_hash=source_hash,scenario_hash=scenario_hash,rows=rows,
        metric_weight_policy=weight_policy,weighted_RMSE=None if not covered_weight else sqrt(fsum(terms)/covered_weight),
        covered_metric_weight=covered_weight,total_metric_weight=total_weight,
        interval_inclusion_fraction=None if not intervals else sum(intervals)/len(intervals),
        interval_inclusion_count=len(intervals),interval_warning='Descriptive interval inclusion, not independent coverage probability or external validation.',
        sex_specific={s:group_summary(s) for s in ('female','male','both')},
        observed_age_decreases=observed_decreases,predicted_age_decreases=predicted_decreases,
        unresolved_tested_population=unresolved_n if counts_complete else None,known_tested_count=known_n,
        population_warning='Tested counts, not national population; survey weights remain metric weights, not people.',
        scientific_acceptance='not_assessed_pending_reviewer_threshold_approval',
        candidate_thresholds={'absolute_error':[.02,.05,.10],'weighted_RMSE':[.025,.05,.10],
            'interval_inclusion':[.8,.9,.95],'unresolved_population_share':[0,.001,.01]},thresholds_approved=False)


def _eligible(row,assay):
    if row['age_upper_exclusive'] is None:return 'open_age_or_unspecified_age_mix_requires_reviewed_weights'
    if row['age_lower']<assay['age_lower']:return 'assay_performance_external_required_for_children'
    if row['assay']!=assay['operator']['name']:return 'assay_generation_or_mixed_assay_mismatch'
    if row['view']!='age':return 'joint_age_distribution_not_reported_for_this_marginal'
    if row['compatibility']['tbi_anchor'] not in {'diagnostic_subgroup_only','diagnostic_us_born_only'}:
        return 'exposure_history_or_population_compatibility_unresolved'
    return None


def _record(row,assay):
    uncertainty={'kind':'unspecified'} if row['lower'] is None else dict(kind='interval',lower=row['lower'],upper=row['upper'],method='published95CI')
    return dict(observation_id=row['record_id'],geography_id=row['geography_id'],geography_level=row['geography_level'],
        age_lower=row['age_lower'],age_upper_exclusive=row['age_upper_exclusive'],sex=row['sex'],
        survey_year=row['survey_year'],survey_time=str(row['survey_year']),numerator=row['numerator'],denominator=row['denominator'],
        prevalence=row['observed_positivity'],count_prevalence_tolerance='1e-50',review=REVIEW,
        source={'source_id':row['source_id'],'source_record_sha256':row['source_record_sha256']},sampling_frame=row['sampling_frame'],
        survey_design=row['survey_design'],uncertainty=uncertainty,assay_operator=assay['operator'],
        weight_source='explicit_uniform_within_reported_finite_band_diagnostic_assumption_not_source_demography',
        strata=[dict(age_lower=row['age_lower'],age_upper_exclusive=row['age_upper_exclusive'],sex=row['sex'],
            age_policy='reviewed_uniform_within_interval',weight='1',age_review=REVIEW)],
        observation_policy='assay_positivity_not_direct_latent_TBI_diagnostic_at_survey_year_start_v1')


def tbi_anchored_validation(inputs,split,history,assay,*,metric_weights,weight_policy,scenario_id,
    method='overall_anchor',forward=None,incidence_target=None):
    verify(inputs,OBSERVATIONS);verify(split,SPLIT)
    require(split['input_hash']==inputs['profile_hash'],'Split belongs to different immutable input.')
    before=canonical(inputs);byid={r['record_id']:r for r in inputs['rows']}
    fits=[byid[i] for i in split['fit_ids']];valid=[byid[i] for i in split['validation_ids']]
    common=dict(direction='A_TBI_anchored_then_incidence_comparison',source_hash=inputs['profile_hash'],split=split,
        scenario_id=scenario_id,assay_scenario=assay,pressure_history_hash=history['profile_hash'],
        timing_policy='survey_year_start; multi_year_surveys_use_explicit_representative_year; not_2024_validation',
        acceptance='internal_diagnostic_only_no_external_signoff')
    reasons=[_eligible(r,assay) for r in fits]
    if any(reasons):return artifact(COMPARISON,status='blocked_anchor_inputs',fitted_scale=None,failures=reasons,**common)
    require(all(r['geography_id']==history['geography_id'] for r in fits+valid),'No silent country transport.')
    observations=observations_profile([_record(r,assay) for r in fits],review=REVIEW)
    require(method in {'overall_anchor','weighted_least_squares'},'Empirical cluster/weighted survey rows do not imply independent binomial sampling.')
    require(method!='overall_anchor' or len(fits)==1,'Single anchor required for exact root.')
    policy=calibration_policy(alpha_upper='.5',observation_weights={r['record_id']:'1' for r in fits},
        residual_tolerance='1',review=REVIEW,method=method)
    try:calibration=calibrate_scale(history,policy,observations=observations)
    except DataValidationError as exc:
        return artifact(COMPARISON,status='infeasible_assay_or_pressure_anchor',fitted_scale=None,failures=[str(exc)],**common)
    alpha=calibration['alpha'];predictions={};latent=[];forward_rows=[]
    for r in valid:
        reason=_eligible(r,assay)
        if reason:predictions[r['record_id']]={'predicted_positivity':None,'unavailable_reason':reason};continue
        result=evaluate_age(history,alpha=alpha,age=r['age_lower'],age_upper_exclusive=r['age_upper_exclusive'],sex=r['sex'],
            time=str(r['survey_year']),age_policy='reviewed_uniform_within_interval',review=REVIEW,include_event_times=False)
        predicted=assay_probability(result['any_tbi'],assay['operator'])
        predictions[r['record_id']]={'predicted_positivity':numeric_text(predicted)}
        latent.append({'record_id':r['record_id'],'model_derived_latent_TBI':result,'not_observed_recency':True})
        if forward is not None:forward_rows.append({'record_id':r['record_id'],'prediction':forward(alpha,r)})
    metrics=comparison_metrics(valid,predictions,weights=metric_weights,weight_policy=weight_policy,
        source_hash=inputs['profile_hash'],scenario_hash=assay['profile_hash'])
    require(before==canonical(inputs),'Empirical validation changed a source observation.')
    return artifact(COMPARISON,status='TBI_anchored_internal_diagnostic',fitted_scale=alpha,calibration=calibration,
        held_out_metrics=metrics,latent_profiles=latent,forward_rows=forward_rows,incidence_target=incidence_target,
        incidence_comparison_status='unavailable_no_population_compatible_first_episode_target' if incidence_target is None else 'validation_only_contextual_population_or_episode_mismatch_not_calibration',
        sources_unchanged=True,**common)


def incidence_anchored_validation(mapping,observation,gate,inputs,validation_ids,assay,*,lower,upper,source,review,provenance_record,
    metric_weights,weight_policy,independence_provenance):
    """Uses the existing3C gate; no empirical gate exception is created."""
    verify(inputs,OBSERVATIONS)
    require(independence_provenance and independence_provenance.get('TBI_observations_used_to_construct_incidence_target') is False,
        'Independent incidence/TBI outcome provenance required; no derived-target leakage.')
    byid={r['record_id']:r for r in inputs['rows']}
    require(validation_ids and len(set(validation_ids))==len(validation_ids) and set(validation_ids)<=byid.keys(),'Explicit independent validation records required.')
    records=[byid[i] for i in validation_ids]
    require(all(r['geography_id']==mapping.history['geography_id'] for r in records),'Validation geography mismatch.')
    result=diagnostic_calibration(mapping,observation,gate,lower=lower,upper=upper,grid_intervals=64,
        source=source,review=review,provenance_record=provenance_record)
    predictions={}
    for r in records:
        reason='blocked_incidence_scientific_or_numerical_gate' if result['fitted_scale'] is None else _eligible(r,assay)
        if reason:predictions[r['record_id']]={'predicted_positivity':None,'unavailable_reason':reason};continue
        age=evaluate_age(mapping.history,alpha=result['fitted_scale'],age=r['age_lower'],age_upper_exclusive=r['age_upper_exclusive'],
            sex=r['sex'],time=str(r['survey_year']),age_policy='reviewed_uniform_within_interval',review=REVIEW,include_event_times=False)
        predictions[r['record_id']]={'predicted_positivity':numeric_text(assay_probability(age['any_tbi'],assay['operator']))}
    metrics=comparison_metrics(records,predictions,weights=metric_weights,weight_policy=weight_policy,
        source_hash=inputs['profile_hash'],scenario_hash=mapping.specification['profile_hash'])
    return artifact(COMPARISON,direction='B_incidence_anchored_independent_TBI_validation',
        status=result['status'],diagnostic_status='diagnostic_only_not_approved_for_production',
        fitted_scale=result['fitted_scale'],incidence_calibration=result,
        validation_source_hash=inputs['profile_hash'],validation_ids=sorted(validation_ids),assay_scenario=assay,
        TBI_observations_used_in_incidence_fit=False,
        independence_provenance=independence_provenance,held_out_metrics=metrics,
        empirical_metrics_status='unavailable_scientific_gate' if result['fitted_scale'] is None else 'requires_matched_validation_history_and_year')


def model_discrepancy(*,source_hash,scenario_hash,entries):
    allowed={'WHO_estimation_uncertainty','observation_incompatibility','migration','imported_disease','recurrence',
        'prior_TB','unobserved_risk_factors','assay_error','historical_incidence_revision','pressure_shape_misspecification',
        'endogenous_transmission','age_sex_aggregation','mortality'}
    require({r['component'] for r in entries}==allowed and len(entries)==len(allowed),'Every discrepancy component must be explicit.')
    require(all(set(r)>={'component','status','direction','magnitude','evidence','required_data'} for r in entries),'Complete discrepancy record required.')
    return artifact(DISCREPANCY,source_hash=source_hash,scenario_hash=scenario_hash,rows=sorted(deepcopy(entries),key=lambda r:r['component']),
        discrepancy_fitted_into_alpha=False,estimated_total_discrepancy=None,
        conclusion='A residual does not identify its causal source; source/scenario/numerical uncertainty kept separate.')


def reviewer_response(component_ids):
    require(len(set(component_ids))==len(component_ids),'Unique review components required.')
    return artifact(RESPONSE,rows=[{'component_id':i,'decision':None,'allowed_decisions':['accept','accept_with_changes','reject'],
        'comments':None,'alternative_source_or_value':None,'reviewer':None,'date':None} for i in sorted(component_ids)],
        status='external_review_not_yet_obtained',automatically_approved=False)
