"""Deterministic baseline stratification. No active disease or intervention operator."""
from collections import defaultdict
from fractions import Fraction
from math import prod

from engine.country_tb.risk_common import (
    artifact,verify,require,metadata,ratio,amount,text_number,canonical,cell_key,clinical_group,matches,
)
from engine.country_tb.risk_factor_catalog import catalogue_records,applicable_effect
from engine.country_tb.risk_prevalence import resolve_risk_prevalence
from engine.country_tb.risk_policies import (
    OVERLAP,ASSOCIATION,MULTIMORBIDITY,INDEPENDENT,JOINT,HIGHEST,
    control_for_cell,feasible_association,phi,
)
from engine.country_tb.population_subgroups import validate_resolution

ALLOCATION='country_tb.baseline_population_risk_allocation.v1'
MIXTURE='country_tb.population_risk_allocation_mixture.v1'
TB_STATES=('any_tbi','never_infected','recent','remote_only','prior_remote_plus_recent')


def _risk_distribution(marginals,effects,overlap,cell):
    """DP on (highest RR, factor count); unknown effects form a null-RR class."""
    ids=sorted(marginals);dp={(Fraction(1),0):Fraction(1)};winner={};patterns=None
    if overlap['method']==INDEPENDENT:
        for factor in ids:
            p=marginals[factor];rr=None if effects[factor]['rr'] is None else ratio(effects[factor]['rr'])
            updated=defaultdict(Fraction)
            for (highest,n),mass in dp.items():
                if p<1:updated[highest,n]+=mass*(1-p)
                if p:updated[None if highest is None or rr is None else max(highest,rr),n+1]+=mass*p
            dp=dict(updated)
        for factor in ids:
            rr=effects[factor]['rr']
            winner[factor]=None if rr is None else marginals[factor]*prod((1-marginals[k] for k in ids if
                effects[k]['rr'] is None or ratio(effects[k]['rr'])>ratio(rr)),start=Fraction(1))
    else:
        selected=[r for r in overlap['rows'] if matches(r,cell)]
        require(len(selected)==1,'Reviewed joint data must cover the whole age/sex cell.')
        patterns=selected[0]['patterns'];dp=defaultdict(Fraction);observed={k:Fraction(0) for k in ids}
        winner={k:Fraction(0) if effects[k]['rr'] is not None else None for k in ids}
        for pattern in patterns:
            present=pattern['factor_ids'];mass=ratio(pattern['probability'],probability=True)
            require(set(present)<=set(ids),'Reviewed joint patterns include factors outside the included scope.')
            for k in present:observed[k]+=mass
            rr=None if any(effects[k]['rr'] is None for k in present) else max([Fraction(1),*[ratio(effects[k]['rr']) for k in present]])
            dp[rr,len(present)]+=mass
            if rr is not None:
                for k in present:
                    if ratio(effects[k]['rr'])==rr:winner[k]+=mass
        require(observed==marginals,'Reviewed joint distribution does not preserve supplied factor marginals.')
        dp=dict(dp)
    require(sum(dp.values())==1,'Risk pattern mass does not reconcile exactly.')
    counts=defaultdict(Fraction)
    for (_,n),mass in dp.items():counts[n]+=mass
    q=1-counts.get(0,Fraction(0))
    return dp,q,counts,winner,patterns


def _state_counts(tbi,any_count,base_count):
    p=tbi['any_tbi'];result={'any_tbi':any_count,'never_infected':base_count-any_count}
    for name in ('recent','remote_only','prior_remote_plus_recent'):
        value=tbi[name]
        result[name]=None if value is None else any_count*value/p if p else Fraction(0)
    return result


def _packed_counts(values):return {k:amount(v) if v is not None else None for k,v in values.items()}


def _metrics(rows):
    total=sum((ratio(r['population']) for r in rows),Fraction(0))
    fields=('tbi_prevalence','any_risk_prevalence','joint_tbi_any_risk','multimorbidity_prevalence','expected_factor_count')
    summary={'population':amount(total)}
    for field in fields:
        covered=sum((ratio(r['population']) for r in rows if r[field] is not None),Fraction(0))
        count=sum((ratio(r['population'])*ratio(r[field]) for r in rows if r[field] is not None),Fraction(0))
        summary[field]={'covered_population':amount(covered),'unresolved_population':amount(total-covered),
            'expected_count':amount(count),'prevalence_covered':amount(count/covered) if covered else None,
            'prevalence_full':amount(count/total) if total and covered==total else None}
    unresolved=sum((ratio(r['population']) for r in rows if not r['complete_expected_stratification']),Fraction(0))
    summary.update(unresolved_population=amount(unresolved),complete_expected_stratification=unresolved==0)
    state_summary={}
    for state in TB_STATES:
        selected=[r for r in rows if r['tbi_state_prevalences'][state] is not None]
        covered=sum((ratio(r['population']) for r in selected),Fraction(0))
        count=sum((ratio(r['population'])*ratio(r['tbi_state_prevalences'][state]) for r in selected),Fraction(0))
        state_summary[state]={'covered_population':amount(covered),'expected_count':amount(count),
            'complete_count':amount(count) if covered==total else None}
    summary['tbi_state_counts']=state_summary
    rr=defaultdict(Fraction)
    for row in rows:
        for entry in row['effective_rr_distribution'] or []:
            key=ratio(entry['effective_rr']) if entry['effective_rr'] is not None else None
            rr[key]+=ratio(entry['expected_population'])
    summary['effective_rr_distribution']=[{'effective_rr':amount(k) if k is not None else None,'expected_population':amount(v)}
        for k,v in sorted(rr.items(),key=lambda x:(x[0] is None,x[0] or 0))]
    summary['factor_marginals']=[]
    factor_ids=sorted({f['factor_id'] for row in rows for f in row['factor_marginals']})
    for factor in factor_ids:
        selected=[(r,f) for r in rows for f in r['factor_marginals'] if f['factor_id']==factor and f['prevalence'] is not None]
        covered=sum((ratio(r['population']) for r,_ in selected),Fraction(0))
        count=sum((ratio(r['population'])*ratio(f['prevalence']) for r,f in selected),Fraction(0))
        summary['factor_marginals'].append({'factor_id':factor,'covered_population':amount(covered),
            'unresolved_or_outside_scope_population':amount(total-covered),'expected_population':amount(count),
            'prevalence_covered':amount(count/covered) if covered else None,
            'prevalence_full':amount(count/total) if covered==total and total else None})
    factor_counts=defaultdict(Fraction);count_covered=Fraction(0)
    for row in rows:
        if row['factor_count_distribution'] is None:continue
        count_covered+=ratio(row['population'])
        for entry in row['factor_count_distribution']:factor_counts[entry['factor_count']]+=ratio(entry['expected_population'])
    summary['factor_count_distribution']={'covered_population':amount(count_covered),
        'unresolved_population':amount(total-count_covered),'rows':[{'factor_count':n,'expected_population':amount(v)} for n,v in sorted(factor_counts.items())]}
    # Overall covariance includes structural age/sex/subgroup effects. Controls
    # and phi are never averaged before allocation.
    eligible=[r for r in rows if r['association'] is not None]
    covered=sum((ratio(r['population']) for r in eligible),Fraction(0))
    if covered:
        p=sum((ratio(r['population'])*ratio(r['tbi_prevalence']) for r in eligible),Fraction(0))/covered
        q=sum((ratio(r['population'])*ratio(r['any_risk_prevalence']) for r in eligible),Fraction(0))/covered
        j=sum((ratio(r['population'])*ratio(r['joint_tbi_any_risk']) for r in eligible),Fraction(0))/covered
        j0=sum((ratio(r['population'])*ratio(r['tbi_prevalence'])*ratio(r['any_risk_prevalence']) for r in eligible),Fraction(0))/covered
        achieved,reason=phi(p,q,j)
        summary['aggregate_association']={'covered_population':amount(covered),'complete_population':covered==total,
            'tbi_prevalence':amount(p),'any_risk_prevalence':amount(q),'joint_prevalence':amount(j),'achieved_phi':achieved,
            'phi_unavailable_reason':reason,'within_stratum_covariance':amount(j-j0),
            'structural_age_sex_subgroup_covariance':amount(j0-p*q),'requested_controls_are_not_averaged':True}
    else:summary['aggregate_association']=None
    return summary


def allocate_baseline_risk(resolved,*,catalog,prevalence,overlap,association,progression_policy,
    source,provenance,review,diagnostic_only=False):
    metadata(source,provenance,review);validate_resolution(resolved)
    factors=catalogue_records(catalog);verify(overlap,OVERLAP);verify(association,ASSOCIATION);verify(progression_policy,MULTIMORBIDITY)
    require(progression_policy['method']==HIGHEST and not progression_policy['calculation_only_sensitivity'],
            'Population allocation uses highest RR centrally; sensitivity combinations are calculation-only.')
    require(progression_policy['effect_use']!='diagnostic_provisional_adult_effects_v1' or diagnostic_only,
            'Provisional effects require explicitly diagnostic allocation.')
    require(prevalence['subgroup_id']==resolved['subgroup_id'] and prevalence['geography_id']==resolved['geography_id'],
            'Risk profile must match the resolved subgroup and current-residence geography.')
    resolved_prevalence,pattern_checks=resolve_risk_prevalence(prevalence,catalog,resolved['rows'])
    rows=[]
    for cell,exposure in zip(resolved['rows'],resolved_prevalence):
        require(cell_key(cell)==exposure['cell_key'],'Age/sex risk input alignment failed.')
        population=ratio(cell['population']);group=clinical_group(cell)
        tbi={k:None if cell[k] is None else ratio(cell[k],probability=True) for k in TB_STATES}
        c,control_provenance=control_for_cell(association,cell)
        effects={k:applicable_effect(factors[k],cell,effect_use=progression_policy['effect_use'],diagnostic_only=diagnostic_only) for k in prevalence['factor_ids']}
        if group=='unresolved_crosses_clinical_age_boundary':
            effects={k:{'rr':None,'state':'unresolved_crosses_clinical_age_boundary'} for k in effects}
        factor_rows=[];marginals={};unknown=[]
        for record in exposure['factors']:
            k=record['factor_id'];value=record['prevalence']
            factor_rows.append({**record,'expected_population':amount(population*ratio(value)) if value is not None else None,
                                'effect':effects[k]})
            if record['state']=='value':marginals[k]=ratio(value)
            elif record['state'] not in {'excluded','not_applicable'}:unknown.append(k)
        row={k:cell[k] for k in ('sex','age_lower','age_upper_exclusive')}
        row.update(subgroup_id=resolved['subgroup_id'],geography_id=resolved['geography_id'],baseline_year=2024,
            population=amount(population),clinical_age_group=group,tbi_prevalence=amount(tbi['any_tbi']) if tbi['any_tbi'] is not None else None,
            tbi_state_prevalences=_packed_counts(tbi),factor_marginals=factor_rows,
            any_risk_prevalence=None,joint_tbi_any_risk=None,multimorbidity_prevalence=None,expected_factor_count=None,
            factor_count_distribution=None,effective_rr_distribution=None,compressed_cells=[],individual_factor_state_counts=[],
            normalized_association_input=str(c.numerator) if c.denominator==1 else association['global_control'] if control_provenance['scope']=='global' else None,
            control_provenance=control_provenance,association=None,highest_rr_contributor_probabilities=[],
            complete_expected_stratification=False,unresolved_factor_ids=unknown,
            source_state_closure_residual=None,warnings=[])
        # Store the requested control itself, not a substituted aggregate value.
        row['normalized_association_input']=text_number(c)
        if unknown:
            row['warnings'].append('missing_unreviewed_mismatched_or_pediatric_prevalence_not_zero_filled')
            rows.append(row);continue
        dp,q,counts,winners,patterns=_risk_distribution(marginals,effects,overlap,cell)
        row.update(any_risk_prevalence=amount(q),multimorbidity_prevalence=amount(sum((m for n,m in counts.items() if n>=2),Fraction(0))),
            expected_factor_count=amount(sum((n*m for n,m in counts.items()),Fraction(0))),
            factor_count_distribution=[{'factor_count':n,'probability':amount(m),'expected_population':amount(population*m)} for n,m in sorted(counts.items())],
            highest_rr_contributor_probabilities=[{'factor_id':k,'probability':amount(v) if v is not None else None,
                'interpretation':'nonexclusive_highest_contributor_probability_all_ties_counted'} for k,v in sorted(winners.items())])
        assoc=None
        if tbi['any_tbi'] is not None:
            assoc=feasible_association(tbi['any_tbi'],q,c);row['association']=assoc;row['joint_tbi_any_risk']=assoc['joint_prevalence']
        rr_mass=defaultdict(Fraction)
        for (rr,n),mass in sorted(dp.items(),key=lambda x:(x[0][1],x[0][0] is None,x[0][0] or 0)):
            rr_mass[rr]+=mass
            states={k:None for k in TB_STATES}
            if assoc is not None:
                j=ratio(assoc['joint_prevalence']);p=tbi['any_tbi']
                conditional=j/q if n and q else (p-j)/(1-q) if not n and q<1 else Fraction(0)
                states=_state_counts(tbi,population*mass*conditional,population*mass)
            row['compressed_cells'].append({'factor_count':n,'any_non_recency_risk_factor':n>0,
                'effective_rr':amount(rr) if rr is not None else None,'probability':amount(mass),
                'expected_population':amount(population*mass),'tbi_state_counts':_packed_counts(states)})
        row['effective_rr_distribution']=[{'effective_rr':amount(rr) if rr is not None else None,
            'probability':amount(mass),'expected_population':amount(population*mass),
            'possible_responsible_factor_ids':sorted(k for k in marginals if marginals[k]>0 and effects[k]['rr'] is not None and ratio(effects[k]['rr'])==rr),
            'tie_handling':'all_possible_tied_contributors_retained_actual_membership_evaluated_by_combination_API'}
            for rr,mass in sorted(rr_mass.items(),key=lambda x:(x[0] is None,x[0] or 0))]
        for k,pj in marginals.items():
            positive={state:None for state in TB_STATES};negative=dict(positive)
            if assoc is not None:
                joint_factor=population*pj*(ratio(assoc['joint_prevalence'])/q if q else 0)
                positive=_state_counts(tbi,joint_factor,population*pj)
                negative={state:None if positive[state] is None else population*tbi[state]-positive[state] for state in TB_STATES}
            row['individual_factor_state_counts'].append({'factor_id':k,
                'present_population':amount(population*pj),'absent_population':amount(population*(1-pj)),
                'present_tbi_state_counts':_packed_counts(positive),'absent_tbi_state_counts':_packed_counts(negative)})
        if assoc is not None:
            for state in TB_STATES:
                if tbi[state] is not None:
                    require(sum((ratio(r['tbi_state_counts'][state]) for r in row['compressed_cells']),Fraction(0))==population*tbi[state],
                            'Allocation changed a source TBI-state marginal.')
            closure=None if tbi['recent'] is None else tbi['any_tbi']-tbi['recent']-tbi['remote_only']
            row['source_state_closure_residual']=amount(closure) if closure is not None else None
            if closure:row['warnings'].append('inherited_source_decimal_state_closure_roundoff_recorded_no_source_values_modified')
            row['complete_expected_stratification']=(tbi['recent'] is not None and not closure and None not in rr_mass and
                group!='unresolved_crosses_clinical_age_boundary')
        rows.append(row)
    summary=_metrics(rows)
    display=[]
    for label,lo,hi in [('0-4',0,5),('5-14',5,15),('15-44',15,45),('45-64',45,65),('65+',65,None)]:
        subset=[r for r in rows if r['age_lower']>=lo and (hi is None or (r['age_upper_exclusive'] is not None and r['age_upper_exclusive']<=hi))]
        display.append({'band':label,'age_lower':lo,'age_upper_exclusive':hi,'summary':_metrics(subset)})
    display_covered=sum((ratio(d['summary']['population']) for d in display),Fraction(0))
    return artifact(ALLOCATION,geography_id=resolved['geography_id'],subgroup_id=resolved['subgroup_id'],baseline_year=2024,
        resolved_baseline=resolved,resolved_baseline_hash=resolved['profile_hash'],catalogue=catalog,prevalence_profile=prevalence,
        overlap_policy=overlap,association_policy=association,multimorbidity_policy=progression_policy,
        rows=sorted(rows,key=canonical),summary=summary,display_bands=display,display_band_unresolved_population=amount(ratio(summary['population'])-display_covered),
        clinical_groups=[{'band':band,'summary':_metrics([r for r in rows if r['clinical_age_group']==band])} for band in ('0-4','5-14','15+')],
        pattern_reconciliation=pattern_checks,diagnostic_only=diagnostic_only or resolved['profile_kind']=='diagnostic',
        ready_as_input_for_future_operator=summary['complete_expected_stratification'] and not diagnostic_only and resolved['usable_for_downstream_simulation'],
        complete_result=summary if summary['complete_expected_stratification'] else None,
        complexity={'method':'DP_over_highest_RR_and_factor_count','maximum_factors':11,'maximum_independent_DP_states':144,
                    'reviewed_joint_maximum_patterns':2048,'no_independent_pattern_enumeration':True},
        stochastic_allocator='deferred_explicit_RNG_individual_generator_contract_documented',
        individual_generator_contract={'explicit_rng_or_seed_required':True,'global_rng_forbidden':True,
            'steps':['sample_unconditional_factor_pattern_from_declared_overlap_policy','derive_any_risk_union',
                'sample_TBI_using_J_over_q_or_p_minus_J_over_1_minus_q','sample_recent_remote_given_TBI_from_unchanged_source',
                'retain_prior_remote_subset_without_double_counting','evaluate_all_present_factors_and_retain_all_highest_RR_ties'],
            'infection_times_require_original_reviewed_distribution_not_a_band_midpoint':True,
            'validation':'compare_realized_marginals_to_exact_expected_counts_with_predeclared_Monte_Carlo_tolerances'},
        conditional_pattern_assumption='same_pattern_given_any_risk_for_TBI_and_no_TBI_no_factor_specific_association_inferred',
        infection_history_policy='source_recency_and_most_recent_event_summary_retained_no_new_events_or_recalibration',
        source=source,provenance=provenance,review=review,forward_operator=None)


def combine_risk_allocations(baseline_mixture,allocations,*,source,provenance,review):
    from engine.country_tb.population_subgroups import MIXTURE as BASE_MIXTURE
    from engine.country_tb.composition_common import verify as verify_base
    metadata(source,provenance,review);verify_base(baseline_mixture,BASE_MIXTURE)
    expected={g['subgroup_id']:g['tbi_profile_hash'] for g in baseline_mixture['subgroups']}
    require(isinstance(allocations,list) and len({a['subgroup_id'] for a in allocations})==len(allocations) and
        {a['subgroup_id'] for a in allocations}==set(expected),'Risk allocations must cover each baseline subgroup exactly once.')
    for a in allocations:
        verify(a,ALLOCATION);require(a['resolved_baseline_hash']==expected[a['subgroup_id']],'Subgroup allocation belongs to a different baseline TBI profile.')
    rows=[r for a in sorted(allocations,key=lambda a:a['subgroup_id']) for r in a['rows']];summary=_metrics(rows)
    require(ratio(summary['population'])==ratio(baseline_mixture['total_population']),'Risk-mixture population reconciliation failed.')
    return artifact(MIXTURE,geography_id=baseline_mixture['geography_id'],baseline_year=2024,
        baseline_mixture=baseline_mixture,baseline_mixture_hash=baseline_mixture['profile_hash'],
        subgroups=sorted(allocations,key=lambda a:a['subgroup_id']),summary=summary,
        complete_result=summary if summary['complete_expected_stratification'] else None,
        unresolved_subgroups=[{'subgroup_id':a['subgroup_id'],'unresolved_population':a['summary']['unresolved_population']}
            for a in sorted(allocations,key=lambda a:a['subgroup_id']) if ratio(a['summary']['unresolved_population'])],
        diagnostic_only=any(a['diagnostic_only'] for a in allocations),
        ready_as_input_for_future_operator=all(a['ready_as_input_for_future_operator'] for a in allocations),
        aggregation='population_counts_after_within_subgroup_age_sex_allocation_no_RR_or_prevalence_premixing',
        source=source,provenance=provenance,review=review,forward_operator=None)
