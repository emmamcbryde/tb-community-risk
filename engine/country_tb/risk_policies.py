"""Overlap, feasible residual association and explicit RR combination rules."""
from copy import deepcopy
from decimal import Decimal
from fractions import Fraction
from math import prod

from engine.country_tb.calibration_common import scientific_context
from engine.country_tb.risk_factor_catalog import FACTOR_IDS
from engine.country_tb.risk_common import (
    artifact,verify,require,metadata,ratio,text_number,amount,canonical,age_bounds,
    validate_disjoint,matches,overlaps,sex_overlaps,
)

OVERLAP='country_tb.risk_factor_overlap_policy.v1'
ASSOCIATION='country_tb.tbi_any_risk_association.v1'
MULTIMORBIDITY='country_tb.multimorbidity_progression_policy.v1'
INDEPENDENT='conditionally_independent_risk_factors_within_age_sex_stratum_v1'
JOINT='reviewed_joint_risk_factor_distribution_v1'
HIGHEST='highest_applicable_reviewed_rr_v1'
MULTIPLICATIVE='multiplicative_rr_sensitivity_v1'
ADDITIVE='additive_excess_rr_sensitivity_v1'


def overlap_policy(*,source,provenance,review,method=INDEPENDENT,joint_rows=None):
    metadata(source,provenance,review)
    require(method in {INDEPENDENT,JOINT},'Unknown factor-overlap policy; pairwise correlations are not inferred.')
    rows=deepcopy(joint_rows or [])
    require(method==JOINT or not rows,'Independent overlap cannot also consume joint-distribution rows.')
    if method==JOINT:
        require(rows,'Reviewed joint overlap needs full pattern distributions.');validate_disjoint(rows)
        for row in rows:
            require(isinstance(row.get('patterns'),list) and row['patterns'],'Full reviewed pattern distribution required.')
            seen=set();mass=Fraction(0)
            for pattern in row['patterns']:
                require(isinstance(pattern,dict) and set(pattern)=={'factor_ids','probability'},'Joint patterns need exactly factor IDs and probability.')
                ids=pattern['factor_ids'];require(isinstance(ids,list) and len(set(ids))==len(ids) and all(isinstance(x,str) for x in ids),'Joint pattern IDs must be unique strings.')
                require(set(ids)<=set(FACTOR_IDS),'Joint patterns cannot contain recency, contact or unknown factors.')
                key=tuple(sorted(ids));require(key not in seen,'Duplicated joint factor pattern.');seen.add(key)
                p=ratio(pattern['probability'],probability=True);mass+=p
                pattern.update(factor_ids=list(key),probability=amount(p))
            require(mass==1,'Reviewed joint pattern probabilities must sum exactly to one.')
            require(len(row['patterns'])<=2048,'Eleven-factor joint contract is bounded at 2048 patterns.')
            row['patterns']=sorted(row['patterns'],key=canonical)
    return artifact(OVERLAP,method=method,rows=sorted(rows,key=canonical),source=source,provenance=provenance,review=review,
        assumptions='independent_factors_given_age_and_sex_before_TBI_coupling' if method==INDEPENDENT else 'complete_reviewed_joint_pattern_table',
        pairwise_correlations_invented=False)


def association_policy(*,global_control,source,provenance,review,overrides=None):
    metadata(source,provenance,review);c=ratio(global_control)
    require(-1<=c<=1,'Normalized residual association control must lie in [-1,1].')
    rows=deepcopy(overrides or []);ids=set()
    for row in rows:
        require(isinstance(row,dict) and {'control_id','scope','sex','control','source','provenance','review'}<=row.keys(),'Incomplete association override.')
        require(row['scope'] in {'sex','age_band','exact_age_sex'} and row['sex'] in {'female','male','both'},'Unknown association override scope.')
        require(row['control_id'] and row['control_id'] not in ids,'Association control IDs must be unique.');ids.add(row['control_id'])
        metadata(row['source'],row['provenance'],row['review'])
        v=ratio(row['control']);require(-1<=v<=1,'Association override outside [-1,1].');row['control']=text_number(v)
        if row['scope']=='sex':
            require(row['sex']!='both' and 'age_lower' not in row and 'age_upper_exclusive' not in row,'Sex-only controls must not supply ages.')
        else:
            require({'age_lower','age_upper_exclusive'}<=row.keys(),'Association age bounds required.');lo,hi=age_bounds(row)
            if row['scope']=='exact_age_sex':require(hi==lo+1,'Exact-age controls need a single-year interval.')
    for i,a in enumerate(rows):
        for b in rows[i+1:]:
            if a['scope']!=b['scope'] or a['sex']!=b['sex']:continue
            require(a['scope']!='sex' and not overlaps(a,b),'Ambiguous association overrides at the same precedence.')
    return artifact(ASSOCIATION,global_control=text_number(c),rows=sorted(rows,key=canonical),
        precedence=['exact_age_sex','exact_age_both','age_band_sex','age_band_both','sex','global'],
        control_interpretation='normalized_feasible_residual_association_not_Pearson_correlation',
        conditioning='within_source_age_sex_subgroup_stratum',source=source,provenance=provenance,review=review,
        factor_composition_policy='preserve_factor_pattern_distribution_given_any_risk_in_both_TBI_groups',
        recency_policy='preserve_recent_remote_and_prior_remote_subset_given_any_TBI')


def control_for_cell(policy,cell):
    verify(policy,ASSOCIATION);candidates=[]
    for row in policy['rows']:
        if row['scope']=='sex':
            if row['sex']==cell['sex']:candidates.append((1,row))
        else:
            if overlaps(row,cell) and sex_overlaps(row['sex'],cell['sex']):
                require(matches(row,cell),'Association boundary cuts a population cell; an aligned age/sex grid is required.')
                priority=(5 if row['scope']=='exact_age_sex' else 3)+(row['sex']!='both')
                candidates.append((priority,row))
    if not candidates:return ratio(policy['global_control']),{'control_id':'global','scope':'global','policy_hash':policy['profile_hash']}
    _,row=max(candidates,key=lambda x:x[0])
    return ratio(row['control']),{'control_id':row['control_id'],'scope':row['scope'],'source':row['source'],
        'provenance':row['provenance'],'review':row['review'],'policy_hash':policy['profile_hash']}


def phi(p,q,joint):
    variance=p*(1-p)*q*(1-q)
    if variance==0:return None,'undefined_degenerate_marginal'
    covariance=joint-p*q
    if covariance*covariance==variance:return ('1' if covariance>0 else '-1'),None
    with scientific_context():
        value=(Decimal(covariance.numerator)/Decimal(covariance.denominator))/(
            Decimal(variance.numerator)/Decimal(variance.denominator)).sqrt()
        return text_number(value),None


def feasible_association(p,q,control):
    p,q,c=ratio(p,probability=True),ratio(q,probability=True),ratio(control)
    require(-1<=c<=1,'Residual association must lie in [-1,1].')
    lower,upper=max(Fraction(0),p+q-1),min(p,q);independent=p*q
    joint=independent+c*((independent-lower) if c<0 else (upper-independent))
    cells={'tbi_and_risk':joint,'tbi_without_risk':p-joint,'risk_without_tbi':q-joint,'neither':1-p-q+joint}
    require(all(x>=0 for x in cells.values()) and sum(cells.values())==1,'Invalid feasible association table; no probability repair.')
    achieved,reason=phi(p,q,joint)
    return {'tbi_prevalence':amount(p),'any_risk_prevalence':amount(q),'lower_bound':amount(lower),'upper_bound':amount(upper),
        'independence_joint':amount(independent),'joint_prevalence':amount(joint),'normalized_association_input':text_number(c),
        'achieved_phi':achieved,'phi_unavailable_reason':reason,'cells':{k:amount(v) for k,v in cells.items()},
        'control_is_achieved_correlation':False}


def multimorbidity_policy(*,source,provenance,review,method=HIGHEST,effect_use='reviewed_effects_only_v1',cap_policy='none',cap=None):
    metadata(source,provenance,review)
    require(method in {HIGHEST,MULTIPLICATIVE,ADDITIVE},'Unknown RR combination policy.')
    require(effect_use in {'reviewed_effects_only_v1','diagnostic_provisional_adult_effects_v1'},'Unknown effect-review policy.')
    require(cap_policy in {'none','explicit_reviewed_cap_v1'},'Explicit cap behavior required.')
    require((cap_policy=='none' and cap is None) or (cap_policy!='none' and cap is not None and ratio(cap)>=1),'Invalid RR cap; no arbitrary implicit cap.')
    require(method!=HIGHEST or cap is None,'Highest-RR central allocation is not capped.')
    return artifact(MULTIMORBIDITY,method=method,effect_use=effect_use,cap_policy=cap_policy,
        cap=text_number(cap) if cap is not None else None,calculation_only_sensitivity=method!=HIGHEST,
        tie_policy='retain_all_tied_highest_contributors_sorted_by_stable_factor_ID',
        recent_infection_included=False,source=source,provenance=provenance,review=review)


def combine_risk_factors(factors_present,effects,policy):
    """Auditable combination contract for a supplied membership pattern; no disease rates."""
    verify(policy,MULTIMORBIDITY)
    require(isinstance(factors_present,list) and len(set(factors_present))==len(factors_present),'Factor memberships must be a unique list.')
    require(set(factors_present)<=effects.keys() and set(factors_present)<=set(FACTOR_IDS),'Unknown factor membership; recency/contact are not ordinary factors.')
    active={k:ratio(effects[k]['rr']) for k in factors_present if effects[k]['rr'] is not None and
        (effects[k]['state']=='reviewed' or (policy['effect_use']=='diagnostic_provisional_adult_effects_v1' and
            effects[k]['state']=='diagnostic_provisional_pending_contemporary_review'))}
    if len(active)!=len(factors_present):
        return {'effective_rr':None,'all_factors_present':sorted(factors_present),'responsible_factors':[],
                'status':'unresolved_effect_review','effect_provenance':{k:effects[k] for k in sorted(factors_present)}}
    method=policy['method'];highest=max([Fraction(1),*active.values()])
    value=highest if method==HIGHEST else prod(active.values(),start=Fraction(1)) if method==MULTIPLICATIVE else 1+sum((r-1 for r in active.values()),Fraction(0))
    if policy['cap'] is not None:value=min(value,ratio(policy['cap']))
    return {'effective_rr':amount(value),'all_factors_present':sorted(factors_present),
        'responsible_factors':sorted(k for k,r in active.items() if r==highest) if method==HIGHEST else sorted(active),
        'tied_highest_contributors':sorted(k for k,r in active.items() if r==highest),
        'status':'diagnostic_provisional' if any('diagnostic' in effects[k]['state'] for k in factors_present) else 'reviewed_or_no_factors',
        'calculation_only_sensitivity':policy['calculation_only_sensitivity'],'cap_policy':policy['cap_policy'],
        'effect_provenance':{k:effects[k] for k in sorted(factors_present)}}
