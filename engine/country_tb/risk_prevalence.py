"""Definition-aware risk prevalence and explicit age/sex allocation policies."""
from copy import deepcopy
from fractions import Fraction

from engine.country_tb.risk_common import (
    artifact,verify,require,metadata,ratio,text_number,amount,validate_disjoint,validate_rows,
    age_bounds,matches,overlaps,sex_overlaps,canonical,digest,cell_key,
)
from engine.country_tb.risk_factor_catalog import catalogue_records

PROFILE='country_tb.risk_factor_prevalence_profile.v1'
MAPPING='country_tb.risk_factor_definition_mapping.v1'
STATES=('value','missing','not_applicable','excluded','not_reviewed','unavailable_due_to_definition_mismatch')
DIRECT='direct_age_sex_values_v1'
CONSTANT='piecewise_constant_within_supplied_band_v1'
UNIFORM='uniform_within_eligible_strata_v1'
PATTERN='reviewed_age_sex_pattern_scaled_to_total_v1'
POLICIES=(DIRECT,CONSTANT,UNIFORM,PATTERN)


def definition_mapping(catalog,*,factor_id,supplied_definition,source,provenance,review,rationale):
    factors=catalogue_records(catalog)
    require(factor_id in factors,'Unknown risk factor; recency/contact are not ordinary factors.')
    metadata(source,provenance,review)
    require(isinstance(supplied_definition,str) and supplied_definition and isinstance(rationale,str) and rationale,
            'A reviewed definition mapping needs its original definition and rationale.')
    return artifact(MAPPING,factor_id=factor_id,supplied_exposure_definition=supplied_definition,
        catalogue_definition=factors[factor_id]['exposure_definition'],catalogue_record_hash=factors[factor_id]['profile_hash'],
        compatibility_status='reviewed_mapping',mapping_policy='explicit_reviewed_definition_transport_v1',
        rationale=rationale,source=source,provenance=provenance,review=review,
        warning='Mapped exposure is an explicit assumption; it does not approve the legacy effect for contemporary use.')


def risk_factor_prevalence_profile(records,*,catalog,factor_ids,subgroup_id,geography_id,
    source,provenance,review):
    factors=catalogue_records(catalog);metadata(source,provenance,review)
    require(isinstance(factor_ids,list) and factor_ids and len(set(factor_ids))==len(factor_ids) and set(factor_ids)<=factors.keys(),
            'An explicit unique catalogue factor scope is required; recency/contact are excluded.')
    require(isinstance(subgroup_id,str) and subgroup_id and isinstance(geography_id,str) and geography_id,'Subgroup/geography identity required.')
    require(isinstance(records,list),'Explicit prevalence records required, including an empty missing-data list.')
    originals=deepcopy(records);rows=[];ids=set()
    for original in records:
        fields={'record_id','factor_id','age_lower','age_upper_exclusive','sex','state','prevalence',
            'allocation_policy','supplied_exposure_definition','definition_mapping','pattern','source','provenance','review','reason'}
        require(isinstance(original,dict) and fields<=original.keys(),'Incomplete risk prevalence contract; missing is explicit, not zero.')
        row=deepcopy(original);validate_rows([row])
        require(isinstance(row['record_id'],str) and row['record_id'] and row['record_id'] not in ids,'Unique risk prevalence record ID required.')
        ids.add(row['record_id']);factor=row['factor_id']
        require(factor in factor_ids,'Prevalence factor is outside the explicit analysis factor scope.')
        metadata(row['source'],row['provenance'],row['review'],required_review=False)
        require(row['state'] in STATES and row['allocation_policy'] in POLICIES,'Unknown prevalence state or allocation policy; no implicit default.')
        p=None if row['prevalence'] is None else ratio(row['prevalence'],probability=True)
        if row['state']=='value':require(p is not None,'A supplied value requires a prevalence; zero is valid.')
        if row['state'] in {'missing','not_applicable','excluded'}:
            require(p is None,'Missing/not-applicable/excluded is not a numerical zero.')
        if row['state'] in {'not_applicable','excluded'}:require(row['reason'],'Explicit exclusion/non-applicability needs a reason.')
        lo,hi=age_bounds(row)
        if row['allocation_policy']==DIRECT:require(hi==lo+1,'Direct age/sex values need a single-year cell; open age is not exact age.')
        if row['allocation_policy'] in {UNIFORM,PATTERN}:
            require(lo==15 and hi is None,'Overall allocation policies refer to the explicitly eligible adult population 15+.')
        expected=factors[factor]['exposure_definition'];mapping=row['definition_mapping']
        status='compatible_exact_definition' if row['supplied_exposure_definition']==expected else 'definition_mismatch'
        if mapping is not None:
            verify(mapping,MAPPING)
            require(mapping['factor_id']==factor and mapping['supplied_exposure_definition']==row['supplied_exposure_definition'] and
                mapping['catalogue_record_hash']==factors[factor]['profile_hash'],'Definition mapping does not match this supplied/catalogue definition.')
            metadata(mapping['source'],mapping['provenance'],mapping['review']);status='reviewed_mapping'
        state=row['state']
        if p is not None and status=='definition_mismatch':state='unavailable_due_to_definition_mismatch'
        if p is not None and row['review']['status']!='reviewed':state='not_reviewed'
        if row['allocation_policy']==PATTERN:
            pattern=row['pattern']
            require(isinstance(pattern,dict) and {'rows','source','provenance','review'}<=pattern.keys(),'A reviewed age/sex pattern is required.')
            metadata(pattern['source'],pattern['provenance'],pattern['review']);validate_disjoint(pattern['rows'])
            for cell in pattern['rows']:
                require(cell['age_lower']>=15 and (row['sex']=='both' or cell['sex']==row['sex']),'Pattern lies outside its eligible adult/sex scope.')
                require('relative_weight' in cell and ratio(cell['relative_weight'])>=0,'Nonnegative pattern weights required.')
                cell['relative_weight']=text_number(cell['relative_weight'])
            pattern['rows']=sorted(pattern['rows'],key=canonical)
        else:require(row['pattern'] is None,'Pattern input requires the explicitly selected pattern allocation policy.')
        row.update(prevalence=text_number(p) if p is not None else None,effective_state=state,
            catalogue_definition=expected,catalogue_record_hash=factors[factor]['profile_hash'],compatibility_status=status,
            mapping_policy='exact_definition_identity_v1' if status=='compatible_exact_definition' else
                'explicit_reviewed_definition_transport_v1' if status=='reviewed_mapping' else 'unavailable_no_implicit_mapping',
            warnings=[mapping['warning']] if mapping is not None else [])
        rows.append(row)
    for factor in factor_ids:validate_disjoint([r for r in rows if r['factor_id']==factor])
    return artifact(PROFILE,subgroup_id=subgroup_id,geography_id=geography_id,baseline_year=2024,
        catalogue_hash=catalog['profile_hash'],factor_ids=sorted(factor_ids),rows=sorted(rows,key=canonical),
        excluded_from_analysis_scope=sorted(set(factors)-set(factor_ids)),
        union_scope='explicit_analysis_scope_less_explicit_user_exclusions_not_all_possible_progression_factors',
        source=source,provenance=provenance,review=review,
        receipt={'source_records':originals,'source_records_sha256':digest(canonical(originals))})


def resolve_risk_prevalence(profile,catalog,cells):
    verify(profile,PROFILE);require(profile['catalogue_hash']==catalog['profile_hash'],'Prevalence catalogue reference mismatch.')
    catalogue=catalogue_records(catalog)
    receipt=profile['receipt'];require(digest(canonical(receipt['source_records']))==receipt['source_records_sha256'],'Risk prevalence source receipt mismatch.')
    rebuilt=risk_factor_prevalence_profile(receipt['source_records'],catalog=catalog,factor_ids=profile['factor_ids'],
        subgroup_id=profile['subgroup_id'],geography_id=profile['geography_id'],source=profile['source'],provenance=profile['provenance'],review=profile['review'])
    require(rebuilt['profile_hash']==profile['profile_hash'],'Risk profile does not reproduce from its original records.')
    scaled={};diagnostics=[]
    for record in profile['rows']:
        if record['allocation_policy']!=PATTERN or record['effective_state']!='value':continue
        selected=[c for c in cells if matches(record,c)]
        pop=sum((ratio(c['population']) for c in selected),Fraction(0));require(pop>0,'Pattern target has no eligible population.')
        values=[]
        for cell in selected:
            rows=[r for r in record['pattern']['rows'] if matches(r,cell)]
            require(len(rows)==1,'Reviewed pattern must cover every eligible age/sex cell, including open ages; missing weights are not zero.')
            values.append((cell,ratio(rows[0]['relative_weight'])))
        total_weight=sum((ratio(c['population'])*w for c,w in values),Fraction(0));target=ratio(record['prevalence'])
        require(total_weight>0 or target==0,'Positive prevalence cannot be allocated from an all-zero pattern.')
        scale=target*pop/total_weight if total_weight else Fraction(0)
        allocated={cell_key(c):w*scale for c,w in values}
        require(all(0<=p<=1 for p in allocated.values()),'Scaled age/sex pattern is infeasible; no prevalence cap or clipping.')
        achieved=sum((ratio(c['population'])*allocated[cell_key(c)] for c,_ in values),Fraction(0))/pop
        require(achieved==target,'Scaled pattern failed exact marginal reconciliation.')
        scaled[record['record_id']]=allocated
        diagnostics.append({'record_id':record['record_id'],'population':amount(pop),'requested_prevalence':amount(target),
                            'achieved_prevalence':amount(achieved),'scale':amount(scale)})
    output=[]
    for cell in cells:
        factors=[]
        for factor in profile['factor_ids']:
            candidates=[r for r in profile['rows'] if r['factor_id']==factor and matches(r,cell)]
            cut=[r for r in profile['rows'] if r['factor_id']==factor and overlaps(r,cell) and sex_overlaps(r['sex'],cell['sex']) and not matches(r,cell)]
            require(not cut,'Prevalence boundary cuts a demographic cell; supply reviewed population age/sex weights or an aligned grid.')
            require(len(candidates)<=1,'Ambiguous factor prevalence allocation.')
            r=candidates[0] if candidates else None
            state=r['effective_state'] if r else ('not_applicable_pending_pediatric_review' if cell['age_lower']<15 else 'missing')
            p=None
            if state=='value':p=scaled[r['record_id']].get(cell_key(cell),Fraction(0)) if r['allocation_policy']==PATTERN else ratio(r['prevalence'])
            factors.append({'factor_id':factor,'state':state,'prevalence':amount(p) if p is not None else None,
                'source_record_id':r['record_id'] if r else None,'allocation_policy':r['allocation_policy'] if r else None,
                'supplied_exposure_definition':r['supplied_exposure_definition'] if r else None,
                'catalogue_definition':catalogue[factor]['exposure_definition'],
                'compatibility_status':r['compatibility_status'] if r else 'unavailable_missing_input',
                'mapping_policy':r['mapping_policy'] if r else None,'review':r['review'] if r else None,
                'warnings':r['warnings'] if r else []})
        output.append({'cell_key':cell_key(cell),'factors':factors})
    return output,diagnostics
