"""Source-preserving direct TBI profiles and explicitly selected hybrid recency."""
from copy import deepcopy
from decimal import Decimal

from engine.country_tb.calibration_common import probability, scientific_context
from engine.country_tb.composition_common import (
    artifact, verify, require, metadata, age_bounds, contains, overlaps, sex_overlaps,
    validate_nonoverlap, counts, uncertainty,
)
from engine.country_tb.data_common import canonical, decimal, numeric_text, digest

DIRECT = 'country_tb.direct_tbi_profile.v1'
UNAVAILABLE_RECENCY = 'recency_unavailable_v1'
CONDITIONAL_RECENCY = 'model_conditional_recency_allocation_v1'
AGE_POLICIES = ('exact_single_year_values_v1','piecewise_constant_within_band_v1',
                'supplied_within_band_distribution_v1')


def _values(row):
    fields=('recent','remote_only','never_infected')
    supplied=[k for k in ('recent','remote_only') if row.get(k) is not None]
    require(not supplied or (len(supplied)==2 and row.get('never_infected') is not None),
            'Direct recent, remote-only and never-infected must all be supplied together.')
    if supplied:
        values={k:probability(row[k]) for k in fields}
        require(sum(values.values())==1,'Direct mutually exclusive infection states must sum exactly to one.')
        total=values['recent']+values['remote_only']
        require(row.get('any_tbi') is None or probability(row['any_tbi'])==total,'Direct total TBI disagrees with recent plus remote-only.')
        prior=None if row.get('prior_remote_plus_recent') is None else probability(row['prior_remote_plus_recent'])
        require(prior is None or prior<=values['recent'],'Prior remote plus recent is a subset of recent, never an additional state.')
        return {**{k:numeric_text(v) for k,v in values.items()},'any_tbi':numeric_text(total),
                'prior_remote_plus_recent':numeric_text(prior) if prior is not None else None}
    require(row.get('any_tbi') is not None,'A direct profile needs total TBI or all three exclusive states.')
    require(row.get('prior_remote_plus_recent') is None,'A remote-plus-recent subset cannot be supplied without direct recent states.')
    total=probability(row['any_tbi'])
    require(row.get('never_infected') is None or probability(row['never_infected'])==1-total,
            'Supplied never infected disagrees with total TBI.')
    return {'any_tbi':numeric_text(total),'never_infected':numeric_text(1-total),
            'recent':None,'remote_only':None,'prior_remote_plus_recent':None}


def direct_tbi_profile(records, *, source, provenance, review, value_target='latent_any_tbi'):
    metadata(source,provenance,review)
    require(value_target=='latent_any_tbi','Direct overrides require latent TBI values; use assay contracts for test positivity anchors.')
    require(isinstance(records,list) and bool(records),'Direct age/sex records required.')
    originals=deepcopy(records)
    rows=[]
    with scientific_context():
        for original in records:
            require(isinstance(original,dict) and {'record_id','age_lower','age_upper_exclusive','sex','age_allocation_policy'}<=original.keys(),
                    'Incomplete direct-profile record.')
            row=deepcopy(original)
            require(isinstance(row['record_id'],str) and row['record_id'],'Stable direct record ID required.')
            lo,hi=age_bounds(row)
            require(row.get('year',2024)==2024,'Direct profile year must be the 2024 baseline.')
            policy=row['age_allocation_policy']
            require(policy in AGE_POLICIES,'Explicit direct-profile age allocation policy required; no interpolation.')
            if policy=='exact_single_year_values_v1':
                require(hi==lo+1,'Single-year values require an exact one-year cell; an open tail is not an exact age.')
            if policy=='supplied_within_band_distribution_v1':
                children=row.get('within_band_values')
                require(isinstance(children,list) and bool(children),'Supplied within-band values required.')
                normalized=[]
                for child in children:
                    require(isinstance(child,dict) and {'age_lower','age_upper_exclusive','sex'}<=child.keys(),
                            'Incomplete supplied within-band age/sex value.')
                    require(contains(row,child) and (row['sex']=='both' or row['sex']==child['sex']),
                            'Within-band values fall outside their parent age/sex band.')
                    normalized.append({**deepcopy(child),**_values(child)})
                validate_nonoverlap(normalized)
                row['within_band_values']=sorted(normalized,key=canonical)
                if row.get('any_tbi') is not None:
                    row['any_tbi']=numeric_text(probability(row['any_tbi']))
            else:
                require('within_band_values' not in row,'Within-band values require the explicitly selected supplied-distribution policy.')
                row.update(_values(row))
            if 'numerator' in row or 'denominator' in row:
                require(row.get('any_tbi') is not None,'A count observation needs its reported total prevalence.')
                counts(row,probability(row['any_tbi']))
            if row.get('any_tbi') is not None and 'uncertainty' in row:
                row['uncertainty']=uncertainty(row['uncertainty'],probability(row['any_tbi']))
            rows.append(row)
        require(len({r['record_id'] for r in rows})==len(rows),'Direct record IDs must be unique.')
        validate_nonoverlap(rows)
    return artifact(DIRECT,rows=sorted(rows,key=canonical),value_target=value_target,source=source,
        provenance=provenance,review=review,source_preservation='all_original_values_and_metadata_retained_unchanged',
        receipt={'source_records':originals,'source_records_sha256':digest(canonical(originals))})


def allocate_direct_profile(profile, demography, *, recency_policy=UNAVAILABLE_RECENCY, model_profile=None):
    verify(profile,DIRECT)
    require(profile['receipt']['source_records_sha256']==digest(canonical(profile['receipt']['source_records'])),
            'Original direct-profile source receipt checksum mismatch.')
    rebuilt=direct_tbi_profile(profile['receipt']['source_records'],source=profile['source'],
        provenance=profile['provenance'],review=profile['review'])
    require(rebuilt['profile_hash']==profile['profile_hash'],'Direct profile does not reproduce from original source records.')
    require(recency_policy in {UNAVAILABLE_RECENCY,CONDITIONAL_RECENCY},'Unknown direct-profile recency policy.')
    if model_profile is not None:
        from engine.country_tb.tbi_baseline import restore_tbi_profile
        restore_tbi_profile(model_profile)
        require(model_profile['geography_id']==demography['geography_id'] and model_profile['baseline_year']==2024,
                'Comparison model must match the analysis geography and baseline year.')
        expected={(r['sex'],r['age_lower'],r['age_upper_exclusive']):decimal(r['population']) for r in demography['rows']}
        actual={(r['sex'],r['age_lower'],r['age_upper_exclusive']):decimal(r['population']) for r in model_profile['rows']}
        require(expected==actual,'Comparison model population/age/sex cells do not match the supplied population.')
    leaves=[]
    for parent in profile['rows']:
        for leaf in parent.get('within_band_values',[parent]):
            leaves.append((parent,leaf))
    require(recency_policy!=CONDITIONAL_RECENCY or model_profile is not None or all(r.get('recent') is not None for _,r in leaves),
            'Explicit model-conditional recency allocation requires a comparison model profile.')
    model={} if model_profile is None else {(r['sex'],r['age_lower'],r['age_upper_exclusive']):r for r in model_profile['rows']}
    output,gaps=[],[]
    used=set()
    with scientific_context():
        for cell in demography['rows']:
            eligible=[(p,r) for p,r in leaves if contains(r,cell) and (r['sex']==cell['sex'] or r['sex']=='both')]
            # A finite demographic count cell is never split without supplied
            # within-cell demographic weights; direct prevalence is not weights.
            crossed=[r for _,r in leaves if overlaps(r,cell) and sex_overlaps(r['sex'],cell['sex']) and not contains(r,cell)]
            require(not crossed,'A direct profile boundary cuts a demographic cell; supply a compatible population age grid.')
            require(len(eligible)<=1,'Overlapping direct allocation for an age/sex cell.')
            row={k:cell[k] for k in ('age_lower','age_upper_exclusive','sex','population')}
            row.update({k:None for k in ('any_tbi','never_infected','recent','remote_only','prior_remote_plus_recent')})
            row.update(most_recent_infection_time=None,recent_time_summary=None,remote_only_time_summary=None)
            if not eligible:
                row.update(status='unresolved_direct_profile_gap',recency_composition='unavailable',source_record_id=None)
                gaps.append({k:row[k] for k in ('age_lower','age_upper_exclusive','sex','population')})
            else:
                parent,leaf=eligible[0]
                used.add(parent['record_id'])
                row.update(_values(leaf),status='direct_values',source_record_id=parent['record_id'],
                    age_allocation_policy=parent['age_allocation_policy'],
                    recency_composition='user_defined' if leaf.get('recent') is not None else 'unavailable')
                if row['recent'] is None and recency_policy==CONDITIONAL_RECENCY:
                    m=model[(cell['sex'],cell['age_lower'],cell['age_upper_exclusive'])]
                    p=decimal(row['any_tbi'])
                    if m['any_tbi'] is None and decimal(cell['population'])==0:
                        row['status']='hybrid_recency_unavailable_empty_model_cell'
                        output.append(row)
                        continue
                    require(m['any_tbi'] is not None,'Conditional recency is unavailable for an unresolved model age cell.')
                    mt=decimal(m['any_tbi'])
                    require(mt>0 or p==0,'Positive user TBI cannot use a zero-probability model conditional distribution.')
                    recent=p*decimal(m['recent'])/mt if mt else Decimal(0)
                    prior=p*decimal(m['prior_remote_plus_recent'])/mt if mt else Decimal(0)
                    row.update(recent=numeric_text(recent),remote_only=numeric_text(p-recent),
                        prior_remote_plus_recent=numeric_text(prior),
                        recency_composition='model_derived',status='hybrid_model_conditional_recency',
                        recency_allocation_label='hybrid_model_derived_not_observed_recent_fraction')
            output.append(row)
        require(used=={r['record_id'] for r in profile['rows']},'Unused direct-profile record outside the population age/sex scope.')
        for parent in profile['rows']:
            if parent['age_allocation_policy']=='supplied_within_band_distribution_v1' and parent.get('any_tbi') is not None:
                selected=[r for r in output if r['source_record_id']==parent['record_id']]
                pop=sum(decimal(r['population']) for r in selected)
                expected_pop=sum(decimal(r['population']) for r in output if contains(parent,r) and
                    (parent['sex']=='both' or parent['sex']==r['sex']))
                require(pop==expected_pop,'Supplied within-band distribution has a population coverage gap.')
                require(all(r['any_tbi'] is not None for r in selected),'Supplied within-band distribution has unresolved cells.')
                achieved=sum(decimal(r['population'])*decimal(r['any_tbi']) for r in selected)
                require(not pop or abs(achieved/pop-decimal(parent['any_tbi']))<=Decimal('1e-30'),
                        'Supplied within-band values disagree with the supplied parent-band prevalence.')
    return sorted(output,key=canonical),gaps
