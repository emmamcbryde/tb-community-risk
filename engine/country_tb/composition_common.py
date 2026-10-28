"""Versioned, offline primitives for baseline source and subgroup contracts."""
from copy import deepcopy
from decimal import Decimal

from engine.country_tb.calibration_common import reviewed, scientific_context, probability
from engine.country_tb.data_common import (
    DataValidationError, canonical, decimal, digest, numeric_text, profile_hash, seal,
)

VERSION = '1.0.0'
IMPLEMENTATION = 'country_tb_baseline_composition_1.0.0'
DEMOGRAPHY = 'country_tb.subgroup_demography.v1'
STATES = ('never_infected', 'recent', 'remote_only', 'any_tbi', 'prior_remote_plus_recent')


def require(condition, message):
    if not condition:
        raise DataValidationError(message)


def metadata(source, provenance, review, *, required_review=True):
    require(isinstance(source, dict) and bool(source.get('source_id')), 'Source with source_id required.')
    require(isinstance(provenance, dict) and bool(provenance), 'Explicit provenance required.')
    reviewed(review, required=required_review)


def artifact(schema_id, **fields):
    return seal({'schema_id': schema_id, 'contract_version': VERSION,
        'implementation_version': IMPLEMENTATION,
        'purpose': 'planning_and_sequencing_not_denying_care', **deepcopy(fields)})


def verify(value, family):
    require(isinstance(value, dict) and value.get('schema_id') == family,
            'Unexpected composition contract family.')
    require(value.get('contract_version') == VERSION and value.get('implementation_version') == IMPLEMENTATION,
            'Explicit composition contract migration required.')
    require(value.get('profile_hash') == profile_hash(value), 'Composition scientific hash mismatch.')
    return value


def identity(geography_id, geography_level, subgroup_id, baseline_year):
    require(isinstance(geography_id, str) and bool(geography_id), 'Explicit geography required.')
    require(geography_level in {'national', 'subnational', 'custom'}, 'Explicit geography level required.')
    require(isinstance(subgroup_id, str) and bool(subgroup_id), 'Stable subgroup ID required.')
    require(type(baseline_year) is int and baseline_year == 2024, 'This milestone requires baseline year 2024.')


def age_bounds(row):
    lo, hi = row.get('age_lower'), row.get('age_upper_exclusive')
    require(type(lo) is int and lo >= 0, 'Age lower bound must be a non-negative integer.')
    require(hi is None or (type(hi) is int and hi > lo), 'Invalid age upper-exclusive bound.')
    return lo, hi


def contains(outer, inner):
    a,b = age_bounds(outer)
    c,d = age_bounds(inner)
    return a <= c and (b is None or (d is not None and d <= b))


def overlaps(left, right):
    a,b = age_bounds(left)
    c,d = age_bounds(right)
    return (b is None or c < b) and (d is None or a < d)


def sex_overlaps(a, b):
    return a == b or 'both' in {a,b}


def validate_nonoverlap(rows):
    for i,left in enumerate(rows):
        require(left.get('sex') in {'female','male','both'}, 'Sex must be female, male or both.')
        age_bounds(left)
        for right in rows[i+1:]:
            require(not (sex_overlaps(left['sex'], right['sex']) and overlaps(left,right)),
                    'Overlapping age/sex records; both-sex rows cannot duplicate separate sexes.')


def uncertainty(value, prevalence):
    require(isinstance(value, dict) and value.get('kind') in {'exact','unspecified','interval'},
            'Explicit uncertainty object required; missing uncertainty is not zero.')
    result = deepcopy(value)
    if value['kind'] == 'interval':
        require({'lower','upper','method'} <= value.keys() and bool(value['method']),
                'Uncertainty interval requires bounds and method.')
        lo, hi = probability(value['lower']), probability(value['upper'])
        require(lo <= prevalence <= hi, 'Prevalence lies outside its supplied uncertainty interval.')
        result.update(lower=numeric_text(lo), upper=numeric_text(hi))
    return result


def counts(row, prevalence):
    require('numerator' in row and 'denominator' in row, 'Counts must be explicit, including jointly null.')
    x,n = row['numerator'], row['denominator']
    require((x is None) == (n is None), 'Numerator and denominator must be jointly known or missing.')
    if x is not None:
        require(type(x) is int and type(n) is int and 0 <= x <= n and n > 0, 'Invalid numerator/denominator.')
        tolerance = decimal(row.get('count_prevalence_tolerance','0'))
        require(abs(Decimal(x)/n-prevalence) <= tolerance or bool(row.get('weighted_count_discrepancy_rationale')),
                'Counts/prevalence disagree without reviewed rounding/weighting rationale.')


def population_demography(rows, *, geography_id, geography_level, population, source, provenance,
                          review, baseline_year=2024):
    """Explicit age/sex partition; never splits counts or invents a terminal age."""
    identity(geography_id, geography_level, 'demography', baseline_year)
    metadata(source, provenance, review)
    require(isinstance(rows,list) and bool(rows), 'Explicit age/sex population rows required.')
    original = deepcopy(rows)
    with scientific_context():
        result = []
        for row in rows:
            require(isinstance(row,dict) and {'age_lower','age_upper_exclusive','sex','population'} <= row.keys(),
                    'Incomplete age/sex population cell.')
            lo,hi = age_bounds(row)
            require(row.get('geography_id',geography_id)==geography_id and
                row.get('geography_level',geography_level)==geography_level and
                row.get('year',row.get('baseline_year',baseline_year))==baseline_year,
                'Population source rows cannot be silently relabelled to another geography or year.')
            result.append({'age_lower':lo, 'age_upper_exclusive':hi, 'sex':row['sex'],
                           'population':numeric_text(decimal(row['population']))})
        validate_nonoverlap(result)
        total = decimal(population)
        require(total > 0 and sum(decimal(r['population']) for r in result) == total,
                'Age/sex population cells must reconcile exactly to the positive denominator.')
        require(not ('both' in {r['sex'] for r in result} and len({r['sex'] for r in result})>1),
                'Both-sex population cannot be combined with separate sex population rows.')
        grids=[{(r['age_lower'],r['age_upper_exclusive']) for r in result if r['sex']==s}
               for s in sorted({r['sex'] for r in result})]
        require(all(g==grids[0] for g in grids),'Sex-specific age grids must align to derive population age margins without rebinning.')
        age = {}
        sex = {}
        for row in result:
            key = (row['age_lower'], row['age_upper_exclusive'])
            age[key] = age.get(key,Decimal(0))+decimal(row['population'])
            sex[row['sex']] = sex.get(row['sex'],Decimal(0))+decimal(row['population'])
        return artifact(DEMOGRAPHY, geography_id=geography_id, geography_level=geography_level,
            baseline_year=baseline_year, population=numeric_text(total),
            rows=sorted(result,key=canonical), source=source, provenance=provenance, review=review,
            age_distribution=[{'age_lower':a,'age_upper_exclusive':b,'population':numeric_text(n),
                'weight':numeric_text(n/total)} for (a,b),n in sorted(age.items(),key=lambda x:x[0][0])],
            sex_distribution=[{'sex':s,'population':numeric_text(n),'weight':numeric_text(n/total)} for s,n in sorted(sex.items())],
            population_scope='supplied_cells_are_the_declared_complete_population_partition',
            open_age_policy='unresolved_unless_direct_values_or_reviewed_within_cell_ages_supplied',
            receipt={'source_rows':original,'source_rows_sha256':digest(canonical(original))})


def validate_demography(value):
    if value.get('schema_id') == 'country_tb.demography.v1':
        from engine.country_tb.calibration_common import verify as verify_m2c
        verify_m2c(value, 'country_tb.demography.v1')
        return value
    verify(value, DEMOGRAPHY)
    require(digest(canonical(value['receipt']['source_rows']))==value['receipt']['source_rows_sha256'],
            'Original population source-row receipt checksum mismatch.')
    rebuilt = population_demography(value['receipt']['source_rows'], geography_id=value['geography_id'],
        geography_level=value['geography_level'], population=value['population'],source=value['source'],
        provenance=value['provenance'],review=value['review'],baseline_year=value['baseline_year'])
    require(rebuilt['profile_hash']==value['profile_hash'], 'Demographic margins do not reproduce from cells.')
    return value


def as_catalytic_demography(value):
    validate_demography(value)
    if value['schema_id'] == 'country_tb.demography.v1':
        return value
    from engine.country_tb.demography import demographic_profile
    rows = [{**r,'source_id':value['source']['source_id'],'source_vintage':value['profile_hash'],
        'geography_id':value['geography_id'],'geography_level':value['geography_level'],'year':2024,
        'estimate_status':'estimate','projection_scenario':'historical','unit':'persons',
        'population_variant':'explicit_user_population_partition'} for r in value['rows']]
    candidate={'snapshot_id':'subgroup-'+value['profile_hash'], 'rows':rows,
        'input_source_hash':value['profile_hash'], 'transformed_file_hash':digest(canonical(rows)),
        'manifest':{'source_id':value['source']['source_id'],'source_category':'user',
            'source_vintage':value['profile_hash'],'review_status':'reviewed','manifest_hash':value['profile_hash'],
            'schema_version':DEMOGRAPHY,'transformation_version':IMPLEMENTATION,
            'transformations':['Explicit subgroup cells adapted to unchanged Milestone 2C demography API'],
            'source_url':value['source'].get('url','urn:reviewed-user-demography'),
            'import_date':'not_applicable_pure_contract','terms':{'status':'user_supplied'},
            'citation':value['source']['source_id'],'projection_boundary':None}}
    return demographic_profile(value['geography_id'],[candidate])


def profile_summary(rows):
    """Component-specific coverage. Unknown subset/recency is never zero-filled."""
    with scientific_context():
        total = sum(decimal(r['population']) for r in rows)
        components = {}
        for field in STATES:
            covered = sum(decimal(r['population']) for r in rows if r.get(field) is not None)
            count = sum((decimal(r['population'])*decimal(r[field]) for r in rows if r.get(field) is not None),Decimal(0))
            components[field] = {'estimated_population':numeric_text(covered),
                'unresolved_population':numeric_text(total-covered), 'estimated_count':numeric_text(count),
                'complete_count':numeric_text(count) if covered==total else None,
                'prevalence_covered':numeric_text(count/covered) if covered else None,
                'prevalence_full':numeric_text(count/total) if total and covered==total else None}
        return {'population':numeric_text(total),'components':components,
            'estimated_population':components['any_tbi']['estimated_population'],
            'unresolved_population':components['any_tbi']['unresolved_population'],
            'unresolved_population_proportion':numeric_text(decimal(components['any_tbi']['unresolved_population'])/total) if total else None,
            'total_tbi_complete':components['any_tbi']['complete_count'] is not None,
            'recency_complete':all(components[k]['complete_count'] is not None for k in ('recent','remote_only','never_infected'))}
