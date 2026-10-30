"""Offline, immutable empirical observations. No clinical parameter defaults.

Source cells and transformed inputs are separate, independently checksummed
files. Internal extraction review never means external scientific validation.
"""
from copy import deepcopy
from decimal import Decimal, localcontext
import json
from pathlib import Path
import re

from engine.country_tb.data_common import (
    canonical,digest,profile_hash,seal,decimal,numeric_text,contained_file,publish_immutable,REPO_ROOT,
)
from engine.country_tb.composition_common import require

IMPLEMENTATION='country_tb_empirical_validation_1.0.0'
ROOT=REPO_ROOT/'data/country_tb_validation/snapshots'
OBSERVATIONS='country_tb.empirical_tbi_validation_inputs.v1'
SPLIT='country_tb.empirical_validation_split.v1'
COMPARISON='country_tb.empirical_tbi_comparison.v1'
PARAMETERS='country_tb.review_parameter_set.v1'
REGISTRY='country_tb.scientific_evidence_registry.v1'
DISCREPANCY='country_tb.model_discrepancy.v1'
RESPONSE='country_tb.scientific_reviewer_response.v1'


def artifact(family,**fields):
    return seal(dict(schema_id=family,contract_version='1.0.0',implementation_version=IMPLEMENTATION,
        production_available=False,externally_validated=False,
        purpose='planning_and_sequencing_not_denying_care',**deepcopy(fields)))


def verify(value,family):
    require(value.get('schema_id')==family and value.get('implementation_version')==IMPLEMENTATION and
        value.get('contract_version')=='1.0.0','Explicit validation contract migration required.')
    require(value.get('profile_hash')==profile_hash(value),'Validation scientific hash mismatch.')
    require(value.get('production_available') is False and value.get('externally_validated') is False,
        'Internal synthesis cannot activate production or claim external validation.')
    return value


def probability(value):
    p=decimal(value);require(p<=1,'Probability outside [0,1]; no clipping.');return p


def normalize_observation(original):
    required={'record_id','source_id','citation','url','access_date','geography_id','geography_level',
        'subgroup_id','participant_pool','membership_group','survey_year','survey_period','age_lower',
        'age_upper_exclusive','sex','view','numerator','denominator','reported_prevalence_percent',
        'uncertainty','assay','threshold','sampling_frame','survey_design','migration',
        'review_status','transformation_notes','source_locator','source_cells','compatibility'}
    require(isinstance(original,dict) and required<=original.keys(),'Incomplete source observation; missing fields must be explicit.')
    r=deepcopy(original)
    for k in ('record_id','source_id','citation','url','geography_id','subgroup_id','participant_pool','membership_group','sampling_frame','source_locator'):
        require(isinstance(r[k],str) and r[k],f'Nonempty {k} required.')
    require(re.fullmatch(r'\d{4}-\d{2}-\d{2}',r['access_date']) is not None,'Access date required.')
    require(r['geography_level'] in {'national','subnational','custom'} and r['sex'] in {'both','female','male','unknown'},'Geography/sex semantics required.')
    require(type(r['survey_year']) is int and r['survey_year']<=2024,'Explicit historical survey year required.')
    lo,hi=r['age_lower'],r['age_upper_exclusive']
    require(type(lo) is int and lo>=0 and (hi is None or type(hi) is int and hi>lo),'Invalid source ages.')
    require(r['view'] in {'age','sex','overall','assay','visa'},'Unknown source partition view.')
    require(r['assay'] in {'QFT-Plus','QFT-GIT','TST','mixed_TST_IGRA','direct_latent_TBI'},'Explicit assay generation required.')
    require(r['threshold'] and isinstance(r['source_cells'],list) and r['source_cells'],'Original cells and threshold must be preserved.')
    require(isinstance(r['migration'],dict) and set(r['migration'])>=
        {'current_residence','birth_geography','previous_residence','exposure_geography','arrival_year','residence_duration'},'Migration fields must be explicit, including unknowns.')
    require(r['review_status']=='internally_extracted_pending_scientific_review','No implicit scientific approval of source transport.')
    n,k=r['denominator'],r['numerator']
    require((n is None)==(k is None),'Numerator and tested denominator jointly known or jointly missing.')
    if n is not None:
        require(type(n) is int and n>0 and type(k) is int and 0<=k<=n,'Invalid numerator/denominator.')
    reported=None if r['reported_prevalence_percent'] is None else decimal(r['reported_prevalence_percent'])/100
    if reported is not None:probability(reported)
    require(reported is not None or n is not None,'No observed value; missing is not zero.')
    with localcontext() as ctx:
        ctx.prec=60
        observed=Decimal(k)/n if n is not None else reported
        # Rounded source percentages may differ, but cannot contradict counts.
        if n is not None and reported is not None:
            require(abs(observed-reported)<=Decimal('.00051'),'Source counts and rounded percentage disagree beyond stated0.1 percentage-point precision.')
        interval=r['uncertainty']
        require(isinstance(interval,dict) and interval.get('kind') in {'reported_95CI','not_reported'},'Explicit source uncertainty status required.')
        if interval['kind']=='reported_95CI':
            a=probability(decimal(interval['lower_percent'])/100);b=probability(decimal(interval['upper_percent'])/100)
            require(a<=observed<=b,'Invalid published uncertainty interval.')
        else:a=b=None
        return dict(record_id=r['record_id'],source_record_sha256=digest(canonical(original)),
            source_id=r['source_id'],geography_id=r['geography_id'],geography_level=r['geography_level'],
            subgroup_id=r['subgroup_id'],participant_pool=r['participant_pool'],membership_group=r['membership_group'],
            survey_year=r['survey_year'],survey_period=r['survey_period'],age_lower=lo,age_upper_exclusive=hi,
            sex=r['sex'],view=r['view'],assay=r['assay'],threshold=r['threshold'],numerator=k,denominator=n,
            observed_positivity=numeric_text(observed),reported_prevalence_percent=r['reported_prevalence_percent'],
            lower=None if a is None else numeric_text(a),upper=None if b is None else numeric_text(b),
            uncertainty_kind=interval['kind'],survey_design=r['survey_design'],sampling_frame=r['sampling_frame'],
            migration=r['migration'],compatibility=r['compatibility'],review_status=r['review_status'],
            transformation='exact_count_ratio' if n is not None else 'published_weighted_percent_divided_by100_no_counts_invented',
            source_preserved=True)


def validation_inputs(records,sources):
    require(isinstance(records,list) and records and isinstance(sources,list) and sources,'Sources and observations required.')
    ids=set();sids=set()
    for s in sources:
        require(set(s)>={'source_id','citation','url','access_date','source_bytes_sha256','licence','extraction'},'Complete source provenance required.')
        require(s['source_id'] not in sids,'Duplicate source ID.');sids.add(s['source_id'])
        require(re.fullmatch('[0-9a-f]{64}',s['source_bytes_sha256']) is not None,'Invalid source checksum.')
        require(s['licence'].get('redistribution')=='permitted_with_attribution' and s['licence'].get('url'),'Explicit redistribution permission required.')
    normalized=[]
    for r in records:
        require(r['record_id'] not in ids and r['source_id'] in sids,'Duplicate observation or missing source.');ids.add(r['record_id'])
        normalized.append(normalize_observation(r))
    return artifact(OBSERVATIONS,sources=sorted(deepcopy(sources),key=lambda s:s['source_id']),
        rows=sorted(normalized,key=lambda r:r['record_id']),source_records=sorted(deepcopy(records),key=lambda r:r['record_id']),
        review_status='internal_transcription_and_consistency_review_only',
        no_clearance_caveat='Assay positivity may persist after treatment or clearance; not proof of viable latent bacilli.')


def freeze_validation_snapshot(records,sources,*,directory=ROOT):
    value=validation_inputs(records,sources)
    files={'source_observations.json':canonical(value['source_records'])+b'\n',
        'sources.json':canonical(value['sources'])+b'\n','model_inputs.json':canonical(value)+b'\n'}
    manifest=seal(dict(schema_id='country_tb.validation_snapshot_manifest.v1',snapshot_id='tbi-validation-'+value['profile_hash'][:12],
        scientific_hash=value['profile_hash'],files={n:digest(b) for n,b in sorted(files.items())},
        immutable=True,licence_policy='source_attribution_retained_no_unlicensed_table_redistribution',
        extraction_version=IMPLEMENTATION))
    files['manifest.json']=canonical(manifest)+b'\n'
    return publish_immutable(directory,manifest['snapshot_id'],files)


def load_validation_snapshot(snapshot_id,*,directory=ROOT):
    require(isinstance(snapshot_id,str) and re.fullmatch('tbi-validation-[0-9a-f]{12}',snapshot_id),'Unsafe snapshot ID.')
    folder=Path(directory)/snapshot_id
    manifest=json.loads(contained_file(folder,'manifest.json').read_text(encoding='utf-8'))
    require(manifest['profile_hash']==profile_hash(manifest) and manifest['snapshot_id']==snapshot_id,'Manifest hash/identity mismatch.')
    require(set(manifest['files'])=={'sources.json','source_observations.json','model_inputs.json'},'Unexpected snapshot file contract.')
    values={}
    for n,h in manifest['files'].items():
        b=contained_file(folder,n).read_bytes();require(digest(b)==h,f'Snapshot checksum mismatch:{n}')
        values[n]=json.loads(b.decode('utf-8'))
    rebuilt=validation_inputs(values['source_observations.json'],values['sources.json'])
    require(canonical(rebuilt)==canonical(values['model_inputs.json']) and rebuilt['profile_hash']==manifest['scientific_hash'],
        'Transformed input does not reproduce from preserved observations.')
    return rebuilt


def load_mortality_validation_snapshot(snapshot_id,*,directory=REPO_ROOT/'data/country_tb_validation/mortality',geography_id):
    require(isinstance(snapshot_id,str) and re.fullmatch('abs-life2022-2024-[0-9a-f]{12}',snapshot_id),'Unknown mortality snapshot contract.')
    folder=Path(directory)/snapshot_id
    manifest=json.loads(contained_file(folder,'manifest.json').read_text(encoding='utf-8'))
    require(manifest['profile_hash']==profile_hash(manifest) and manifest['snapshot_id']==snapshot_id,'Mortality manifest mismatch.')
    require(set(manifest['files'])=={'mortality.json'},'Unexpected mortality payload.')
    b=contained_file(folder,'mortality.json').read_bytes();require(digest(b)==manifest['files']['mortality.json'],'Mortality byte checksum mismatch.')
    value=json.loads(b.decode('utf-8'))
    require(value['profile_hash']==profile_hash(value)==manifest['scientific_hash'],'Mortality scientific hash mismatch.')
    require(value['metadata']['geography_id']==geography_id,'No cross-country mortality fallback.')
    require(value['metadata']['adoption_status']=='validation_snapshot_only_not_a_model_default_or_approved_competing_hazard',
        'Mortality source adoption cannot be inferred from successful extraction.')
    require(value['production_available'] is False and value['externally_validated'] is False,'Mortality validation is internal only.')
    raw={(r['age'],r['sex']):r for r in value['source_rows']};seen=set()
    require(len(raw)==202 and set(raw)=={(a,s) for a in range(101) for s in ('female','male')},'Incomplete source age/sex table.')
    for row in value['rows']:
        key=(row['age_lower'],row['sex']);require(key not in seen,'Duplicate mortality interval.');seen.add(key)
        require(row['age_upper_exclusive']==row['age_lower']+1 and key[0]<100,'No invented terminal-age expansion.')
        require(decimal(row['survival_probability_one_year'])==1-decimal(raw[key]['qx']),'Survival does not reproduce as1-qx.')
    require(seen=={(a,s) for a in range(100) for s in ('female','male')},'Incomplete finite age/sex survival.')
    return value


def _overlap(a,b):
    if a['participant_pool']!=b['participant_pool']:return False
    if a['membership_group']!=b['membership_group'] and 'all' not in {a['membership_group'],b['membership_group']}:return False
    if a['sex']!=b['sex'] and 'both' not in {a['sex'],b['sex']}:return False
    return (a['age_upper_exclusive'] is None or b['age_lower']<a['age_upper_exclusive']) and (b['age_upper_exclusive'] is None or a['age_lower']<b['age_upper_exclusive'])


def validation_split(inputs,fit_ids,validation_ids,*,kind,rationale):
    verify(inputs,OBSERVATIONS);fit_ids=sorted(fit_ids);validation_ids=sorted(validation_ids)
    require(kind in {'held_out_age_records','cross_cohort_transport_comparison','independent_survey'} and rationale,'Explicit split design/rationale required.')
    require(fit_ids and validation_ids and len(set(fit_ids))==len(fit_ids) and len(set(validation_ids))==len(validation_ids),'Nonempty unique split IDs required.')
    require(not set(fit_ids)&set(validation_ids),'Fit/validation data leakage.')
    byid={r['record_id']:r for r in inputs['rows']};require(set(fit_ids+validation_ids)<=byid.keys(),'Unknown split record.')
    if kind=='independent_survey':
        require(all(byid[i]['participant_pool']!=byid[j]['participant_pool'] for i in fit_ids for j in validation_ids),
            'Disjoint age records in one survey are not an independent second survey.')
    for i in fit_ids:
        for j in validation_ids:require(not _overlap(byid[i],byid[j]),'Overlapping participants/parent totals/paired assays cannot be held-out validation.')
    for collection in (fit_ids,validation_ids):
        for k,i in enumerate(collection):
            for j in collection[k+1:]:require(not _overlap(byid[i],byid[j]),'Overlapping fit/validation records would double-count a cohort.')
    return artifact(SPLIT,input_hash=inputs['profile_hash'],fit_ids=fit_ids,validation_ids=validation_ids,kind=kind,rationale=rationale,
        independence_claim=kind=='independent_survey',cluster_caveat='Disjoint age records can remain correlated within sampled households/clusters.')
