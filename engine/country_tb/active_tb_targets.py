"""Active-disease observations for a future, deliberately unavailable operator."""
from copy import deepcopy

from engine.country_tb.calibration_common import reviewed, scientific_context
from engine.country_tb.composition_common import artifact, metadata, require, age_bounds,verify
from engine.country_tb.data_common import decimal, numeric_text, canonical, digest

TARGETS = 'country_tb.active_tb_inference_targets.v1'
TYPES = ('incident_disease','notifications','baseline_prevalent_disease',
         'screen_detected_prevalent_disease','imported_active_disease','mixed_or_unclear')
BLOCKED = 'unavailable_pending_progression_and_risk_forward_operator'
DEPENDENCIES = ('recent_and_remote_progression','age_specific_progression','risk_factor_prevalence',
    'risk_factor_effects','competing_mortality','migration_and_imported_disease',
    'incidence_ascertainment','dynamic_transmission_where_applicable')


def active_tb_targets(records, *, source, provenance, review):
    metadata(source,provenance,review,required_review=False)
    require(isinstance(records,list) and bool(records),'Explicit active-TB observations required.')
    rows,ids = [],set()
    with scientific_context():
        for original in records:
            required={'observation_id','active_tb_target_type','rate','count','unit','denominator',
                'denominator_unit','observation_period','calendar_year','age_lower','age_upper_exclusive',
                'sex','geography_id','geography_level','uncertainty','source','source_class','ascertainment',
                'review','provenance','observation_model'}
            require(isinstance(original,dict) and required <= original.keys(),'Incomplete active-TB target contract.')
            row=deepcopy(original)
            require(isinstance(row['observation_id'],str) and row['observation_id'] and row['observation_id'] not in ids,
                    'Active-TB observation IDs must be unique and nonempty.')
            ids.add(row['observation_id'])
            kind=row['active_tb_target_type']
            require(kind in TYPES,'Unknown active-TB observation type.')
            age_bounds(row)
            require(row['sex'] in {'female','male','both'},'Explicit target sex required.')
            require(row['geography_id'] and row['geography_level'] in {'national','subnational','custom'},'Explicit target geography and level required.')
            metadata(row['source'],row['provenance'],row['review'],required_review=False)
            require(row['source_class'] in {'WHO_estimate','notification','registry','user_supplied','other'},
                    'Distinguish WHO estimates from notifications and other sources.')
            require((row['rate'] is None) != (row['count'] is None),'Supply exactly one active-TB rate or count.')
            require(decimal(row['denominator'])>0 and row['denominator_unit'] in {'persons','person_years'},
                    'A positive denominator and denominator unit are mandatory.')
            period=row['observation_period']
            require(isinstance(period,dict) and {'start_time','end_time'} <= period.keys(),'Observation period required, including counts.')
            start,end=decimal(period['start_time']),decimal(period['end_time'])
            year=row['calendar_year']
            require(type(year) is int and 0<year<=2024 and year<=start<year+1 and start<=end<=2025,
                    'Calendar year and observation period are incompatible.')
            prevalent=kind in {'baseline_prevalent_disease','screen_detected_prevalent_disease'}
            require(end>start or prevalent,'Incident, notification and imported counts/rates need a positive observation period.')
            if row['count'] is not None:
                require(row['unit']=='cases','Active-TB count unit must be cases.')
                row['count']=numeric_text(decimal(row['count']))
            else:
                allowed={'cases_per_100000_persons','cases_per_person'} if prevalent else {
                    'cases_per_100000_person_year','cases_per_person_year'}
                require(row['unit'] in allowed,'Observation units do not match incident versus prevalent disease.')
                row['rate']=numeric_text(decimal(row['rate']))
            u=row['uncertainty']
            require(isinstance(u,dict) and u.get('kind') in {'exact','unspecified','interval'},'Explicit target uncertainty required.')
            if u['kind']=='interval':
                require({'lower','upper','method'} <= u.keys() and bool(u['method']),'Target uncertainty bounds/method required.')
                value=decimal(row['rate'] if row['rate'] is not None else row['count'])
                require(decimal(u['lower'])<=value<=decimal(u['upper']),'Invalid active-TB uncertainty bounds.')
            a=row['ascertainment']
            require(isinstance(a,dict) and a.get('status') in {'reviewed','not_reviewed'} and a.get('description'),
                    'Ascertainment/review status and description required.')
            if row['observation_model'] is not None:
                require(isinstance(row['observation_model'],dict) and {'source','review','provenance'}<=row['observation_model'].keys(),
                        'Observation model needs explicit source, review and provenance.')
                metadata(row['observation_model']['source'],row['observation_model']['provenance'],row['observation_model']['review'])
            eligible=(kind=='incident_disease' and row['review']['status']=='reviewed' and a['status']=='reviewed' and review['status']=='reviewed')
            row.update(denominator=numeric_text(row['denominator']),observation_period={**period,
                'start_time':numeric_text(start),'end_time':numeric_text(end),'duration_years':numeric_text(end-start)},
                eligible_for_future_pressure_calibration=eligible,
                eligibility_reason='reviewed_incident_disease_only_operator_still_unavailable' if eligible else
                    'not_a_reviewed_ascertainment_assessed_incident_disease_target',
                provenance_hash=digest(canonical({'source':row['source'],'provenance':row['provenance'],'review':row['review']})),
                receipt={'source_record':deepcopy(original),'source_record_sha256':digest(canonical(original))})
            rows.append(row)
    return artifact(TARGETS,rows=sorted(rows,key=lambda r:r['observation_id']),source=source,provenance=provenance,
        review=review,status=BLOCKED,required_forward_operator_components=list(DEPENDENCIES),
        infection_pressure_inference=None,tbi_estimate=None,
        interpretation='observation_contract_only_no_disease_to_infection_conversion')


def validate_active_tb_targets(value):
    verify(value,TARGETS)
    originals=[]
    for row in value['rows']:
        receipt=row['receipt']
        require(digest(canonical(receipt['source_record']))==receipt['source_record_sha256'],
                'Active-TB source receipt checksum mismatch.')
        originals.append(receipt['source_record'])
    rebuilt=active_tb_targets(originals,source=value['source'],provenance=value['provenance'],review=value['review'])
    require(rebuilt['profile_hash']==value['profile_hash'],'Active-TB target contract does not reproduce from original observations.')
    return value
