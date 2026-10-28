"""Audited adult candidate values, separated from approval of effect use."""
from copy import deepcopy
from engine.country_tb.risk_common import (
    artifact,verify,require,metadata,ratio,text_number,age_bounds,contains,validate_disjoint,canonical,digest,
)

CATALOG='country_tb.progression_risk_factor_catalog.v1'
RECORD='country_tb.progression_risk_factor_record.v1'
SOURCE={'source_id':'Menzies_2008_Table_2','url':'https://www.tstin3d.com/thinking3D_2008.pdf',
    'citation':'Menzies et al. Int J Tuberc Lung Dis 2008;12:498-505, Table 2, page 500',
    'pubmed':'https://pubmed.ncbi.nlm.nih.gov/18419884/','local_audit':'docs/country_subnational_model/tstin3d_evidence.md'}
SOURCE['audit_commit']='5f56715b4f338338abcdbe04faabbf5a46a7471a'
CANDIDATES=(
 ('cxr_granuloma','CXR granuloma','Granulomatous abnormality on chest radiography','2','24-25','Imaging definition and past treatment need review.'),
 ('heavy_smoking','Smoking more than one pack/day','Cigarette smoking more than one pack/day','3','26','Current smoking is not this intensity definition.'),
 ('underweight','Underweight below 90% ideal weight','Body weight less than 90% of ideal body weight','3','27','BMI below 18.5 is not an automatic substitute.'),
 ('diabetes','Diabetes','Diabetes mellitus; record type, glycaemic control and duration when available','3.6','28-31','Type, control, age and contemporary adjusted effects need review.'),
 ('tnf_inhibitor','TNF inhibitor treatment','TNF inhibitor treatment; record agent, dose and duration; legacy estimate originated with infliximab','4','32','Other agents and modern screening/prophylaxis eras need review.'),
 ('cxr_fibronodular','Fibronodular CXR disease','Fibronodular chest radiographic disease consistent with past TB','6','25,33-34','Lower legacy estimate selected; imaging, active disease exclusion and prior treatment matter.'),
 ('head_neck_carcinoma','Head or neck carcinoma','Carcinoma of the head or neck; record site and treatment status','16','37','All cancer is not an equivalent exposure.'),
 ('haemodialysis','Haemodialysis','Renal failure requiring maintenance haemodialysis','25.3','38-41','All chronic kidney disease is not equivalent.'),
 ('silicosis','Silicosis','Diagnosed pulmonary silicosis','30','42-43','Occupational exposure, severity and HIV confounding need review.'),
 ('hiv','HIV','HIV infection; record stage, CD4 and antiretroviral treatment when known','50','44-47','Legacy all-stage value; treated/untreated HIV and contemporary effects must not be conflated.'),
 ('transplant_immunosuppression','Transplantation with immunosuppression','Transplantation requiring immunosuppressive treatment','74','48-51','All immunosuppression is not equivalent; organ, regimen and era matter.'),
)
FACTOR_IDS=tuple(row[0] for row in CANDIDATES)


def progression_risk_factor_catalog(*,reviewed_effects=None):
    """Overrides add reviewed age/sex-specific effects without altering the inventory."""
    overrides=deepcopy(reviewed_effects or [])
    for row in overrides:
        require(isinstance(row,dict) and {'factor_id','age_lower','age_upper_exclusive','sex','rr','source','evidence_vintage','review','provenance'}<=row.keys(),
                'Reviewed effects need factor, age/sex, RR, evidence and provenance.')
        require(row['factor_id'] in FACTOR_IDS,'Unknown factor; recency/contact are excluded.')
        metadata(row['source'],row['provenance'],row['review']);age_bounds(row)
        require(ratio(row['rr'])>=1,'This catalogue requires a progression risk RR of at least one.')
        row['rr']=text_number(row['rr'])
    for factor in FACTOR_IDS:validate_disjoint([r for r in overrides if r['factor_id']==factor])
    rows=[]
    for factor,label,definition,rr,refs,notes in CANDIDATES:
        rows.append(artifact(RECORD,factor_id=factor,display_label=label,definition_id=factor+':adult_legacy_definition_v1',
            exposure_definition=definition,candidate_effect=rr,effect_type='candidate_relative_risk_not_causal_HR',
            outcome='active_TB_disease_progression_or_reactivation_after_infection_candidate_not_acquisition_pressure',
            source={**SOURCE,'table_reference_chain':refs},evidence_vintage='2008_summary_of_older_studies',
            age_applicability={'age_lower':15,'age_upper_exclusive':None,'rule':'adult_evidence_only_no_automatic_pediatric_use'},
            sex_applicability=['female','male','both'],
            review={'status':'not_reviewed','source_id':'contemporary_effect_review_pending','rationale':'Audited legacy transcription does not approve contemporary progression use.'},
            transcription_review={'status':'reviewed','source_id':'M3A_legacy_inventory_audit','rationale':'Requested values checked against Menzies 2008 Table 2; no contemporary effect appraisal.'},
            evidence_status='provisional_pending_contemporary_scientific_review',definition_compatibility_notes=notes,
            reviewed_effects=sorted([r for r in overrides if r['factor_id']==factor],key=canonical)))
    return artifact(CATALOG,rows=sorted(rows,key=lambda r:r['factor_id']),source=SOURCE,
        excluded_inputs=['recent_infection','close_contact','aids_duplicate_of_hiv'],
        forbidden_uses=['baseline_progression_rate','lifetime_to_age80','TST_diagnostic_interpretation','BCG_or_NTM_adjustment','Styblo_conversion','raw_additive_TSTin3D_scoring'],
        evidence_status='provisional_inventory_not_an_approved_progression_operator')


def catalogue_records(catalog):
    verify(catalog,CATALOG)
    require({r['factor_id'] for r in catalog['rows']}==set(FACTOR_IDS) and len(catalog['rows'])==11,'Catalogue must retain the eleven audited factors.')
    for row in catalog['rows']:verify(row,RECORD)
    return {r['factor_id']:r for r in catalog['rows']}


def applicable_effect(record,cell,*,effect_use,diagnostic_only):
    selected=[r for r in record['reviewed_effects'] if contains(r,cell) and (r['sex']=='both' or r['sex']==cell['sex'])]
    require(len(selected)<=1,'Ambiguous reviewed effect applicability.')
    if selected:
        return {'rr':selected[0]['rr'],'state':'reviewed','source':selected[0]['source'],'review':selected[0]['review'],
                'provenance':selected[0]['provenance'],'catalogue_record_hash':record['profile_hash']}
    if cell['age_lower']<15:
        return {'rr':None,'state':'not_applicable_pending_pediatric_review','catalogue_record_hash':record['profile_hash']}
    if effect_use=='diagnostic_provisional_adult_effects_v1':
        require(diagnostic_only,'Provisional legacy effects may only be used in explicitly diagnostic allocation.')
        return {'rr':record['candidate_effect'],'state':'diagnostic_provisional_pending_contemporary_review',
                'source':record['source'],'review':record['review'],'catalogue_record_hash':record['profile_hash']}
    return {'rr':None,'state':'not_reviewed','catalogue_record_hash':record['profile_hash']}
