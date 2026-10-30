"""Offline source extraction. Requires licensed primary XML, never fetches data.

Example: python scripts/freeze_country_tb_validation_3d.py --primary ...
Outputs are NEW immutable validation inputs, not production defaults.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
from xml.etree import ElementTree as ET

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from engine.country_tb.validation_data import freeze_validation_snapshot

DATE='2026-10-30'
LICENSES={
    'PMC10446657':('CC BY4.0','https://creativecommons.org/licenses/by/4.0/'),
    'PMC4633161':('Public domain','https://creativecommons.org/publicdomain/zero/1.0/'),
    'PMC9423895':('US Government public domain; CDC attribution','https://wwwnc.cdc.gov/eid/page/copyright-and-disclaimers'),
}


def text(node):return ' '.join(''.join(node.itertext()).split())


def primary(folder,pmc):
    b=(folder/(pmc+'.xml')).read_bytes();root=ET.fromstring(b)
    doi=next(x.text for x in root.findall('.//article-id') if x.get('pub-id-type')=='doi')
    title=text(root.find('.//article-title'))
    authors=[]
    for a in root.findall('.//article-meta/contrib-group/contrib/name'):
        authors.append(text(a.find('surname'))+' '+text(a.find('given-names')))
    journal=text(root.find('.//journal-title'))
    label,url=LICENSES[pmc]
    permission=' '.join(text(x) for x in root.findall('.//permissions'))
    assert ('CC-BY 4.0' in permission or 'public domain' in permission.lower()),'Unexpected source rights'
    source=dict(source_id=pmc,citation='; '.join(authors)+'. '+title+'. '+journal+'. doi:'+doi,
        url='https://pmc.ncbi.nlm.nih.gov/articles/'+pmc+'/',doi=doi,access_date=DATE,
        source_bytes_sha256=hashlib.sha256(b).hexdigest(),
        licence={'label':label,'url':url,'redistribution':'permitted_with_attribution','source_permission_statement':permission},
        extraction='Exact selected XML table cells; no digitization; source observations separate from transforms')
    return root,source


def table(root,identifier):
    t=next(t for t in root.findall('.//table-wrap') if t.get('id')==identifier)
    return [[text(c) for c in r if c.tag in {'td','th'}] for r in t.findall('.//tr')]


def base(source,*,identifier,geo,subgroup,pool,membership,year,period,lo,hi,sex,view,assay,threshold,frame,design,cells,locator,migration,compatibility):
    return dict(record_id=identifier,source_id=source['source_id'],citation=source['citation'],url=source['url'],access_date=DATE,
        geography_id=geo,geography_level='national' if geo=='who:USA' else 'custom',subgroup_id=subgroup,
        participant_pool=pool,membership_group=membership,survey_year=year,survey_period=period,
        age_lower=lo,age_upper_exclusive=hi,sex=sex,view=view,assay=assay,threshold=threshold,
        sampling_frame=frame,survey_design=design,migration=migration,
        numerator=None,denominator=None,reported_prevalence_percent=None,uncertainty={'kind':'not_reported'},
        review_status='internally_extracted_pending_scientific_review',transformation_notes='Source cells unchanged; normalization is separately stored.',
        source_locator=locator,source_cells=cells,compatibility=compatibility)


def extract(folder):
    records=[];sources=[]
    root,s=primary(folder,'PMC10446657');sources.append(s)
    rows=table(root,'i2220-8372-13-3-83-t03');allowed={'Overall':(0,None,'both','overall'),
        'Female':(0,None,'female','sex'),'Male':(0,None,'male','sex')}
    labels=['0–<10','10–<20','20–<30','30–<40','40–<50','50–<60','60–<70','⩾70']
    allowed.update({label:(i*10,None if i==7 else (i+1)*10,'both','age') for i,label in enumerate(labels)})
    for row in rows:
        if not row or row[0] not in allowed:continue
        lo,hi,sex,view=allowed[row[0]]
        for year,offset in [(2020,1),(2021,5)]:
            for j,assay,threshold in [(1,'QFT-Plus','>=0.35 IU/mL, manufacturer positive algorithm'),
                (2,'TST','>=5mm, 5TU BulBio PPD'),(3,'TST','>=10mm, 5TU BulBio PPD')]:
                value=row[offset+j];m=re.fullmatch(r'([\d,]+)\s*\(([\d.]+)\)',value);assert m,value
                r=base(s,identifier=f'vnm_hhc_{year}_{view}_{sex}_{lo}_{hi}_{j}',geo='who:VNM',subgroup='household_contacts_'+str(year),
                    pool='vnm_campaign_'+str(year),membership='all',year=year,period=str(year),lo=lo,hi=hi,sex=sex,view=view,
                    assay=assay,threshold=threshold,frame='Convenience QFT subset of household contacts in routine community campaigns; different districts by year',
                    design={'kind':'convenience_household_cluster_sample','design_effect':None,'indeterminate':'excluded from paired analysis'},
                    cells=row,locator='Table3/'+row[0]+'/'+str(year)+'/'+threshold,
                    migration={'current_residence':'Viet Nam campaign district','birth_geography':None,'previous_residence':None,
                        'exposure_geography':None,'arrival_year':None,'residence_duration':None},
                    compatibility={'nationally_representative':False,'tbi_anchor':'diagnostic_subgroup_only',
                        'incidence_target':'no_matched_first_episode_followup','sex_age_joint_weights':'not_reported'})
                r.update(numerator=int(m[1].replace(',','')),denominator=int(row[offset].replace(',','')),reported_prevalence_percent=m[2]);records.append(r)
    sources[-1]['reconciliation']={'2020':{'paired_total':1330,'known_age':1327,'missing_age':3},
        '2021':{'paired_total':1158,'known_age':1157,'missing_age':1,'known_sex':1124,'missing_sex':34},
        'warning':'Age,sex,overall and paired assay views overlap. They must never be pooled as independent observations.'}
    root,s=primary(folder,'PMC4633161');sources.append(s)
    rows=table(root,'pone.0140881.t004')
    allowed={'Total':(6,None,'both','overall'),'Female':(6,None,'female','sex'),'Male':(6,None,'male','sex'),
        '6–14':(6,15,'both','age'),'15–24':(15,25,'both','age'),'25–44':(25,45,'both','age'),
        '45–64':(45,65,'both','age'),'≥65':(65,None,'both','age')}
    for row in rows:
        if not row or row[0] not in allowed:continue
        lo,hi,sex,view=allowed[row[0]]
        for col,group in [(1,'all'),(2,'us_born'),(3,'foreign_born')]:
            m=re.search(r'([\d.]+)\s*\(([\d.]+)[–-]([\d.]+)\)',row[col]);assert m,row
            r=base(s,identifier=f'usa_nhanes_2012_{group}_{view}_{sex}_{lo}_{hi}',geo='who:USA',subgroup=group,
                pool='nhanes_2011_2012',membership=group,year=2012,period='2011–2012',lo=lo,hi=hi,sex=sex,view=view,
                assay='QFT-GIT',threshold='Manufacturer QFT-GIT positive algorithm, >=0.35 IU/mL',
                frame='US civilian noninstitutionalized population aged6+; survey-weighted estimates',
                design={'kind':'complex_multistage_weighted_survey','sampling_variance':'published design-based95CI','counts_from_prevalence':'forbidden'},
                cells=row,locator='Table4/'+row[0]+'/'+group,
                migration={'current_residence':'United States','birth_geography':'United States' if group=='us_born' else None,
                    'previous_residence':None,'exposure_geography':'United States proxy assumption only' if group=='us_born' else None,
                    'arrival_year':None,'residence_duration':None},
                compatibility={'nationally_representative':True,'tbi_anchor':'diagnostic_us_born_only' if group=='us_born' else 'migration_history_required',
                    'incidence_target':'WHO_total_episode_not_first_episode','sex_age_joint_weights':'not_reported'})
            r.update(reported_prevalence_percent=m[1],uncertainty={'kind':'reported_95CI','lower_percent':m[2],'upper_percent':m[3],
                'method':'published design-based CI; source reliability flags retained in source_cells'});records.append(r)
    sources[-1]['reconciliation']={'warning':'Weighted prevalence is not k/n. No age-specific counts or population weights invented; source reliability flags retained.'}
    root,s=primary(folder,'PMC9423895');sources.append(s)
    rows=table(root,'T2');allowed={'Total':(0,15,'both','overall'), 'F':(0,15,'female','sex'),'M':(0,15,'male','sex'),
        '<4':(0,4,'both','age'),'5–9':(5,10,'both','age'),'10–14':(10,15,'both','age')}
    for row in rows:
        if not row or row[0] not in allowed:continue
        lo,hi,sex,view=allowed[row[0]];m=re.fullmatch(r'([\d,]+)\s*\(([\d.]+)\)',row[4]);assert m,row
        r=base(s,identifier=f'aus_visa_2016_{view}_{sex}_{lo}_{hi}',geo='who:AUS',subgroup='premigration_child_visa_applications',
            pool='aus_visas_2015_2017',membership='all',year=2016,period='November2015–June2017',lo=lo,hi=hi,sex=sex,view=view,
            assay='mixed_TST_IGRA',threshold='As reported by panel physicians; mixed assay-specific criteria, not a common positivity rule',
            frame='Targeted premigration child visa applications, not distinct resident persons; repeat applicants possible',
            design={'kind':'administrative_selected_screening','BCG':'not_available','unit':'visa_application'},cells=row,locator='Table2/'+row[0],
            migration={'current_residence':'outside_Australia_premigration','birth_geography':None,'previous_residence':None,
                'exposure_geography':None,'arrival_year':None,'residence_duration':None},
            compatibility={'nationally_representative':False,'tbi_anchor':'unavailable_mixed_assays_and_migration',
                'incidence_target':'not_a_resident_incidence_risk_set','published_age_label_ambiguity':row[0]=='<4'})
        r.update(numerator=int(m[1].replace(',','')),denominator=int(row[3].replace(',','')),reported_prevalence_percent=m[2]);records.append(r)
    sources[-1]['reconciliation']={'eligible_applications':134759,'tested_applications':45060,'valid_application_results':44841,
        'positive_applications':1513,'test_count_not_person_count':45345,
        'age_label_warning':'<4 is preserved literally; age4 assignment ambiguous. It is not silently changed to0–4.',
        'assay_warning':'Age-specific combined TST/IGRA rows cannot use a single assay observation operator.'}
    return records,sources


def main():
    p=argparse.ArgumentParser();p.add_argument('--primary',type=Path,required=True);p.add_argument('--output',type=Path,default=None);a=p.parse_args()
    records,sources=extract(a.primary)
    result=freeze_validation_snapshot(records,sources,**({} if a.output is None else {'directory':a.output}))
    print(json.dumps({'snapshot':str(result),'records':len(records),'sources':len(sources)}))


if __name__=='__main__':main()
