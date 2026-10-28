"""Invented offline M3A diagnostics; no country estimates or disease outcomes."""
import argparse
from pathlib import Path
import time

from engine.country_tb.composition_common import population_demography
from engine.country_tb.data_common import canonical,csv_bytes,digest,publish_immutable,REPO_ROOT,DataValidationError
from engine.country_tb.risk_common import artifact,ratio,text_number
from engine.country_tb.risk_factor_catalog import progression_risk_factor_catalog,catalogue_records,applicable_effect
from engine.country_tb.risk_prevalence import risk_factor_prevalence_profile,UNIFORM
from engine.country_tb.risk_policies import (
    overlap_policy,association_policy,multimorbidity_policy,feasible_association,combine_risk_factors,
    MULTIPLICATIVE,ADDITIVE,
)
from engine.country_tb.risk_allocation import allocate_baseline_risk,combine_risk_allocations
from engine.country_tb.baseline_tbi_source import resolve_baseline_tbi
from scripts.preview_country_tb_sources import (
    REVIEW,SOURCE,PROVENANCE,GEO,overall_selection,incidence_selection,make_subgroup,make_mixture,
)


def toy_demography(population=1000,*,children=False):
    ages=(3,10,20) if children else (20,60)
    per=ratio(population)/(2*len(ages))
    rows=[{'age_lower':age,'age_upper_exclusive':age+1,'sex':sex,'population':text_number(per)}
          for sex in ('female','male') for age in ages]
    rows += [{'age_lower':100,'age_upper_exclusive':None,'sex':sex,'population':'0'} for sex in ('female','male')]
    return population_demography(rows,geography_id=GEO,geography_level='national',population=population,
        source=SOURCE,provenance=PROVENANCE,review=REVIEW)


def prevalence_record(catalog,factor,p,*,identifier=None,sex='both',age_lower=15,age_upper_exclusive=None,
    policy=UNIFORM,state='value',definition=None,mapping=None,pattern=None,review=None,reason=None):
    entry=catalogue_records(catalog)[factor]
    return {'record_id':identifier or factor+':'+sex,'factor_id':factor,'state':state,'prevalence':p,
        'age_lower':age_lower,'age_upper_exclusive':age_upper_exclusive,'sex':sex,'allocation_policy':policy,
        'supplied_exposure_definition':entry['exposure_definition'] if definition is None else definition,
        'definition_mapping':mapping,'pattern':pattern,'source':SOURCE,'provenance':PROVENANCE,
        'review':review or REVIEW,'reason':reason}


def allocate_example(resolved,values,*,control='0',catalog=None,records=None):
    catalog=catalog or progression_risk_factor_catalog()
    records=records if records is not None else [prevalence_record(catalog,k,v) for k,v in sorted(values.items())]
    scope=sorted(values)
    profile=risk_factor_prevalence_profile(records,catalog=catalog,factor_ids=scope,
        subgroup_id=resolved['subgroup_id'],geography_id=resolved['geography_id'],source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    overlap=overlap_policy(source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    association=association_policy(global_control=control,source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    policy=multimorbidity_policy(source=SOURCE,provenance=PROVENANCE,review=REVIEW,
        effect_use='diagnostic_provisional_adult_effects_v1')
    return allocate_baseline_risk(resolved,catalog=catalog,prevalence=profile,overlap=overlap,association=association,
        progression_policy=policy,source=SOURCE,provenance=PROVENANCE,review=REVIEW,diagnostic_only=True)


def diagnostics():
    catalog=progression_risk_factor_catalog();groups=[];allocations=[]
    specifications=[('australian_born',8000,'.04',{'diabetes':'.08','heavy_smoking':'.10','haemodialysis':'.01'},'0'),
        ('recent_arrivals',1000,'.30',{'diabetes':'.04','heavy_smoking':'.20','haemodialysis':'.005'},'.5'),
        ('long_term_migrants',1000,'.18',{'diabetes':'.15','heavy_smoking':'.03','haemodialysis':'.02'},'-.5')]
    for name,n,p,values,c in specifications:
        demo=toy_demography(n);group=make_subgroup(name,demo,overall_selection(name,demo,p))
        groups.append(group);allocations.append(allocate_example(group['tbi_profile'],values,control=c,catalog=catalog))
    baseline_mixture=make_mixture(groups)
    australian=combine_risk_allocations(baseline_mixture,allocations,source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    demo=toy_demography();resolved=resolve_baseline_tbi(overall_selection('synthetic_vnm_sexes',demo,'.2'),demo)
    smoking=[prevalence_record(catalog,'heavy_smoking',p,sex=s) for s,p in [('female','.01'),('male','.20')]]
    vnm=allocate_example(resolved,{'heavy_smoking':None},records=smoking,catalog=catalog)
    mismatch=allocate_example(resolved,{'heavy_smoking':None},records=[prevalence_record(catalog,'heavy_smoking','.20',definition='Current smoking')],catalog=catalog)
    multimorbidity=allocate_example(resolved,{'diabetes':'.20','heavy_smoking':'.10','haemodialysis':'.05'},catalog=catalog)
    effects={k:applicable_effect(v,{'age_lower':20,'age_upper_exclusive':21,'sex':'female'},
        effect_use='diagnostic_provisional_adult_effects_v1',diagnostic_only=True) for k,v in catalogue_records(catalog).items()}
    policies=[multimorbidity_policy(source=SOURCE,provenance=PROVENANCE,review=REVIEW,method=m,
        effect_use='diagnostic_provisional_adult_effects_v1') for m in ('highest_applicable_reviewed_rr_v1',MULTIPLICATIVE,ADDITIVE)]
    combinations=[{'policy':p,'combination':combine_risk_factors(['diabetes','heavy_smoking','haemodialysis'],effects,p)} for p in policies]
    young=toy_demography(1200,children=True)
    children=allocate_example(resolve_baseline_tbi(overall_selection('children',young,'.2'),young),
        {'diabetes':'.1','heavy_smoking':'.2','haemodialysis':'.01'},catalog=catalog)
    blocked_born=make_subgroup('australian_born',groups[0]['demography'],incidence_selection('australian_born',groups[0]['demography']))
    blocked_alloc=allocate_example(blocked_born['tbi_profile'],specifications[0][3],catalog=catalog)
    incomplete=combine_risk_allocations(make_mixture([blocked_born,*groups[1:]]),[blocked_alloc,*allocations[1:]],
        source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    return artifact('country_tb.risk_allocation_diagnostic_set.v1',diagnostic_only=True,
        warning='INVENTED VALUES ONLY; NOT AUSTRALIAN OR VIET NAM ESTIMATES',
        examples={'Australian_mixture':australian,'Viet_Nam_synthetic_heavy_smoking_by_sex':vnm,
            'current_smoking_definition_mismatch':mismatch,'multimorbidity':multimorbidity,
            'children_adult_effects_unavailable':children,'incomplete_incidence_selected_mixture':incomplete},
        association_examples=[feasible_association('.3','.2',c) for c in ('-1','0','.5','1')],
        combination_examples=combinations,catalogue=catalog,source=SOURCE,provenance=PROVENANCE,review=REVIEW)


def summary_rows(value):
    rows=[]
    for name,allocation in sorted(value['examples'].items()):
        s=allocation['summary']
        rows.append({'example':name,'population':s['population']['value'],'unresolved_population':s['unresolved_population']['value'],
            'TBI_full':s['tbi_prevalence']['prevalence_full']['value'] if s['tbi_prevalence']['prevalence_full'] is not None else None,
            'any_risk_full':s['any_risk_prevalence']['prevalence_full']['value'] if s['any_risk_prevalence']['prevalence_full'] is not None else None,
            'joint_TBI_risk_full':s['joint_tbi_any_risk']['prevalence_full']['value'] if s['joint_tbi_any_risk']['prevalence_full'] is not None else None,
            'diagnostic_only':True})
    return rows


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=REPO_ROOT/'abm/output/country_tb_m3a_preview')
    args=parser.parse_args(argv);output=args.output.resolve()
    if not output.is_relative_to((REPO_ROOT/'abm/output').resolve()):raise DataValidationError('Developer output must stay under abm/output.')
    start=time.perf_counter();value=diagnostics();rows=summary_rows(value)
    files={'diagnostics.json':canonical(value)+b'\n','summary.csv':csv_bytes(rows,tuple(rows[0])),
        'README.md':b'# Invented risk-allocation diagnostics\n\nNOT COUNTRY ESTIMATES. Provisional adult RRs are diagnostic only. No disease progression or intervention is calculated.\n'}
    files['manifest.json']=canonical(artifact('country_tb.risk_allocation_export_receipt.v1',
        profile_hash_reference=value['profile_hash'],files={k:digest(v) for k,v in sorted(files.items())}))+b'\n'
    print(publish_immutable(output,'risk-diagnostics-'+value['profile_hash'][:12],files))
    for row in rows:print(row)
    print(f'Runtime seconds: {time.perf_counter()-start:.3f}')


if __name__=='__main__':main()
