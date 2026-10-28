"""Six invented Milestone 2D scenarios; offline developer output, no country estimates."""
from copy import deepcopy
from decimal import Decimal
from pathlib import Path
import argparse
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from engine.country_tb.composition_common import artifact,population_demography,as_catalytic_demography
from engine.country_tb.calibration_common import scientific_context
from engine.country_tb.data_common import canonical,csv_bytes,digest,numeric_text,publish_immutable,DataValidationError
from engine.country_tb.pressure_shape import supplied_shape,apply_backcast
from engine.country_tb.catalytic_calibration import calibration_policy,calibrate_scale
from engine.country_tb.tbi_baseline import baseline_tbi_profile
from engine.country_tb.tbi_observations import assay_operator
from engine.country_tb.direct_tbi_profile import direct_tbi_profile,UNAVAILABLE_RECENCY,CONDITIONAL_RECENCY
from engine.country_tb.active_tb_targets import active_tb_targets
from engine.country_tb.baseline_tbi_source import (
    baseline_tbi_source,resolve_baseline_tbi,interface_specification,INCIDENCE,OBSERVED,OVERALL,DIRECT_TOTAL,
)
from engine.country_tb.population_subgroups import population_subgroup,population_mixture

REVIEW={'status':'reviewed','source_id':'synthetic:M2D',
        'rationale':'Invented software-validation scenario; not an empirical review or country estimate.'}
SOURCE={'source_id':'synthetic:M2D','citation':'Invented values only; not Australian estimates'}
PROVENANCE={'source_id':'synthetic:M2D','diagnostic_only':True,'real_evidence_values_adopted':False}
GEO='who:AUS'


def synthetic_demography(population=1000,tail_population=20,*,geography_id=GEO,geography_level='national'):
    with scientific_context():
        total,tail=Decimal(str(population)),Decimal(str(tail_population))
        rows=[{'age_lower':a,'age_upper_exclusive':a+1 if a<100 else None,'sex':s,
            'population':numeric_text((total-tail)/200 if a<100 else tail/2)}
            for s in ('female','male') for a in range(101)]
    return population_demography(rows,geography_id=geography_id,geography_level=geography_level,
        population=population,source=SOURCE,provenance=PROVENANCE,review=REVIEW)


def synthetic_history(geography_id=GEO):
    q=supplied_shape([{'start':'2000','end':'2025','shape':'1'}],geography_id=geography_id,
        review=REVIEW,provenance=PROVENANCE,normalization_time='2024')
    return apply_backcast(q,start_time='1900')


def direct_operator():
    return assay_operator(assay_class='direct_latent_tbi',name='Invented direct latent-TBI observation',
        positivity_threshold='not_applicable_direct_latent_TBI',sensitivity='1',specificity='1',review=REVIEW,provenance=PROVENANCE)


def overall_selection(subgroup_id,demography,prevalence,*,profile_kind='diagnostic',history=None):
    identifier='overall:'+subgroup_id
    policy=calibration_policy(alpha_upper='2',observation_weights={identifier:'1'},residual_tolerance='.000001',review=REVIEW)
    anchor={'observation_id':identifier,'population_or_subgroup_id':subgroup_id,'population_denominator':demography['population'],
        'prevalence':prevalence,'numerator':None,'denominator':None,'year':2024,'source':SOURCE,'provenance':PROVENANCE,
        'review':REVIEW,'uncertainty':{'kind':'exact'},'assay_operator':direct_operator(),
        'coverage_policy':'apply_anchor_to_model_covered_population_v1'}
    return baseline_tbi_source(mode=OBSERVED,geography_id=demography['geography_id'],geography_level=demography['geography_level'],
        subgroup_id=subgroup_id,source=SOURCE,provenance=PROVENANCE,review=REVIEW,profile_kind=profile_kind,
        configuration={'submode':OVERALL,'baseline_time':'2024.5','recent_window_years':'5',
            'pressure_history':history or synthetic_history(demography['geography_id']),'calibration_policy':policy,'anchor':anchor})


def direct_selection(subgroup_id,demography,records,*,recency_policy=UNAVAILABLE_RECENCY,comparison=None,profile_kind='diagnostic'):
    profile=direct_tbi_profile(records,source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    return baseline_tbi_source(mode=OBSERVED,geography_id=demography['geography_id'],geography_level=demography['geography_level'],
        subgroup_id=subgroup_id,source=SOURCE,provenance=PROVENANCE,review=REVIEW,profile_kind=profile_kind,
        configuration={'submode':DIRECT_TOTAL,'baseline_time':'2024.5','recent_window_years':'5',
            'direct_profile':profile,'recency_policy':recency_policy,'comparison_profile':comparison})


def incidence_selection(subgroup_id,demography,*,target_type='incident_disease'):
    target={'observation_id':'synthetic:active-disease','active_tb_target_type':target_type,'rate':'10','count':None,
        'unit':'cases_per_100000_person_year','denominator':demography['population'],'denominator_unit':'persons',
        'observation_period':{'start_time':'2024','end_time':'2025'},'calendar_year':2024,'age_lower':0,
        'age_upper_exclusive':None,'sex':'both','geography_id':demography['geography_id'],'geography_level':demography['geography_level'],
        'uncertainty':{'kind':'exact'},'source':SOURCE,'source_class':'user_supplied',
        'ascertainment':{'status':'reviewed','description':'Invented observation for blocked-path testing'},
        'review':REVIEW,'provenance':PROVENANCE,'observation_model':None}
    targets=active_tb_targets([target],source=SOURCE,provenance=PROVENANCE,review=REVIEW)
    return baseline_tbi_source(mode=INCIDENCE,geography_id=demography['geography_id'],geography_level=demography['geography_level'],
        subgroup_id=subgroup_id,source=SOURCE,provenance=PROVENANCE,review=REVIEW,profile_kind='diagnostic',
        configuration={'targets':targets,'pressure_history':None,'forward_operator':None})


def comparison_model(demography,*,reviewed_tail=True):
    h=synthetic_history(demography['geography_id'])
    p=calibration_policy(alpha_upper='2',observation_weights={},residual_tolerance='0',review=REVIEW)
    fit=calibrate_scale(h,p,supplied_alpha='.01',scale_provenance=PROVENANCE)
    tails={s+':100':{'review':REVIEW,'provenance':PROVENANCE,'ages':[{'age':'101','weight':'1'}]}
           for s in ('female','male')} if reviewed_tail else None
    return baseline_tbi_profile(as_catalytic_demography(demography),fit,baseline_time='2024.5',tail_distributions=tails)


def membership(subgroup_id,person_ids=None):
    return {'partition_id':'synthetic:Australian-residents','category_id':subgroup_id,
        'definition':'Invented mutually exclusive resident partition: '+subgroup_id,'source':SOURCE,
        'provenance':PROVENANCE,'review':REVIEW,'person_ids':person_ids}


def make_subgroup(subgroup_id,demography,selection,**metadata):
    return population_subgroup(subgroup_id=subgroup_id,label=subgroup_id.replace('_',' '),
        description='Invented diagnostic subgroup; not an Australian population estimate',
        current_residence_geography=demography['geography_id'],geography_level=demography['geography_level'],
        population=demography['population'],demography=demography,tbi_source=selection,membership=membership(subgroup_id),
        source=SOURCE,provenance=PROVENANCE,review=REVIEW,**metadata)


def make_mixture(groups,*,mixture_id='synthetic:Australian-mixture',total_population=None,**kwargs):
    policy={'mode':'reviewed_disjoint_partition_v1','partition_id':'synthetic:Australian-residents',
        'subgroup_ids':[g['subgroup_id'] for g in groups],'source':SOURCE,'provenance':PROVENANCE,'review':REVIEW}
    return population_mixture(groups,mixture_id=mixture_id,geography_id=GEO,geography_level='national',
        total_population=total_population if total_population is not None else numeric_text(sum(Decimal(g['population']) for g in groups)),
        membership_policy=policy,source=SOURCE,provenance=PROVENANCE,review=REVIEW,**kwargs)


def synthetic_scenarios():
    demo=synthetic_demography()
    one=resolve_baseline_tbi(overall_selection('single_population',demo,'.20'),demo)
    groups=[]
    for name,pop,prevalence,origin,arrival,duration in [
        ('australian_born',8000,'.04','who:AUS',None,None),
        ('recent_arrivals',1000,'.30','who:VNM',2023,'1.5'),
        ('long_term_migrants',1000,'.18','who:IND',2010,'14.5')]:
        d=synthetic_demography(pop,0)
        groups.append(make_subgroup(name,d,overall_selection(name,d,prevalence),
            country_or_region_of_birth={'geography_id':origin,'geography_level':'national'},
            exposure_history_geography=None,arrival_year=arrival,duration_in_current_setting_years=duration))
    mixture=make_mixture(groups)
    records=[{'record_id':name,'age_lower':lo,'age_upper_exclusive':hi,'sex':'both',
        'age_allocation_policy':'piecewise_constant_within_band_v1','any_tbi':value}
        for name,lo,hi,value in [('child',0,15,'0.1000'),('adult',15,65,'0.2000'),('older',65,None,'0.3000')]]
    model=comparison_model(demo)
    direct=resolve_baseline_tbi(direct_selection('override',demo,records,comparison=model),demo)
    hybrid=resolve_baseline_tbi(direct_selection('override_hybrid',demo,records,
        recency_policy=CONDITIONAL_RECENCY,comparison=model),demo)
    blocked=resolve_baseline_tbi(incidence_selection('blocked',demo),demo)
    blocked_born=make_subgroup('australian_born',groups[0]['demography'],incidence_selection('australian_born',groups[0]['demography']))
    incomplete=make_mixture([blocked_born,*groups[1:]],mixture_id='synthetic:incomplete-mixture')
    return artifact('country_tb.baseline_tbi_synthetic_scenario_set.v1',diagnostic_only=True,
        warning='INVENTED VALUES ONLY; NOT AUSTRALIAN OR OTHER COUNTRY ESTIMATES',
        scenarios={'1_overall_anchor':one,'2_Australian_mixture':mixture,'3_complete_direct_override':direct,
            '4_total_only_default':direct,'4_explicit_conditional_recency':hybrid,
            '5_incidence_mode_blocked':blocked,'6_incomplete_mixture':incomplete},
        interface=interface_specification(),source=SOURCE,provenance=PROVENANCE,review=REVIEW)


def summary_rows(scenarios):
    rows=[]
    for name,value in sorted(scenarios['scenarios'].items()):
        summary=value['summary'];tbi=summary['components']['any_tbi']
        rows.append({'scenario':name,'population':summary['population'],'estimated_population':summary['estimated_population'],
            'unresolved_population':summary['unresolved_population'],'TBI_prevalence_covered':tbi['prevalence_covered'],
            'TBI_prevalence_full':tbi['prevalence_full'],'recent_prevalence_full':summary['components']['recent']['prevalence_full'],
            'status':value['status'],'diagnostic_only':True,'usable_for_downstream_simulation':value['usable_for_downstream_simulation']})
    return rows


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'abm/output/country_tb_m2d_preview')
    args=parser.parse_args(argv);output=args.output.resolve()
    if not output.is_relative_to((ROOT/'abm/output').resolve()):raise DataValidationError('Developer output must stay under abm/output.')
    start=time.perf_counter();scenarios=synthetic_scenarios();rows=summary_rows(scenarios)
    files={'synthetic_scenarios.json':canonical(scenarios)+b'\n','interface_specification.json':canonical(scenarios['interface'])+b'\n',
        'summary.csv':csv_bytes(rows,tuple(rows[0])),'README.md':b'# Invented Milestone 2D diagnostics\n\nNOT AUSTRALIAN OR OTHER COUNTRY ESTIMATES. No empirical input values adopted.\nIncidence inference remains deliberately blocked. No production UI or engine is connected.\n'}
    files['manifest.json']=canonical(artifact('country_tb.baseline_tbi_export_receipt.v1',
        scenario_hash=scenarios['profile_hash'],files={k:digest(v) for k,v in sorted(files.items())}))+b'\n'
    path=publish_immutable(output,'synthetic-sources-'+scenarios['profile_hash'][:12],files)
    print(path)
    for row in rows:print(row)
    print(f'Runtime seconds: {time.perf_counter()-start:.3f}')
    return 0


if __name__=='__main__':raise SystemExit(main())
