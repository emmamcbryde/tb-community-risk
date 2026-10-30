"""Proposed evidence/parameter package. No central configuration is approved."""
from copy import deepcopy

from engine.country_tb.validation_data import artifact,REGISTRY,PARAMETERS,require
from engine.country_tb.data_common import seal
from engine.country_tb.progression_natural_history import curve_specification,natural_history
from engine.country_tb.tbi_observations import assay_operator
from engine.country_tb.risk_factor_catalog import CANDIDATES
from engine.country_tb.empirical_validation import REVIEW

ADULT={'source_id':'Menzies_NA_2021','citation':'Menzies NA et al. Epidemiology2021;32:70–78. doi:10.1097/EDE.0000000000001271',
    'url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC7707158/'}
USPSTF={'source_id':'USPSTF_Jonas_2023_assay_review','citation':'Jonas DE et al. Screening for Latent Tuberculosis Infection in Adults: Updated Evidence Report and Systematic Review for the USPSTF. JAMA2023. doi:10.1001/jama.2023.3954',
    'url':'https://www.uspreventiveservicestaskforce.org/uspstf/document/RecommendationStatementFinal/latent-tuberculosis-infection-screening'}


def progression_sets():
    result=[]
    for name,risks,tail in [('conservative',['.032','.058','.064','.070'],'.00032'),
        ('central_candidate',['.038','.066','.072','.079'],'.00038'),
        ('higher_progression',['.045','.077','.083','.089'],'.00045')]:
        curve=curve_specification([{'time_years':t,'cumulative_risk':p} for t,p in [('0','0'),*zip(['1','5','10','25'],risks)]],
            tail='fixed_remote_hazard_v1',tail_parameters={'hazard_per_year':tail},source=ADULT,
            provenance={'scope':'diagnostic_source_marginal_bounds_not_joint_confidence_band','tail':'constant_extension_is_unapproved_modelling_assumption'},
            review=REVIEW,status='diagnostic_adult_surrogate',risk_interpretation='diagnostic_cumulative_risk_surrogate_pending_cause_specific_review_v1')
        history=natural_history({'0-4':None,'5-14':None,'15+':curve},source=ADULT,
            provenance={'internal_synthesis_only':True,'no_automatic_pediatric_or_sex_effect':True},review=REVIEW)
        result.append(artifact(PARAMETERS,parameter_set_id=name,natural_history=history,
            cumulative_risk_knots=curve['knots'],segment_hazards=curve['segments'],interpolation=curve['interpolation'],
            central_adult_curve=None,central_remote_tail=None,pediatric_curves={'0-4':None,'5-14':None},
            diagnostic_remote_hazard=tail,review_status='proposed_sensitivity_only_external_review_required',
            applicable_population='focal_modelled_adult_infected_at20; no_country_transport_approval',
            evidence_strength='published_model_outputs_with_historical_cohort_inputs; net_hazard_transport_unresolved',
            limitations=['2year_value_is_interpolation_not_an_observed_knot','within_year1_constant_hazard_is_assumed',
                'source_marginal_UI_endpoints_are_not_a_joint_curve_interval','Sloot_incident_contact5year2.4percent_has_a_different_estimand',
                'older_age_and_late_tail_not_identified','competing_risk_CIF_not_automatically_net_risk',
                'ordered_cumulative_risk_knots_do_not_order_all_segment_hazards_or_remote_incidence'],
            mortality_default=None,risk_factor_default=None,migration_default=None,recurrence_adjustment=None,
            source=ADULT,scientific_approval='pending',selected_by_default=False))
    return result


def review_parameter_set(parameter_set_id):
    require(parameter_set_id in {'conservative','central_candidate','higher_progression'},
        'Explicit review-set selection required; no default or unknown-set fallback.')
    return next(r for r in progression_sets() if r['parameter_set_id']==parameter_set_id)


def assay_scenario(product,*,variant='central_candidate'):
    require(product in {'QFT-Plus','QFT-GIT'},'No interchangeability or unreviewed mixed/TST operator.')
    se,sp=('.89','.98') if product=='QFT-Plus' else ('.81','.99')
    changes={'central_candidate':{},'lower_sensitivity':{'se':'.875' if product=='QFT-Plus' else '.79'},
        'higher_sensitivity':{'se':'.942' if product=='QFT-Plus' else '.84'},
        'lower_specificity':{'sp':'.955' if product=='QFT-Plus' else '.98'},
        'higher_specificity':{'sp':'.989' if product=='QFT-Plus' else '.995'}}
    require(variant in changes,'Unknown fixed assay scenario; no fitting of performance.')
    se=changes[variant].get('se',se);sp=changes[variant].get('sp',sp)
    provenance={'source':USPSTF,'one_way_bounds_source':{'source_id':'Oh_2021_QFT_review','url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC8423471/'},
        'bound_attribution':'QFT-Plus bounds:Oh2021; QFT-GIT sensitivity/lower-specificity:USPSTF2023; upper-specificity:Oh2021. Bounds from different reviews are separate one-way scenarios.',
        'sensitivity':'active_TB_reference_populations_not_a_gold_standard_for_latent_viability',
        'specificity':'low_risk_reference_populations; country_BCG_immunosuppression_transport_not_approved',
        'scenario':'fixed_source_bound_sensitivity_not_joint_assay_CI; no_simultaneous_alpha_performance_fit'}
    op=assay_operator(assay_class='IGRA',name=product,positivity_threshold='manufacturer_positive_algorithm_0.35IU_per_mL',
        sensitivity=se,specificity=sp,review=REVIEW,provenance=provenance)
    return seal(dict(schema_id='country_tb.assay_validation_scenario.v1',operator=op,age_lower=15,variant=variant,
        source=USPSTF,provenance=provenance,review_status='internal_candidate_pending_population_specific_scientific_review',
        pediatric_performance=None,production_available=False,externally_validated=False))


def evidence_registry():
    # Numeric effect measures remain separate; none updates the frozen3A catalogue.
    definitions={r[0]:r for r in CANDIDATES}
    data=[
        ('cxr_granuloma','incompatible_definition','CDC2020_LTBI_guide','https://stacks.cdc.gov/view/cdc/118719/cdc_118719_DS1.pdf',
            'Fully calcified discrete granulomas have no increased progression risk per guidance; broad granuloma RR2 is not interchangeable with fibrosis.',None,None,None),
        ('heavy_smoking','incompatible_definition','JGlobHealth2026_doi10.7189/jogh.16.04079','https://pmc.ncbi.nlm.nih.gov/articles/PMC12964331/',
            '29 Chinese studies,any smoking versus nonsmoking; not >one pack/day or conditional on TBI. OR and HR pooled separately.',
            'adjusted_or_study_specific_HR',2.39,[1.28,4.45]),
        ('underweight','incompatible_definition','Sinha_et_al_CID2024_79_233','https://pmc.ncbi.nlm.nih.gov/articles/PMC11259213/',
            'Prospective Indian household contacts: BMI-defined undernutrition differs from<90%ideal weight; disease risk distinguished from acquisition. No definition mapping approved.',None,None,None),
        ('diabetes','sensitivity_only','Franco_et_al_Cochrane2024_CD016013','https://pmc.ncbi.nlm.nih.gov/articles/PMC11342417/',
            '48 cohorts>61million,median5year follow-up. AdjustedHR1.90,RR1.60 andOR1.61 are distinct. HR<10year1.52(1.47–1.57); not conditional-TBI HR.',
            'adjusted_HR',1.90,[1.51,2.40]),
        ('tnf_inhibitor','sensitivity_only','Korean_AS_biologics_cohort2025','https://pmc.ncbi.nlm.nih.gov/articles/PMC12107047/',
            '2010–2021 claims;43643 adults,agent-specific screening/treatment era. Non-exposure13.96 versus infliximab33.62 and adalimumab27.79/10000PY. No universal HR4.',None,None,None),
        ('cxr_fibronodular','insufficient_evidence','Zhu_et_al_NatureCommun2025_doi10.1038/s41467-025-64834-9','https://pmc.ncbi.nlm.nih.gov/articles/PMC12603200/',
            '183808 elderly Chinese screened2020,followed2022. Stable,uncertain,active-suggestive findings differ; not a common untreated LTBI fibronodular class. Active-suggestive HR16.7 must not become scar RR.',None,None,None),
        ('head_neck_carcinoma','sensitivity_only','Lai_et_al_2017_head_neck_cohort','https://pmc.ncbi.nlm.nih.gov/articles/PMC5671858/',
            'Taiwan retrospective matched cancer cohort; pulmonary TB4.70 versus1.64/1000PY; conditional TBI status/treatment-era transport not established.',
            'adjusted_HR',2.90,[2.11,3.99]),
        ('haemodialysis','diagnostic_legacy','Dobler_McDonald_Marks2011_and_CKD_reviews2023','https://pmc.ncbi.nlm.nih.gov/articles/PMC3246480/',
            'Australian dialysis2001–2006,14506patients,mean3.8years,37TB; general-population adjusted RR7.8.2023CKDwithoutfailure review is not dialysis; prevalence meta-analysis is not progression.',
            'adjusted_disease_RR',7.8,[3.3,18.7]),
        ('silicosis','sensitivity_only','Ehrlich_et_al_BMC_PublicHealth2021','https://pubmed.ncbi.nlm.nih.gov/34016067/',
            'Eight silicosis studies; severity,ongoing exposure,HIV and populations differ; no universal conditional-TBI HR30.',
            'pooled_disease_RR',4.01,[2.88,5.58]),
        ('hiv','insufficient_evidence','Gatechompol_et_al_BMC_InfectDis2022_doi10.1186/s12879-022-07332-3','https://pmc.ncbi.nlm.nih.gov/articles/PMC8988401/',
            'Thai adults onART,2636people,24229PY,median7.6years. Death-competing risk; TB4.66/1000PY;1/2/5/7year CIF0.7/1.7/3.3/4.3% afterART,not infection. CD4/ART/TPT-specific input required.',None,None,None),
        ('transplant_immunosuppression','sensitivity_only','TBnet_JInfect2026_92_106668_doi10.1016/j.jinf.2025.106668','https://pubmed.ncbi.nlm.nih.gov/41429410/',
            '5805adult European recipients,medianage51;33.8%infection-tested,10.3%TPT;68.1/100000PY,SIR6.2. Donor/acquired/prevalent disease and organ/regimen differ; not TBI-conditionalHR74.',
            'standardized_incidence_ratio',6.2,[3.9,9.3]),
    ]
    rows=[]
    for fid,status,citation,url,notes,measure,value,interval in data:
        original=definitions[fid]
        rows.append(dict(factor_id=fid,exact_legacy_definition=original[2],legacy_RR=original[3],status=status,
            citation=citation,url=url,age_range='adult evidence; no automatic pediatric application',
            treatment_status='study-specific; ART,TPT,immunosuppressive_regimen_or_cancer_therapy_not_uniform',
            TBI_ascertainment='not_a_common_confirmed_latent_infection_denominator',reference_population='see extraction_notes; not assumed TBI-positive',
            effect_measure=measure,adjusted_estimate=value,uncertainty=interval,follow_up=None,absolute_event_rate=None,
            adjustment_variables=None,unavailable_fields_reason='Full effect-specific primary extraction not established for a transportable progression estimand; numerical use withheld.',
            extraction_notes=notes,applicable_as_progression_hazard_multiplier=False,
            proposed_use='definition_review_or_explicit_sensitivity_only_no_catalogue_override',
            contemporary_relevance='targeted_search2026-10-30; older cohort retained where no matching modern conditional estimate identified',
            external_scientific_review='pending'))
    details={
        'diabetes':{'follow_up':'median5years,IQR1.5–10,range1–16.9','age_range':'all studies adults;nine also children/adolescents,not a pediatric conditional effect',
            'reference_population':'people without diabetes in general-population/cohort settings','adjustment_variables':'study-specific multivariable models; no common covariate set',
            'additional_effect_measures':{'RR':1.60,'RR_95CI':[1.42,1.80],'OR':1.61,'OR_95CI':[1.27,2.04],'HR_prediction_interval':[.83,4.39]}},
        'tnf_inhibitor':{'follow_up':'2010–2021claims,variable patient/agent exposure','absolute_event_rate':{'unit':'cases_per10000PY','non_exposure':13.96,'adalimumab':27.79,'etanercept':14.28,'golimumab':21.19,'infliximab':33.62},
            'treatment_status':'agent-specific biologic exposure; background screening/prophylaxis not untreated infection'},
        'head_neck_carcinoma':{'absolute_event_rate':{'unit':'cases_per1000PY','cancer':4.70,'no_cancer':1.64},'reference_population':'matched people without cancer'},
        'haemodialysis':{'follow_up':'mean3.8years,Australia2001–2006','absolute_event_rate':{'unit':'cases_per100000PY','dialysis':66.8,'general_population':5.7},
            'reference_population':'Australian general population','adjustment_variables':['age','sex','Indigenous_status','TB_incidence_country_of_birth']},
        'hiv':{'age_range':'>=18years;median32','follow_up':'median7.6years,IQR1.9–15.7','treatment_status':'onART; incidentTB>3months afterART start',
            'absolute_event_rate':{'unit':'cases_per1000PY','rate':4.66,'95CI':[3.87,5.60]},'reference_population':'no common HIV-negative TBI-positive comparator'},
        'transplant_immunosuppression':{'age_range':'adult,median51years','treatment_status':'organ-specific immunosuppression;10.3%TPT;33.8%infection-tested',
            'absolute_event_rate':{'unit':'cases_per100000PY','rate':68.1},'reference_population':'age/background standardized general population,not TBI-positive'}
    }
    for row in rows:
        row['estimate']=row['adjusted_estimate']
        if row['factor_id'] in {'heavy_smoking','silicosis','transplant_immunosuppression'}:row['adjusted_estimate']=None
        row.update(details.get(row['factor_id'],{}))
    contextual=[
        ('adult_curve',ADULT['url'],'included_sensitivity_only','Focal modelled infection at20; historical15/10/15year source cohorts; four risks are not a validated national cause-specific curve.'),
        ('adult_incident_contacts','https://pubmed.ncbi.nlm.nih.gov/25265362/','included_estimand_comparison','Sloot2014: five-year incident TB2.4%(1.2–4.7) in untreated infected contacts; contact ascertainment not known infection time. Do not merge with prevalent+incident risks.'),
        ('adult_contact_timing','https://pmc.ncbi.nlm.nih.gov/articles/PMC6534268/','included_context','Prospective contact follow-up with registry matching; early disease and baseline screening must be separated.'),
        ('pediatric_IPD','https://pmc.ncbi.nlm.nih.gov/articles/PMC7289654/','included_no_curve_adopted','Martinez2020:19%under5 risk combines prevalent/incident outcomes;83%cases<90days; contact,TPT,BCG and age definitions prevent untreated last-infection curve transport.'),
        ('pediatric_birth_cohort','https://pmc.ncbi.nlm.nih.gov/articles/PMC11579303/','included_no_curve_adopted','daCosta2024: repeatedTST,BCG,selectiveTPT,ongoing community exposure to10years; no complete untreated5–14curve.'),
        ('remote_tail','https://pubmed.ncbi.nlm.nih.gov/33891908/','included_external_required','Dale2021: late risk persists,declines; unknown-exposure cohorts<100/100000PY beyond10years; age/reinfection limit inference. No universal beyond25year tail.'),
        ('vietnam_HHC','https://pmc.ncbi.nlm.nih.gov/articles/PMC10446657/','included_licensed_validation','Exact age/sex/assay counts,CCBY4; convenience contacts and different campaign districts by year,not national. TST excludes<5/pregnancy.'),
        ('vietnam_CaMau','https://pubmed.ncbi.nlm.nih.gov/29471900/','candidate_data_required','Marks2018,cluster-random community sample,1319adults; acquire exact age table and rights before snapshot adoption.'),
        ('vietnam_two_provinces','https://pmc.ncbi.nlm.nih.gov/articles/PMC9923314/','context_only_not_redistributed','Vo2023 mixed HHC,vulnerable and healthcare groups; sourceCCBY-NC. No unreviewed commercial redistribution or population transport.'),
        ('vietnam_active_prevalence','https://pmc.ncbi.nlm.nih.gov/articles/PMC7179905/','excluded_as_TBI_anchor','Active disease prevalence is not TBI or incident disease.'),
        ('australia_arrival_children','https://pmc.ncbi.nlm.nih.gov/articles/PMC9423895/','included_descriptive_validation','Publicdomain test counts; visa applications not persons,selected overseas screening,mixed assays and age<4 label ambiguity.'),
        ('australia_refugee_clinic','https://pubmed.ncbi.nlm.nih.gov/30565764/','candidate_data_required','28/235 tested children,2014–2015 clinic; selective testing and age/assay/rights extraction incomplete.'),
        ('australia_born_migrant_model','https://doi.org/10.1183/13993003.01218-2018','excluded_as_independent_empirical_validation','Dale2018 is incidence/history-derived TBI modelling; would be circular validation. Useful migration structure only.'),
        ('US_NHANES','https://pmc.ncbi.nlm.nih.gov/articles/PMC4633161/','included_licensed_validation','Probability survey,weighted QFT-GIT prevalence/CI; publicdomain. Nativity retained; no counts reconstructed.'),
        ('US_NHANES_other_analyses','https://pmc.ncbi.nlm.nih.gov/articles/PMC5003328/','excluded_duplicate_cohort','SameNHANES2011–2012; agreement evidence,not independent second validation dataset.'),
        ('Korea_intermediate_setting','https://doi.org/10.3904/kjim.2015.095','candidate_data_required','Prospective adult contacts; separate from2016national survey. Primary age/count/assay extraction and representativeness review needed before adoption.'),
        ('US_prior_treatment','https://wwwnc.cdc.gov/eid/article/25/10/19-0439_article','included_discrepancy','Prior treatment changes eligibility; positive assays do not identify untreated viable TBI.'),
        ('assay_review',USPSTF['url'],'included_candidate_only','No latent gold standard; product,age,BCG and immunosuppression transport require review.'),
        ('mortality_ABS','https://www.abs.gov.au/statistics/people/population/life-expectancy/2022-2024','included_validation_snapshot_only','CCBY4,period age/sex survival; all-cause death scope and100+ unresolved for country adoption.'),
        ('mortality_UN','https://population.un.org/wpp/','candidate_data_required','WPP2024 country/sex complete/abridged tables,CCBY3IGO; licensing/coverage/vintage frozen extraction required forVNM.'),
        ('mortality_WHO','https://www.who.int/publications/m/item/who-methods-and-data-sources-for-life-tables-2000-2021','candidate_data_required','Country/sex life tables; period/model inputs and source licensing require snapshot-specific review.'),
        ('WHO_observation','https://data.who.int/dashboards/tuberculosis/routine-tb-surveillance-indicators','included_compatibility_only','New/recurrent episodes; failure re-registration differs. No unsupported recurrence/import subtraction.'),
    ]
    return artifact(REGISTRY,search_date='2026-10-30',protocol='docs/country_subnational_model/milestone3d_protocol.md',
        search_queries=['Vietnam latent tuberculosis infection population survey age IGRA prevalence QuantiFERON community age table',
            'NHANES2011 2012 latent tuberculosis infection age group prevalence table IGRA United States',
            'Australia refugee latent tuberculosis IGRA age prevalence cohort table open access',
            'Prevalence of latent tuberculous infection adults general population Ca Mau Viet Nam',
            'adult tuberculosis time since infection Sloot2014 contact cohort first five years2024 risk systematic review',
            'late reactivation tuberculosis systematic review Dale2021 25years long term age hazard',
            'diabetes tuberculosis risk systematic review meta analysis2024 adjusted hazard ratio cohort',
            'solid organ transplantation tuberculosis cohort2024 2023 incidence latent infection risk hazard',
            'silicosis tuberculosis meta analysis2021 relative risk4.01',
            'TNF inhibitors tuberculosis2024 cohort incidence Korea',
            'HIV tuberculosis ART incidence cohort2023 CD4',
            'calcified granuloma tuberculosis risk guideline',
            'head and neck cancer tuberculosis hazard2023 cohort',
            'dialysis tuberculosis meta-analysis2024 latent',
            'undernutrition tuberculosis incidence2024 cohort RATIONS household contacts body mass index',
            'USPSTF2023 evidence review latent tuberculosis QFTPlus sensitivity89 specificity98 QFTGIT81 99',
            'ABS2022-2024 life expectancy life tables xlsx copyright creative commons',
            'WHO tuberculosis data definitions new relapse recurrent treatment failure incidence estimates2025'],
        scope='focused_internal_synthesis_not_exhaustive_systematic_review_or_external_signoff',
        factor_rows=rows,rows=[dict(evidence_id=i,url=u,decision=d,reason=r) for i,u,d,r in sorted(contextual)],
        no_independent_cohort_double_counting=True,review_status='internal_consolidation_external_judgement_pending',
        RR74_warning={'one_year_probability':.9431205575309726,'interpretation':'scientific_warning; no universal cap or central legacy multiplier approved'},
        pending_components=['adult_net_hazard_transport','0-4','5-14','remote_tail','risk_effects','mortality_TB_death_scope',
            'recurrence_importation','migration','pediatric_and_BCG_assay_performance','validation_acceptance_thresholds'])
