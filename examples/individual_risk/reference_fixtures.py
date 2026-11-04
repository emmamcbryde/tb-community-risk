"""FICTIONAL reference people. Historical outputs are not clinical recommendations."""
from copy import deepcopy

def fictional_reference_people():
    invented={"citation":"invented_exposure_input_for_software_reference_only","review_status":"synthetic_not_clinical"}
    base=dict(record_id="young_adult",current_age=25,sex="female",tst_result="positive",tst_mm=15,
        bcg=False,historical_ari_probability={"state":"known_present","value":.001,"source":invented},
        historical_ntm_sensitization_prevalence={"state":"known_present","value":.1,"source":invented},
        historical_exposure_policy="constant_full_age_exposure_confirmed",historical_contact_type="none",close_contact=False,
        recent_conversion=False,time_since_infection=10,diabetes=False,underweight=False,smoking=False,
        silicosis=False,haemodialysis=False,hiv=False,hiv_treatment=False,transplantation=False,head_neck_carcinoma=False,
        tnf_inhibitor=False,cxr="normal",active_tb_at_baseline=False,previous_active_tb=False,previous_preventive_treatment=False,
        horizon_years=5,igra_result="not_done",pretest_tbi_probability=.2,igra_product="SYNTHETIC_ASSAY")
    variants=[{},dict(record_id="older_adult",current_age=70,time_since_infection=50),
        dict(record_id="bcg_infancy",bcg=True,bcg_age=.5),dict(record_id="bcg_after_infancy",bcg=True,bcg_age=6),
        dict(record_id="negative_igra_only",tst_result="not_done",tst_mm=None,igra_result="negative"),
        dict(record_id="positive_igra_only",tst_result="not_done",tst_mm=None,igra_result="positive"),
        dict(record_id="close_contact",historical_contact_type="close",close_contact=True,time_since_infection=0),
        dict(record_id="recent_conversion",recent_conversion=True,time_since_infection=0),
        dict(record_id="diabetes",diabetes=True),dict(record_id="underweight",underweight=True,percent_ideal_body_weight=85),
        dict(record_id="heavy_smoking",smoking=True,smoking_packs_per_day=1.5),dict(record_id="silicosis",silicosis=True),
        dict(record_id="haemodialysis",haemodialysis=True),dict(record_id="hiv",hiv=True),
        dict(record_id="transplant",transplantation=True,immunosuppression=True),dict(record_id="fibronodular",cxr="fibronodular"),
        dict(record_id="granuloma",cxr="granuloma"),dict(record_id="head_neck_carcinoma",head_neck_carcinoma=True),
        dict(record_id="infliximab",tnf_inhibitor=True,tnf_inhibitor_product="infliximab"),
        dict(record_id="multiple_factors",diabetes=True,smoking=True,smoking_packs_per_day=1.5),
        dict(record_id="unknown_bcg",bcg=None),dict(record_id="unknown_diabetes",diabetes=None),
        dict(record_id="child",current_age=8,time_since_infection=1),dict(record_id="near80",current_age=79),
        dict(record_id="above80",current_age=85),dict(record_id="tst10_14",tst_mm=12)]
    return [dict(deepcopy(base),**change) for change in variants]
