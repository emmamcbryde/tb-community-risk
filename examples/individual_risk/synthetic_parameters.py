"""ALL NUMBERS INVENTED. Algorithm demonstration only, NEVER clinical parameters."""
from tb_individual_risk import ModelParameters,Parameter,PopulationPrior
from tb_individual_risk.parameters import COMPONENTS


def synthetic_parameters():
    label="synthetic_not_clinical"
    source="invented_M1_arithmetic_demonstration_no_empirical_claim"
    def parameter(name,config):return Parameter(name,label,source,config)
    effects={name:dict(field=name,rr=rr,definition=definition,age_min=15,age_max=150,status=label,source=source,effect_measure="invented_RR_as_hazard_approximation")
             for name,rr,definition in [("diabetes",2.0,"Invented diabetes state"),("smoking",3.0,"Invented explicitly defined smoking exposure"),("haemodialysis",5.0,"Invented haemodialysis state"),("hiv",8.0,"Invented untreated HIV state"),("transplantation",12.0,"Invented transplantation state")]}
    effects["hiv"]["required_context"]={"hiv_treatment":False}
    configs={
        "infection_probability":{"mode":"bayes_assays","assays":["igra"]},
        "tst_interpretation":{"sensitivity":.8,"specificity":.85,"positive_threshold_mm":10},
        "igra_interpretation":{"sensitivity":.8,"specificity":.9,"product":"SYNTHETIC_ASSAY"},
        "baseline_progression_hazard":{"age_curves":[dict(age_min=15,age_max=150,knots=[[0,0],[2,.04],[10,.08],[30,.10]],status=label,source=source)],"required_history":["active_tb_at_baseline","previous_active_tb","previous_preventive_treatment"]},
        "recent_infection_progression":{"clock":"most_recent_infection","recent_window_years":2},
        "remote_infection_progression":{"tail_policy":"external_constant_hazard","annual_hazard":.001},
        "age_effects":{"policy":"none_explicit"},
        "radiographic_effects":{"factors":{"granuloma":dict(field="cxr",present_value="granuloma",rr=2.0,definition="Invented radiographic granuloma",age_min=15,age_max=150,status=label,source=source,effect_measure="invented_RR_as_hazard_approximation"),"fibronodular":dict(field="cxr",present_value="fibronodular",rr=4.0,definition="Invented fibronodular finding",age_min=15,age_max=150,status=label,source=source,effect_measure="invented_RR_as_hazard_approximation")}},
        "medical_risk_effects":{"factors":effects,"explicit_scope_exclusions":{name:"Not modelled in this invented demonstration; no clinical inference" for name in ("underweight","silicosis","severe_renal_disease","tnf_inhibitor","immunosuppression","head_neck_carcinoma")}},
        "competing_mortality":{"mode":"none_explicit"},"prediction_horizons":{"selected_years":5,"maximum_years":20},
    }
    distributions={"diabetes":[dict(value=False,probability=.8),dict(value=True,probability=.2)],
                   "smoking":[dict(value=False,probability=.7),dict(value=True,probability=.3)],
                   "pretest_tbi_probability":[dict(value=.1,probability=.5),dict(value=.3,probability=.5)],
                   "time_since_infection":[dict(value=1,probability=.4),dict(value=10,probability=.6)]}
    prior=PopulationPrior(source=source,status=label,distributions=distributions,
                          bounds={k:[r["value"] for r in rows] for k,rows in distributions.items()},dependence_policy="explicit_independence")
    return ModelParameters(components={n:parameter(n,configs[n]) for n in COMPONENTS},population_prior=prior,classification=label)


def fictional_people():
    base=dict(current_age=40,sex="female",current_country="FICTIONAL",country_of_birth="FICTIONAL",
        pretest_tbi_probability=.2,igra_result="positive",igra_product="SYNTHETIC_ASSAY",time_since_infection=1,
        diabetes=False,smoking=False,haemodialysis=False,hiv=False,hiv_treatment=False,transplantation=False,cxr="normal",
        active_tb_at_baseline=False,previous_active_tb=False,previous_preventive_treatment=False,
        horizon_years=5,observation_date="2026-11-03")
    variants=[{},dict(current_age=25,sex="male",smoking=True),dict(current_age=70,diabetes=True,time_since_infection=20),
        dict(current_age=55,smoking=True,diabetes=True),dict(current_age=60,haemodialysis=True),
        dict(current_age=30,igra_result="negative",pretest_tbi_probability=.05),dict(current_age=19,recent_conversion=True,close_contact=True,time_since_infection=.2),
        dict(current_age=50,diabetes=None,smoking=None),dict(current_age=35,pretest_tbi_probability=None,time_since_infection=None),
        dict(current_age=80,transplantation=True,hiv=True),dict(current_age=8),dict(current_age=45,previous_active_tb=True)]
    return [dict(base,**changes,record_id="fictional_%02d"%(i+1)) for i,changes in enumerate(variants)]
