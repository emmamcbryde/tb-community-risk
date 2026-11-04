"""Publication-only historical arithmetic. No website code and no clinical defaults."""
from copy import deepcopy
from dataclasses import asdict
import math
from .contracts import PersonRecord,ValidationError,finite
from .infection_probability import Unavailable
from .provenance import seal,source_safe

FORMULATION="menzies_2008_compatibility_v1"
BOUNDARY="Baseline infection only: no future acquisition, reinfection, transmission, migration, recurrence or changing risk factors."
FACTOR_FIELDS=("diabetes","underweight","smoking","silicosis","haemodialysis","hiv","transplantation","head_neck_carcinoma","tnf_inhibitor")


def _known(record,name):
    obs=record.fields[name]
    if not obs.known:raise Unavailable(name+": "+obs.state+"; no historical blank-field default applied")
    return obs.value


def prevalence_from_ari(age,ari,contact_prevalence):
    age=finite(age,"age");ari=finite(ari,"ARI probability",maximum=1);c=finite(contact_prevalence,"contact prevalence",maximum=1)
    if age==0:return c
    return 1.0 if ari==1 else -math.expm1(math.log1p(-c)+age*math.log1p(-ari)) if c<1 else 1.0


def ppv_from_components(true_positive,fp_bcg,fp_ntm):
    tp=finite(true_positive,"TP",maximum=1);bcg=finite(fp_bcg,"FP BCG",maximum=1);ntm=finite(fp_ntm,"FP NTM",maximum=1)
    if tp+bcg+ntm==0:raise Unavailable("Zero TP+FP denominator; no invented PPV")
    return tp/(tp+bcg+ntm)


def tst_ppv(record,package):
    v=record.values();age=finite(_known(record,"current_age"),"age")
    if age<18:raise Unavailable("Adult-only publication; conservative software scope excludes under-18s, no pediatric extrapolation")
    if _known(record,"tst_result")!="positive" or float(_known(record,"tst_mm"))<10:
        raise Unavailable("Publication implements positive TST >=10mm; negative, not-done or 5-9mm interpretation unavailable")
    _known(record,"historical_exposure_policy")
    # Immigration correction is not specified completely; never silently truncate exposure.
    if v.get("age_at_immigration") not in (None,0) or v.get("immigration_year") is not None:
        raise Unavailable("Historical immigration/exposure mapping unavailable; no residence-based substitution")
    contact=_known(record,"historical_contact_type")
    if v.get("close_contact") is not None and v["close_contact"] != (contact=="close"):
        raise Unavailable("Close-contact flag conflicts with explicit historical contact type")
    if v.get("historical_ari_probability") is not None and v.get("historical_smear_positive_incidence_per_100k") is not None:
        raise Unavailable("Select exactly one supplied historical ARI or incidence input")
    if v.get("historical_ari_probability") is not None:ari=float(v["historical_ari_probability"]);ari_method="supplied_historical_ARI_probability"
    elif v.get("historical_smear_positive_incidence_per_100k") is not None:
        ari=float(v["historical_smear_positive_incidence_per_100k"])/package.value("styblo_incidence_per_percent_ari")/100
        ari=finite(ari,"historical ARI probability",maximum=1);ari_method="historical_incidence_div49_percent_then_div100"
    else:raise Unavailable("Original country/state ARI lookup absent; supply a provenance-bearing historical input")
    c=0.0 if contact=="none" else package.value("close_contact_prevalence" if contact=="close" else "casual_contact_prevalence")
    tp=prevalence_from_ari(age,ari,c)
    if _known(record,"bcg") is True:
        # Table1 aggregates >=10mm; methods describe size-specific differences not provided numerically.
        _known(record,"bcg_age")
        raise Unavailable("BCG aggregate .085/.418 documented, but exact 10-14/>=15mm mapping unresolved; no current-website value imported")
    if float(v["tst_mm"])>=15:fp_ntm=0.0;ntm_method="published_ge15_NTM_crossreaction_rule"
    else:
        fp_ntm=package.value("ntm_cross_reactivity")*float(_known(record,"historical_ntm_sensitization_prevalence"));ntm_method=".027_times_supplied_NTM_sensitization"
    return dict(probability=ppv_from_components(tp,0.0,fp_ntm),p_ltbi=tp,ari=ari,contact_term=c,fp_bcg=0.0,fp_ntm=fp_ntm,
                method="published_PPVTST_equation",ari_method=ari_method,ntm_method=ntm_method)


def factor_summary(record,package):
    v=record.values();present={}
    for field in FACTOR_FIELDS:
        if _known(record,field) is not True:continue
        if field=="smoking" and float(_known(record,"smoking_packs_per_day"))<=1:
            raise Unavailable("Current smoking is not the published >1 pack/day definition")
        if field=="underweight" and float(_known(record,"percent_ideal_body_weight"))>=90:
            raise Unavailable("Underweight input incompatible with <90% ideal body weight")
        if field=="tnf_inhibitor" and str(_known(record,"tnf_inhibitor_product")).lower()!="infliximab":
            raise Unavailable("Table2 infliximab effect not automatically applied to other TNF inhibitors")
        if field=="transplantation" and _known(record,"immunosuppression") is not True:
            raise Unavailable("Published transplant effect requires immunosuppressant therapy")
        present[field]=package.value("rr_"+field)
    cxr=_known(record,"cxr")
    if cxr=="other":raise Unavailable("Unclassified radiograph has no verified historical effect")
    if cxr!="normal":present[cxr]=package.value("rr_"+cxr)
    if v.get("severe_renal_disease") is True and v.get("haemodialysis") is not True:
        raise Unavailable("Severe renal disease is not automatically haemodialysis")
    if v.get("immunosuppression") is True and v.get("transplantation") is not True:
        raise Unavailable("Unspecified immunosuppression is not the transplant exposure definition")
    known_custom=[k for k,o in record.fields.items() if k.startswith("factor:") and o.known and o.value is True]
    if known_custom:raise Unavailable("Unverified historical custom effects: "+", ".join(known_custom))
    contact=_known(record,"historical_contact_type");tau=v.get("time_since_infection")
    conversion=_known(record,"recent_conversion")
    if tau is not None and float(tau)==package.value("early_period_years"):
        raise Unavailable("Exactly two-year recency boundary differs between Table2 wording and footnote; no undocumented convention")
    if conversion and tau is not None and float(tau)>=package.value("early_period_years"):
        raise Unavailable("Confirmed recent conversion conflicts with supplied infection duration")
    recent=contact=="close" or conversion or tau is not None and float(tau)<package.value("early_period_years")
    all_effects=dict(present)
    if recent:all_effects["recent_infection_or_close_contact"]=package.value("early_progression_rr")
    # Table3 suggests a sum for one pair; that does not verify a general rule.
    effective=next(iter(all_effects.values())) if len(all_effects)==1 else 1.0 if not all_effects else None
    return dict(individual_rr=present,all_effects_including_historical_recency=all_effects,recent=bool(recent),effective_rr=effective,
                combination_policy="historical_single_factor_only_general_rule_unresolved",factors_present=sorted(all_effects))


def linear_rate_projection(ppv,annual_conditional_rate,horizon):
    p=finite(ppv,"PPV",maximum=1);rate=finite(annual_conditional_rate,"annual conditional rate");t=finite(horizon,"horizon")
    conditional_raw=rate*t;absolute_raw=p*conditional_raw
    return dict(conditional_raw=conditional_raw,absolute_raw=absolute_raw,
                conditional_probability=conditional_raw if 0<=conditional_raw<=1 else None,
                absolute_probability=absolute_raw if 0<=conditional_raw<=1 and 0<=absolute_raw<=1 else None)


def predict_reference(person,package):
    base=dict(schema_id="individual_tb.prediction.v1",classification="historical_only",clinical_use_approved=False,
        parameter_hash=package.record()["scientific_hash"],parameter_record=package.record(),
        combination_policy="historical_single_factor_only_general_rule_unresolved",
        central_absolute_risk=None,lower_missing_scenario=None,upper_missing_scenario=None,
        lower_conditional_missing_scenario=None,upper_conditional_missing_scenario=None,
        tbi_probability=None,conditional_progression_risk=None,effective_rr=None,expected_effective_rr=None,
        annual_risk=None,two_year_risk=None,five_year_risk=None,historical_cumulative_risk_to_age80=None,
        principal_risk_drivers=[],warnings=["Historical publication reference only; not a contemporary clinical recommendation.",BOUNDARY],
        unresolved_assumptions=[],validation_errors=[],historical_arithmetic={},horizon_years=None)
    try:
        r=PersonRecord.from_dict(person.original if isinstance(person,PersonRecord) else person)
    except (ValidationError,TypeError,ValueError) as exc:
        base.update(record_id=person.get("record_id") if isinstance(person,dict) and isinstance(person.get("record_id"),(str,int)) else None,status="invalid",original_record=source_safe(person),validation_errors=[exc.as_dict() if isinstance(exc,ValidationError) else {"field":"record","message":str(exc)}])
        return seal(base)
    v=r.values();base.update(record_id=r.record_id,identifier_state=r.identifier_state,original_record=source_safe(r.original),
        input_states={k:asdict(o) for k,o in r.fields.items()},missing_critical_variables=[],input_completeness={"known_fields":sum(o.known for o in r.fields.values()),"total_fields":len(r.fields)})
    issues=[];infection=None;factors=None
    try:infection=tst_ppv(r,package);base["tbi_probability"]=infection["probability"];base["infection_calculation"]=infection
    except (Unavailable,ValidationError) as exc:issues.append(str(exc))
    try:factors=factor_summary(r,package);base["factor_calculation"]=factors;base["effective_rr"]=base["expected_effective_rr"]=factors["effective_rr"]
    except (Unavailable,ValidationError) as exc:issues.append(str(exc))
    for field in ("active_tb_at_baseline","previous_active_tb","previous_preventive_treatment"):
        if not r.fields[field].known:issues.append(field+": "+r.fields[field].state+"; first-event/history scope unresolved")
        elif v[field] is not False:issues.append(field+": recurrence/treatment/current-disease model not provided")
    if v.get("igra_result") not in (None,"not_done"):
        base["warnings"].append("IGRA is not implemented in the 2008 paper; any available result here is TST-only, not combined assay evidence.")
    horizon=v.get("horizon_years");age=v.get("current_age")
    if horizon is not None:base["horizon_years"]=float(horizon)
    else:issues.append("horizon_years unknown; no undocumented historical default")
    if factors and factors["effective_rr"] is None:issues.append("General historical combination rule unverified; concurrent factors not forced into a modern policy")
    if factors and factors["effective_rr"] is not None:
        rate=package.value("baseline_annual_rate")*factors["effective_rr"]
        base["historical_arithmetic"]["annual_conditional_rate"]=rate
        base["principal_risk_drivers"]=[dict(factor_id=k,reported_rr=value) for k,value in factors["all_effects_including_historical_recency"].items()]
        if factors["effective_rr"]>=20:base["warnings"].append("Extreme historical RR; large linear risks are a scientific warning, not a numerical validation.")
        if infection:
            ppv=infection["probability"];base["annual_risk"]=ppv*rate
            def project(t):
                if factors["recent"]:
                    tau=v.get("time_since_infection")
                    if tau is None:return None,"Recent infection duration is unknown; remaining early period unresolved"
                    remaining=max(0,package.value("early_period_years")-float(tau))
                    raw=ppv*(rate*min(t,remaining)+package.value("baseline_annual_rate")*max(0,t-remaining))
                    base["historical_arithmetic"]["footnote_schedule_raw_"+str(t)]=raw
                    if t>remaining:return None,"Early-period/lifetime integration conflicts with published Table3; raw schedule shown, not adopted"
                result=linear_rate_projection(ppv,rate,t)
                base["historical_arithmetic"]["linear_projection_"+str(t)]=result
                return result,None if result["absolute_probability"] is not None else "Historical linear probability outside [0,1]; raw arithmetic retained without clipping"
            for key,t in (("two_year_risk",2.0),("five_year_risk",5.0)):
                projected,error=project(t)
                if projected:base[key]=projected["absolute_probability"]
                if error:base["unresolved_assumptions"].append(key+": "+error)
            if age is not None:
                lifetime=package.value("assumed_lifetime_age")-float(age)
                base["historical_arithmetic"]["cumulative_to_age80_literal_raw"]=ppv*rate*lifetime
                if lifetime<0:issues.append("Age above80: negative historical lifetime duration; no clipping to zero")
                elif factors["recent"]:
                    project(lifetime);base["unresolved_assumptions"].append("cumulative_age80: equation, footnote and Table3 do not reconcile")
                else:
                    result=linear_rate_projection(ppv,rate,lifetime);base["historical_cumulative_risk_to_age80"]=result["absolute_probability"]
                    if result["absolute_probability"] is None:base["warnings"].append("Linear lifetime arithmetic exceeds probability bounds; result unavailable, raw value retained")
            if horizon is not None and age is not None and float(age)<=80 and float(horizon)<=80-float(age):
                projected,error=project(float(horizon))
                if projected:base["central_absolute_risk"]=projected["absolute_probability"];base["conditional_progression_risk"]=projected["conditional_probability"]
                if error:issues.append(error)
            elif horizon is not None:issues.append("Prediction horizon outside the historical age80 lifetime scope")
    if issues:base["central_absolute_risk"]=base["conditional_progression_risk"]=None
    if any(message.startswith(("active_tb_at_baseline:","previous_active_tb:","previous_preventive_treatment:")) for message in issues):
        for key in ("annual_risk","two_year_risk","five_year_risk","historical_cumulative_risk_to_age80"):base[key]=None
    base["unresolved_assumptions"]+=issues;base["missing_critical_variables"]=sorted({message.split(":")[0] for message in issues})
    base["warnings"]+=issues;base["warnings"]=sorted(set(base["warnings"]))
    base["status"]="calculated_historical_reference" if base["central_absolute_risk"] is not None else "unavailable"
    base["warning_status"]="warning";base["assumptions"]=[BOUNDARY,"No undocumented historical missing-field defaults.","Linear arithmetic retained; no exponential replacement."]
    return seal(base)
