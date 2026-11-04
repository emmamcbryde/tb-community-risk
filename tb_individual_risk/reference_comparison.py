"""Publication reproduction stages and controlled arithmetic comparisons, not validation."""
import hashlib
import json
import math
from pathlib import Path
from .provenance import seal
from .menzies_2008 import linear_rate_projection
from .risk_factors import combine_rr
from .contracts import ValidationError,finite

PUBLISHED_EXAMPLES_SHA256="b3ac7586a9c1086aaeaa94b80e27f20c06b40524f465ed5ec216de23e6bb6e43"


def published_examples():
    raw=(Path(__file__).parent/"reference_data/published_examples_2008.json").read_bytes()
    if hashlib.sha256(raw).hexdigest()!=PUBLISHED_EXAMPLES_SHA256:raise ValidationError("examples","Published fixture checksum mismatch")
    return json.loads(raw)


def _resolution(percent_text):
    digits=len(percent_text.split(".")[1]) if "." in percent_text else 0
    return 10**(-digits)/100


def _comparison(published,reproduced,printed_text):
    difference=None if reproduced is None else abs(reproduced-published)
    return dict(published_value=published,reproduced_value=reproduced,absolute_difference=difference,
        relative_difference=None if difference is None or published==0 else difference/abs(published),
        declared_tolerance=_resolution(printed_text),tolerance_policy="one_printed_resolution_unit_partial_stage_only_rounding_rule_unverified",
        within_declared_tolerance=None if difference is None else difference<=_resolution(printed_text)+1e-15)


def reproduce_published_examples(package):
    source=published_examples();rows=[]
    for entry in source["table3"]:
        row=dict(entry);p=float(entry["published_ppv_percent"])/100;published=float(entry["published_lifetime_percent"])/100
        # This is a conditional progression-stage check using published PPV as an INPUT.
        # It is never an independent reconstruction of the PPV or a whole-algorithm match.
        if entry["group"] in {"healthy","diabetes"}:
            rr=1 if entry["group"]=="healthy" else package.value("rr_diabetes")
            reproduced=p*package.value("baseline_annual_rate")*rr*(package.value("assumed_lifetime_age")-entry["age"])
            reason="Original ARI/NTM/BCG mapping absent; using printed PPV only to test the progression equation"
        else:
            reproduced=None
            reason="General factor-combination rule unverified" if entry["group"]=="granuloma_underweight" else "Published age80 equation, early-period footnote and contact Table3 do not reconcile"
        row.update(_comparison(published,reproduced,entry["published_lifetime_percent"]));row.update(full_algorithm_reproduced=False,ppv_reproduced=None,comparison_scope="partial_progression_only_using_published_PPV",unresolved_implementation_detail=reason)
        if entry["group"]=="granuloma_underweight":
            row["explicit_modern_combination_math_only"]={policy:dict(rr=combine_rr({"granuloma":2,"underweight":3},policy)["effective_rr"],lifetime_raw=p*.001*combine_rr({"granuloma":2,"underweight":3},policy)["effective_rr"]*45) for policy in ("highest_applicable_rr","additive_excess_rr","multiplicative_rr")}
            row["pair_specific_implied_rr_not_adopted"]=published/(p*.001*45)
        if entry["group"]=="recent_contact":
            row["literal_initial_annual_times45_raw"]=p*.001*15*45
            row["hypothetical_tau0_footnote_schedule_raw"]=p*(.001*15*2+.001*43)
            row["footnote_schedule_difference"]=abs(row["hypothetical_tau0_footnote_schedule_raw"]-published)
        rows.append(row)
    figure=source["figure2"];p=float(figure["published_ppv_percent"])/100
    for name,computed in (("annual",p*.001*package.value("rr_diabetes")),("lifetime",p*.001*package.value("rr_diabetes")*(80-figure["age"]))):
        printed=figure["published_"+name+"_percent"]
        rows.append(dict(figure,example_id=figure["example_id"]+"_"+name,**_comparison(float(printed)/100,computed,printed),full_algorithm_reproduced=False,ppv_reproduced=None,comparison_scope="partial_progression_only_using_published_PPV",unresolved_implementation_detail="Original ARI/NTM mapping and rounding rule unavailable; original UI combined never-BCG with unknown"))
    return seal(dict(schema_id="individual_tb.menzies_reproduction_report.v1",source=source["source_citation"],parameter_package=package.package_id,parameter_package_hash=package.record()["scientific_hash"],source_pdf_sha256=source["source_pdf_sha256"],fictional=True,rows=rows,full_algorithm_matches=0,partial_comparisons=sum(r["reproduced_value"] is not None for r in rows),partial_within_printed_resolution=sum(r.get("within_declared_tolerance") is True for r in rows),not_clinical_validation=True))


def formulation_comparisons():
    rows=[]
    # Predeclared mathematical cases. Same annual number interpreted as a hazard is
    # a calculation-only sensitivity, not a transfer into contemporary parameters.
    for name,age,p,rate,m,t in [("low_short",35,.5,.001,1,1),("low_long",35,.5,.001,1,45),("moderate",35,.5,.01,1,20),("high_rr_short",35,1,.001,74,1),("high_rr_long",35,1,.001,74,45),("near80",79,1,.001,1,1),("above80_literal",85,1,.001,1,-5)]:
        raw=p*rate*m*t;exponential=None if t<0 else p*(-math.expm1(-rate*m*t))
        rows.append(dict(scenario=name,age=age,p_tbi=p,annual_number=rate,effective_rr=m,horizon=t,historical_linear_raw=raw,
            linear_probability=raw if t>=0 and 0<=rate*m*t<=1 else None,exponential_math_sensitivity=exponential,
            difference=None if exponential is None else raw-exponential,interpretation="controlled_arithmetic_only_not_adopted_natural_history"))
    for policy in ("highest_applicable_rr","additive_excess_rr","multiplicative_rr"):
        m=combine_rr({"diabetes":3.6,"smoking":3},policy)["effective_rr"]
        rows.append(dict(scenario="multiple_"+policy,age=35,p_tbi=1,annual_number=.001,effective_rr=m,horizon=45,historical_linear_raw=.001*m*45,
            linear_probability=None,exponential_math_sensitivity=-math.expm1(-.001*m*45),interpretation="historical_combination_unverified_modern_policy_arithmetic_only"))
    rows.append(dict(scenario="recent_tau0_five_year",historical_literal_initial_rate_raw=.015*5,footnote_schedule_raw=.015*2+.001*3,
        exponential_schedule_math_sensitivity=-math.expm1(-(.015*2+.001*3)),published_integration_unresolved=True))
    return seal(dict(schema_id="individual_tb.formulation_comparison.v1",rows=rows,not_clinical=True,contemporary_adoption=False))
