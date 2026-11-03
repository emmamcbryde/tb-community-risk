"""First-event conditional progression through cumulative hazards."""
import math
from .contracts import finite, ValidationError
from .infection_probability import require_component, Unavailable
from .parameters import STATUSES


def validate_knots(knots):
    if not knots or tuple(knots[0]) != (0, 0):
        raise ValidationError("knots", "Cumulative-hazard curve must start at (0,0).")
    points = [(finite(t,"knot time"),finite(h,"cumulative hazard")) for t,h in knots]
    if len(points) < 2 or any(t1 <= t0 or h1 < h0 for (t0,h0),(t1,h1) in zip(points,points[1:])):
        raise ValidationError("knots", "Times strictly increasing and cumulative hazards nondecreasing required.")
    return points


def cumulative_hazard(t, knots, tail_policy="external_required", tail_hazard=None):
    t = finite(t, "duration")
    points = validate_knots(knots)
    for (a,ha),(b,hb) in zip(points,points[1:]):
        if t <= b:
            return ha + (t-a)*(hb-ha)/(b-a)
    if tail_policy == "constant_final_segment":
        a,ha = points[-2];b,hb=points[-1];rate=(hb-ha)/(b-a)
    elif tail_policy == "external_constant_hazard" and tail_hazard is not None:
        rate=finite(tail_hazard,"remote hazard")
    else:
        raise Unavailable("Remote tail beyond evidence window is unresolved, not zero hazard.")
    return points[-1][1] + rate*(t-points[-1][0])


def hazard_segments(tau, horizon, knots, tail_policy, tail_hazard=None):
    points=validate_knots(knots)
    end=tau+horizon
    if not math.isfinite(end):
        raise Unavailable("Prediction interval exceeds finite numeric range.")
    # Force explicit coverage check even when no interior knots occur.
    cumulative_hazard(end,points,tail_policy,tail_hazard)
    boundaries=sorted({tau,end}|{t for t,_ in points if tau<t<end})
    segments=[]
    for a,b in zip(boundaries,boundaries[1:]):
        rate=None
        for (x,hx),(y,hy) in zip(points,points[1:]):
            if x <= a < y:
                rate=(hy-hx)/(y-x);break
        if rate is None:
            rate=(points[-1][1]-points[-2][1])/(points[-1][0]-points[-2][0]) if tail_policy=="constant_final_segment" else finite(tail_hazard,"remote hazard")
        if not math.isfinite(rate):raise Unavailable("Segment hazard exceeds finite numeric range; no silent approximation")
        segments.append((b-a,rate))
    return segments


def hazard_to_risk(hazard, multiplier=1):
    h=finite(hazard,"cumulative hazard");m=finite(multiplier,"multiplier",minimum=1)
    return -math.expm1(-m*h)


def horizon_for(values,parameters):
    cfg=require_component(parameters.component("prediction_horizons"))
    selected=values.get("horizon_years")
    if selected is None:selected=cfg.get("selected_years")
    if selected is None:raise Unavailable("Prediction horizon external_parameter_required")
    t=finite(selected,"horizon_years")
    if "maximum_years" not in cfg or t>finite(cfg["maximum_years"],"maximum horizon"):
        raise Unavailable("Requested horizon outside explicitly supported prediction horizons")
    return t


def conditional_progression(values,parameters,multiplier):
    horizon=horizon_for(values,parameters)
    cfg=require_component(parameters.component("baseline_progression_hazard"))
    for name in cfg.get("required_history", []):
        if values.get(name) is None:raise Unavailable(name+": history unknown")
        if values[name] is not False:
            raise Unavailable(name+": history-specific natural history external_parameter_required; ordinary first-event model unavailable")
    age=values.get("current_age");tau=values.get("time_since_infection")
    if age is None or tau is None:raise Unavailable("Age and time since most recent infection require values, priors or bounds")
    age=finite(age,"current_age");tau=finite(tau,"time_since_infection",maximum=age)
    curves=[c for c in cfg.get("age_curves",[]) if float(c["age_min"])<=age<float(c["age_max"])]
    if len(curves)!=1:raise Unavailable("No unique age-appropriate progression curve; adult curve not borrowed for children")
    curve=curves[0]
    if curve.get("status") not in STATUSES or curve.get("status")=="external_parameter_required" or not curve.get("source"):
        raise Unavailable("Age-specific curve external_parameter_required")
    recent=require_component(parameters.component("recent_infection_progression"))
    if recent.get("clock")!="most_recent_infection":raise Unavailable("Most-recent-infection progression clock must be explicitly selected")
    window=finite(recent.get("recent_window_years"),"recent window")
    remote=require_component(parameters.component("remote_infection_progression"))
    tail=remote.get("tail_policy","external_required")
    segments=hazard_segments(tau,horizon,curve["knots"],tail,remote.get("annual_hazard"))
    age_cfg=require_component(parameters.component("age_effects"))
    if age_cfg.get("policy")=="none_explicit":age_multiplier=1.0
    elif age_cfg.get("policy")=="age_band_multiplier":
        bands=[b for b in age_cfg.get("bands",[]) if float(b["age_min"])<=age<float(b["age_max"])]
        if len(bands)!=1:raise Unavailable("Age effect lacks unique applicability")
        age_multiplier=finite(bands[0]["multiplier"],"age multiplier",minimum=1)
    else:raise Unavailable("Age effects external_parameter_required")
    m=multiplier*age_multiplier
    if not math.isfinite(m):raise Unavailable("Combined age/effective multiplier overflow; no cap")
    mortality=require_component(parameters.component("competing_mortality"))
    mode=mortality.get("mode")
    if mode=="none_explicit":mu=0.0
    elif mode=="constant_external_hazard":mu=finite(mortality.get("annual_hazard"),"death hazard")
    else:raise Unavailable("Competing mortality model external_parameter_required")
    disease=death=0.0;survival=1.0
    for duration,rate in segments:
        h=m*rate;total=h+mu
        if math.isinf(h):raise Unavailable("Scaled segment hazard exceeds finite numeric range; no cap or approximation")
        if math.isinf(total):
            scale=max(h,mu);denominator=h/scale+mu/scale
            event=-math.expm1(-(scale*duration)*denominator)
            disease+=survival*(h/scale)/denominator*event
            death+=survival*(mu/scale)/denominator*event
            survival*=1-event
            continue
        exit_probability=-math.expm1(-total*duration)
        if total:
            disease+=survival*(h/total)*exit_probability
            death+=survival*(mu/total)*exit_probability
        survival*=math.exp(-total*duration)
    return {"probability":disease,"competing_death_probability":death,"surviving_disease_free_probability":survival,
            "baseline_cumulative_hazard_interval":math.fsum(d*r for d,r in segments),
            "age_multiplier":age_multiplier,"effective_hazard_multiplier":m,"horizon_years":horizon,
            "time_since_most_recent_infection":tau,"baseline_recency":"recent" if tau<window else "remote",
            "mortality_policy":mode,"curve_source":curve["source"],
            "warnings":["Competing mortality not modelled; long-horizon disease risk can be biased upward."] if mu==0 else []}
