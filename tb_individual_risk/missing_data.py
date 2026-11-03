"""Exact finite marginalisation under explicit population priors."""
from itertools import product
import math
from .contracts import Observation, ValidationError, NUMBER_FIELDS
from .infection_probability import Unavailable
from .parameters import EXTERNAL
from .provenance import canonical


def validated_value(name,value):
    obs=Observation.parse(value,name)
    if not obs.known:raise ValidationError(name,"Prior support states must be known values.")
    return obs.value


def equivalent(name,a,b):
    return float(a)==float(b) if name in NUMBER_FIELDS else a==b


def scenarios(person,required,prior,*,bounds=False,budget=8192):
    values=person.values()
    unknown=sorted(name for name in required if values.get(name) is None)
    eligible=[name for name in unknown if person.fields.get(name,Observation()).state=="unknown"]
    weights=[(1.0,{},values)]
    unresolved=set(unknown)
    if prior.status==EXTERNAL:return weights,unresolved
    if prior.dependence_policy=="explicit_joint":
        rows=prior.joint_bound_states if bounds else prior.joint_states
        if not rows:return weights,unresolved
        joint_fields=set(rows[0]["values"])
        compatible=[]
        for row in rows:
            assignment={k:validated_value(k,v) for k,v in row["values"].items()}
            if any(values.get(k) is not None and not equivalent(k,values[k],v) for k,v in assignment.items()):continue
            if any(k in unknown and k not in eligible for k in assignment):continue
            updates={k:v for k,v in assignment.items() if k in eligible}
            w=1.0 if bounds else float(row["probability"])
            if w>0:compatible.append((w,updates,dict(values,**updates)))
        if not compatible:raise Unavailable("No joint prior/scenario support is compatible with known inputs.")
        if len(compatible)>budget:raise Unavailable("Explicit scenario budget exceeded; no rows dropped.")
        total=math.fsum(w for w,_,_ in compatible)
        weights=[(w/total,a,v) for w,a,v in compatible]
        unresolved-=joint_fields & set(eligible)
        return weights,unresolved
    supports=[]
    for name in eligible:
        rows=prior.bounds.get(name) if bounds else prior.distributions.get(name)
        if not rows:continue
        if bounds:
            support=[(validated_value(name,v),1.0) for v in rows]
        else:
            support=[(validated_value(name,r["value"]),float(r["probability"])) for r in rows if float(r["probability"])>0]
        if len({canonical(v) for v,_ in support})!=len(support):
            raise ValidationError(name,"Duplicate prior support states.")
        supports.append((name,support));unresolved.discard(name)
    count=math.prod(len(s) for _,s in supports)
    if count>budget:raise Unavailable("Explicit scenario budget exceeded; no independence approximation or row truncation applied.")
    weights=[]
    for combination in product(*(support for _,support in supports)):
        updates={name:v for (name,_),(v,_) in zip(supports,combination)}
        weight=math.prod(w for _,w in combination)
        if weight>0:weights.append((weight,updates,dict(values,**updates)))
    total=math.fsum(w for w,_,_ in weights)
    return [(w/total,a,v) for w,a,v in weights],unresolved
