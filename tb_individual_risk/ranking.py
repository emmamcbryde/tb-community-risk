"""Deterministic risk ranking with visible incompleteness and overlap."""
from collections import Counter
from copy import deepcopy
import math
from .contracts import ValidationError,finite
from .line_list import records
from .provenance import seal,canonical


def rank(model,line_list,*,by="unconditional_absolute_risk",likely_infected_threshold=None,horizon_years=None):
    if by not in {"unconditional_absolute_risk","conditional_progression_risk","expected_cases_per_1000"}:
        raise ValidationError("ranking","Unknown ranking metric.")
    if by=="conditional_progression_risk":
        if likely_infected_threshold is None:raise ValidationError("ranking","Select an explicit likely-infected threshold for conditional ranking.")
        threshold=finite(likely_infected_threshold,"likely infected threshold",maximum=1)
    source=records(line_list)
    if horizon_years is not None:
        finite(horizon_years,"ranking horizon")
        source=[dict(row,horizon_years=horizon_years) if isinstance(row,dict) else row for row in source]
    results=[model.predict(row) for row in source]
    lower_key,upper_key=("lower_conditional_missing_scenario","upper_conditional_missing_scenario") if by=="conditional_progression_risk" else ("lower_missing_scenario","upper_missing_scenario")
    counts=Counter(r.get("record_id") for r in results if r.get("record_id") is not None)
    horizons={r.get("horizon_years") for r in results if r.get("central_absolute_risk") is not None}
    incomparable=len(horizons)>1
    for i,row in enumerate(results):
        row["source_row_index"]=i
        row["ranking_metric"]=by
        row["warnings"]=list(row.get("warnings",[]))
        if row.get("record_id") is None:row["warnings"].append("Person identifier is unknown; source_row_index is an audit locator, not an invented person ID.")
        if row.get("record_id") is not None and counts[row["record_id"]]>1:
            row["status"]="invalid";row["central_absolute_risk"]=None
            row.setdefault("validation_errors",[]).append({"field":"record_id","code":"duplicate_identifier","message":"Duplicate identifiers retained as invalid; resolve membership."})
        score=row.get("central_absolute_risk")
        if by=="conditional_progression_risk":
            score=row.get("conditional_progression_risk") if row.get("tbi_probability") is not None and row["tbi_probability"]>=threshold else None
        if incomparable:
            score=None;row["warnings"].append("Mixed prediction horizons are not ranked. Select a common horizon explicitly.")
        if row.get("status")=="invalid":score=None
        row["ranking_score"]=score
        row["expected_cases_per_1000"]=None if row.get("central_absolute_risk") is None else 1000*row["central_absolute_risk"]
        row["rank"]=None
        if by=="expected_cases_per_1000" and score is not None:row["ranking_score"]=score*1000
        row["overlapping_missing_information_ranges"]=False
        row["ranking_uncertain"]=row.get(lower_key) is None or row.get(upper_key) is None
    # For each interval, overlap exists with a preceding prefix maximum or next lower bound.
    intervals=sorted((r[lower_key],r[upper_key],i) for i,r in enumerate(results)
        if r.get(lower_key) is not None and r.get(upper_key) is not None)
    prefix=-math.inf
    for j,(lo,hi,i) in enumerate(intervals):
        overlap=prefix>=lo or (j+1<len(intervals) and intervals[j+1][0]<=hi)
        results[i]["overlapping_missing_information_ranges"]=overlap
        results[i]["ranking_uncertain"]|=overlap
        prefix=max(prefix,hi)
    results.sort(key=lambda r:(r["ranking_score"] is None,-r["ranking_score"] if r["ranking_score"] is not None else 0,canonical(r.get("record_id")),r["source_row_index"]))
    previous=object();assigned=None
    for position,row in enumerate(results,1):
        score=row["ranking_score"]
        if score is not None:
            if score!=previous:assigned=position
            row["rank"]=assigned;previous=score
        if row["ranking_uncertain"]:row["warnings"].append("Ranking is uncertain or missing-information ranges overlap; adjacent ranks need not be meaningfully different.")
        row["warning_status"]="warning" if row["warnings"] else "none"
        row.pop("scientific_hash",None)
    return [seal(r) for r in results]
