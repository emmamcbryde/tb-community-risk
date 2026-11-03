"""Canonical, finite JSON and content-addressed scientific records. No I/O/network."""
from copy import deepcopy
from dataclasses import asdict, is_dataclass
import hashlib
import json

VERSION = "individual_tb_risk_1.0.0"

def plain(value):
    if is_dataclass(value):
        return plain(asdict(value))
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    return value

def canonical(value):
    return json.dumps(plain(value), sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False, allow_nan=False)

def scientific_hash(value):
    return hashlib.sha256(canonical(value).encode("utf-8")).hexdigest()

def seal(value):
    result = deepcopy(plain(value))
    result.pop("scientific_hash", None)
    result["scientific_hash"] = scientific_hash(result)
    return result


def source_safe(value):
    """Preserve non-JSON input anomalies as explicit source markers, never scores."""
    import math
    from numbers import Real
    if isinstance(value, Real) and not isinstance(value, bool) and not math.isfinite(float(value)):
        return {"original_nonfinite_number": repr(value)}
    if isinstance(value, dict): return {str(k):source_safe(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)): return [source_safe(v) for v in value]
    if value is None or isinstance(value,(str,bool,int,float)):return value
    if isinstance(value,Real):return float(value)
    return {"original_python_type":type(value).__name__,"representation":repr(value)}
