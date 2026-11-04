"""Immutable, hash-pinned reference packages. No network and no implicit adoption."""
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from .contracts import ValidationError
from .provenance import canonical,seal

PINS = {'menzies_2008_published_reference_v1': '36642676423aec352af2fb448effcfe245a803afd32a17558186ba915acbed22', 'synthetic_demonstration_v1': '715ce4ffe46c542b64456a25574ae8cd5381525de249c27db6ee2622e194ab43', 'contemporary_parameters_external_required_v1': 'a72c446521d6b28639d4fa502d66af2b93ca5664ac7180d1572bef480dce99ed'}

@dataclass(frozen=True)
class ParameterPackage:
    _json: str

    def __post_init__(self):
        value=json.loads(self._json)
        if value.get("schema_id")!="individual_tb.parameter_package.v1" or seal(value)["scientific_hash"]!=value.get("scientific_hash"):
            raise ValidationError("package","Invalid package contract/hash")
        package_id=value.get("package_id")
        if package_id not in PINS:raise ValidationError("package","Unknown immutable package")
        official=(Path(__file__).parent/"reference_data"/(package_id+".json")).read_bytes()
        if hashlib.sha256(official).hexdigest()!=PINS[package_id] or canonical(value)!=canonical(json.loads(official)):
            raise ValidationError("package","Immutable package changed; create a new version instead")

    @property
    def package_id(self):return json.loads(self._json)["package_id"]

    def record(self):return json.loads(self._json)

    def value(self,name):
        entries=[r for r in self.record()["records"] if r["name"]==name]
        if len(entries)!=1:raise ValidationError("parameter", "No unique package parameter: "+name)
        return entries[0]["value"]


def load_parameter_package(name):
    if name not in PINS:raise ValidationError("package","Unknown package; no fallback")
    raw=(Path(__file__).parent/"reference_data"/(name+".json")).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=PINS[name]:raise ValidationError("package","Immutable package checksum mismatch")
    return ParameterPackage(raw.decode("utf-8"))
