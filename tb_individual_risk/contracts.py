"""Versioned person records with explicit observation states."""
from copy import deepcopy
from dataclasses import dataclass, field
from datetime import date
import math
from numbers import Real

PERSON_SCHEMA = "individual_tb.person_record.v1"
STATES = {"unknown", "known_absent", "known_present", "not_applicable", "excluded", "invalid"}
BOOLEAN_FIELDS = {
    "close_contact", "recent_exposure", "recent_conversion", "bcg",
    "previous_active_tb", "previous_preventive_treatment", "active_tb_at_baseline",
    "hiv", "hiv_treatment", "diabetes", "smoking", "underweight", "silicosis",
    "haemodialysis", "severe_renal_disease", "tnf_inhibitor", "transplantation",
    "immunosuppression", "head_neck_carcinoma",
}
NUMBER_FIELDS = {"current_age", "age_at_immigration", "immigration_year", "birth_year",
                 "time_since_infection", "tst_mm", "bcg_age", "horizon_years",
                 "tbi_probability", "pretest_tbi_probability"}
TEXT_FIELDS = {"country_of_birth", "current_country", "igra_product", "observation_date"}
CATEGORIES = {
    "sex": {"female", "male", "other"},
    "tst_result": {"positive", "negative", "not_done"},
    "igra_result": {"positive", "negative", "indeterminate", "borderline", "not_done"},
    "cxr": {"normal", "granuloma", "fibronodular", "other"},
}
FIELDS = BOOLEAN_FIELDS | NUMBER_FIELDS | TEXT_FIELDS | set(CATEGORIES)

class ValidationError(ValueError):
    def __init__(self, field, message):
        self.field = field
        super().__init__(message)
    def as_dict(self):
        return {"field": self.field, "message": str(self), "code": "invalid_record"}


def finite(value, field_name, minimum=0, maximum=None):
    if isinstance(value, bool):
        raise ValidationError(field_name, "Boolean is not a numeric value.")
    try:
        number = float(value)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValidationError(field_name, "Expected a finite number.") from exc
    if not math.isfinite(number) or number < minimum or (maximum is not None and number > maximum):
        raise ValidationError(field_name, "Number outside permitted finite range; no clipping.")
    return number


@dataclass(frozen=True)
class Observation:
    state: str = "unknown"
    value: object = None
    source: object = None

    @property
    def known(self):
        return self.state in {"known_present", "known_absent"}

    @classmethod
    def parse(cls, raw, name):
        if isinstance(raw, cls):
            raw = dict(state=raw.state, value=raw.value, source=raw.source)
        if isinstance(raw, dict):
            if set(raw) - {"state", "value", "source"}:
                raise ValidationError(name, "Unknown observation keys.")
            state = raw.get("state", "unknown")
            value = raw.get("value")
            source = deepcopy(raw.get("source"))
        elif raw is None or (isinstance(raw, Real) and math.isnan(float(raw))):
            state, value, source = "unknown", None, None
        elif isinstance(raw, str) and raw in STATES:
            state, value, source = raw, None, None
        elif isinstance(raw, bool):
            state, value, source = ("known_present" if raw else "known_absent"), raw, None
        else:
            state, value, source = "known_present", raw, None
        if state not in STATES:
            raise ValidationError(name, "Unknown observation state.")
        if state not in {"known_present", "known_absent"}:
            if value is not None:
                raise ValidationError(name, "Non-known states cannot contain a value.")
            return cls(state, None, source)
        if name in BOOLEAN_FIELDS or name.startswith("factor:"):
            expected = state == "known_present"
            if value is None:
                value = expected
            if not isinstance(value, bool) or value != expected:
                raise ValidationError(name, "Boolean value contradicts present/absent state.")
        elif state == "known_absent":
            if name in {"tst_result", "igra_result"} and value in {None, "not_done"}:
                value = "not_done"
            else:
                raise ValidationError(name, "Use an explicit known value for non-boolean fields.")
        elif name in NUMBER_FIELDS:
            limit = 1 if name in {"tbi_probability", "pretest_tbi_probability"} else None
            number = finite(value, name, maximum=limit)
            if name in {"birth_year", "immigration_year"} and number != int(number):
                raise ValidationError(name, "Calendar year must be integral.")
            # Preserve supplied numeric strings without display rounding.
        elif name in CATEGORIES:
            if not isinstance(value, str) or value not in CATEGORIES[name]:
                raise ValidationError(name, "Invalid categorical value.")
        elif name in TEXT_FIELDS:
            if not isinstance(value, str) or not value.strip():
                raise ValidationError(name, "Expected nonempty text.")
            if name == "observation_date":
                try:
                    date.fromisoformat(value)
                except ValueError as exc:
                    raise ValidationError(name, "Use ISO observation date YYYY-MM-DD.") from exc
        return cls(state, value, source)


@dataclass(frozen=True)
class PersonRecord:
    record_id: object
    fields: dict
    original: dict = field(repr=False)
    schema_id: str = PERSON_SCHEMA
    identifier_state: str = "known_present"

    @classmethod
    def from_dict(cls, raw):
        if not isinstance(raw, dict):
            raise ValidationError("record", "Expected a person dictionary.")
        if raw.get("schema_id", PERSON_SCHEMA) != PERSON_SCHEMA:
            raise ValidationError("schema_id", "Unsupported person contract.")
        if set(raw) - FIELDS - {"record_id", "factors", "schema_id"}:
            raise ValidationError("record", "Unexpected fields: " + ", ".join(sorted(set(raw) - FIELDS - {"record_id", "factors", "schema_id"})))
        identifier = raw.get("record_id")
        identifier_state="known_present"
        if isinstance(identifier,dict):
            if set(identifier)-{"state","value","source"}:raise ValidationError("record_id","Unknown identifier observation keys")
            identifier_state=identifier.get("state","unknown");value=identifier.get("value")
            if identifier_state not in STATES or identifier_state=="invalid":raise ValidationError("record_id","Invalid identifier state")
            if identifier_state!="known_present" and value is not None:raise ValidationError("record_id","Non-present identifier state has a value")
            identifier=value
        elif identifier is None or identifier=="":identifier_state="unknown";identifier=None
        if identifier_state=="known_present" and identifier is None:raise ValidationError("record_id","Known-present identifier requires a value")
        if identifier is not None and (isinstance(identifier,bool) or not isinstance(identifier,(str,int)) or isinstance(identifier,str) and not identifier.strip()):
            raise ValidationError("record_id","Use a string/integer identifier or explicit unknown state; no generated replacement")
        fields = {name: Observation.parse(raw.get(name), name) for name in sorted(FIELDS)}
        factors = raw.get("factors", {})
        if not isinstance(factors, dict):
            raise ValidationError("factors", "Expected factor-ID mapping.")
        for name, value in factors.items():
            if not isinstance(name, str) or not name or name in BOOLEAN_FIELDS:
                raise ValidationError("factors", "Custom factor IDs must be distinct from built-in fields.")
            fields["factor:" + name] = Observation.parse(value, "factor:" + name)
        for name, value in fields.items():
            if value.state == "invalid":
                raise ValidationError(name, "Explicit invalid observation cannot enter prediction.")
        def number(name):
            v = fields[name]
            return float(v.value) if v.known else None
        age, immigration, tau = number("current_age"), number("age_at_immigration"), number("time_since_infection")
        if age is not None and any(x is not None and x > age for x in (immigration, tau, number("bcg_age"))):
            raise ValidationError("current_age", "Immigration/vaccination age or infection duration exceeds current age.")
        observed = fields["observation_date"]
        if observed.known:
            year = date.fromisoformat(observed.value).year
            for name in ("birth_year", "immigration_year"):
                if number(name) is not None and number(name) > year:
                    raise ValidationError(name, "Year is after observation date.")
            if number("birth_year") is not None and number("immigration_year") is not None and number("immigration_year") < number("birth_year"):
                raise ValidationError("immigration_year", "Immigration precedes birth.")
            if age is not None and number("birth_year") is not None and not -1 < year - number("birth_year") - age <= 1:
                raise ValidationError("birth_year", "Age and birth year are inconsistent with observation year.")
        return cls(identifier, fields, deepcopy(raw),identifier_state=identifier_state)

    def values(self):
        return {k: v.value if v.known else None for k, v in self.fields.items()}
