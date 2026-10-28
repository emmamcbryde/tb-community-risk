"""Shared primitives for isolated reinfection-capable catalytic calibration."""
from __future__ import annotations

from copy import deepcopy
from decimal import Context, Decimal, ROUND_HALF_EVEN, localcontext

from engine.country_tb.data_common import DataValidationError, decimal, numeric_text, profile_hash, seal

CONTRACT_VERSION = "1.0.0"
IMPLEMENTATION_VERSION = "country_tb_catalytic_calibration_1.0.0"
EXCLUSIONS = (
    "clearance_or_reversion", "migration",
    "infection_selective_mortality", "active_tb_exits", "endogenous_transmission",
)
STATE_FIELDS = ("never_infected", "recent", "remote_only")
SEMANTICS = "nonhomogeneous_poisson_repeated_events_no_protection_no_clearance"
REINFECTION_POLICY = "same_pressure_regardless_of_previous_infection"
CLOCK_POLICY = "most_recent_reinfection_resets_clock_no_active_tb_progression_link"
QUANTUM = Decimal("1e-40")


def scientific_context():
    """Fresh fixed context: caller precision, rounding and flags cannot leak in."""
    return localcontext(Context(prec=80, rounding=ROUND_HALF_EVEN))


def reviewed(value, *, required=True):
    if not isinstance(value, dict) or value.get("status") not in {"reviewed", "not_reviewed"}:
        raise DataValidationError("Explicit reviewed/not_reviewed metadata required.")
    if required and value["status"] != "reviewed":
        raise DataValidationError("This operation requires reviewed scientific inputs.")
    if not value.get("source_id") or not value.get("rationale"):
        raise DataValidationError("Review metadata requires source_id and rationale.")
    return deepcopy(value)


def verify(value, schema_id=None):
    if not isinstance(value, dict) or value.get("profile_hash") != profile_hash(value):
        raise DataValidationError("Scientific artifact hash mismatch.")
    if schema_id is not None and value.get("schema_id") != schema_id:
        raise DataValidationError("Unexpected scientific artifact family.")
    if schema_id is not None:
        expected = "1.1.0" if schema_id in {
            "country_tb.incidence_history.v1", "country_tb.demography.v1",
            "country_tb.baseline_data_profile.v1",
        } else CONTRACT_VERSION
        if value.get("contract_version") != expected:
            raise DataValidationError("Explicit contract migration required; design proposals are not executable.")
        if expected == CONTRACT_VERSION and (
            value.get("exclusions") != list(EXCLUSIONS)
            or value.get("infection_state_semantics") != SEMANTICS
            or value.get("disease_forward_operator") is not None
        ):
            raise DataValidationError("Incompatible model semantics/exclusions.")
    return value


def artifact(schema_id, **fields):
    return seal({
        "schema_id": schema_id, "contract_version": CONTRACT_VERSION,
        "implementation_version": IMPLEMENTATION_VERSION,
        "exclusions": list(EXCLUSIONS), "infection_state_semantics": SEMANTICS,
        "disease_forward_operator": None,
        "numerical_policy": {"working_decimal_precision": 80, "state_decimal_places": 40, "rounding": "ROUND_HALF_EVEN",
                             "probability_repair": "forbidden"},
        "interpretation": "conditional_reinfection_capable_catalytic_diagnostic",
        "reinfection_policy": REINFECTION_POLICY,
        "progression_clock_policy": CLOCK_POLICY,
        **deepcopy(fields),
    })


def probability(value):
    result = decimal(value)
    if result > 1:
        raise DataValidationError("Probability must be within [0,1]; no clipping.")
    return result


def states_from_survivors(never, no_recent, prior_remote_plus_recent=0):
    with scientific_context():
        never, no_recent = probability(never), probability(no_recent)
        if no_recent < never:
            raise DataValidationError("Inconsistent cumulative hazards/survivors.")
        never = never.quantize(QUANTUM)
        no_recent = no_recent.quantize(QUANTUM)
        recent, remote = 1 - no_recent, no_recent - never
        return {
            "never_infected": numeric_text(never),
            "recent": numeric_text(recent),
            "remote_only": numeric_text(remote),
            "prior_remote_plus_recent": numeric_text(probability(prior_remote_plus_recent).quantize(QUANTUM)),
            "any_tbi": numeric_text(1 - never),
        }


def validate_states(row):
    with scientific_context():
        values = [probability(row[key]) for key in STATE_FIELDS]
        if sum(values) != 1 or probability(row["any_tbi"]) != values[1] + values[2]:
            raise DataValidationError("Infection states do not reconcile exactly.")
        if not 0 <= probability(row["prior_remote_plus_recent"]) <= values[1]:
            raise DataValidationError("Prior remote plus recent must be a subset of recent.")


def finite_integer(value, *, minimum=0):
    if type(value) is not int or value < minimum:
        raise DataValidationError("Expected an explicit integer.")
    return value
