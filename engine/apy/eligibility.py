from __future__ import annotations

from typing import Any


RESTRICTED_ELIGIBILITY_ERROR = (
    "Restricted eligibility is configured but no eligibility-selection rule is implemented."
)


def resolve_eligibility(config: dict[str, Any]) -> dict[str, Any]:
    population = float(config.get("N"))
    scenario = config.get("scenario") if isinstance(config.get("scenario"), dict) else {}
    eligible = scenario.get("eligible") if isinstance(scenario.get("eligible"), dict) else {}
    number = eligible.get("number")
    proportion = eligible.get("proportion")
    explicit_all = False
    if number is None and proportion is None:
        resolved = population
        proportion_value = 1.0
        explicit_all = True
    elif number is not None:
        resolved = float(number)
        proportion_value = resolved / population if population else 0.0
    else:
        proportion_value = float(proportion)
        resolved = population * proportion_value
    if resolved < 0 or resolved - population > 1e-9:
        raise ValueError("eligible population must be between 0 and N.")
    if abs(proportion_value - 1.0) > 1e-9:
        raise ValueError(RESTRICTED_ELIGIBILITY_ERROR)
    return {
        "number": resolved,
        "proportion": proportion_value,
        "selectionRule": "all_population" if explicit_all or abs(proportion_value - 1.0) <= 1e-9 else None,
    }


def screening_coverage_of_population(config: dict[str, Any], eligible_number: float) -> float:
    population = float(config.get("N"))
    coverage = float(config.get("screenCoverage"))
    requested = eligible_number * coverage
    return min(requested, eligible_number) / population if population else 0.0
