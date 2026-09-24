"""Demonstration working-default population profile.

The demonstration profile reuses the currently validated engine parameter values
so that the general application is runnable out of the box. These values are not
representative of any particular country and must be reviewed and replaced with
locally applicable evidence.
"""

from __future__ import annotations

from dataclasses import replace
from functools import lru_cache
from typing import Any

from engine.profiles.population_profile import (
    AgeBand,
    EffectMeasure,
    IncidenceData,
    Location,
    LocationKind,
    PopulationProfile,
    Provenance,
    RiskFactor,
    TrendSpec,
    bundled_value,
)


DEMONSTRATION_PROFILE_ID = "demonstration-working-defaults"
DEMONSTRATION_PROFILE_LABEL = "Demonstration working defaults"
DEMONSTRATION_PROFILE_VERSION = "demonstration_working_defaults_v1"
GENERAL_DEFAULT_POPULATION = 10_000
GENERAL_DEFAULT_STOCHASTIC_SIMULATIONS = 1_000
DEMONSTRATION_WARNING = (
    "These are demonstration working defaults, not evidence for any particular country. "
    "Review and replace the epidemiological, risk-factor, intervention and cost assumptions "
    "with locally applicable evidence before using results for planning."
)
ENGINE_EFFECT_APPLICATION = "progression_hazard_multiplier"
DEMONSTRATION_SOURCE = "Bundled demonstration input table"

# (profile id, engine key, label, age-specific prevalence key in engine parameters)
ENGINE_RISK_FACTORS: tuple[tuple[str, str, str, str], ...] = (
    ("contact", "contact", "Contact with an infectious TB case", "contactPrevByAge"),
    ("smoking", "smoking", "Current smoking", "smokingPrevByAge"),
    ("cannabis", "MJ", "Cannabis use", "mjPrevByAge"),
    ("renal", "renal", "Renal impairment", "renalPrevByAge"),
    ("diabetes", "diabetes", "Diabetes", "diabetesPrevByAge"),
    ("chronic_lung_disease", "cld", "Chronic lung disease", "cldPrevByAge"),
    ("alcohol_drugs", "alcohol", "Harmful alcohol or drug use", "alcoholPrevByAge"),
)
ENGINE_KEYS = tuple(item[1] for item in ENGINE_RISK_FACTORS)
BROAD_AGE_BANDS = (("0-4 years", 0, 4), ("5-14 years", 5, 14), ("15+ years", 15, None))
DEFAULT_LTBI_CALIBRATION_TARGET = 47 / 624


def build_demonstration_profile() -> PopulationProfile:
    """Return a fresh demonstration profile (no user overrides)."""
    return _cached_demonstration_profile()


@lru_cache(maxsize=1)
def _cached_demonstration_profile() -> PopulationProfile:
    pars = _engine_default_parameters()
    pop_frac = [float(value) for value in pars["popFrac"]]
    age_bands = tuple(
        AgeBand(
            label=label,
            lower_age=lower,
            upper_age=upper,
            proportion=bundled_value(fraction, "proportion", "Bundled demonstration age distribution"),
        )
        for (label, lower, upper), fraction in zip(BROAD_AGE_BANDS, pop_frac)
    )
    risk_factors = []
    for factor_id, engine_key, label, prev_key in ENGINE_RISK_FACTORS:
        by_age = [float(value) for value in pars[prev_key]]
        prevalence = sum(value * weight for value, weight in zip(by_age, pop_frac))
        risk_factors.append(
            RiskFactor(
                risk_factor_id=factor_id,
                label=label,
                enabled=True,
                prevalence=bundled_value(
                    prevalence,
                    "proportion",
                    DEMONSTRATION_SOURCE,
                    notes="Population-weighted summary of age-specific demonstration prevalence.",
                ),
                effect_estimate=bundled_value(
                    float(pars["disOR"][engine_key]),
                    "ratio",
                    DEMONSTRATION_SOURCE,
                ),
                effect_measure=EffectMeasure.OR,
                engine_key=engine_key,
                engine_application=ENGINE_EFFECT_APPLICATION,
                evidence_source=DEMONSTRATION_SOURCE,
            )
        )
    return PopulationProfile(
        profile_id=DEMONSTRATION_PROFILE_ID,
        name=DEMONSTRATION_PROFILE_LABEL,
        location=Location(name=DEMONSTRATION_PROFILE_LABEL, kind=LocationKind.DEMONSTRATION),
        population_size=bundled_value(float(GENERAL_DEFAULT_POPULATION), "people", "General application default"),
        age_distribution=age_bands,
        age_distribution_source="Bundled demonstration age distribution",
        ltbi_prevalence=bundled_value(
            DEFAULT_LTBI_CALIBRATION_TARGET,
            "proportion",
            "Bundled demonstration calibration target",
        ),
        incidence=IncidenceData(
            source="No incidence series linked",
            snapshot_id=None,
            provenance=Provenance.BUNDLED,
        ),
        trend=TrendSpec(method="not_estimated"),
        risk_factors=tuple(risk_factors),
        data_vintage=DEMONSTRATION_PROFILE_VERSION,
        notes=DEMONSTRATION_WARNING,
        demonstration=True,
    )


def demonstration_profile_without_risk_factors() -> PopulationProfile:
    """Demonstration profile with no risk-factor stratification."""
    return replace(
        build_demonstration_profile(),
        profile_id=f"{DEMONSTRATION_PROFILE_ID}-no-risk-factors",
        risk_factors=(),
    )


def _engine_default_parameters() -> dict[str, Any]:
    from engine.apy.data import load_parameters_from_config
    from engine.apy.working_defaults import build_unified_working_default_preset

    return load_parameters_from_config(build_unified_working_default_preset()["config"])
