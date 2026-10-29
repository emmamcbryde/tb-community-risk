"""Versioned forward-only primitives; no inverse calibration or application hooks."""
from copy import deepcopy
from decimal import Decimal
from math import isfinite

from engine.country_tb.data_common import canonical,seal,profile_hash,decimal,numeric_text,DataValidationError
from engine.country_tb.composition_common import require,metadata
from engine.country_tb.calibration_common import scientific_context

IMPLEMENTATION='country_tb_progression_forward_1.0.0'
EXCLUSIONS=['baseline_active_TB','recurrent_TB_after_treatment','relapse','treatment_failure',
    'imported_active_TB_after_baseline','migration','endogenous_transmission','interventions','economics',
    'historical_TB_depletion','infection_selective_historical_mortality']


def number(value,*,probability=False):
    d=decimal(value,nonnegative=False)
    v=float(d)
    require(isfinite(v),'Number is outside binary64 numerical range.')
    if probability:require(0<=d<=1,'Probability outside [0,1]; no clipping.')
    return v


def num(value):return numeric_text(value)


def artifact(family,**fields):
    return seal({'schema_id':family,'contract_version':'1.0.0','implementation_version':IMPLEMENTATION,
        'purpose':'planning_and_sequencing_not_denying_care',**deepcopy(fields)})


def verify(value,family):
    require(isinstance(value,dict) and value.get('schema_id')==family and value.get('contract_version')=='1.0.0'
        and value.get('implementation_version')==IMPLEMENTATION,'Unknown forward contract/version; migration required.')
    require(value.get('profile_hash')==profile_hash(value),'Forward scientific hash mismatch.')
    return value


# Positive nodes/weights for deterministic 8-point Gauss-Legendre integration.
GAUSS=((.1834346424956498,.3626837833783620),(.5255324099163290,.3137066458778873),
       (.7966664774136267,.2223810344533745),(.9602898564975363,.1012285362903763))


def quadrature(lo,hi):
    center=(lo+hi)/2;half=(hi-lo)/2
    return sorted((center+sign*half*x,half*w) for x,w in GAUSS for sign in (-1,1))
