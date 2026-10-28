"""Exact expected-value arithmetic for isolated baseline risk allocation."""
from copy import deepcopy
from decimal import Decimal
from fractions import Fraction

from engine.country_tb.calibration_common import scientific_context
from engine.country_tb.composition_common import require,metadata,age_bounds,contains,overlaps,sex_overlaps
from engine.country_tb.data_common import decimal,numeric_text,seal,profile_hash,canonical,digest

VERSION='1.0.0'
IMPLEMENTATION='country_tb_risk_allocation_1.0.0'


def ratio(value, *, probability=False):
    if isinstance(value,Fraction):
        result=value
    elif isinstance(value,dict) and set(value)=={'value','exact'}:
        exact=value['exact']
        require(isinstance(exact,dict) and set(exact)=={'numerator','denominator'},'Exact rational fields required.')
        require(isinstance(exact['numerator'],str) and isinstance(exact['denominator'],str),'Rational integers must be strings.')
        try:result=Fraction(int(exact['numerator']),int(exact['denominator']))
        except (ValueError,ZeroDivisionError):require(False,'Invalid rational number.')
        require(amount(result)==value,'Decimal display does not match exact rational value.')
    else:
        result=Fraction(decimal(value,nonnegative=False))
    if probability:require(0<=result<=1,'Prevalence must lie in [0,1]; no clipping or percentage conversion.')
    return result


def text_number(value):
    f=ratio(value)
    # Preserve every supplied decimal digit, independently of caller context.
    # Non-terminating derived values always also retain an exact rational pair.
    denominator=f.denominator;a=b=0
    while denominator%2==0:denominator//=2;a+=1
    while denominator%5==0:denominator//=5;b+=1
    if denominator==1:
        scale=max(a,b);integer=abs(f.numerator)*2**(scale-a)*5**(scale-b)
        digits=str(integer).rjust(scale+1,'0')
        text=(digits[:-scale]+'.'+digits[-scale:]).rstrip('0').rstrip('.') if scale else digits
        return ('-' if f<0 else '')+text
    with scientific_context():return numeric_text(Decimal(f.numerator)/Decimal(f.denominator))


def amount(value):
    f=value if isinstance(value,Fraction) else ratio(value)
    with scientific_context():display=numeric_text(Decimal(f.numerator)/Decimal(f.denominator))
    return {'value':display,'exact':{'numerator':str(f.numerator),'denominator':str(f.denominator)}}


def artifact(family, **fields):
    return seal({'schema_id':family,'contract_version':VERSION,'implementation_version':IMPLEMENTATION,
        'purpose':'planning_and_sequencing_not_denying_care',
        'numerical_policy':'exact_rational_probabilities_and_expected_counts_decimal80_display',
        **deepcopy(fields)})


def verify(value,family):
    require(isinstance(value,dict) and value.get('schema_id')==family,'Unexpected risk-allocation contract family.')
    require(value.get('contract_version')==VERSION and value.get('implementation_version')==IMPLEMENTATION,'Explicit risk contract migration required.')
    require(value.get('profile_hash')==profile_hash(value),'Risk-allocation scientific hash mismatch.')
    return value


def cell_key(row):return (row['sex'],row['age_lower'],row['age_upper_exclusive'])


def clinical_group(row):
    lo,hi=age_bounds(row)
    if hi is not None and hi<=5:return '0-4'
    if lo>=5 and hi is not None and hi<=15:return '5-14'
    if lo>=15:return '15+'
    return 'unresolved_crosses_clinical_age_boundary'


def matches(record,cell):
    return contains(record,cell) and (record['sex']=='both' or record['sex']==cell['sex'])


def validate_rows(rows):
    require(isinstance(rows,list),'An explicit record list is required.')
    for row in rows:
        require(isinstance(row,dict) and {'age_lower','age_upper_exclusive','sex'}<=row.keys(),'Age/sex bounds must be explicit.')
        age_bounds(row);require(row['sex'] in ('female','male','both'),'Unknown sex.')


def validate_disjoint(rows):
    validate_rows(rows)
    for i,a in enumerate(rows):
        for b in rows[i+1:]:
            require(not(sex_overlaps(a['sex'],b['sex']) and overlaps(a,b)),'Overlapping records at the same input scope.')
