"""Pure configuration, source-preserving CSV boundary and stale-result rules."""
from copy import deepcopy
from decimal import Decimal,localcontext
import csv
import io
import json

from engine.country_tb.data_common import canonical,digest,seal,DataValidationError,decimal,numeric_text
from country_model_review_ui import APP_VERSION

PAGES = ('Set up population','Baseline TBI','Risk factors','Natural history and scientific gates',
         'Run diagnostic','Results and validation','Evidence and unresolved decisions')
SOURCE_MODES = ('Observed or user-specified TBI prevalence','Estimate TBI from active-TB incidence')
GATE_ITEMS = ('Adult progression','Pediatric progression or reviewed exclusion','Remote-infection tail',
    'Reinfection/reset policy','Mortality source','Risk-factor effect policy','Assay observation model',
    'Migration and historical exposure','Imported disease','Prior TB and recurrence',
    'Compatible active-TB observation','Incidence observation period','Unresolved-population threshold',
    'Numerical refinement','Root uniqueness and identifiability')
GATE_STATES = ('Resolved and reviewed','Selected for sensitivity only','Unresolved','Incompatible','Blocked')
TBI_COLUMNS = ('record_id','subgroup_id','age_lower','age_upper_exclusive','sex','any_tbi','recent','remote_only',
    'never_infected','prior_remote_plus_recent','allocation_policy','source','review_status','geography_id','baseline_year')
POP_COLUMNS = ('age_lower','age_upper_exclusive','sex','population','source','review_status','geography_id','baseline_year')
RISK_COLUMNS = ('record_id','subgroup_id','factor_id','age_lower','age_upper_exclusive','sex','state','prevalence',
    'allocation_policy','supplied_exposure_definition','source','review_status','reason')
ASSOCIATION_COLUMNS = ('control_id','scope','subgroup_id','sex','age_lower','age_upper_exclusive','control','source','review_status')


def require(condition,message):
    if not condition:raise DataValidationError(message)


def default_config():
    return dict(schema_id='country_tb.review_ui_configuration.v1',application_version=APP_VERSION,
        population=dict(country='AUS',scope='National',custom_id='custom:review_population',label='Review population',
            source='Recommended frozen source',baseline_year=2024,size='',scale_source=False,national_proxy=False,
            csv='',subgroup_csv='',subgroups=[],subgroup_proportions=False,partition_confirmed=False),
        tbi=dict(mode=SOURCE_MODES[0],input_kind='Overall prevalence',overall_percent='',derive=False,
            source='Reviewer entry — not scientific approval',review_status='not_reviewed',csv='',hybrid_recency=False,
            conditional_clock=False,pressure='WHO historical shape',proxy_confirmed=False,assay='Direct latent TBI',
            sensitivity_percent='',specificity_percent='',assay_source='',threshold='',recent_percent='',remote_percent='',prior_percent=''),
        risk=dict(records=[],no_factors=False,legacy_effects=False,csv=''),association=dict(global_control='0',csv=''),
        natural=dict(adult='External required',tail='No extrapolation',remote_hazard='',children='External required',
            mortality='No mortality — diagnostic only',reset='Most recent reinfection resets clock',
            horizon='Instantaneous baseline',compare_scenarios=False,active_baseline='Unresolved',
            future='No future infection — baseline cohort only'),
        diagnostic=dict(enabled=False,bundle='',lower='0',upper='.2',grid_intervals=64),
        display=dict(sexes='Separate female and male'),example=None)


def material_hash(config):
    value=deepcopy(config);value.pop('display',None)
    return digest(canonical(value))


def review_binding_hash(config):
    value=deepcopy(config);value['diagnostic']['bundle']=''
    return material_hash(value)


def initialize(session):
    if 'review_config' not in session:session['review_config']=default_config()
    if 'review_timings' not in session:session['review_timings']={}


def stale(session):
    result=session.get('review_result')
    return result is not None and result['configuration_hash']!=material_hash(session['review_config'])


def invalidate_exports(session):
    if session.get('review_exports_hash')!=material_hash(session['review_config']):
        session.pop('review_exports',None);session.pop('review_exports_hash',None)


def percent_input(value,*,required=True):
    if str(value).strip()=='':
        require(not required,'Enter a prevalence; missing is not zero.');return None
    n=decimal(str(value).strip());require(n<=100,'Percentage must be between0 and100; no clipping.')
    with localcontext() as ctx:
        ctx.prec=max(80,len(n.as_tuple().digits)+4)
        return numeric_text(n/100)


def probability(value,*,optional=False):
    if value in ('',None):
        require(optional,'A probability is missing; missing is not zero.');return None
    p=decimal(value);require(p<=1,'CSV probabilities must lie in[0,1], not percentages; no clipping.')
    return numeric_text(p)


def parse_csv(payload,kind,*,geography_id=None,subgroup_ids=None):
    require(kind in {'population','subgroup_population','tbi','risk','association'},'Unknown upload kind.')
    if isinstance(payload,bytes):
        require(len(payload)<=2_000_000,'CSV exceeds the2MB review limit.')
        payload=payload.decode('utf-8-sig')
    require(isinstance(payload,str) and len(payload)<=2_000_000,'UTF-8 CSV required, at most2MB.')
    reader=csv.DictReader(io.StringIO(payload.lstrip('\ufeff')))
    cols={'population':POP_COLUMNS,'subgroup_population':('subgroup_id',*POP_COLUMNS),'tbi':TBI_COLUMNS,'risk':RISK_COLUMNS,'association':ASSOCIATION_COLUMNS}[kind]
    require(reader.fieldnames is not None and len(set(reader.fieldnames))==len(reader.fieldnames) and
        set(reader.fieldnames)==set(cols),'CSV must contain exactly the documented '+kind+' columns.')
    rows=list(reader);require(rows and len(rows)<=5000,'Supply1–5000 rows; no empty upload.')
    originals=deepcopy(rows);ids=set()
    for r in rows:
        require(None not in r and all(v is not None for v in r.values()),'Malformed CSV row.')
        require(r['source'].strip() and r['review_status'] in {'reviewed','not_reviewed','diagnostic'},'Every row needs source and explicit review status.')
        if 'record_id' in r:
            require(r['record_id'] and r['record_id'] not in ids,'Unique nonempty record IDs required.');ids.add(r['record_id'])
        if 'subgroup_id' in r and subgroup_ids is not None:
            require(r['subgroup_id'] in subgroup_ids or kind in {'risk','association'} and r['subgroup_id']=='all',
                'CSV subgroup is outside the selected partition.')
        if 'geography_id' in r:
            require(r['geography_id'] and (geography_id is None or r['geography_id']==geography_id),'CSV geography mismatch; no national-to-subgroup relabelling.')
            require(r['baseline_year']=='2024','Only the frozen2024 baseline is supported.')
            r['baseline_year']=2024
        require(r['sex'] in {'both','female','male'},'Sex must be both,female or male.')
        if kind!='association' or r['scope']!='sex':
            require(r['age_lower'].isdigit(),'Nonnegative integer age_lower required.')
            r['age_lower']=int(r['age_lower']);hi=r['age_upper_exclusive']
            require(hi=='' or hi.isdigit() and int(hi)>r['age_lower'],'Invalid age interval; leave open age_upper_exclusive empty.')
            r['age_upper_exclusive']=None if hi=='' else int(hi)
        else:r.pop('age_lower');r.pop('age_upper_exclusive')
        if kind in {'population','subgroup_population'}:
            r['population']=numeric_text(decimal(r['population']));require(r['sex']!='both','Population CSV must retain separate female/male counts.')
        elif kind=='tbi':
            for k in ('any_tbi','recent','remote_only','never_infected','prior_remote_plus_recent'):r[k]=probability(r[k],optional=True)
            require(r['any_tbi'] is not None or all(r[k] is not None for k in ('recent','remote_only','never_infected')),'TBI total or all exclusive states required.')
        elif kind=='risk':r['prevalence']=probability(r['prevalence'],optional=True)
        else:
            c=decimal(r['control'],nonnegative=False);require(-1<=c<=1,'Association control must lie in[-1,1].');r['control']=numeric_text(c)
    return dict(rows=rows,source_rows=originals,source_sha256=digest(payload.encode('utf-8')),kind=kind)


def csv_text(rows,columns):
    s=io.StringIO(newline='');w=csv.DictWriter(s,fieldnames=columns,extrasaction='ignore',lineterminator='\n')
    w.writeheader();w.writerows({k:'' if r.get(k) is None else r.get(k,'') for k in columns} for r in rows)
    return s.getvalue()


def working_example(*,mixture=False):
    c=default_config();c['example']='Explicitly restored synthetic review example — not Australian estimates'
    c['population'].update(scope='Custom / subnational',custom_id='custom:synthetic_review',label='Synthetic adult review cohort',
        source='User CSV',size='1000',csv=csv_text([dict(age_lower=a,age_upper_exclusive=a+1,sex=s,population='250',
            source='Invented working example',review_status='diagnostic',geography_id='custom:synthetic_review',baseline_year=2024)
            for s in ('female','male') for a in (20,60)],POP_COLUMNS))
    c['tbi'].update(overall_percent='20',derive=True,pressure='Explicit constant-shape diagnostic',proxy_confirmed=True,
        source='Invented working example — no empirical TBI adopted')
    c['risk'].update(no_factors=True)
    c['natural'].update(adult='Central-candidate four-point sensitivity',tail='Fixed remote-hazard sensitivity',remote_hazard='.00038',
        active_baseline='Assume no baseline active TB — diagnostic only')
    if mixture:
        c['population'].update(subgroup_proportions=True,partition_confirmed=True,subgroups=[
            dict(subgroup_id=i,label=l,population=n,tbi_percent=p,source_mode='User supplied',exposure_geography=g,
                arrival_year='',prior_tb_prevalence='',imported_active_contribution='',recurrence_contribution='')
            for i,l,n,p,g in [('locally_born','Locally born — invented','700','5','who:AUS'),
                ('recent_arrivals','Recent arrivals — invented','200','25','who:VNM'),
                ('long_term','Longer-term migrants — invented','100','15','who:IND')]])
    return c


def restore(session,config):
    session['review_config']=deepcopy(config)
    # Widget values from other pages cannot override an explicit restore.
    for k in list(session):
        if k.startswith('edit_'):session.pop(k,None)
    session.pop('review_population',None);session.pop('review_exports',None)
