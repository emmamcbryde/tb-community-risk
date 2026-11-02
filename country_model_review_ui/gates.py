"""Fail-closed UI boundary to the unchanged3C diagnostic inversion."""
from copy import deepcopy
from math import floor,fsum
import json
from engine.country_tb.data_common import seal,canonical,digest,DataValidationError
from engine.country_tb.incidence_mapping import IncidenceMapping
from engine.country_tb.risk_policies import association_policy,control_for_cell
from engine.country_tb.incidence_inversion import inversion_gate,diagnostic_calibration,assess_identifiability,PUBLIC_STATUS
from country_model_review_ui.state import GATE_ITEMS,require,review_binding_hash,material_hash,parse_csv
from country_model_review_ui.model import SOURCE,REVIEW,PROVENANCE,geography,META,catalog,_groups,scaled_rows

DMETA=dict(source=SOURCE,review=REVIEW,provenance_record=PROVENANCE)


def gate_panel(config):
    n=config['natural'];t=config['tbi'];rows=[]
    for label in GATE_ITEMS:
        state='Unresolved';reason='Reviewed inputs and source-bound decisions required.'
        if label=='Adult progression' and n['adult']!='External required':
            state='Selected for sensitivity only';reason='Four-point surrogate is not an approved country curve.'
        elif label=='Pediatric progression or reviewed exclusion':
            reason='Child curves external-required; exclusion needs an explicit denominator,target and threshold review.'
        elif label=='Remote-infection tail':
            if n['tail']!='No extrapolation':state='Selected for sensitivity only';reason='Tail extension is an explicit sensitivity,not a central approval.'
            else:reason='No extrapolation means unresolved progression beyond evidence support,not zero hazard.'
        elif label=='Reinfection/reset policy':state='Resolved and reviewed';reason='Retained approved central policy: most recent reinfection resets the clock.'
        elif label=='Mortality source':
            state='Incompatible' if n['mortality']!='No mortality — diagnostic only' and config['population']['country']!='AUS' else 'Selected for sensitivity only'
            reason='ABS mortality cannot apply to this country.' if state=='Incompatible' else 'Selected mortality is a diagnostic input; country adoption and TB-death scope need review.'
        elif label=='Risk-factor effect policy':
            state='Selected for sensitivity only' if config['risk']['legacy_effects'] or config['risk']['no_factors'] else 'Unresolved'
            reason='Highest applicable RR centrally. Legacy effects and RR1 reference are not approved country risks.'
        elif label=='Assay observation model':
            state='Selected for sensitivity only';reason='Direct latent input or fixed supplied assay operator; transport review is separate.'
        elif label in {'Imported disease','Prior TB and recurrence','Migration and historical exposure'}:
            reason='No imported/recurrent/migration fraction is invented or automatically subtracted.'
        elif label in {'Numerical refinement','Root uniqueness and identifiability'}:
            reason='Bounded native3C scan/refinement must pass after scientific pre-fit review; no root is accepted before these checks.'
        rows.append(dict(item=label,status=state,reason=reason))
    if config['population']['country']=='PRK':
        next(r for r in rows if r['item']=='Compatible active-TB observation').update(status='Blocked',reason='Usable WHO incidence unavailable; no target is invented.')
    return seal(dict(schema_id='country_tb.review_ui_gate_panel.v1',rows=rows,
        enabled=config['diagnostic']['enabled'],prefit_allowed=False,production_available=False,public_status=PUBLIC_STATUS,
        review_binding_hash=review_binding_hash(config)))


def prepare_diagnostic(config,context):
    require(context['population_definition_hash']==digest(canonical(config['population'])),'Population context is stale.')
    require(config['diagnostic']['enabled'],'Incidence diagnostic is disabled by default; no source fallback.')
    require(context['observed_incidence'] is not None and context['observed_incidence']['central'] is not None,
        'Country incidence is unavailable,including PRK; no invented target or inherited calibration.')
    text=config['diagnostic']['bundle'];require(text and len(text)<=2_000_000,'Upload an explicit reviewed diagnostic JSON bundle,at most2MB.')
    try:bundle=json.loads(text)
    except (ValueError,TypeError) as exc:raise DataValidationError('Invalid review-bundle JSON.') from exc
    required={'review_binding_hash','component_reviews','mapping_inputs','observation','unresolved_threshold_approval',
        'scientific_review','concentration_policy','numerical_protocol','assay_decision'}
    require(isinstance(bundle,dict) and required<=bundle.keys(),'Incomplete diagnostic review bundle.')
    require(bundle['review_binding_hash']==review_binding_hash(config),'Review receipt is stale or belongs to different material inputs.')
    for item in GATE_ITEMS:
        r=bundle['component_reviews'].get(item,{})
        require(r.get('status')=='Resolved and reviewed' and all(isinstance(r.get(k),str) and r[k].strip() for k in ('source_id','reviewer','rationale')),
            'Unresolved scientific gate: '+item)
    numerical=bundle['numerical_protocol']
    require(numerical.get('status')=='reviewed' and numerical.get('source_id') and
        numerical.get('policy')=='retained_3C_bounded_grid_all_roots_and_refinement_v1',
        'An explicit reviewed numerical/refinement protocol is required; tolerances are not relaxed in the UI.')
    require(bundle['assay_decision'].get('source_id') and bundle['assay_decision'].get('status')=='reviewed','Explicit assay/latent observation decision required.')
    args=deepcopy(bundle['mapping_inputs'])
    require(args.get('geography_id')==geography(config),'Diagnostic population geography differs from current review inputs.')
    require(args.get('measure')=='instantaneous_baseline_per_100000_person_years',
        'This UI exposes bounded instantaneous inversion only; annual/prospective targets cannot be silently relabelled.')
    require(len(args.get('strata',[]))<=300,'Diagnostic review limit is300explicit strata.')
    require(args.get('population_kind')!='synthetic_closed_population' or config['example'] is not None,
        'A synthetic protocol cannot approve an empirical country population.')
    # Verify every support count against source age/sex cells. Discrete ages
    # inside a source band need their explicit reviewed within-cell policy.
    mapping_counts={};expected={};scope_lo=args.get('age_lower',0);scope_hi=args.get('age_upper_exclusive')
    for row in context['demography']['rows']:
        overlaps=row['age_upper_exclusive'] is None or row['age_upper_exclusive']>scope_lo
        overlaps=overlaps and (scope_hi is None or row['age_lower']<scope_hi)
        fits=row['age_lower']>=scope_lo and (scope_hi is None or row['age_upper_exclusive'] is not None and row['age_upper_exclusive']<=scope_hi)
        require(not overlaps or fits,'Diagnostic scope boundary cuts a source age cell; reviewed age disaggregation is required.')
        if row['age_lower']<scope_lo or scope_hi is not None and (row['age_upper_exclusive'] is None or row['age_upper_exclusive']>scope_hi):continue
        expected[row['sex'],row['age_lower'],row['age_upper_exclusive']]=float(row['population'])
    for row in args['strata']:
        age=row.get('age');require(age is not None,'Unresolved open ages cannot enter this diagnostic as exact ages.')
        candidates=[key for key in expected if key[0]==row['sex'] and key[1]<=float(age) and (key[2] is None or float(age)<key[2])]
        require(len(candidates)==1,'Diagnostic age support does not match source cells.')
        require(row.get('age_support_policy')=='reviewed_within_cell_age_distribution' and row.get('age_distribution_review'),
            'Do not convert demographic age cells to invented birthdays; supply reviewed within-cell support.')
        key=candidates[0];mapping_counts.setdefault(key,[]).append(float(row['population']))
    group_ids={g['subgroup_id'] for g in config['population']['subgroups']} or {'main'}
    require({r['subgroup_id'] for r in args['strata']}==group_ids,'Diagnostic subgroup IDs must match the selected population definition.')
    groups=_groups(config,context);group_cells={};achieved={}
    uploaded_pop=parse_csv(config['population']['subgroup_csv'],'subgroup_population',geography_id=geography(config),subgroup_ids=group_ids)['rows'] if config['population']['subgroup_csv'] else None
    for group in groups:
        if uploaded_pop:cells=[r for r in uploaded_pop if r['subgroup_id']==group['subgroup_id']]
        elif len(groups)==1:cells=context['demography']['rows']
        else:
            require(config['population']['subgroup_proportions'],'Explicit subgroup demographic distributions required for diagnostic inversion.')
            cells,_=scaled_rows(context['demography']['rows'],group['population'])
        for cell in cells:
            key=(cell['sex'],cell['age_lower'],cell['age_upper_exclusive'])
            if key in expected:group_cells[group['subgroup_id'],*key]=float(cell['population'])
    for row in args['strata']:
        candidates=[key for key in group_cells if key[0]==row['subgroup_id'] and key[1]==row['sex'] and key[2]<=float(row['age']) and (key[3] is None or float(row['age'])<key[3])]
        require(len(candidates)==1,'Diagnostic subgroup/age/sex support differs from the declared distribution.')
        achieved.setdefault(candidates[0],[]).append(float(row['population']))
    require(set(achieved)==set(group_cells) and all(abs(fsum(achieved[k])-v)<=max(1e-8,v*1e-12) for k,v in group_cells.items()),
        'Diagnostic subgroup age/sex counts do not reconcile with selected demographic inputs.')
    risk=config['risk'];raw=parse_csv(risk['csv'],'risk',subgroup_ids=group_ids)['rows'] if risk['csv'] else risk['records']
    require(risk['no_factors'] or raw,'Risk prevalences are missing; an empty uploaded risk vector is not an RR1 selection.')
    overrides=parse_csv(config['association']['csv'],'association',subgroup_ids=group_ids)['rows'] if config['association']['csv'] else []
    for row in args['strata']:
        if risk['no_factors']:
            require(not raw and not row['factors'],'Uploaded factors conflict with the explicitly selected RR1 reference.')
        else:
            eligible=[r for r in raw if r.get('subgroup_id') in {row['subgroup_id'],'all'} and r['sex'] in {row['sex'],'both'} and
                r['age_lower']<=float(row['age']) and (r['age_upper_exclusive'] is None or float(row['age'])<r['age_upper_exclusive'])]
            require(eligible,'Uploaded stratum has no selected risk-prevalence coverage.')
            definitions={r['factor_id']:r['exposure_definition'] for r in catalog()['rows']}
            require(all(r['factor_id'] in definitions and r['supplied_exposure_definition']==definitions[r['factor_id']] for r in eligible),
                'Risk exposure-definition mismatch remains unavailable; a mapping label cannot silently repair it.')
            require(all(r['state'] in {'value','excluded','not_applicable'} for r in eligible),'Missing/unreviewed risk prevalence cannot authorize inversion.')
            expected_factors={r['factor_id']:r['prevalence'] for r in eligible if r['state']=='value'}
            require(len({r['factor_id'] for r in eligible})==len(eligible),'Overlapping selected risk inputs.')
            require(set(expected_factors)=={r['factor_id'] for r in row['factors']} and
                all(abs(float(f['prevalence'])-float(expected_factors[f['factor_id']]))<=1e-12 for f in row['factors']),
                'Uploaded mapping risk marginals differ from the selected review inputs.')
        controls=[{k:v for k,v in {**r,'source':{'source_id':r['source']},'review':REVIEW,'provenance':PROVENANCE}.items()
            if k not in {'subgroup_id','review_status'}} for r in overrides if r['subgroup_id'] in {row['subgroup_id'],'all'}]
        policy=association_policy(global_control=config['association']['global_control'],overrides=controls,**META)
        cell=next(dict(sex=k[0],age_lower=k[1],age_upper_exclusive=k[2]) for k in expected if k[0]==row['sex'] and k[1]<=float(row['age']) and (k[2] is None or float(row['age'])<k[2]))
        requested,_=control_for_cell(policy,cell)
        require(abs(float(row['association_control'])-float(requested))<=1e-12,'Uploaded association differs from the selected within-stratum control.')
    require(set(mapping_counts)==set(expected) and all(abs(fsum(mapping_counts[k])-v)<=max(1e-8,v*1e-12) for k,v in expected.items()),
        'Diagnostic age/sex counts do not reconcile with the selected population scope.')
    mapping=IncidenceMapping(**args);obs=bundle['observation']
    gate=inversion_gate(mapping,obs,unresolved_threshold_approval=bundle['unresolved_threshold_approval'],
        scientific_review=bundle['scientific_review'],concentration_policy=bundle['concentration_policy'],**DMETA)
    require(gate['prefit_allowed'],'Native3C gate blocks inference: '+', '.join(gate['failed_criteria']))
    return mapping,obs,gate,bundle


def audit_diagnostic(config,context):
    mapping,obs,gate,bundle=prepare_diagnostic(config,context)
    grid=config['diagnostic']['grid_intervals'];require(type(grid) is int and 64<=grid<=512,'Bounded grid64–512required.')
    target=obs['calibration_target']
    report=assess_identifiability(mapping.incidence,lower=config['diagnostic']['lower'],upper=config['diagnostic']['upper'],
        target=target,grid_intervals=grid,scenario_hash=mapping.specification['profile_hash'],**DMETA)
    valid=report['identifiable_over_search_interval'] and report['numerical_convergence'] and report['target_attainable']
    difference=None
    tolerance=bundle['numerical_protocol'].get('absolute_refinement_tolerance')
    require(tolerance is not None and 0<float(tolerance)<=max(1e-6,abs(float(target))*1e-5),
        'A positive reviewed refinement tolerance no looser than native3C is required.')
    if valid:
        root=report['all_roots'][0]['alpha'];difference=abs(mapping.incidence(root)-mapping.refined().incidence(root))
        valid=difference<=float(tolerance)
    return seal(dict(schema_id='country_tb.review_ui_numerical_gate.v1',configuration_hash=material_hash(config),
        passed=bool(valid),identifiability=report,integration_refinement_difference=difference,
        reviewed_allowed_difference=tolerance,scientific_prefit_gate=gate,production_available=False,
        fitting_accepted=False,classification='gate_audit_only_no_TBI_profile_adopted'))


def execute_diagnostic(config,context,*,numerical_clearance=None):
    require(numerical_clearance is not None and numerical_clearance.get('configuration_hash')==material_hash(config) and numerical_clearance.get('passed') is True,
        'Root uniqueness and numerical refinement gates are unresolved. Check the diagnostic gates first.')
    require(numerical_clearance.get('profile_hash')==seal({k:v for k,v in numerical_clearance.items() if k!='profile_hash'})['profile_hash'],
        'Numerical gate receipt hash mismatch.')
    # Receipts identify inputs; they cannot authorize an unchecked or hand-resealed map.
    verified=audit_diagnostic(config,context);require(verified['passed'],'Current numerical gates failed; no inversion may run.')
    mapping,obs,gate,bundle=prepare_diagnostic(config,context)
    grid=config['diagnostic']['grid_intervals'];require(type(grid) is int and 64<=grid<=512,'Bounded grid64–512required.')
    result=diagnostic_calibration(mapping,obs,gate,lower=config['diagnostic']['lower'],upper=config['diagnostic']['upper'],grid_intervals=grid,**DMETA)
    return seal(dict(schema_id='country_tb.review_ui_diagnostic_run.v1',configuration_hash=material_hash(config),configuration=deepcopy(config),
        context=context,diagnostic=result,review_bundle=bundle,gate=gate,production_available=False,externally_validated=False,
        numerical_gate=verified,
        public_status=PUBLIC_STATUS,status=result['status'],classification='diagnostic_only_not_approved_for_production',
        declared_population=context['demography']['population'],diagnostic_population=mapping.population,
        excluded_population=float(context['demography']['population'])-mapping.population,
        warning='A successful mathematical diagnostic is not national scientific approval or an empirical TBI estimate.'))
