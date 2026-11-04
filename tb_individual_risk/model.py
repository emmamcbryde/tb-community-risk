"""Public individual prediction orchestration; no UI or network imports."""
from copy import deepcopy
from dataclasses import replace
import math
from .contracts import PersonRecord,ValidationError,FIELDS,BOOLEAN_FIELDS,Observation
from .parameters import ModelParameters,POLICIES,EXTERNAL
from .infection_probability import infection_probability,required_fields as infection_fields,Unavailable
from .progression import conditional_progression,horizon_for
from .risk_factors import risk_factors,factor_definitions,combine_rr
from .missing_data import scenarios,validated_value
from .provenance import seal,canonical,VERSION,source_safe


class IndividualTBRiskModel:
    def __init__(self,parameters=None,*,formulation="hazard_based_progression_v1"):
        from .packages import ParameterPackage
        if formulation not in {"hazard_based_progression_v1","menzies_2008_compatibility_v1"}:
            raise ValidationError("formulation","Unknown formulation; no fallback")
        self.formulation=formulation
        self._package=None
        if isinstance(parameters,dict) and parameters.get("schema_id")=="individual_tb.parameter_package.v1":
            parameters=ParameterPackage(canonical(parameters))
        if isinstance(parameters,ParameterPackage):
            package=parameters.record()
            if formulation not in package["allowed_formulations"]:
                raise ValidationError("formulation","Parameter package and formulation incompatible; no silent historical/modern mixing")
            self._package=parameters
            if formulation=="menzies_2008_compatibility_v1":
                self._parameters=None;self._parameter_record=package;return
            parameters=package["model_parameters"]
        elif formulation=="menzies_2008_compatibility_v1":
            raise ValidationError("parameters","Compatibility requires the immutable menzies_2008_published_reference_v1 package")
        parameters=ModelParameters() if parameters is None else parameters
        if isinstance(parameters,dict):parameters=ModelParameters.from_dict(parameters)
        if not isinstance(parameters,ModelParameters):raise TypeError("Supply versioned ModelParameters.")
        self._parameters=deepcopy(parameters)
        # Roundtrip validates even mutable nested dataclass payloads before use.
        self._parameters=ModelParameters.from_dict(self._parameters.record())
        prior=self._parameters.population_prior
        for name,rows in prior.distributions.items():
            if name not in FIELDS and not name.startswith("factor:"):raise ValidationError(name,"Unknown prior field.")
            for row in rows:validated_value(name,row["value"])
        for name,values in prior.bounds.items():
            if name not in FIELDS and not name.startswith("factor:"):raise ValidationError(name,"Unknown bound field.")
            for v in values:validated_value(name,v)
        for row in (*prior.joint_states,*prior.joint_bound_states):
            for name,v in row["values"].items():
                if name not in FIELDS and not name.startswith("factor:"):raise ValidationError(name,"Unknown joint field.")
                validated_value(name,v)
        for _,spec in factor_definitions(self._parameters).values():
            if spec.get("field") not in FIELDS and not str(spec.get("field","")).startswith("factor:"):
                raise ValidationError("effect field","Unsupported factor field.")
        self._parameter_record=self._parameters.record()

    @property
    def parameters(self):return self._package if self.formulation=="menzies_2008_compatibility_v1" else deepcopy(self._parameters)

    def _required(self,person):
        p=self._parameters
        required=infection_fields(p,person.values())|{"current_age","time_since_infection"}
        required|={s["field"] for _,s in factor_definitions(p).values()}
        required|={name for _,spec in factor_definitions(p).values() if person.values().get(spec["field"]) in (None,spec.get("present_value",True)) for name in spec.get("required_context",{})}
        required|=set(p.component("baseline_progression_hazard").config.get("required_history",[]))
        return required

    def _evaluate(self,values):
        p=self._parameters
        output={"values":values,"infection":None,"progression":None,"factors":None,"absolute_risk":None,"errors":[],"warnings":[]}
        try:output["infection"]=infection_probability(values,p)
        except (Unavailable,ValidationError) as exc:output["errors"].append(str(exc))
        try:
            factors=risk_factors(values,p);output["factors"]=factors
            output["progression"]=conditional_progression(values,p,factors["effective_rr"])
            output["warnings"]+=factors["warnings"]+output["progression"]["warnings"]
        except (Unavailable,ValidationError) as exc:output["errors"].append(str(exc))
        if output["infection"] is not None:
            probability=output["infection"]["probability"]
            if probability==0:
                output["absolute_risk"]=0.0
                if output["progression"] is None:output["warnings"].append("Known zero infection gives zero baseline-infection progression risk; conditional progression remains unavailable. No future infection model.")
            elif output["progression"] is not None:output["absolute_risk"]=probability*output["progression"]["probability"]
        return output

    def _predict_hazard(self,person):
        try:
            record=person if isinstance(person,PersonRecord) else PersonRecord.from_dict(person)
            # Revalidate caller-supplied records, including original explicit states.
            record=PersonRecord.from_dict(record.original)
        except (ValidationError,TypeError,ValueError) as exc:
            error=exc.as_dict() if isinstance(exc,ValidationError) else {"field":"record","message":str(exc),"code":"invalid_record"}
            return seal(dict(schema_id="individual_tb.prediction.v1",record_id=person.get("record_id") if isinstance(person,dict) and isinstance(person.get("record_id"),(str,int)) and not isinstance(person.get("record_id"),bool) else None,
                status="invalid",validation_errors=[error],central_absolute_risk=None,lower_missing_scenario=None,upper_missing_scenario=None,
                original_record=source_safe(person) if isinstance(person,dict) else repr(person),combination_policy=self._parameters.combination_policy,
                clinical_use_approved=False,classification=self._parameters.classification,parameter_hash=self._parameter_record["scientific_hash"],parameter_record=deepcopy(self._parameter_record)))
        required=self._required(record);prior=self._parameters.population_prior
        outputs={};unresolved={};warnings=[]
        for label,bounds in (("central",False),("bounds",True)):
            try:
                supports,unresolved[label]=scenarios(record,required,prior,bounds=bounds,budget=self._parameters.max_scenarios)
                outputs[label]=[(w,assign,self._evaluate(values)) for w,assign,values in supports]
            except (Unavailable,ValidationError) as exc:
                outputs[label]=[];unresolved[label]=set(required);warnings.append(str(exc))
        central_rows=outputs["central"];bound_rows=outputs["bounds"]
        def expectation(key):
            vals=[(w,row.get(key)) for w,_,row in central_rows]
            return math.fsum(w*v for w,v in vals)/math.fsum(w for w,_ in vals) if vals and all(v is not None for _,v in vals) else None
        central=expectation("absolute_risk")
        infection_values=[(w,row["infection"]["probability"] if row["infection"] else None) for w,_,row in central_rows]
        p_tbi=math.fsum(w*v for w,v in infection_values)/math.fsum(w for w,_ in infection_values) if infection_values and all(v is not None for _,v in infection_values) else None
        conditional=central/p_tbi if central is not None and p_tbi is not None and p_tbi>0 else None
        available_bounds=bool(bound_rows) and all(row["absolute_risk"] is not None for _,_,row in bound_rows)
        lower=min(row["absolute_risk"] for _,_,row in bound_rows) if available_bounds else None
        upper=max(row["absolute_risk"] for _,_,row in bound_rows) if available_bounds else None
        conditional_bound_values=[row["progression"]["probability"] for _,_,row in bound_rows if row["progression"] is not None and row["infection"] is not None and row["infection"]["probability"]>0] if available_bounds else []
        lower_conditional=min(conditional_bound_values) if conditional_bound_values else None
        upper_conditional=max(conditional_bound_values) if conditional_bound_values else None
        rr_values=[(w,row["factors"]["effective_rr"] if row["factors"] else None) for w,_,row in central_rows]
        expected_rr=math.fsum(w*v for w,v in rr_values) if rr_values and all(v is not None for _,v in rr_values) else None
        effective_rr=expected_rr if expected_rr is not None and len({v for _,v in rr_values})==1 else None
        contribution={}
        for name in sorted(required):
            obs=record.fields.get(name,Observation())
            groups={}
            for weight,assignment,row in central_rows:
                if name not in assignment:continue
                value=assignment[name];key=canonical(value)
                g=groups.setdefault(key,{"value":value,"weight":0.0,"weighted_risk":0.0,"available":True})
                g["weight"]+=weight
                if row["absolute_risk"] is None:g["available"]=False
                else:g["weighted_risk"]+=weight*row["absolute_risk"]
            for g in groups.values():
                g["conditional_mean_absolute_risk"]=g.pop("weighted_risk")/g["weight"] if g.pop("available") else None
            contribution[name]={"input_state":obs.state,"supplied_value":obs.value,"source":obs.source,
                "resolution":"known" if obs.known else ("marginalised_explicit_prior" if groups else "unresolved"),
                "conditional_state_contributions":list(groups.values())}
        drivers={}
        for weight,_,row in central_rows:
            if not row["factors"] or row["absolute_risk"] is None or not row["progression"]:continue
            for factor_id in row["factors"]["factors_present"]:
                other={k:v for k,v in row["factors"]["individual_rr"].items() if k!=factor_id}
                rr=combine_rr(other,self._parameters.combination_policy)["effective_rr"]
                counterfactual=conditional_progression(row["values"],self._parameters,rr)["probability"]*row["infection"]["probability"]
                d=drivers.setdefault(factor_id,{"factor_id":factor_id,"probability_present":0.0,"absolute_risk_difference_if_absent":0.0,"highest_multiplier_contributor_probability":0.0})
                d["probability_present"]+=weight
                d["absolute_risk_difference_if_absent"]+=weight*(row["absolute_risk"]-counterfactual)
                if factor_id in row["factors"]["highest_contributors"]:d["highest_multiplier_contributor_probability"]+=weight
        errors=sorted({e for rows in outputs.values() for _,_,row in rows for e in row["errors"]})
        warnings+=errors+[w for rows in outputs.values() for _,_,row in rows for w in row["warnings"]]
        if self._parameters.classification=="synthetic_not_clinical":warnings.append("INVENTED SYNTHETIC PARAMETERS: NOT FOR CLINICAL USE.")
        if self._parameters.classification!="reviewed":warnings.append("No reviewed clinical parameter set is supplied.")
        if central is not None and lower is not None and not lower-1e-14<=central<=upper+1e-14:
            warnings.append("Central estimate falls outside the prespecified missing-information scenarios; these are not confidence bounds.")
        known=sum(obs.known for obs in record.fields.values());critical_known=sum(record.fields.get(n,Observation()).known for n in required)
        used_fields=required|{"horizon_years"}
        infection_cfg=self._parameters.component("infection_probability").config
        if "tst" in infection_cfg.get("assays",[]):used_fields.add("tst_mm")
        unused=sorted(n for n,o in record.fields.items() if n not in used_fields and o.known)
        ignored_evidence=[n for n in ("close_contact","recent_exposure","recent_conversion","bcg","hiv_treatment","tst_result","igra_result","tst_mm") if n in unused and record.fields[n].value not in (False,"not_done")]
        if ignored_evidence:warnings.append("Known evidence not modelled individually by selected parameter scope: "+", ".join(ignored_evidence)+". Supply an appropriate evidence-conditioned infection probability/history model; no automatic multiplier.")
        modeled={spec["field"] for _,spec in factor_definitions(self._parameters).values()}
        exclusions=self._parameters.component("medical_risk_effects").config.get("explicit_scope_exclusions",{})
        active_unmodelled=[n for n,o in record.fields.items() if (n in BOOLEAN_FIELDS or n.startswith("factor:")) and o.known and o.value is True and n not in modeled and n not in {"bcg","close_contact","recent_exposure","recent_conversion","hiv_treatment","previous_active_tb","previous_preventive_treatment","active_tb_at_baseline"}]
        unexplained=[n for n in active_unmodelled if n not in exclusions]
        if unexplained:
            warnings.append("Known progression findings lack an effect or explicit scope exclusion: "+", ".join(unexplained))
            central=lower=upper=conditional=None
            lower_conditional=upper_conditional=None
        if any(record.fields[n].known and record.fields[n].value is True for n in ("previous_active_tb","previous_preventive_treatment","active_tb_at_baseline")):
            # This milestone does not supply recurrence/treatment/prevalent disease policies.
            warnings.append("Prior TB, preventive treatment or active TB requires an external history-specific model.")
            central=lower=upper=conditional=None
            lower_conditional=upper_conditional=None
        try:horizon=horizon_for(record.values(),self._parameters)
        except (Unavailable,ValidationError):horizon=None
        status=("calculated" if critical_known==len(required) else "calculated_with_missing_information") if central is not None else ("bounded_only" if lower is not None else "unavailable")
        return seal(dict(schema_id="individual_tb.prediction.v1",record_id=record.record_id,identifier_state=record.identifier_state,status=status,validation_errors=[],
            central_absolute_risk=central,lower_missing_scenario=lower,upper_missing_scenario=upper,
            lower_conditional_missing_scenario=lower_conditional,upper_conditional_missing_scenario=upper_conditional,
            tbi_probability=p_tbi,conditional_progression_risk=conditional,effective_rr=effective_rr,expected_effective_rr=expected_rr,
            combination_policy=self._parameters.combination_policy,principal_risk_drivers=sorted(drivers.values(),key=lambda d:(-d["absolute_risk_difference_if_absent"],d["factor_id"])),
            missing_critical_variables=sorted(n for n in required if not record.fields.get(n,Observation()).known),
            unresolved_fields=sorted(unresolved["central"]),
            input_completeness={"known_fields":known,"total_fields":len(record.fields),"fraction_known":known/len(record.fields),"critical_known":critical_known,"critical_total":len(required)},
            horizon_years=horizon,input_contributions=contribution,known_inputs_not_used_by_selected_model=unused,
            explicit_scope_exclusions=exclusions,assumptions=["Baseline infection only; no future infection or reinfection generated.","First disease event only; no prior-disease depletion or recurrence model.","Multipliers applied to cumulative hazard, never directly to probability.","Missing-information extrema are over prespecified states, not sampling confidence intervals.","Prior dependence policy: "+prior.dependence_policy],
            warnings=sorted(set(warnings)),warning_status="warning" if warnings else "none",original_record=source_safe(record.original),
            scenarios=[{"weight":w,"marginalised_inputs":a,"infection":r["infection"],"progression":r["progression"],"factors":r["factors"],"absolute_risk":r["absolute_risk"],"errors":r["errors"]} for w,a,r in central_rows],
            bound_scenarios=[{"inputs":a,"absolute_risk":r["absolute_risk"],"errors":r["errors"]} for _,a,r in bound_rows],
            prior_hash=seal(self._parameters.population_prior)["scientific_hash"],parameter_hash=self._parameter_record["scientific_hash"],parameter_record=deepcopy(self._parameter_record),
            implementation_version=VERSION,classification=self._parameters.classification,clinical_use_approved=False))

    def predict(self,person):
        if self.formulation=="menzies_2008_compatibility_v1":
            from .menzies_2008 import predict_reference
            result=predict_reference(person,self._package)
        else:result=self._predict_hazard(person)
        package=self._package.record() if self._package is not None else None
        result.setdefault("unresolved_assumptions",list(result.get("unresolved_fields",[])))
        result["implementation_version"]=VERSION
        result.update(formulation=self.formulation,
            parameter_package=package["package_id"] if package else "custom_parameters_"+self._parameter_record["scientific_hash"][:16],
            parameter_package_hash=package["scientific_hash"] if package else self._parameter_record["scientific_hash"],
            source_version=package["source_version"] if package else VERSION)
        if package and self.formulation=="hazard_based_progression_v1":
            result["parameter_package_overrides"]={} if self._parameter_record==package["model_parameters"] else {"combination_policy":self._parameters.combination_policy}
        return seal(result)

    def rank(self,line_list,**kwargs):
        from .ranking import rank
        return rank(self,line_list,**kwargs)

    def compare_policies(self,line_list,policies=None,**kwargs):
        from .line_list import records
        if self.formulation=="menzies_2008_compatibility_v1":raise ValidationError("combination_policy","Historical general combination rule unverified; use separate calculation-only policy summaries")
        rows=records(line_list)
        chosen=list(policies) if policies is not None else ["highest_applicable_rr","additive_excess_rr","multiplicative_rr"]
        if not chosen or len(chosen)!=len(set(chosen)) or any(p not in POLICIES for p in chosen):
            raise ValidationError("policies","Select distinct supported policies.")
        result={}
        for policy in chosen:
            model=IndividualTBRiskModel(replace(self._parameters,combination_policy=policy),formulation=self.formulation)
            model._package=self._package
            result[policy]=model.rank(rows,**kwargs)
        return result
