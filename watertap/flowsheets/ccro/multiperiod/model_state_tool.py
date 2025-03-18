from pyomo.common.collections import ComponentMap
from pyomo.environ import Var, Constraint, Param, Expression, value, Block, Objective
import idaes.core.util.scaling as iscale
import numpy as np
import json

__author__ = "Alexander V. Dudchenko"

""" This tool is designed to help store model states and transfer them 
other models, as well as save model results to jason file in similar
manner as the ParameterSweep tool"""


class ModelState:
    def get_model_state(self, model):
        """Get model state of passed in model and store it in a dict"""
        self.variableStatesDict = {}
        for v in model.component_data_objects(Var):
            self.variableStatesDict[str(v)] = {
                "value": v.value,
                "lb": v.lb,
                "ub": v.ub,
                "fixed": v.fixed,
            }
        self.paramStatesDict = {}
        for v in model.component_data_objects(Param):
            self.paramStatesDict[str(v)] = {"value": v.value}
        self.constraintStatesDict = {}
        for v in model.component_data_objects(Constraint):
            self.constraintStatesDict[str(v)] = v.active
        self.variableScalingStateDict = {}
        self.constraintScalingStateDict = {}
        for v in model.component_data_objects(Var):
            scale = iscale.get_scaling_factor(v)
            self.variableScalingStateDict[str(v)] = scale
        for v in model.component_data_objects(Constraint):
            scaling_value = iscale.get_constraint_transform_applied_scaling_factor(v)
            self.constraintScalingStateDict[str(v)] = scaling_value
        self.old_model = model.name

    def set_model_state(self, model):
        """Set model state of a stored model into new model"""
        for v, val in self.variableStatesDict.items():

            v = self.update_model_reference(v, self.old_model)
            try:
                model.find_component(v).setlb(val["lb"])
                model.find_component(v).setub(val["ub"])
                model.find_component(v).set_value(val["value"])
                if val["fixed"]:
                    model.find_component(v).fix()
                else:
                    model.find_component(v).unfix()
            except:
                print("failed setitng", v)
        for v, active in self.constraintStatesDict.items():
            v = self.update_model_reference(v, self.old_model)
            try:
                if active:
                    model.find_component(v).activate()
                else:
                    model.find_component(v).deactivate()
            except:
                print("Model state tool failed to copy over", v)

        for v, scale in self.variableScalingStateDict.items():
            v = self.update_model_reference(v, self.old_model)
            if scale is None:
                iscale.unset_scaling_factor(model.find_component(v))
            else:
                iscale.set_scaling_factor(model.find_component(v), scale)
        for v, scale in self.variableScalingStateDict.items():
            v = self.update_model_reference(v, self.old_model)
            if scale is None:
                iscale.unset_scaling_factor(model.find_component(v))
            else:
                iscale.set_scaling_factor(model.find_component(v), scale)
        for v, scale in self.constraintScalingStateDict.items():
            v = self.update_model_reference(v, self.old_model)
            if model.find_component(v) is not None:
                if scale is None:
                    iscale.constraint_scaling_transform_undo(model.find_component(v))
                else:
                    iscale.constraint_scaling_transform(model.find_component(v), scale)
        for v in model.component_data_objects(Expression):
            v = self.update_model_reference(v, self.old_model)
            try:
                value(v)
            except:
                pass

    def update_model_reference(self, var, old_model):
        """this will update variable name with new model name"""
        if isinstance(var, str) == False:
            var = var.name
        return var.replace(old_model + ".", "")

    def add_expected_outputs(
        self, model, scenario="unspecified_scenario", number_of_cases=1
    ):
        """This will create a dict of output keys based on model structure so
        we can populate it with desired outputs

        Keywords:
            model -- concrete model object to store results for
            scenario -- scenarios to store, should be string or list, this is useful
            if different options or scenarios are being ran, for example when one changes
            costs for pumps and wants to sweep across different water recoveries, the scenario is pump cost type and
            cases are water recovery.
            number_of_cases -- this specifies how many simulations are ran for a given scenario
        """
        if hasattr(self, "output_dict") == False:
            self.output_dict = {}
        if scenario is not None:
            if isinstance(scenario, str):
                scenario = [scenario]
            for s in scenario:
                self.output_dict[s] = {"outputs": {}}

        self.default_scenario = scenario[0]
        # No outputs are specified, so every Var, Expression, and Objective on the model should be saved
        for key, items in self.output_dict.items():
            for pyo_obj in model.component_data_objects(
                (Var, Expression, Objective, Param), active=True
            ):
                items["outputs"][pyo_obj.name] = self._create_component_output_skeleton(
                    pyo_obj, number_of_cases
                )

    def _create_component_output_skeleton(self, component, number_of_cases=1):
        comp_dict = {}
        comp_dict["value"] = np.zeros(number_of_cases, dtype=float).tolist()

        if hasattr(component, "lb"):
            comp_dict["lower bound"] = component.lb
        if hasattr(component, "ub"):
            comp_dict["upper bound"] = component.ub
        if hasattr(component, "get_units"):
            unit_obj = component.get_units()
            if unit_obj is not None:
                comp_dict["units"] = component.get_units().name
            else:
                comp_dict["units"] = "None"
        comp_dict["full_name"] = component.name
        return comp_dict

    def update_stored_outputs(self, model, scenario=None, case_number=0):
        """This will update the data dictionary with model results for specified scenario and case number

        Keywords:
            model -- concrete model object to store results for
            scenario -- scenarios to store, should be string or list, this is useful
            if different options or scenarios are being ran, for example when one changes
            costs for pumps and wants to sweep across different water recoveries, the scenario is pump cost type and
            cases are water recovery.
            case_number -- this specifies simulations number
        """
        if scenario is None:
            scenario = self.default_scenario
        print(f"Saving to {scenario}, case number {case_number}")
        for var_name, specs in self.output_dict[scenario]["outputs"].items():
            pyo_obj = model.find_component(specs["full_name"])
            # incase value is not initialized or can't be evaluated
            # typical case, is a var is created, but not initialized or touched, such is 0 index vars in 1D RO
            try:
                self.output_dict[scenario]["outputs"][var_name]["value"][
                    case_number
                ] = value(pyo_obj)
            except ValueError:
                pass

    def save_model_outputs_to_json(self, json_result_file_name):
        """saves data in a json
        Keyword:
        json_result_file_name -- file name and location should not include .json"""
        with open(f"{json_result_file_name}.json", "w", encoding="utf-8") as f:
            json.dump(self.output_dict, f, ensure_ascii=False, indent=4)
