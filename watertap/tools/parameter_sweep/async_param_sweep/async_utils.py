import numpy as np

from pyomo.common.collections import ComponentSet, ComponentMap
import pyomo.environ as pyo
import copy


def _create_component_output_skeleton(component, value):
    comp_dict = {}
    comp_dict["value"] = value
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

    # Add information to this output that WILL NOT be written as part
    # of the file saving step.
    # comp_dict["_pyo_obj"] = component

    return comp_dict


def _create_local_output_skeleton(model, sweep_params, outputs, run_successful):
    output_dict = {}
    output_dict["sweep_params"] = {}
    output_dict["outputs"] = {}
    output_dict["solve_successful"] = run_successful
    sweep_param_objs = ComponentSet()

    # Store the inputs
    for key, item in sweep_params.items():
        var = model.find_component(key)
        sweep_param_objs.add(var)
        output_dict["sweep_params"][var.name] = _create_component_output_skeleton(
            var, item
        )

    if outputs is None:
        # No outputs are specified, so every Var, Expression, and Objective on the model should be saved
        for pyo_obj in model.component_data_objects(
            (pyo.Var, pyo.Param, pyo.Expression, pyo.Objective), active=True
        ):
            # Only need to save this variable if it isn't one of the value in sweep_params
            try:
                if run_successful:
                    val = pyo.value(pyo_obj)
                else:
                    val = np.nan
                output_dict["outputs"][
                    pyo_obj.name
                ] = _create_component_output_skeleton(pyo_obj, val)
            except ValueError:
                print("Valueerror in getitng vals")

    else:
        # Save only the outputs specified in the outputs dictionary
        for short_name, object_name in outputs.items():
            pyo_obj = model.find_component(object_name)
            if run_successful:
                val = pyo.value(pyo_obj)
            else:
                val = np.nan
            output_dict["outputs"][short_name] = _create_component_output_skeleton(
                pyo_obj, val
            )

    return output_dict


def _merge_copied_dicts(dict_array):
    num_samples = len(dict_array)
    local_copy = copy.deepcopy(dict_array[0])
    local_copy.pop("order_index", None)
    global_output_dict = copy.deepcopy(local_copy)
    global_output_dict["solve_successful"] = np.zeros(num_samples)
    for cn, case_dict in enumerate(dict_array):
        case_number = case_dict["order_index"]
        case_dict.pop("order_index", None)

        for key in global_output_dict["sweep_params"].keys():
            try:
                if (
                    isinstance(
                        global_output_dict["sweep_params"][key]["value"], np.ndarray
                    )
                    == False
                ):
                    global_output_dict["sweep_params"][key]["value"] = np.zeros(
                        num_samples
                    )

                global_output_dict["sweep_params"][key]["value"][
                    case_number
                ] = case_dict["sweep_params"][key]["value"]
            except KeyError:
                print("keyerror global dict", key)
        for key in global_output_dict["outputs"].keys():
            try:
                if (
                    isinstance(global_output_dict["outputs"][key]["value"], np.ndarray)
                    == False
                ):
                    global_output_dict["outputs"][key]["value"] = np.zeros(num_samples)

                global_output_dict["outputs"][key]["value"][case_number] = case_dict[
                    "outputs"
                ][key]["value"]
            except KeyError:
                print("keyerror global dict", key)
        global_output_dict["solve_successful"][case_number] = case_dict[
            "solve_successful"
        ]
    return global_output_dict


def _convert_sweep_params_to_dict(params, values):
    param_dict = []
    for i in range(len(values)):
        param_dict.append({})
        for k, item in enumerate(params.values()):
            # print(k, item)
            param_dict[-1][item.pyomo_object.name] = values[i, :][k]
        param_dict[-1]["order_index"] = i
    # print("param_dict", param_dict)
    return param_dict


def _convert_outputs_to_dict(outputs):
    out_dict = {}
    for key, obj in outputs.items():
        out_dict[key] = obj.name
    return out_dict


def _update_model_values_from_dict(m, values):
    set_vals = {}
    for k, value in values.items():
        # print(k, value)
        # param = item.pyomo_object
        param = m.find_component(str(k))
        # print(type(param))
        # set_vals[str(item.pyomo_object)] = values[k]
        if param.is_variable_type():
            # Fix the single value to values[k]
            param.fix(value)

        elif param.is_parameter_type():
            # Fix the single value to values[k]
            param.set_value(value)
