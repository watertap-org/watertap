import multiprocessing

import pyomo.environ as pyo
import copy
import numpy as np
import os

import platform
import time
from pyomo.common.collections import ComponentSet, ComponentMap


from idaes.core.util import to_json, from_json


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
        # 3if "order_index" in case_dict:
        case_number = case_dict["order_index"]
        case_dict.pop("order_index", None)
        # else:
        #    case_number = cn

        # print("case_number", case_dict)
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


def do_mp_sweep(self, model, sweep_params, outputs, local_values):
    optimize_function = self.config.optimize_function
    optimize_kwargs = self.config.optimize_kwargs
    reinitialize_before_sweep = self.config.reinitialize_before_sweep
    reinitialize_function = self.config.reinitialize_function
    reinitialize_kwargs = self.config.reinitialize_kwargs
    build_function = self.config.custom_do_param_sweep_kwargs.get("build_function")
    build_kwargs = self.config.custom_do_param_sweep_kwargs.get("build_kwargs")
    probe_function = self.config.probe_function
    html_notice = None  # self.config.html_noticem will be implememnted for UI later

    outputs = _convert_outputs_to_dict(outputs)

    run_dict = _convert_sweep_params_to_dict(sweep_params, local_values)
    num_workers = self.config.custom_do_param_sweep_kwargs.get("num_workers")
    if num_workers is None:
        num_workers = multiprocessing.cpu_count()
    if self.config.custom_do_param_sweep_kwargs.get("load_form_json"):
        json_string = to_json(model, return_json_string=True)
    else:
        json_string = None
    run_queue = multiprocessing.Queue()
    return_queue = multiprocessing.Queue()
    for param in run_dict:
        # print(param)
        run_queue.put(param)
    actors = []
    for cpu in range(num_workers):
        actors.append(
            multiprocessing.Process(
                target=setupMCActors,
                args=(
                    run_queue,
                    return_queue,
                    build_function,
                    build_kwargs,
                    reinitialize_before_sweep,
                    reinitialize_function,
                    reinitialize_kwargs,
                    optimize_function,
                    optimize_kwargs,
                    probe_function,
                    outputs,
                    html_notice,
                    json_string,
                    cpu,
                ),
            )
        )
        actors[-1].start()
    for cpu in range(num_workers):
        actors[cpu].join()
    results = []
    while True:
        if return_queue.empty():
            break
        else:
            result = return_queue.get()
            results.append(result)

    output = _merge_copied_dicts(list(results))

    return output


def setupMCActors(
    queue,
    return_queue,
    build_function,
    build_kwargs,
    reinitialize_before_sweep,
    reinitialize_function,
    reinitialize_kwargs,
    optimize_function,
    optimize_kwargs,
    probe_function,
    outputs,
    html_notice,
    json_string,
    process_id,
):
    param_actor = paramActor(
        build_function,
        build_kwargs,
        reinitialize_before_sweep,
        reinitialize_function,
        reinitialize_kwargs,
        optimize_function,
        optimize_kwargs,
        probe_function,
        outputs,
        json_string,
    )
    while True:
        if queue.empty():
            break
        else:
            sweep_params = queue.get()

            result = param_actor.solve_problem(sweep_params)
            return_queue.put(result)


class paramActor:
    def __init__(
        self,
        build_function,
        build_kwargs,
        reinitialize_before_sweep,
        reinitialize_function,
        reinitialize_kwargs,
        optimize_function,
        optimize_kwargs,
        probe_function,
        outputs,
        json_string,
    ):
        self.build_function = build_function
        self.build_kwargs = build_kwargs
        self.reinitialize_before_sweep = reinitialize_before_sweep
        self.reinitialize_function = reinitialize_function
        self.reinitialize_kwargs = reinitialize_kwargs
        self.optimize_function = optimize_function
        self.optimize_kwargs = optimize_kwargs
        self.outputs = outputs
        self.probe_function = probe_function
        self.model_init = False
        self.model = None
        self.json_string = json_string
        self.build_model()

    def build_model(self):
        del self.model

        self.model = self.build_function(**self.build_kwargs)
        if self.json_string is not None:
            from_json(self.model, sd=None, fname=None, s=self.json_string)
            self.model_init = True

    def init_model(self):
        if self.json_string is not None:
            from_json(self.model, sd=None, fname=None, s=self.json_string)
        self.reinitialize_function(self.model, **self.reinitialize_kwargs)
        # self.model_state = modelStateStorage(self.model)
        # self.clean_state = to_json(self.model, return_json_string=True)
        self.model_init = True

    def get_probe_result(self, m, sweep_params):
        if self.probe_function is None:
            return True
        else:
            try:
                probe_result = self.probe_function(m, sweep_params)
            except:
                print("----------------------------------")
                print("error in probe function!!!!!!!!!!!")
                probe_result = True
            return probe_result

    def solve_problem(self, sweep_params):
        order_index = sweep_params["order_index"]
        # print("order_index", order_index)
        sweep_params.pop("order_index", None)

        run_successful = False

        _update_model_values_from_dict(self.model, sweep_params)

        if self.get_probe_result(self.model, self.reinitialize_kwargs):
            for i in ["Try #0", "Try #1"]:
                if (
                    self.reinitialize_function is not None
                    and self.model_init == False
                    or self.reinitialize_before_sweep
                ):
                    try:
                        self.build_model()
                        self.init_model()
                        self.model_init = True
                    except:
                        self.model_init = False
                        run_successful = False
                        init_okay = False

                _update_model_values_from_dict(self.model, sweep_params)

                if self.model_init and self.get_probe_result(
                    self.model, self.reinitialize_kwargs
                ):
                    try:
                        results = self.optimize_function(
                            self.model, **self.optimize_kwargs
                        )
                        pyo.assert_optimal_termination(results)
                        run_successful = True
                        break
                    except:
                        run_successful = False
                        if (
                            self.reinitialize_before_sweep
                            or self.reinitialize_function is None
                        ):
                            break
                        else:
                            self.model_init = False
                else:
                    print("Model failed to init, skipping condition!")
                    break
        output_dict = _create_local_output_skeleton(
            self.model, sweep_params, self.outputs, run_successful
        )
        output_dict["order_index"] = order_index
        print("RUN Complete, solution optimal: {}".format(run_successful))
        print(
            "Param {} \n Build option {} \n Init option {}".format(
                sweep_params, self.build_kwargs, self.reinitialize_kwargs
            )
        )
        return output_dict
