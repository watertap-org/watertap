import multiprocessing

import pyomo.environ as pyo
import copy
import numpy as np
import os

import platform
import time
from pyomo.common.collections import ComponentSet, ComponentMap


from idaes.core.util import to_json, from_json

import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

from watertap.tools.analysis_tools.model_state_tool import modelStateStorage

import watertap.tools.analysis_tools.step_optimize_tool as stepTool
from watertap.tools.parameter_sweep.async_param_sweep.async_utils import (
    _create_component_output_skeleton,
    _create_local_output_skeleton,
    _merge_copied_dicts,
    _convert_sweep_params_to_dict,
    _convert_outputs_to_dict,
    _update_model_values_from_dict,
)


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
        model_state = modelStateStorage(model)
        json_string = model_state.get_dict_state()
        # json_string = to_json(model, return_json_string=True)

    else:
        json_string = None
    run_queue = multiprocessing.Queue()
    return_queue = multiprocessing.Queue()
    for param in run_dict:
        # print(param)
        run_queue.put(param)
    actors = []
    if len(run_dict) < num_workers:
        num_workers = len(run_dict)
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
    results = []
    while len(results) < len(run_dict):
        if return_queue.empty() == False:
            result = return_queue.get()
            results.append(result)

    for cpu in range(num_workers):
        actors[cpu].join()

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
            # from_json(self.model, sd=None, fname=None, s=self.json_string)

            self.model_init = True
            self.model.state_after_box_solve = modelStateStorage(
                self.model, restore_dict=self.json_string
            )

    def init_model(self):
        if self.json_string is not None:
            self.build_model()

        self.reinitialize_function(self.model, **self.reinitialize_kwargs)
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

    def step_optimize(
        self, m, solver=None, try_final_first=False, check_termination=False
    ):
        # fsTools.solve(m)
        solver = get_solver()
        results = stepTool.step_optimize(
            m,
            solver,
            self.optimize_function,
            self.model.state_after_box_solve,
            steps=3,
            try_final_first=True,
            re_steps=2,
        )
        return results

    def solve_problem(self, sweep_params):
        order_index = sweep_params["order_index"]
        sweep_params.pop("order_index", None)

        run_successful = False

        _update_model_values_from_dict(self.model, sweep_params)
        # results = self.step_optimize(self.model)
        if self.get_probe_result(self.model, self.reinitialize_kwargs):
            for i in ["Try #0", "Try #1"]:
                if (
                    self.reinitialize_function is not None
                    and self.model_init == False
                    or self.reinitialize_before_sweep
                ):
                    try:
                        self.build_model()
                        if self.reinitialize_function is not None:
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
                        results = self.step_optimize(
                            self.model
                        )  # , **self.optimize_kwargs)
                        pyo.assert_optimal_termination(results)
                        run_successful = True
                        break
                    except:
                        run_successful = False
                        self.model_init = False
                        if (
                            self.reinitialize_before_sweep == False
                            and self.reinitialize_function is None
                        ):
                            break

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
