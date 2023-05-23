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


from watertap.tools.parameter_sweep.async_param_sweep.async_utils import (
    _create_component_output_skeleton,
    _create_local_output_skeleton,
    _merge_copied_dicts,
    _convert_sweep_params_to_dict,
    _convert_outputs_to_dict,
    _update_model_values_from_dict,
)

from watertap.tools.parameter_sweep.async_param_sweep.async_sweep_logic

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
