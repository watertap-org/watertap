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


from watertap.tools.parameter_sweep.async_param_sweep.async_param_sweep_func import (
    paramActor,
)


def do_mp_sweep(param_options, run_dict, num_workers):
    # creat queses for work load
    if num_workers is None:
        num_workers = int(multiprocessing.cpu_count())
    if len(run_dict) < num_workers:
        num_workers = len(run_dict)
    run_queue = multiprocessing.Queue()
    return_queue = multiprocessing.Queue()
    # load parameters into queue
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
                    param_options,
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

    return list(results)


def setupMCActors(
    queue,
    return_queue,
    param_options,
    process_id,
):
    param_actor = paramActor(param_options)
    while True:
        if queue.empty():
            break
        else:
            sweep_params = queue.get()

            result = param_actor.solve_problem(sweep_params)
            return_queue.put(result)
