import ray
from ray.util import ActorPool
import pyomo.environ as pyo
import copy
import numpy as np
import os
import idaes.logger as idaeslog
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
from watertap.tools.parameter_sweep.async_param_sweep.async_param_sweep_func import (
    paramActor,
)


def setup_ray_cluster():
    if ray.is_initialized() == False:
        try:
            print(os.environ["ip_head"], os.environ["redis_password"])

            ray.init(
                include_dashboard=False,
                address=os.environ["ip_head"],
                _redis_password=os.environ["redis_password"],
            )
            print("Nodes in the Ray cluster:")
            print(ray.nodes())

            print(ray.cluster_resources())
        except KeyError:
            print("Did not find ray cluster address, running in local mode")
            if ray.is_initialized() == False:
                # ray.shutdown()
                ray.init(include_dashboard=False)
    else:
        if platform.system() != "Linux":
            # restart ray on windows machine to deal with memoery issues
            if ray.is_initialized():
                ray.shutdown()
            ray.init(include_dashboard=False)
        print("ray already running")


def do_rayio_sweep(param_options, run_dict, num_workers):
    # start ray cluster
    setup_ray_cluster()
    # ray.init(include_dashboard=False)
    if num_workers is None:
        num_workers = int(ray.cluster_resources()["CPU"])
    if len(run_dict) < num_workers:
        num_workers = len(run_dict)
    actors = setupRayActors(
        param_options,
        num_workers,
    )

    actor_pool = ActorPool(actors)
    results = actor_pool.map_unordered(
        lambda actor, var: actor.solve_problem.remote(var), run_dict
    )
    del actors

    return list(results)


def setupRayActors(
    param_options,
    num_workers,
):
    actors = []
    for cpu in range(num_workers):
        actors.append(paramActorRay.remote(param_options))
    time.sleep(1)
    return actors


@ray.remote(num_cpus=1, scheduling_strategy="SPREAD")
class paramActorRay(paramActor):
    pass
