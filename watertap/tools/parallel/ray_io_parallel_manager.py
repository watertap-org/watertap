#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import numpy
import logging

from watertap.tools.parallel.results import LocalResults
from watertap.tools.parallel.parallel_manager import (
    parallelActor,
    ParallelManager,
)
from pyomo.common.dependencies import attempt_import

ray, ray_available = attempt_import("ray", defer_check=False)
import os
import platform


__author__ = "Alexander V. Dudchenko (SLAC)"

_log = logging.getLogger(__name__)


class RayIoParallelManager(ParallelManager):
    def __init__(self, number_of_subprocesses=1, **kwargs):
        self.max_number_of_subprocesses = number_of_subprocesses

        # this will be updated when child processes are kicked off
        self.actual_number_of_subprocesses = None

        # Future -> (process number, parameters). Used to keep track of the process number and parameters for
        # all in-progress futures
        self.running_futures = dict()
        # TODO: this should be deprciated once max resource option is avaiable
        self.cluster_mode = False
        if ray_available:
            # constrain the number of child processes to the number of unique values to be run
            self.setup_ray_cluster()

    def is_root_process(self):
        return True

    def get_rank(self):
        return self.ROOT_PROCESS_RANK

    def number_of_worker_processes(self):
        if self.actual_number_of_subprocesses is None:
            return self.max_number_of_subprocesses
        return self.actual_number_of_subprocesses

    def sync_with_peers(self):
        pass

    def sync_array_with_peers(self, data):
        pass

    def sync_pyobject_with_peers(self, obj):
        return obj

    def combine_data_with_peers(self, data):
        return [data]

    def sum_values_and_sync(self, sendbuf, recvbuf):
        recvbuf[0][:] = sendbuf[0][:]

    def gather_arrays_to_root(self, sendbuf, recvbuf_spec):
        receive_arr = recvbuf_spec[0]
        receive_sizes = recvbuf_spec[1]

        assert len(receive_arr) == sum(
            receive_sizes
        ), "Gathering arrays to root cannot be done with mismatched sizes"
        receive_arr[:] = sendbuf[:]

    def scatter(
        self,
        do_build,
        do_build_kwargs,
        do_execute,
        all_parameters,
    ):
        # over ride max_number of subprocess if user setus up cluster mode
        # by adding "ip_head" and  "redis_password" to their ENVS
        # and starting ray cluster before running parameter sweep
        if self.cluster_mode:
            self.actual_number_of_subprocesses = min(
                int(ray.cluster_resources()["CPU"]), len(all_parameters)
            )
        else:
            self.actual_number_of_subprocesses = min(
                self.max_number_of_subprocesses, len(all_parameters)
            )
        # split the parameters prameters for async run
        self.expected_samples = len(all_parameters)
        divided_parameters = numpy.array_split(all_parameters, self.expected_samples)
        # create queues, run queue will be used to store paramters we want to run

        actors = []
        paramActor = create_paramActor_class()
        # start ray actirs
        for cpu in range(self.actual_number_of_subprocesses):
            actors.append(
                paramActor.remote(
                    do_build,
                    do_build_kwargs,
                    do_execute,
                    divided_parameters[0],
                )
            )
        # create actor pool for load balancing
        actor_pool = ray.util.ActorPool(actors)
        # run in async.

        # load intoshared memory space
        run_vars_ray = ray.put(divided_parameters)
        self.results = actor_pool.map_unordered(
            lambda actor, var: actor.excute_with_order.remote(var, run_vars_ray),
            numpy.arange(self.expected_samples),
        )

    def gather(self):
        results = []
        for i, values, result in list(self.results):
            results.append(LocalResults(i, values, result))
        # sort the results by the process number to keep a deterministic ordering
        # results.sort(key=lambda result: result.process_number)

        results.sort(key=lambda result: result.process_number)
        self.shut_down_ray()
        return results

    def results_from_local_tree(self, results):
        return results

    # will need clean up
    def setup_ray_cluster(self):
        if ray.is_initialized() == False:
            try:
                """This will try to connect ot existing ray
                cluster, typical usage is for a cluster
                where ray needs to be started as head, with additional
                workers started on each node
                the parllel manaager will connect to this cluster
                useing ip_head address of head node, and its password
                if these are not found it will reveret to local mode."""
                _log.info(
                    "Connected to IP: {}, redis password: {}".format(
                        os.environ["ip_head"], os.environ["redis_password"]
                    )
                )

                ray.init(
                    include_dashboard=False,
                    address=os.environ["ip_head"],
                    _redis_password=os.environ["redis_password"],
                )
                _log.info("Nodes in ray cluster {}".format(ray.nodes()))
                _log.info("Resources, {}".format(ray.cluster_resources()))
                self.cluster_mode = True
            except KeyError:
                _log.info("Did not find ray cluster address, running in local mode")
                if ray.is_initialized() == False:
                    # ray.shutdown()
                    ray.init(include_dashboard=False)
                self.cluster_mode = False
        else:
            if platform.system() != "Linux" and self.cluster_mode == False:
                _log.info("Restarting ray")
                # restart ray on windows machine to deal with memoery issues
                if ray.is_initialized():
                    ray.shutdown()
                ray.init(include_dashboard=False)
                self.cluster_mode = False

    def shut_down_ray(self):
        """used to shutdown ray instance after run only in local mode
        if running on cluster or with head, we assume external head script will shut it down
        as all we do here is connect to it"""
        if self.cluster_mode == False:
            if ray.is_initialized():
                ray.shutdown()


def create_paramActor_class():
    @ray.remote(num_cpus=1)
    class paramActor(parallelActor):
        # this lets us track the order in execution
        def excute_with_order(self, order_index, local_parameters):
            local_parameters = local_parameters[order_index]
            return (
                order_index,
                local_parameters,
                self.execute(local_parameters),
            )

    return paramActor
