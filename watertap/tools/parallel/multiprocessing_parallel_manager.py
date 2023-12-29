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

import logging
import multiprocessing
from queue import Empty as EmptyQueue

import numpy

from watertap.tools.parallel.results import LocalResults
from watertap.tools.parallel.parallel_manager import (
    parallelActor,
    ParallelManager,
)


_logger = logging.getLogger(__name__)


class MultiprocessingParallelManager(ParallelManager):
    def __init__(
        self,
        number_of_subprocesses=1,
        **kwargs,
    ):
        self.max_number_of_subprocesses = number_of_subprocesses

        # this will be updated when child processes are kicked off
        self.actual_number_of_subprocesses = None

        # Future -> (process number, parameters). Used to keep track of the process number and parameters for
        # all in-progress futures
        self.running_futures = dict()

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
        recvbuf[:] = sendbuf

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
        # constrain the number of child processes to the number of unique values to be run
        self.actual_number_of_subprocesses = min(
            self.max_number_of_subprocesses, len(all_parameters)
        )

        # split the parameters prameters for async run
        self.expected_samples = len(all_parameters)
        divided_parameters = numpy.array_split(all_parameters, self.expected_samples)
        # create queues, run queue will be used to store paramters we want to run
        # and return_queue is used to store results
        self.run_queue = multiprocessing.Queue()
        self.return_queue = multiprocessing.Queue()
        for i, param in enumerate(divided_parameters):
            # print(param)
            self.run_queue.put([i, param])
        # setup multiprocessing actors
        self.actors = []

        for cpu in range(self.actual_number_of_subprocesses):
            self.actors.append(
                multiprocessing.Process(
                    target=multiProcessingActor,
                    args=(
                        self.run_queue,
                        self.return_queue,
                        do_build,
                        do_build_kwargs,
                        do_execute,
                        divided_parameters[0],
                    ),
                )
            )
            self.actors[-1].start()

    def gather(self):
        results = []
        # collect result from the actors
        while len(results) < self.expected_samples:
            try:
                i, values, result = self.return_queue.get()
                results.append(LocalResults(i, values, result))
            except EmptyQueue:
                break
        self._shut_down()
        # sort the results by the process number to keep a deterministic ordering
        results.sort(key=lambda result: result.process_number)
        return results

    def _shut_down(self):
        n_shut_down = 0
        for process in self.actors:
            _logger.debug("Attempting to shut down %s", process)
            process.kill()
            n_shut_down += 1
        self.actors.clear()
        _logger.debug("Shut down %d processes", n_shut_down)

    def results_from_local_tree(self, results):
        return results


# This function is used for running the actors in multprocessing
def multiProcessingActor(
    queue,
    return_queue,
    do_build,
    do_build_kwargs,
    do_execute,
    local_parameters,
):
    actor = parallelActor(do_build, do_build_kwargs, do_execute, local_parameters)
    while True:
        try:
            msg = queue.get()
        except EmptyQueue:
            return
        i, local_parameters = msg
        result = actor.execute(local_parameters)
        return_queue.put([i, local_parameters, result])
