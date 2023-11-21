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

import functools
import logging
import multiprocessing

import numpy

from watertap.tools.parallel.results import LocalResults
from watertap.tools.parallel.parallel_manager import (
    build_and_execute,
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
        # constrain the number of child processes to the number of unique values to be run
        self.actual_number_of_subprocesses = min(
            self.max_number_of_subprocesses, len(all_parameters)
        )

        # split the parameters prameters for async run
        self.expected_samples = len(all_parameters)
        divided_parameters = numpy.array_split(all_parameters, self.expected_samples)

        build_and_execute_this_run = functools.partial(
            build_and_execute, do_build, do_build_kwargs, do_execute
        )

        actor = functools.partial(multiProcessingActor, build_and_execute_this_run)

        self._pool = multiprocessing.Pool(self.actual_number_of_subprocesses)
        self._results = self._pool.map_async(actor, enumerate(divided_parameters))
        self._pool.close()

    def gather(self):
        results = self._results.get()
        # sort the results by the process number to keep a deterministic ordering
        results.sort(key=lambda result: result.process_number)
        return results

    def results_from_local_tree(self, results):
        return results


def multiProcessingActor(
    build_and_execute_this_run,
    i_local_parameters,
):
    i, local_parameters = i_local_parameters
    return LocalResults(
        i, local_parameters, build_and_execute_this_run(local_parameters)
    )
