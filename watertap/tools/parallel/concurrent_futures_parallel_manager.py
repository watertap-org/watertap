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

from concurrent import futures

from watertap.tools.parallel.results import LocalResults
from watertap.tools.parallel.parallel_manager import ParallelManager, run_sweep


class ConcurrentFuturesParallelManager(ParallelManager):
    def __init__(self, number_of_subprocesses=1, **kwargs):
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

    def sync_data_with_peers(self, data):
        pass

    def scatter(
        self,
        param_sweep_instance,
        common_sweep_args,
        rebuild_common_sweep_args_fn,
        rebuild_common_sweep_args_kwargs,
        all_parameter_combos,
    ):

        # constrain the number of child processes to the number of unique parameter combinations to be run
        self.actual_number_of_subprocesses = min(
            self.max_number_of_subprocesses, len(all_parameter_combos)
        )

        # split the parameters into chunks, one for each child process
        divided_combo_array = numpy.array_split(
            all_parameter_combos, self.actual_number_of_subprocesses
        )

        # create an executor and kick off the child processes that will perform the computation
        self.executor = futures.ProcessPoolExecutor(
            max_workers=self.actual_number_of_subprocesses
        )

        for i in range(self.actual_number_of_subprocesses):
            local_params = divided_combo_array[i]
            # save the mapping of future -> (process number, params that it's running)
            self.running_futures[
                self.executor.submit(
                    run_sweep,
                    param_sweep_instance,
                    common_sweep_args,
                    rebuild_common_sweep_args_fn,
                    rebuild_common_sweep_args_kwargs,
                    local_params,
                )
            ] = (i, local_params)

    def gather(self):
        results = []
        try:
            execution_results = futures.wait(self.running_futures.keys())
            for future in execution_results.done:
                process_number, params = self.running_futures[future]
                results.append(LocalResults(process_number, params, future.result()))

            if len(execution_results.not_done) > 0:
                print(
                    f"{len(execution_results.not_done)} out of {len(self.running_futures.keys())} total subprocesses did not finish and provide results"
                )
        finally:
            self.executor.shutdown()

        # sort the results by the process number to keep a deterministic ordering
        results.sort(key=lambda result: result.process_number)
        return results

    def results_from_local_tree(self, results):
        return results
