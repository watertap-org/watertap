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

from watertap.tools.parallel.results import LocalResults
from watertap.tools.parallel.parallel_manager import ParallelManager, run_sweep


class MPIParallelManager(ParallelManager):
    def __init__(self, mpi4py_lib, **kwargs):
        self.mpi4py = mpi4py_lib
        self.comm = self.mpi4py.MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.num_procs = self.comm.Get_size()

        self.results = None

    def is_root_process(self):
        return self.rank == self.ROOT_PROCESS_RANK

    def get_rank(self):
        return self.comm.Get_rank()

    def number_of_worker_processes(self):
        return self.num_procs

    def sync_with_peers(self):
        self.comm.Barrier()

    def sync_data_with_peers(self, data):
        """
        Broadcast the array to all processes. this call acts as a synchronization point
        when run by multiple peer mpi processes.
        """
        if self.num_procs > 1:
            self.comm.Bcast(data, root=self.ROOT_PROCESS_RANK)

    def scatter(
        self,
        param_sweep_instance,
        common_sweep_args,
        rebuild_common_sweep_args_fn,
        rebuild_common_sweep_args_kwargs,
        all_parameter_combos,
    ):
        # Split the total list of combinations into NUM_PROCS chunks,
        # one per each of the MPI ranks
        divided_combo_array = numpy.array_split(all_parameter_combos, self.num_procs)

        # The current process's portion of the total workload
        local_combo_array = divided_combo_array[self.rank]

        self.results = LocalResults(
            self.get_rank(),
            local_combo_array,
            run_sweep(
                param_sweep_instance,
                common_sweep_args,
                rebuild_common_sweep_args_fn,
                rebuild_common_sweep_args_kwargs,
                local_combo_array,
            ),
        )

    def gather(self):
        return self.comm.allgather(self.results)

    def results_from_local_tree(self, results):
        return [
            next(
                (
                    result
                    for result in results
                    if result.process_number == self.get_rank()
                )
            )
        ]
