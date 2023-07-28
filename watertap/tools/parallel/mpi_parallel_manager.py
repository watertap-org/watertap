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
from watertap.tools.parallel.parallel_manager import build_and_execute, ParallelManager


class MPIParallelManager(ParallelManager):
    def __init__(self, MPI, **kwargs):
        self.MPI = MPI
        self.comm = self.MPI.COMM_WORLD
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

    def sync_array_with_peers(self, data):
        """
        Broadcast the array to all processes. this call acts as a synchronization point
        when run by multiple peer mpi processes.
        """
        if self.num_procs > 1:
            self.comm.Bcast(data, root=self.ROOT_PROCESS_RANK)

    def sync_pyobject_with_peers(self, obj):
        """
        Broadcast the object to all processes from the root. this call acts as a synchronization point
        when run by multiple peer mpi processes.
        """
        return self.comm.bcast(obj, root=self.ROOT_PROCESS_RANK)

    def combine_data_with_peers(self, data):
        return self.comm.allgather(data)

    def gather_arrays_to_root(self, sendbuf, recvbuf_spec):
        self.comm.Gatherv(sendbuf, recvbuf_spec, root=self.ROOT_PROCESS_RANK)

    def sum_values_and_sync(self, sendbuf, recvbuf):
        self.comm.Allreduce(sendbuf, recvbuf)

    def scatter(
        self,
        do_build,
        do_build_kwargs,
        do_execute,
        all_parameters,
    ):

        # Split the total list of values into NUM_PROCS chunks,
        # one per each of the MPI ranks
        divided_parameters = numpy.array_split(all_parameters, self.num_procs)

        # The current process's portion of the total workload
        local_parameters = divided_parameters[self.rank]

        results = build_and_execute(
            do_build, do_build_kwargs, do_execute, local_parameters
        )
        self.results = LocalResults(
            self.get_rank(),
            local_parameters,
            results,
        )

    def gather(self):
        return self.comm.allgather(self.results)

    def results_from_local_tree(self, results):
        return [
            result for result in results if result.process_number == self.get_rank()
        ]
