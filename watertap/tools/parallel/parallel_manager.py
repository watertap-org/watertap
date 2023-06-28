from abc import abstractmethod, ABC
import os

import numpy

from pyomo.common.dependencies import attempt_import

mpi4py, mpi4py_available = attempt_import("mpi4py")


class ParallelManager(ABC):

    ROOT_PROCESS_RANK = 0

    @classmethod
    def create_parallel_manager(cls, parallel_manager_class=None):
        """
        Create and return an instance of a ParallelManager, based on the libraries available in the
        runtime environment.

        Allows an optional python class to be passed in as parallel_manager_class. If so, this class
        is instantiated and returned rather than checking the local environment.
        """
        if parallel_manager_class is not None:
            return parallel_manager_class()

        if ParallelManager.has_mpi_peer_processes():
            return MPIParallelManager(mpi4py)

        return SingleProcessParallelManager()

    @classmethod
    def has_mpi_peer_processes(cls):
        """
        Returns whether the process was run as part of an MPI group with > 1 processes.
        """
        return mpi4py_available and mpi4py.MPI.COMM_WORLD.Get_size() > 1

    @abstractmethod
    def is_root_process(self):
        """
        Return whether the current process should be considered as the root of the parallel group it
        is a part of.
        """
        raise NotImplementedError

    @abstractmethod
    def get_rank(self):
        """
        Get the process's rank number in its parallel peer group.
        """
        raise NotImplementedError

    @abstractmethod
    def number_of_processes(self):
        """
        Return how many processes are part of the parallel peer group.
        """
        raise NotImplementedError

    @abstractmethod
    def sync_point(self):
        """
        Implement a synchronization point. All processes must reach this point before any continue.
        """
        raise NotImplementedError

    @abstractmethod
    def sync_data(self, data):
        """
        Synchronize a piece of data with all peer processes. The data parameter is either:
        - (if root) the source of truth that will be broadcast to all processes
        - (if not root) an empty buffer that will hold the data sent from root once sync_data returns
        """
        raise NotImplementedError

    def scatter(self, all_parameter_combos, fn):
        """
        Scatter a function out to a set of child processes, to be run
        for a list of parameters.
        - all_parameter_combos is a list, where each item represents the
        parameters for a single run.
        - fn is a function to be run, and should accept a single parameter.
        It will be called for each item in the all_parameter_combos list.
        """
        raise NotImplementedError

    def gather(self):
        """
        Gather the results of the computation that was kicked off via a previous
        scatter.
        - returns: an instance of GatherResults
        """
        raise NotImplementedError


class LocalResults:
    """
    Class representing the results of one process's run of an optimization routine.
    """

    def __init__(self, parameters, results):
        self.parameters = parameters
        self.results = results


class GatherResults:
    """
    Class representing all of the results of a sweep.
    - local_results is an instance of LocalResults representing the results
    from this process.
    - all_results is a list containing one instance of LocalResults for
    each process, including the local one.
    """

    def __init__(self, local_results, all_results):
        self.local_results = local_results
        self.all_results = all_results


class SingleProcessParallelManager(ParallelManager):
    def __init__(self):
        self.results = None

    def is_root_process(self):
        return True

    def get_rank(self):
        return self.ROOT_PROCESS_RANK

    def number_of_processes(self):
        return 1

    def sync_point(self):
        pass

    def sync_data(self, data):
        pass

    def scatter(self, all_parameter_combos, fn):
        self.results = LocalResults(all_parameter_combos, fn(all_parameter_combos))

    def gather(self):
        return GatherResults(self.results, [self.results])


class MPIParallelManager(ParallelManager):
    def __init__(self, mpi4py_lib):
        self.mpi4py = mpi4py_lib
        self.comm = self.mpi4py.MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.num_procs = self.comm.Get_size()

        self.results = None

    def is_root_process(self):
        return self.rank == self.ROOT_PROCESS_RANK

    def get_rank(self):
        return self.comm.Get_rank()

    def number_of_processes(self):
        return self.num_procs

    def sync_point(self):
        self.comm.Barrier()

    def sync_data(self, data):
        """
        Broadcast the array to all processes. this call acts as a synchronization point
        when run by multiple peer mpi processes.
        """
        if self.num_procs > 1:
            self.comm.Bcast(data, root=self.ROOT_PROCESS_RANK)

    def scatter(self, all_parameter_combos, fn):
        # Split the total list of combinations into NUM_PROCS chunks,
        # one per each of the MPI ranks
        divided_combo_array = numpy.array_split(all_parameter_combos, self.num_procs)

        # The current process's portion of the total workload
        local_combo_array = divided_combo_array[self.rank]

        self.results = LocalResults(local_combo_array, fn(local_combo_array))

    def gather(self):
        all_results = self.comm.allgather(self.results)
        return GatherResults(self.results, all_results)
