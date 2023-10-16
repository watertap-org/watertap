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

from abc import abstractmethod, ABC


class ParallelManager(ABC):
    ROOT_PROCESS_RANK = 0

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
    def number_of_worker_processes(self):
        """
        Return how many total processes are running local simulations.
        """
        raise NotImplementedError

    @abstractmethod
    def sync_with_peers(self):
        """
        Implement a synchronization point with any peer processes.
        """
        raise NotImplementedError

    @abstractmethod
    def sync_array_with_peers(self, data):
        """
        Synchronize an array with all peers processes. The data parameter is either:
        - (if root) the source of truth that will be broadcast to all processes
        - (if not root) an empty buffer that will hold the data sent from root once the function returns
        """
        raise NotImplementedError

    @abstractmethod
    def sync_pyobject_with_peers(self, obj):
        """
        Synchronize a python object with all peer processes. The obj parameter is either:
        - (if root) the source of truth that will be broadcast to all processes
        - (if not root) ignored

        Different from sync_array_with_peers in that it returns the synced object rather than
        using an existing buffer.
        """
        raise NotImplementedError

    @abstractmethod
    def combine_data_with_peers(self, data):
        """
        Combine the data from each peer into a list. The return value will be a list of all
        the data elements provided by each process that calls this function.

        With multiple processes, this must:
        - act as a synchronization point
        - return the list in the same order on each process
        """
        raise NotImplementedError

    @abstractmethod
    def gather_arrays_to_root(self, sendbuf, recvbuf_spec):
        """
        Gather the data in the send buffer to the root process. Parameters are:
        - sendbuf: the data to be sent from this process
        - recvbuf_spec: a tuple of (receive buffer, list of sizes) where the list of sizes is
        how much data will come from each process. Ignored if not the root process.
        """
        raise NotImplementedError

    @abstractmethod
    def sum_values_and_sync(self, sendbuf, recvbuf):
        """
        Sum a list of values from each process and sync the result to all processes. Parameters:
        - sendbuf: an array of values the local process is contributing to the sum
        - recvbuf: an array of a single value that will hold the summed result (same on
        each process) when the function returns.

        Note: this is a specific case of a global reduce (with sum() as the reducing function).
        If more than sum() is needed in the future, this function should be extended to receive
        an Op argument.
        """
        raise NotImplementedError

    @abstractmethod
    def scatter(
        self,
        do_build,
        do_build_kwargs,
        do_execute,
        all_parameters,
    ):
        """
        Scatter the specified execution out, as defined by the implementation's parallelism, for
        a list of parameters.
        Args:
        - do_build: a function that builds the arguments necessary for the execution function. expected to
        return a list that will be exploded and passed into the do_execute function as arguments.
        - do_build_kwargs: a dictionary of keyword arguments for the do_build function
        - do_execute: the execution function. expected to take in the list of local parameters
        as its first argument. any arguments after that should match the list returned by the
        do_build function.
        - all_parameters: a list of all parameters to be run. included so that
        different implementations of the parallel manager can make decisions about splitting and
        parallelization.
        """
        raise NotImplementedError

    @abstractmethod
    def gather(self):
        """
        Gather the results of the computation that was kicked off via a previous scatter.
        Returns:
        - a list of LocalResults, representing the results for each process. each result will
        be the return value of the do_execute function from scatter() for one process.
        """
        raise NotImplementedError

    @abstractmethod
    def results_from_local_tree(self, results):
        """
        Given a list of LocalResults objects, return a sublist of the ones the current process
        is responsible for.
        """
        raise NotImplementedError


# TODO this probably should be owned by parameer sweep as its PS specific


def build_and_execute(do_build, do_build_kwargs, do_execute, local_parameters):
    """
    Entrypoint for implementations of the parallel manager to use for running the
    build and execute functions. Defined at the top level so that it's picklable.

    For a description of the first three arguments, see the scatter() function above.
    The fourth argument is the list of local parameters that should be run by this process.
    """
    pa = parallelActor(do_build, do_build_kwargs, do_execute, local_parameters)
    results = pa.execute(local_parameters)
    # execute_args = do_build(**do_build_kwargs)
    # results = do_execute(local_parameters, *execute_args)
    return results


# TODO this probably should be owned by parameer sweep as its PS specific
class parallelActor:
    def __init__(self, do_build, do_build_kwargs, do_execute, local_parameters):
        self.do_build = do_build
        self.do_build_kwargs = do_build_kwargs
        self.do_execute = do_execute
        self.local_parameters = local_parameters
        self.build_model()

    def build_model(self):
        self.exec_args = self.do_build(**self.do_build_kwargs)

    def execute(self, local_parameters):
        result = self.do_execute(local_parameters, *self.exec_args)
        return result
