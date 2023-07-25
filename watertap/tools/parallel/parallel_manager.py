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
    def sync_data_with_peers(self, data):
        """
        Synchronize a piece of data with all peers processes. The data parameter is either:
        - (if root) the source of truth that will be broadcast to all processes
        - (if not root) an empty buffer that will hold the data sent from root once the function returns
        """
        raise NotImplementedError

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

    def gather(self):
        """
        Gather the results of the computation that was kicked off via a previous scatter.
        Returns:
        - a list of LocalResults, representing the results for each process. each result will
        be the return value of the do_execute function from scatter() for one process.
        """
        raise NotImplementedError

    def results_from_local_tree(self, results):
        """
        Given a list of LocalResults objects, return a sublist of the ones the current process
        is responsible for.
        """
        raise NotImplementedError


def build_and_execute(do_build, do_build_kwargs, do_execute, local_parameters):
    """
    Entrypoint for implementations of the parallel manager to use for running the
    build and execute functions. Defined at the top level so that it's picklable.

    For a description of the first three arguments, see the scatter() function above.
    The fourth argument is the list of local parameters that should be run by this process.
    """
    execute_args = do_build(**do_build_kwargs)
    results = do_execute(local_parameters, *execute_args)
    return results
