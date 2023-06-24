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

    @classmethod
    def remove_unpicklable_state(cls, parameter_sweep_instance):
        """
        Remove and return any state from the ParameterSweep object that cannot be
        pickled, to make the instance picklable. Needed in order to use the
        ConcurrentFuturesParallelManager.
        """
        saved_state = {
            "parallel_manager": parameter_sweep_instance.parallel_manager,
            "comm": parameter_sweep_instance.comm,
            "writer": parameter_sweep_instance.writer,
        }

        parameter_sweep_instance.parallel_manager = None
        parameter_sweep_instance.comm = None
        parameter_sweep_instance.writer = None
        return saved_state

    @classmethod
    def restore_unpicklable_state(cls, parameter_sweep_instance, state):
        """
        Restore a collection of saved state that was removed in order to pickle
        the ParameterSweep object.
        """
        parameter_sweep_instance.parallel_manager = state.get("parallel_manager", None)
        parameter_sweep_instance.comm = state.get("comm", None)
        parameter_sweep_instance.writer = state.get("writer", None)

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
        param_sweep_instance,
        common_sweep_args,
        rebuild_common_sweep_args_fn,
        rebuild_common_sweep_args_kwargs,
        all_parameter_combos,
    ):
        """
        Scatter a function out to a set of child processes, to be run for a list of parameters.
        - param_sweep_instance: a reference to the ParameterSweep object
        - common_params: a list of the parameters that every processes' invocation of the sweep should receive.
        - all_parameter_combos: a list, where each item represents the parameters for a single run.
        """
        raise NotImplementedError

    def gather(self):
        """
        Gather the results of the computation that was kicked off via a previous scatter.
        Returns:
        - a list of LocalResults, representing the results for each process
        """
        raise NotImplementedError

    def results_from_local_tree(self, results):
        """
        Given a list of LocalResults objects, return a sublist of the ones the current process
        is responsible for.
        """
        raise NotImplementedError


def return_arguments_as_list(*args):
    return args


def run_sweep(
    param_sweep_instance,
    common_sweep_args,
    rebuild_common_sweep_args_fn,
    rebuild_common_sweep_args_kwargs,
    local_combo_array,
):
    """
    Run the parameter sweep object's sweep function. Implemented as a top-level function in
    this module so that it is picklable for the ConcurrentFuturesParallelManager.
    Parameters are:
    - common_sweep_args: the parameters that all processes' invocations of the sweep function need
    - rebuild_common_sweep_args_fn: an optional function that will be used to rebuild the common sweep args.
    - rebuild_common_sweep_args_kwargs: kwargs for the rebuild_common_sweep_args_fn.
    - local_combo_array: the list of combinations to be run on the current processes
    """

    # if we've been directed to rebuild the sweep args, do so
    if rebuild_common_sweep_args_fn is not None:
        if rebuild_common_sweep_args_kwargs is None:
            rebuild_common_sweep_args_kwargs = dict()
        common_sweep_args = rebuild_common_sweep_args_fn(
            **rebuild_common_sweep_args_kwargs
        )

    if param_sweep_instance.config.custom_do_param_sweep is not None:
        return param_sweep_instance.custom_do_param_sweep(
            *common_sweep_args, local_combo_array
        )

    return param_sweep_instance._do_param_sweep(*common_sweep_args, local_combo_array)
