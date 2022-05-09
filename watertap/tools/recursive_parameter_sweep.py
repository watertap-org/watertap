###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################
import numpy as np
import warnings
from enum import Enum
from watertap.tools.parameter_sweep import (
    _aggregate_results_arr,
    _build_combinations,
    _create_local_output_skeleton,
    _create_global_output,
    _default_optimize,
    _divide_combinations,
    _do_param_sweep,
    _init_mpi,
    _process_sweep_params,
    _process_results_filename,
    _save_results,
)

# ================================================================


def _filter_recursive_solves(model, sweep_params, outputs, recursive_local_dict, comm):

    # Figure out how many filtered solves did this rank actually do
    filter_counter = 0
    for case, content in recursive_local_dict.items():
        filter_counter += sum(content["solve_successful"])

    # Now that we have all of the local output dictionaries, we need to construct
    # a consolidated dictionary of successful solves.
    local_filtered_dict = _create_local_output_skeleton(
        model, sweep_params, outputs, filter_counter
    )
    local_filtered_dict["solve_successful"] = []

    # Populate local_successful_outputs
    offset = 0
    for case_number, content in recursive_local_dict.items():
        # Filter all of the sucessful solves
        optimal_indices = [
            idx for idx, success in enumerate(content["solve_successful"]) if success
        ]
        n_successful_solves = len(optimal_indices)
        stop = offset + n_successful_solves

        for key, item in content.items():
            if key != "solve_successful":
                for subkey, subitem in item.items():
                    local_filtered_dict[key][subkey]["value"][offset:stop] = subitem[
                        "value"
                    ][optimal_indices]

        # Place the solve status
        local_filtered_dict["solve_successful"].extend(
            [content["solve_successful"][i] for i in optimal_indices]
        )

        offset += n_successful_solves

    return local_filtered_dict, filter_counter


# ================================================================


def _aggregate_filtered_input_arr(
    global_filtered_dict, req_num_samples, comm, rank, num_procs
):

    global_filtered_values = np.zeros(
        (req_num_samples, len(global_filtered_dict["sweep_params"])), dtype=np.float64
    )

    if rank == 0:
        for i, (key, item) in enumerate(global_filtered_dict["sweep_params"].items()):
            global_filtered_values[:, i] = item["value"][:req_num_samples]

    if num_procs > 1:  # pragma: no cover
        comm.Bcast(global_filtered_values, root=0)

    return global_filtered_values


# ================================================================


def _aggregate_filtered_results(
    local_filtered_dict, req_num_samples, comm, rank, num_procs
):

    global_filtered_dict = _create_global_output(
        local_filtered_dict, req_num_samples, comm, rank, num_procs
    )
    global_filtered_results = _aggregate_results_arr(
        global_filtered_dict, req_num_samples, comm, rank, num_procs
    )
    global_filtered_values = _aggregate_filtered_input_arr(
        global_filtered_dict, req_num_samples, comm, rank, num_procs
    )

    return global_filtered_dict, global_filtered_results, global_filtered_values


# ================================================================


def recursive_parameter_sweep(
    model,
    sweep_params,
    outputs=None,
    csv_results_file_name=None,
    h5_results_file_name=None,
    optimize_function=_default_optimize,
    optimize_kwargs=None,
    reinitialize_function=None,
    reinitialize_kwargs=None,
    reinitialize_before_sweep=False,
    mpi_comm=None,
    debugging_data_dir=None,
    interpolate_nan_outputs=False,
    req_num_samples=None,
    seed=None,
):

    # Get an MPI communicator
    comm, rank, num_procs = _init_mpi(mpi_comm)

    # Convert sweep_params to LinearSamples
    sweep_params, sampling_type = _process_sweep_params(sweep_params)

    # Set the seed before sampling
    np.random.seed(seed)

    # Set up optimize_kwargs
    if optimize_kwargs is None:
        optimize_kwargs = dict()
    # Set up reinitialize_kwargs
    if reinitialize_kwargs is None:
        reinitialize_kwargs = dict()

    n_samples_remaining = req_num_samples
    num_total_samples = req_num_samples

    local_output_collection = {}
    loop_ctr = 0
    while n_samples_remaining > 0 and loop_ctr < 10:
        # Enumerate/Sample the parameter space
        global_values = _build_combinations(
            sweep_params, sampling_type, num_total_samples, comm, rank, num_procs
        )

        # divide the workload between processors
        local_values = _divide_combinations(global_values, rank, num_procs)
        local_num_cases = np.shape(local_values)[0]
        if loop_ctr == 0:
            true_local_num_cases = local_num_cases

        local_output_collection[loop_ctr] = _do_param_sweep(
            model,
            sweep_params,
            outputs,
            local_values,
            optimize_function,
            optimize_kwargs,
            reinitialize_function,
            reinitialize_kwargs,
            reinitialize_before_sweep,
            comm,
        )

        # Get the number of successful solves on this proc (sum of boolean flags)
        success_count = sum(local_output_collection[loop_ctr]["solve_successful"])
        failure_count = local_num_cases - success_count

        # Get the global number of successful solves and update the number of remaining samples
        if num_procs > 1:  # pragma: no cover
            global_success_count = np.zeros(1, dtype=np.float64)
            global_failure_count = np.zeros(1, dtype=np.float64)
            comm.Allreduce(
                np.array(success_count, dtype=np.float64), global_success_count
            )
            comm.Allreduce(
                np.array(failure_count, dtype=np.float64), global_failure_count
            )
        else:
            global_success_count = success_count
            global_failure_count = failure_count

        success_prob = global_success_count / (
            global_failure_count + global_success_count
        )

        if success_prob < 0.1:
            warnings.warn(
                f"Success rate of solves = {100.0*success_prob}%, consider adjusting sweep limits."
            )

        n_samples_remaining -= global_success_count

        # The total number of samples to generate at the next iteration is a multiple of the total remaining samples
        scale_factor = 2.0 / max(success_prob, 0.10)
        num_total_samples = int(np.ceil(scale_factor * n_samples_remaining))
        loop_ctr += 1

    # Now that we have all of the local output dictionaries, we need to construct
    # a consolidated dictionary based on a filter, e.g., optimal solves.
    local_filtered_dict, local_n_successful = _filter_recursive_solves(
        model, sweep_params, outputs, local_output_collection, comm
    )

    # if we are debugging
    if debugging_data_dir is not None:
        local_filtered_values = np.zeros(
            (local_n_successful, len(local_filtered_dict["sweep_params"])),
            dtype=np.float64,
        )
        for i, (key, item) in enumerate(local_filtered_dict["sweep_params"].items()):
            local_filtered_values[:, i] = item["value"][:]
    else:
        local_filtered_values = None

    # Not that we have all of the successful outputs in a consolidated dictionary locally,
    # we can now construct a global dictionary of successful solves.
    (
        global_filtered_dict,
        global_filtered_results,
        global_filtered_values,
    ) = _aggregate_filtered_results(
        local_filtered_dict, req_num_samples, comm, rank, num_procs
    )

    # Now we can save this
    if num_procs > 1:  # pragma: no cover
        comm.Barrier()

    global_save_data = _save_results(
        sweep_params,
        local_filtered_values,
        global_filtered_values,
        local_filtered_dict,
        global_filtered_dict,
        global_filtered_results,
        csv_results_file_name,
        h5_results_file_name,
        debugging_data_dir,
        comm,
        rank,
        num_procs,
        interpolate_nan_outputs,
    )

    return global_save_data
