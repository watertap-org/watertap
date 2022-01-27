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
import pyomo.environ as pyo
import sys
import os
import itertools
import warnings
import copy, pprint
import h5py

from scipy.interpolate import griddata
from enum import Enum, auto
from abc import abstractmethod, ABC
from idaes.core.util import get_solver

from idaes.surrogate.pysmo import sampling
from idaes.core.util.model_statistics import (variables_in_activated_equalities_set,
    expressions_set, total_objectives_set)
from pyomo.core.base.block import TraversalStrategy
from watertap.tools.parameter_sweep import *

np.set_printoptions(linewidth=200)

# ================================================================

def _force_exception(ctr):

    random_number = np.random.rand()
    # if (ctr % 2) == 0:
    #     return 1/0
    # else:
    #     return 1
    if random_number < 0.5:
        return 1/0
    else:
        return 1

# ================================================================

def _filter_recursive_solves(model, sweep_params, recursive_local_dict, true_local_num_cases, comm, filter_keyword="optimal"):

    pyomo_termination_condition = TerminationConditionMapping()
    try:
        assert filter_keyword in pyomo_termination_condition.mapping.keys()
    except:
        warnings.warn("Invalid filtering option specified. Filtering optimal values")
        filter_keyword = "optimal"

    # Figure out how many filtered solves did this rank actually do
    filter_counter = 0
    for case, content in recursive_local_dict.items():
        filter_counter += content["solve_status"].count(filter_keyword)

    # Now that we have all of the local output dictionaries, we need to construct
    # a consolidated dictionary of successful solves.
    local_filtered_dict = _create_local_output_skeleton(model, sweep_params, filter_counter,
                                                          variable_type="unfixed")
    local_filtered_dict["solve_status"] = []

    # Populate local_successful_outputs
    offset = 0
    for case_number, content in recursive_local_dict.items():
        # Filter all of the sucessful solves
        optimal_indices = [idx for idx, status in enumerate(content["solve_status"]) if status == "optimal"]
        n_successful_solves = len(optimal_indices)
        stop = offset+n_successful_solves

        for key, item in content.items():
            if key != "solve_status":
                for subkey, subitem in item.items():
                    local_filtered_dict[key][subkey]['value'][offset:stop] = subitem['value'][optimal_indices]

        # Place the solve status
        local_filtered_dict["solve_status"].extend([content["solve_status"][i] for i in optimal_indices])

        offset += n_successful_solves

    return local_filtered_dict

# ================================================================

def recursive_parameter_sweep(model, sweep_params, outputs, results_dir=None, results_fname=None,
        optimize_function=_default_optimize, optimize_kwargs=None, reinitialize_function=None,
        reinitialize_kwargs=None, reinitialize_before_sweep=False, mpi_comm=None, debugging_data_dir=None,
        interpolate_nan_outputs=False, req_num_samples=None, seed=None):

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

    n_samples_remaining = copy.deepcopy(req_num_samples)
    num_total_samples = copy.deepcopy(req_num_samples)

    local_output_collection = {}
    loop_ctr = 0
    while n_samples_remaining > 0 :
        # Enumerate/Sample the parameter space
        global_values = _build_combinations(sweep_params, sampling_type, num_total_samples, comm, rank, num_procs)

        # divide the workload between processors
        local_values = _divide_combinations(global_values, rank, num_procs)
        local_num_cases = np.shape(local_values)[0]
        if loop_ctr == 0:
            true_local_num_cases = local_num_cases

        _, local_output_collection[loop_ctr], fail_counter = _do_param_sweep(model, sweep_params, outputs, local_values,
            optimize_function, optimize_kwargs, reinitialize_function, reinitialize_kwargs, reinitialize_before_sweep, comm)

        n_successful_solves = local_num_cases - fail_counter
        n_successful_list = comm.allgather(n_successful_solves)
        n_samples_remaining -= sum(n_successful_list)
        num_total_samples = int(np.ceil(2 * n_samples_remaining))
        loop_ctr += 1

    # Now that we have all of the local output dictionaries, we need to construct
    # a consolidated dictionary based on a filter, e.g., optimal solves.
    local_filtered_dict = _filter_recursive_solves(model, sweep_params, local_output_collection,
        true_local_num_cases, comm, "optimal")

    # Not that we have all of the successful outputs in a consolidated dictionary locally,
    # we can now construct a global dictionary of successful solves.
    global_filtered_dict = _create_global_output(local_filtered_dict, req_num_samples, comm)


    # Now we can save this
    comm.Barrier()
    if rank == 0:
        _write_outputs(global_filtered_dict, results_dir, fname_no_extension=results_fname,
            write_h5=True, write_txt=True, txt_options="keys")

    return

# ================================================================


def parameter_sweep_deprecated(model, sweep_params, outputs, results_file=None, optimize_function=_default_optimize,
        optimize_kwargs=None, reinitialize_function=None, reinitialize_kwargs=None,
        reinitialize_before_sweep=False, mpi_comm=None, debugging_data_dir=None,
        interpolate_nan_outputs=False, num_samples=None, seed=None):

    '''
    This function offers a general way to perform repeated optimizations
    of a model for the purposes of exploring a parameter space while
    monitoring multiple outputs.
    If provided, writes single CSV file to ``results_file`` with all inputs and resulting outputs.

    Arguments:

        model : A Pyomo ConcreteModel containing a watertap flowsheet, for best
                results it should be initialized before being passed to this
                function.

        sweep_params: A dictionary containing the values to vary with the format
                      ``sweep_params['Short/Pretty-print Name'] =
                      (model.fs.variable_or_param[index], lower_limit, upper_limit, num_samples)``.
                      A uniform number of samples ``num_samples`` will be take between
                      the ``lower_limit`` and ``upper_limit``.

        outputs : A dictionary containing "short names" as keys and and Pyomo objects
                  on ``model`` whose values to report as values. E.g.,
                  ``outputs['Short/Pretty-print Name'] = model.fs.variable_or_expression_to_report``.

        results_file (optional) : The path and file name where the results are to be saved;
                                   subdirectories will be created as needed.

        optimize_function (optional) : A user-defined function to perform the optimization of flowsheet
                                       ``model`` and loads the results back into ``model``. The first
                                       argument of this function is ``model``\. The default uses the
                                       default IDAES solver, raising an exception if the termination
                                       condition is not optimal.

        optimize_kwargs (optional) : Dictionary of kwargs to pass into every call to
                                     ``optimize_function``. The first arg will always be ``model``,
                                     e.g., ``optimize_function(model, **optimize_kwargs)``. The default
                                     uses no kwargs.

        reinitialize_function (optional) : A user-defined function to perform the re-initialize the
                                           flowsheet ``model`` if the first call to ``optimize_function``
                                           fails for any reason. After ``reinitialize_function``, the
                                           parameter sweep tool will immediately call
                                           ``optimize_function`` again.

        reinitialize_kwargs (optional) : Dictionary or kwargs to pass into every call to
                                         ``reinitialize_function``. The first arg will always be
                                         ``model``, e.g.,
                                         ``reinitialize_function(model, **reinitialize_kwargs)``.
                                         The default uses no kwargs.

        reinitialize_before_sweep (optional): Boolean option to reinitialize the flow sheet model before
                                              every parameter sweep realization. The default is False.
                                              Note the parameter sweep model will try to reinitialize the
                                              solve regardless of the option if the run fails.

        mpi_comm (optional) : User-provided MPI communicator for parallel parameter sweeps.
                              If None COMM_WORLD will be used. The default is sufficient for most
                              users.

        debugging_data_dir (optional) : Save results on a per-process basis for parallel debugging
                                        purposes. If None no `debugging` data will be saved.

        interpolate_nan_outputs (optional) : When the parameter sweep has finished, interior values
                                             of np.nan will be replaced with a value obtained via
                                             a linear interpolation of their surrounding valid neighbors.
                                             If true, a second output file with the extension "_clean"
                                             will be saved alongside the raw (un-interpolated) values.

        num_samples (optional) : If the user is using sampling techniques rather than a linear grid
                                 of values, they need to set the number of samples

        seed (optional) : If the user is using a random sampling technique, this sets the seed

    Returns:

        save_data : A list were the first N columns are the values of the parameters passed
                    by ``sweep_params`` and the remaining columns are the values of the
                    simulation identified by the ``outputs`` argument.
    '''

    # Get an MPI communicator
    comm, rank, num_procs = _init_mpi(mpi_comm)

    # Convert sweep_params to LinearSamples
    sweep_params, sampling_type = _process_sweep_params(sweep_params)

    # Set the seed before sampling
    np.random.seed(seed)

    # Enumerate/Sample the parameter space
    global_values = _build_combinations(sweep_params, sampling_type, num_samples, comm, rank, num_procs)

    # divide the workload between processors
    local_values = _divide_combinations(global_values, rank, num_procs)

    # Initialize space to hold results
    local_num_cases = np.shape(local_values)[0]
    local_results = np.zeros((local_num_cases, len(outputs)))

    # Set up optimize_kwargs
    if optimize_kwargs is None:
        optimize_kwargs = dict()
    # Set up reinitialize_kwargs
    if reinitialize_kwargs is None:
        reinitialize_kwargs = dict()

    # Create the output skeleton for storing detailed data
    local_output_dict = _create_local_output_skeleton(model, sweep_params, local_num_cases,
                                                      variable_type=output_variable_type)

    local_solve_status_list = []
    fail_counter = 0

    # ================================================================
    # Run all optimization cases
    # ================================================================

    for k in range(local_num_cases):
        # Update the model values with a single combination from the parameter space
        _update_model_values(model, sweep_params, local_values[k, :])

        if reinitialize_before_sweep:
            # Forced reinitialization of the flowsheet if enabled
            reinitialize_function(model, **reinitialize_kwargs)

        try:
            # Simulate/optimize with this set of parameters
            results = optimize_function(model, **optimize_kwargs)
            pyo.assert_optimal_termination(results)
            local_solve_status_list.append(results.solver.termination_condition.name) # We will store status as a string

            # store the values of the optimization
            _update_local_output_dict(model, sweep_params, k, local_values[k, :], local_output_dict)

        except:
            # If the run is infeasible, report nan
            local_results[k, :] = np.nan
            previous_run_failed = True

            fail_counter += 1

        else:
            # If the simulation suceeds, report stats
            local_results[k, :] = [pyo.value(outcome) for outcome in outputs.values()]
            previous_run_failed = False

        if previous_run_failed and (reinitialize_function is not None):
            # We choose to re-initialize the model at this point
            try:
                reinitialize_function(model, **reinitialize_kwargs)
                optimize_function(model, **optimize_kwargs)
            except:
                # do we raise an error here?
                # nothing to do
                pass
            else:
                local_results[k, :] = [pyo.value(outcome) for outcome in outputs.values()]

    # ================================================================
    # Save results
    # ================================================================

    global_results = _aggregate_results(local_results, global_values, comm, num_procs)

    local_output_dict["solve_status"] = local_solve_status_list
    global_output_dict = _create_global_output(local_output_dict, num_samples, comm)

    # Make a directory for saved outputs
    if rank == 0:
        if results_file is not None:
            dirname = os.path.dirname(results_file)
            if dirname != '':
                os.makedirs(dirname, exist_ok=True)

        if debugging_data_dir is not None:
            os.makedirs(debugging_data_dir, exist_ok=True)

    if num_procs > 1:
        comm.Barrier()

    # Write a header string for all data files
    data_header = ','.join(itertools.chain(sweep_params,outputs))

    if debugging_data_dir is not None:
        # Create the local filename and data
        fname = os.path.join(debugging_data_dir, f'local_results_{rank:03}.csv')
        local_save_data = np.hstack((local_values, local_results))

        # Save the local data
        np.savetxt(fname, local_save_data, header=data_header, delimiter=', ', fmt='%.6e')

    # Create the global filename and data
    global_save_data = np.hstack((global_values, global_results))

    if rank == 0 and results_file is not None:
        # Save the global data
        np.savetxt(results_file, global_save_data, header=data_header, delimiter=',', fmt='%.6e')

        # Save the data of output dictionary
        _write_outputs(global_output_dict, dirname, txt_options="keys")
        # _read_output_h5("./output/output_dict.h5")


        if interpolate_nan_outputs:
            global_results_clean = _interp_nan_values(global_values, global_results)
            global_save_data_clean = np.hstack((global_values, global_results_clean))

            head, tail = os.path.split(results_file)

            if head == '':
                interp_file = 'interpolated_%s' % (tail)
            else:
                interp_file = '%s/interpolated_%s' % (head, tail)

            np.savetxt(interp_file, global_save_data_clean, header=data_header, delimiter=',', fmt='%.6e')

    return global_save_data

# ================================================================
