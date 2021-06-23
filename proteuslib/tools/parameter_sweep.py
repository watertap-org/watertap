###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################
import numpy as np
import pyomo.environ as pyo
import sys
import os
import itertools

from idaes.core.util import get_solver

# ================================================================

def _init_mpi(mpi_comm=None):

    if mpi_comm is None:
        try:
            from mpi4py import MPI

        except:
            return None, 0, 1

        else:
            mpi_comm = MPI.COMM_WORLD

    return mpi_comm, mpi_comm.Get_rank(), mpi_comm.Get_size()

# ================================================================

def _build_and_divide_combinations(d, rank, num_procs):

    param_values = []

    for k, v in d.items():
        # Build a vector of discrete values for this parameter
        # and record the parameter's name
        # v[1] = start, v[2] = stop, v[3] = resolution (number of elements)
        p = np.linspace(v[1], v[2], v[3])
        param_values.append(p)

    num_var_params = len(param_values)

    # Form an array with every possible combination of parameter values
    full_combo_array = np.array(np.meshgrid(*param_values))
    full_combo_array = full_combo_array.T.reshape(-1, num_var_params)

    # Split the total list of combinations into NUM_PROCS chunks,
    # one per each of the MPI ranks
    divided_combo_array = np.array_split(full_combo_array, num_procs, axis=0)

    # Return only this rank's portion of the total workload
    local_combo_array = divided_combo_array[rank]

    return local_combo_array, full_combo_array

# ================================================================

def _update_model_values(m, param_dict, values):

    for k, item in enumerate(param_dict.values()):

        param = item[0]

        # Fix the single value to values[k]
        if param.is_variable_type():
            param.fix(values[k])

        elif param.is_parameter_type():
            # Fix the single value to values[k]
            param.set_value(values[k])

        else:
            raise RuntimeError(f"Unrecognized Pyomo object {param}")

# ================================================================

def _aggregate_results(local_results, global_values, comm, num_procs):

    if num_procs > 1:
        local_results = local_results.astype(np.float64)

        global_results = np.zeros((np.shape(global_values)[0], np.shape(local_results)[1]), dtype=np.float64)

        # Collect the number of result values to be sent from each process
        send_counts = np.zeros(num_procs, dtype=np.int64)
        comm.Gather(np.int64(np.size(local_results)), send_counts, root=0)

        # Collect the global results results onto rank 0
        comm.Gatherv(local_results, (global_results, send_counts), root=0)

    else:
        global_results = np.copy(local_results)

    return global_results

# ================================================================

def _default_optimize(model, options=None, tee=False):
    '''
    Default optimization function used in parameter_sweep.
    Optimizes ``model`` using the IDAES default solver.
    Raises a RuntimeError if the TerminationCondition is not optimal

    Arguments:

        model : A Pyomo ConcreteModel to optimize

        options (optional) : Solver options to pass into idaes.core.util.get_solver.
                             Default is None
        tee (options) : To display the solver log. Default it False

    '''
    solver = get_solver(options=options)
    results = solver.solve(m, tee=tee)

    if results.solver.termination_condition != pyo.TerminationCondition.optimal:
        raise RuntimeError("The solver failed to converge to an optimal solution. "
                           "This suggests that the user provided infeasible inputs "
                           "or that the model is poorly scaled.")

# ================================================================

def parameter_sweep(model, sweep_params, outputs, results_file, optimize_function=_default_optimize,
        optimize_kwargs=None, reinitialize_function=None, reinitialize_kwargs=None,
        mpi_comm=None, debugging_data_dir=None):

    '''
    This function offers a general way to perform repeated optimizations
    of a model for the purposes of exploring a parameter space while
    monitoring multiple outputs. 
    Writes single CSV file to ``results_file`` with all inputs and resulting outputs.

    Arguments:

        model : A Pyomo ConcreteModel containing a proteuslib flowsheet, for best 
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

        results_file : The path and file name where the results are to be saved; subdirectories
                       will be created as needed.

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

        mpi_comm (optional) : User-provided MPI communicator for parallel parameter sweeps.
                              If None COMM_WORLD will be used. The default is sufficient for most
                              users.

        debugging_data_dir (optional) : Save results on a per-process basis for parallel debugging
                                        purposes. If None no `debugging` data will be saved.

    Returns:
        None
    '''

    # Get an MPI communicator
    comm, rank, num_procs = _init_mpi(mpi_comm)

    # Enumerate all possibilities and divide the workload between processors
    local_values, global_values = _build_and_divide_combinations(sweep_params, rank, num_procs)

    # Initialize space to hold results
    local_num_cases = np.shape(local_values)[0]
    local_results = np.zeros((local_num_cases, len(outputs)))

    # Set up optimize_kwargs
    if optimize_kwargs is None:
        optimize_kwargs = dict()
    # Set up reinitialize_kwargs
    if reinitialize_kwargs is None:
        reinitialize_kwargs = dict()

    # ================================================================
    # Run all optimization cases
    # ================================================================

    for k in range(local_num_cases):
        # Update the model values with a single combination from the parameter space
        _update_model_values(model, sweep_params, local_values[k, :])

        try:
            # Simulate/optimize with this set of parameters
            optimize_function(model, **optimize_kwargs)

        except:
            # If the run is infeasible, report nan
            local_results[k, :] = np.nan
            previous_run_failed = True

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

    # Make a directory for saved outputs
    if rank == 0:
        dirname = os.path.dirname(results_file)
        if dirname != '':
            os.makedirs(dirname, exist_ok=True)
        if debugging_data_dir is not None:
            os.makedirs(debugging_data_dir, exist_ok=True)

    if num_procs > 1:
        comm.Barrier()

    # Write a header string for all data files
    data_header = ', '.join(itertools.chain(sweep_params,outputs))

    if debugging_data_dir is not None:
        # Create the local filename and data
        fname = os.path.join(debugging_data_dir, f'local_results_{rank:03}.csv')
        save_data = np.hstack((local_values, local_results))

        # Save the local data
        np.savetxt(fname, save_data, header=data_header, delimiter=', ', fmt='%.6e')

    if rank == 0:
        # Create the global filename and data
        save_data = np.hstack((global_values, global_results))

        # Save the global data
        np.savetxt(results_file, save_data, header=data_header, delimiter=', ', fmt='%.6e')

# ================================================================
