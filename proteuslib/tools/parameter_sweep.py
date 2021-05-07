import numpy as np
import pyomo.core.base as pyobase
import sys
import os

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

    for k, item in enumerate(param_dict.items()):

        key = item[0]
        param = item[1][0]

        if isinstance(param, pyobase.var._GeneralVarData):
            if param.is_indexed():
                # If an indexed variable, fix all values to values[k]
                for p in param:
                    p.fix(values[k])
            else:
                # Otherwise, fix the single value to values[k]
                param.fix(values[k])

        elif isinstance(param, pyobase.param._ParamData):
            if param.is_indexed():
                # If an indexed param, set all values to values[k]
                for p in param:
                    p.value(values[k])
            else:
                # Otherwise, set the single value to values[k]
                param.value(values[k])

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

def parameter_sweep(m, sweep_params, outputs, objective=None, output_dir='output', mpi_comm=None,
    num_stages=2, optimization=None, display_metrics=None, save_debugging_data=False):

    '''
    This function offers a general way to perform repeated optimizations
    of a model for the purposes of exploring a parameter space while
    monitoring multiple outputs:
        m : A Pyomo model containing a proteuslib flowsheet, for best results
            it should be initialized before being passed to this function
        sweep_params: A dictionary containing the values to vary with the 
                      format sweep_params['Short/Pretty-print Name'] =
                      (path.to.model.variable, lower_limit, upper_limit, num_samples)
        outputs : A list of strings indicating which elements from the potential
                  output dictionary to monitor
        objective : A string indicating the objective function to use
        output_dir : The directory to save all output files
        mpi_comm : User-provided MPI communicator
        num_stages : The number of stages to use (specific to the NStage model)
        optimization : A user-defined function to perform the optimization of flowsheet m
        display_metrics : Au ser-defined function to calculate outputs of flowsheet m
        save_debugging_data : Optionally, save results on a per-process basis for 
                              parallel debugging purposes
    Returns:
        None,
        Writes global CSV file with all inputs and resulting outputs
    '''

    # Get an MPI communicator
    comm, rank, num_procs = _init_mpi(mpi_comm)

    # Enumerate all possibilities and divide the workload between processors
    local_values, global_values = _build_and_divide_combinations(sweep_params, rank, num_procs)

    # Initialize space to hold results
    local_num_cases = np.shape(local_values)[0]
    local_results = np.zeros((local_num_cases, len(outputs)))

    # ================================================================
    # Run all optimization cases
    # ================================================================

    for k in range(local_num_cases):
        # Update the model values with a single combination from the parameter space
        _update_model_values(m, sweep_params, local_values[k, :])

        try:
            # Simulate/optimize with this set of parameters
            m = optimization(m, objective, N=num_stages)

        except:
            # If the run is infeasible, report nan
            local_results[k, :] = np.nan
            previous_run_failed = True

        else:
            # If the simulation suceeds, report stats
            local_results[k, :] = display_metrics(m, N=num_stages, outputs=outputs)
            previous_run_failed = False

        if previous_run_failed:
            # We might choose to re-initialize the model at this point
            pass

    # ================================================================
    # Save results
    # ================================================================

    global_results = _aggregate_results(local_results, global_values, comm, num_procs)

    # Make a directory for saved outputs
    if rank == 0:
        os.makedirs('%s' % (output_dir), exist_ok=True)

    if num_procs > 1:
        comm.Barrier()

    # Write a header string for all data files
    data_header = ''
    for k, v in sweep_params.items():
        data_header += '%s, ' % (k)
    data_header += ', '.join(outputs)

    if save_debugging_data:
        # Create the local filename and data
        fname = '%s/local_results_%03d.csv' % (output_dir, rank)
        save_data = np.hstack((local_values, local_results))

        # Save the local data
        np.savetxt(fname, save_data, header=data_header, delimiter=', ', fmt='%.6e')

    if rank == 0:
        # Create the global filename and data
        fname = '%s/global_results.csv' % (output_dir)
        save_data = np.hstack((global_values, global_results))

        # Save the global data
        np.savetxt(fname, save_data, header=data_header, delimiter=', ', fmt='%.6e')

# ================================================================