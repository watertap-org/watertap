import numpy as np
import pyomo.core.base as pyobase
import sys
import os

# ================================================================

def build_and_divide_combinations(d, rank, num_procs):

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

def update_model_values(m, param_dict=None, values=None):

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

def run_param_sweep(m, sweep_params, outputs, output_dir='output', mpi_comm=None,
    num_stages=2, optimization=None, display_metrics=None):

    # Get an MPI communicator
    if mpi_comm is None:
        try:
            from mpi4py import MPI

        except:
            print('The parallel manager functions require a version')
            print('of mpi4py to be installed in this environment.')
            print('Defaulting to rank = 0, num_procs = 1.')

            comm = None
            rank = 0
            num_procs = 1

        else:
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            num_procs = comm.Get_size()

    else:
        print('Passing MPI communicators is not currently supported')
        comm = None
        rank = 0
        num_procs = 1

    # Make a directory for saved outputs
    os.makedirs('%s' % (output_dir), exist_ok=True)

    # Enumerate all possibilities and divide the workload between processors
    local_values, global_values = build_and_divide_combinations(sweep_params, rank, num_procs)

    # Initialize space to hold results
    local_num_cases = np.int64(np.shape(local_values)[0])
    local_results = np.zeros((local_num_cases, len(outputs)), dtype=np.float64)

    global_num_cases = np.int64(np.shape(global_values)[0])
    global_results = np.zeros((global_num_cases, len(outputs)), dtype=np.float64)

    # ================================================================
    # Run all optimization cases
    # ================================================================

    for k in range(local_num_cases):
        # Update the model values with a single combination from the parameter space
        update_model_values(m, param_dict=sweep_params, values=local_values[k, :])

        try:
            # Simulate/optimize with this set of parameters
            m = optimization(m, 'LCOW', N=num_stages)

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

    # Save the local results to a file (helpful for parallel debugging)
    fname = '%s/local_results_%03d.csv' % (output_dir, rank)
    save_data = np.hstack((local_values, local_results))

    data_header = ''
    for k, v in sweep_params.items():
        data_header += '%s, ' % (k)
    data_header += ', '.join(outputs)

    np.savetxt(fname, save_data, header=data_header, delimiter=', ', fmt='%.6e')

    # Collect the number of result values to be sent from each process
    send_counts = np.zeros(num_procs, dtype=np.int64)
    comm.Gather(np.int64(np.size(local_results)), send_counts, root=0)

    # Collect the global results results onto rank 0
    comm.Gatherv(local_results, (global_results, send_counts), root=0)

    if rank == 0:
        # Save the global results to a file
        fname = '%s/%d_stage.csv' % (output_dir, num_stages)
        save_data = np.hstack((global_values, global_results))
        np.savetxt(fname, save_data, header=data_header, delimiter=', ', fmt='%.6e')

# ================================================================