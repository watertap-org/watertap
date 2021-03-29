import numpy as np
import pyomo.core.base as pyobase
import sys
import os

from sim_LSRRO_Nstage import optimization, display_metrics, display_state

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

        # OLD CODE SAVED AS EXAMPLE OF WHERE/HOW VALUES WERE CHANGED
        # TO GENERATE LSRRO RESULTS FOR NAWI WEBINAR.
        #
        # if 'recovery' in key:
        #     product_recovery = value
        # elif 'flow_mass_comp' in key:
        #     feed_mass_frac_NaCl = value
        #     # m.fs.M1.feed.flow_mass_comp[0, 'NaCl'].fix(feed_mass_frac_NaCl)
        #     # m.fs.M1.feed.flow_mass_comp[0, 'H2O'].fix(1-feed_mass_frac_NaCl)
        #     # m.fs.M1.upstream.flow_mass_comp[0, 'NaCl'].fix(feed_mass_frac_NaCl)
        #     # m.fs.M1.upstream.flow_mass_comp[0, 'H2O'].fix(1-feed_mass_frac_NaCl)
        #     m.fs.P1.inlet.flow_mass_comp[0, 'NaCl'].fix(feed_mass_frac_NaCl)
        #     m.fs.P1.inlet.flow_mass_comp[0, 'H2O'].fix(1-feed_mass_frac_NaCl)
        # elif "ele_cost" in key:
        #     m.fs.costing_param.electricity_cost.fix(value)
        # elif "mem_cost" in key:
        #     m.fs.costing_param.mem_cost.fix(value)
        # elif "water_perm" in key:
        #     for i in range(N):
        #         # These could be indexed stage[2], etc.
        #         stage = getattr(m.fs,"Stage"+repr(i+1))
        #         stage.A.fix(value)

# ================================================================

def run_param_sweep(m, num_stages, sweep_params, output_dir='output', mpi_comm=None):

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

    # Initialize space to hold local results
    num_cases = np.int64(np.shape(local_values)[0])
    local_EC = np.zeros(num_cases, dtype=np.float64)
    local_LCOW = np.zeros(num_cases, dtype=np.float64)

    for k in range(num_cases):
        # Update the model values with a single combination from the parameter space
        update_model_values(m, param_dict=sweep_params, values=local_values[k, :])

        try:
            # Simulate/optimize with this set of parameters
            # print(param_keys, local_values[k, :])
            m = optimization(m, 'LCOW', N=num_stages)

        except:
            # If the run is infeasible, report nan
            EC = np.nan
            LCOW = np.nan
            previous_run_failed = True

        else:
            # If the simulation suceeds, report stats
            # states = display_state(m, N=num_stages)
            EC, LCOW = display_metrics(m, N=num_stages)
            previous_run_failed = False

        local_EC[k] = EC
        local_LCOW[k] = LCOW

        if previous_run_failed:
            pass
            # We might choose to re-initialize the model at this point

    # Save the local results to a file (helpful for parallel debugging)
    fname = '%s/local_results_%03d.csv' % (output_dir, rank)
    save_data = np.hstack((local_values, local_EC[:, None], local_LCOW[:, None]))

    header = ''
    for k, v in sweep_params.items():
        header += '%s, ' % (k)
    header += 'EC, LCOW'

    np.savetxt(fname, save_data, header=header, delimiter=', ', fmt='%.6e')

    # ================================================================
    # Save results
    # ================================================================

    # Collect all results onto rank 0
    global_EC = np.zeros(np.shape(global_values)[0], dtype=np.float64)
    global_LCOW = np.zeros(np.shape(global_values)[0], dtype=np.float64)

    send_counts = np.zeros(num_procs, dtype=np.int64)
    comm.Gather(num_cases, send_counts, root=0)

    comm.Gatherv(local_EC, (global_EC, send_counts), root=0)
    comm.Gatherv(local_LCOW, (global_LCOW, send_counts), root=0)

    if rank == 0:
        # Save the global results to a file
        fname = '%s/%d_stage.csv' % (output_dir, num_stages)
        save_data = np.hstack((global_values, global_EC[:, None], global_LCOW[:, None]))
        np.savetxt(fname, save_data, header=header, delimiter=', ', fmt='%.6e')


# ================================================================

# def set_nested_attr(m, full_attr_path, value):
    
#     parent = m
#     attr_list = full_attr_path.split('.')
#     stripped_list = []
#     nn = len(attr_list)
    
#     for attr in attr_list:
#         if '[' in attr:
#             attr_name = attr.split('[')[0]
#             idx = int(attr.split('[')[1].split(']')[0])
#         else:
#             attr_name = attr
#             idx = None
            
#         stripped_list.append([attr_name, idx])
    
#     for k, item in enumerate(stripped_list):
#         attr_name = item[0]
#         idx = item[1]
#         child = getattr(parent, attr_name)

#         if type(child) == list and idx is None:
#             raise ValueError('The "%s" attribute is a list and must have an element specified with an integer index' % attr_name)
        
#         if k < nn-1:
#             if idx is not None:
#                 child = child[idx]
    
#             parent = child

#         else:
#             if idx is not None:
#                 child[idx] = value
#             else:
#                 setattr(parent, attr_name, value)

#     return m

# ================================================================

# def set_nested_attr_eval(m, key, value):


#     temp = eval('%s' % (key))
#     print(type(temp))
#     print(temp)
#     if type(temp) == var:
#         print('yep')
#     exit()


#     if 'recovery' in key:
#         product_recovery = value
#     elif 'flow_mass_comp' in key:
#         feed_mass_frac_NaCl = value
#         # m.fs.M1.feed.flow_mass_comp[0, 'NaCl'].fix(feed_mass_frac_NaCl)
#         # m.fs.M1.feed.flow_mass_comp[0, 'H2O'].fix(1-feed_mass_frac_NaCl)
#         # m.fs.M1.upstream.flow_mass_comp[0, 'NaCl'].fix(feed_mass_frac_NaCl)
#         # m.fs.M1.upstream.flow_mass_comp[0, 'H2O'].fix(1-feed_mass_frac_NaCl)
#         m.fs.P1.inlet.flow_mass_comp[0, 'NaCl'].fix(feed_mass_frac_NaCl)
#         m.fs.P1.inlet.flow_mass_comp[0, 'H2O'].fix(1-feed_mass_frac_NaCl)
#     elif "ele_cost" in key:
#         m.fs.costing_param.electricity_cost.fix(value)
#     elif "mem_cost" in key:
#         m.fs.costing_param.mem_cost.fix(value)
#     elif "water_perm" in key:
#         for i in range(N):
#             stage = getattr(m.fs,"Stage"+repr(i+1))
#             stage.A.fix(value)
