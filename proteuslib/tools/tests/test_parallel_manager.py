import pytest
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from proteuslib.tools.parallel_manager import build_and_divide_combinations
# from sim_LSRRO_2stage import build_model, simulate, optimization, display_metrics, display_state
from sim_LSRRO_Nstage import build_model, simulate, optimization, display_metrics, display_state

try:
    from mpi4py import MPI

except:
    print('The parallel manager functions require a version \n \
           of mpi4py be installed in this environemnt. \n \
           Defaulting to rank = 0, num_procs = 1.')

    comm = None
    rank = 0
    num_procs = 1

else:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    num_procs = comm.Get_size()

# ================================================================
# Set parameter sweep options
# ================================================================

sweep_params = dict()
# sweep_params['fs.P1.control_volume.properties_out[0].pressure'] = [55e5, 65e5, 4]
# sweep_params['fs.P2.control_volume.properties_out[0].pressure'] = [65e5, 75e5, 4]
sweep_params['m.fs.M1.feed.flow_mass_comp'] = [0.025, 0.05, 3]
sweep_params['m.fs.recovery'] = [0.3, 0.7, 3]
# sweep_params['m.fs.recovery'] = [0.3, 0.6, 3]

param_keys, values, global_values = build_and_divide_combinations(sweep_params, rank, num_procs)

num_cases = np.int64(np.shape(values)[0])
results_local_EC = np.zeros(num_cases, dtype=np.float64)
results_local_LCOW = np.zeros(num_cases, dtype=np.float64)

fp = open('local_results_%03d.txt' % (rank), 'w')

# ================================================================
# Run all cases
# ================================================================

num_stages = 2

for k in range(num_cases):
    if k == 0:
        # Build the model only once, the first time
        m = build_model(N=num_stages)
        m = simulate(m, N=num_stages);

    try:
        # Simulate/optimize with this set of parameters
        fp.write('%.3e, %.3e, ' % (values[k, 0], values[k, 1]))
        m = optimization(m, 'LCOW', params=param_keys, values=values[k, :], N=num_stages)

    except:
        # If the run is infeasible, report nan
        EC = np.nan
        LCOW = np.nan
        
    else:
        # If the simulation suceeds, report stats
        # states = display_state(m, N=num_stages)
        EC, LCOW = display_metrics(m, N=num_stages)

    fp.write('%.3e, %.3e\n' % (EC, LCOW))

    results_local_EC[k] = EC
    results_local_LCOW[k] = LCOW

fp.close()

# ================================================================
# Print and visualize results
# ================================================================

# Collect all results onto rank 0
results_global_EC = np.zeros(np.shape(global_values)[0], dtype=np.float64)
results_global_LCOW = np.zeros(np.shape(global_values)[0], dtype=np.float64)

send_counts = np.zeros(num_procs, dtype=np.int64)
comm.Gather(num_cases, send_counts, root=0)

comm.Gatherv(results_local_EC, (results_global_EC, send_counts), root=0)
comm.Gatherv(results_local_LCOW, (results_global_LCOW, send_counts), root=0)

if rank == 0:

    data_x = 'm.fs.M1.feed.flow_mass_comp'
    data_y = 'm.fs.recovery'

    fp = open('raw_data2.csv', 'w')
    fp.write('%s, %s, %s, %s\n' % (data_x, data_y, 'EC', 'LCOW'))

    for k in range(np.size(results_global_EC)):
        fp.write('%f, %f, %f, %f\n' % (global_values[k, 0], global_values[k, 1], results_global_EC[k], results_global_LCOW[k]))

    fp.close()

    # print(results_global)


    # id_x = param_keys.index(data_x)
    # id_y = param_keys.index(data_y)

    # ln_x = sweep_params[data_x]
    # ln_y = sweep_params[data_y]

    # x_mesh, y_mesh = np.meshgrid(np.linspace(ln_x[0], ln_x[1], ln_x[2]), 
    #                              np.linspace(ln_y[0], ln_y[1], ln_y[2]))

    # points = global_values[:, (id_x, id_y)]
    # xi = (x_mesh, y_mesh)

    # c_data = griddata(points, results_global, xi)

    # print(c_data)

    # plt.contourf(x_mesh, y_mesh, c_data)
    # plt.plot(points[:, 0], points[:, 1], '.', color='k')
    # plt.title(save_val)
    # plt.xlabel(data_x)
    # plt.ylabel(data_y)
    # plt.colorbar()
    # plt.savefig('contour.png')
    # plt.savefig('contour.pdf')

print('global', global_values)
print('values', values)

