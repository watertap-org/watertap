import pytest
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from proteuslib.tools.parallel_manager import build_and_divide_combinations
from sim_LSRRO_2stage import build_model, simulate, optimization, display_metrics, display_state

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
sweep_params['m.fs.M1.feed.flow_mass_comp'] = [0.025, 0.15, 4]
sweep_params['m.fs.recovery'] = [0.3, 0.7, 4]
# sweep_params['m.fs.recovery'] = [0.3, 0.6, 3]

save_val = 'EC'

param_keys, values, global_values = build_and_divide_combinations(sweep_params, rank, num_procs)

num_cases = np.int64(np.shape(values)[0])
results_local = np.zeros(num_cases, dtype=np.float64)

# ================================================================
# Run all cases
# ================================================================

for k in range(num_cases):
    if k == 0:
        # Build the model only once, the first time
        m = build_model()
        m = simulate(m);

    try:
        # Simulate/optimize with this set of parameters
        m = optimization(m, save_val, param_keys, values[k, :])

    except:
        # If the run is infeasible, report nan
        sim_val = np.nan
        
    else:
        # If the simulation suceeds, report stats
        states = display_state(m)
        metrics = display_metrics(m)

        # sim_val = states['Mixed prod'][0]
        sim_val = metrics[save_val]

    results_local[k] = sim_val

# ================================================================
# Print and visualize results
# ================================================================

# Collect all results onto rank 0
results_global = np.zeros(np.shape(global_values)[0], dtype=np.float64)

send_counts = np.zeros(num_procs, dtype=np.int64)
comm.Gather(num_cases, send_counts, root=0)

comm.Gatherv(results_local, (results_global, send_counts), root=0)

if rank == 0:

    data_x = 'm.fs.M1.feed.flow_mass_comp'
    data_y = 'm.fs.recovery'

    fp = open('raw_data.csv', 'w')
    fp.write('%s, %s, %s\n' % (data_x, data_y, save_val))

    for k in range(np.size(results_global)):
        fp.write('%f, %f, %f\n' % (global_values[k, 0], global_values[k, 1], results_global[k]))

    fp.close()

    print(results_global)


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

