import numpy as np
from scipy.interpolate import griddata
from mpi4py import MPI
from idaes.core.util import get_solver

from proteuslib.tools.parameter_sweep import LinearSample, parameter_sweep
from proteuslib.flowsheets.full_treatment_train.flowsheet_components.flowsheet_limited import (solve_optimization, optimize, set_up_optimization)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Set up the solver
kwargs_flowsheet = {
    'has_bypass': True, 'has_desal_feed': False, 'is_twostage': True,
    'NF_type': 'ZO', 'NF_base': 'ion',
    'RO_type': '0D', 'RO_base': 'TDS', 'RO_level': 'simple'}
# solve_flowsheet_limited_NF(**kwargs_flowsheet)
# m, pass_fail = solve_optimization(system_recovery=0.78, max_conc_factor=3, **kwargs_flowsheet)
m = solve_optimization(system_recovery=0.5, max_conc_factor=3, **kwargs_flowsheet)

# Define the sampling type and ranges for three different variables
sweep_params = {}
sweep_params['Max Conc Factor'] = LinearSample(m.fs.max_conc_factor_target, 2.5, 3.5, 3)
sweep_params['System Recovery'] = LinearSample(m.fs.system_recovery_target, 0.5, 0.95, 24)
#sweep_params['Max RO Pressure'] = LinearSample(m.fs.max_allowable_pressure, 65e5, 300e5, 20)
sweep_params['Max RO Pressure'] = LinearSample(m.fs.max_allowable_pressure, 300e5, 65e5, 15)

# Define the list of output variables
outputs = {}
outputs['LCOW'] = m.fs.costing.LCOW
#outputs['EC'] = m.fs.costing.EC
outputs['NF Area'] = m.fs.NF.area
outputs['NF Flowrate'] = m.fs.NF.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']
outputs['Operating Pressure'] = m.fs.pump_RO.control_volume.properties_out[0].pressure
outputs['RO Recovery'] = m.fs.RO.recovery_mass_phase_comp[0, 'Liq', 'H2O']
outputs['NF Recovery'] = m.fs.NF.mass_transfer_phase_comp[0.0, 'Liq', 'H2O']
outputs['Bypass Fraction'] = m.fs.splitter.split_fraction[0, 'bypass']
#outputs['Required Reinit'] = m.fs.required_reinit

# Define the filename
output_filename = 'output/sweep_flowsheet_limited_%s_%s_%s.csv' % (kwargs_flowsheet['RO_type'], kwargs_flowsheet['RO_base'], kwargs_flowsheet['RO_level'])

global_results = parameter_sweep(m, sweep_params, outputs, output_filename, 
optimize_function=optimize, debugging_data_dir='output/local')#,
#reinitialize_function=solve_optimization, reinitialize_kwargs=kwargs_flowsheet)

if rank == 0:
    print(global_results)
    global_results_clean = np.copy(global_results)

    nc = len(sweep_params)

    # Build a mask of non-nan points
    mask = np.isfinite(global_results[:, nc])

    # Known points and their values
    x0 = global_results[mask, 0:nc]
    xi = global_results[:, 0:nc]

    for k in range(len(outputs)):
        y0 = global_results[mask, nc + k]
        yi = griddata(x0, y0, xi, method='linear', rescale=True)
        global_results_clean[~mask, nc + k] = yi[~mask]

    with open('%s' % (output_filename), 'r') as fp:
        header_line = fp.readline().split('#')[1].strip()

    np.savetxt('%s_clean.csv' % (output_filename.split('.')[0]), global_results_clean, delimiter=',', header=header_line)
