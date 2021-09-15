import numpy as np
from scipy.interpolate import griddata
from idaes.core.util import get_solver

from proteuslib.tools.parameter_sweep import _init_mpi, LinearSample, parameter_sweep
from proteuslib.flowsheets.full_treatment_train.example_flowsheets.flowsheet_mvp import (solve_optimization, optimize, set_up_optimization)

comm, rank, num_procs = _init_mpi()

#================================
# Set up the solver
#================================
kwargs_flowsheet = {
    'has_bypass': True, 'has_desal_feed': False, 'is_twostage': True, 'has_ERD': True,
    'NF_type': 'ZO', 'NF_base': 'ion',
    'RO_type': '0D', 'RO_base': 'TDS', 'RO_level': 'detailed'}
m = solve_optimization(system_recovery=0.80, **kwargs_flowsheet)

#================================
# Define the sampling type and ranges for three different variables
#================================
sweep_params = {}
sweep_params['System Recovery'] = LinearSample(m.fs.system_recovery_target, 0.5, 0.95, 4)
sweep_params['Max RO Pressure'] = LinearSample(m.fs.max_allowable_pressure, 300e5, 65e5, 4)

#================================
# Define the list of output variables
#================================
outputs = {}
outputs['LCOW'] = m.fs.costing.LCOW

# NF Values
outputs['NF Area'] = m.fs.NF.area
outputs['NF Flowrate'] = m.fs.NF.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']
outputs['NF Recovery'] = m.fs.NF.mass_transfer_phase_comp[0.0, 'Liq', 'H2O']

# RO-1 Values
outputs['RO-1 Area'] = m.fs.RO.area
outputs['RO-1 Recovery'] = m.fs.RO.recovery_mass_phase_comp[0, 'Liq', 'H2O']
outputs['RO-1 Operating Pressure'] = m.fs.pump_RO.control_volume.properties_out[0].pressure

# RO-2 Values
outputs['RO-2 Area'] = m.fs.RO2.area
outputs['RO-2 Recovery'] = m.fs.RO2.recovery_mass_phase_comp[0, 'Liq', 'H2O']
outputs['RO-2 Operating Pressure'] = m.fs.pump_RO2.control_volume.properties_out[0].pressure

outputs['Bypass Fraction'] = m.fs.splitter.split_fraction[0, 'bypass']

# Collect capital and operating-cost outputs
for unit_name in ['NF', 'pump_NF', 'RO', 'RO2', 'pump_RO', 'pump_RO2', 'ERD']:
    unit = getattr(m.fs, unit_name)
    
    cap_cost = unit.costing.capital_cost.value
    op_cost = unit.costing.operating_cost.value

    outputs['capital_cost_%s' % (unit_name)] = unit.costing.capital_cost
    outputs['operating_cost_%s' % (unit_name)] = unit.costing.operating_cost

outputs['capital_cost_total'] = m.fs.costing.capital_cost_total
outputs['operating_cost_total'] = m.fs.costing.operating_cost_total

#================================
# Start the parameter sweep
#================================
output_filename = 'output/sweep_mvp_results_%s_%s_%s.csv' % (kwargs_flowsheet['RO_type'], kwargs_flowsheet['RO_base'], kwargs_flowsheet['RO_level'])

global_results = parameter_sweep(m, sweep_params, outputs, output_filename, 
optimize_function=optimize, debugging_data_dir='output/local')#,

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
