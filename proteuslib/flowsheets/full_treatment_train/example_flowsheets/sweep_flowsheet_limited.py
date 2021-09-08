from idaes.core.util import get_solver
from proteuslib.tools.parameter_sweep import LinearSample, parameter_sweep

from proteuslib.flowsheets.full_treatment_train.example_flowsheets.flowsheet_limited import (solve_optimization, optimize)

# Set up the solver
kwargs_flowsheet = {
    'has_bypass': True, 'has_desal_feed': False, 'is_twostage': True,
    'NF_type': 'ZO', 'NF_base': 'ion',
    'RO_type': '0D', 'RO_base': 'TDS', 'RO_level': 'detailed'}
# solve_flowsheet_limited_NF(**kwargs_flowsheet)
m, pass_fail = solve_optimization(system_recovery=0.78, max_conc_factor=3, **kwargs_flowsheet)

# Define the sampling type and ranges for three different variables
sweep_params = {}
sweep_params['System Recovery'] = LinearSample(m.fs.system_recovery_target, 0.5, 0.95, 20)
sweep_params['Max RO Pressure'] = LinearSample(m.fs.max_allowable_pressure, 65e5, 300e5, 20)
sweep_params['Max Conc Factor'] = LinearSample(m.fs.max_conc_factor_target, 2.5, 3.5, 3)

# Define the list of output variables
outputs = {}
outputs['NF Area'] = m.fs.NF.area
outputs['NF Flowrate'] = m.fs.NF.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']
outputs['Operating Pressure'] = m.fs.pump_RO.control_volume.properties_out[0].pressure
outputs['RO Recovery'] = m.fs.RO.recovery_mass_phase_comp[0, 'Liq', 'H2O']
outputs['NF Recovery'] = m.fs.NF.mass_transfer_phase_comp[0.0, 'Liq', 'H2O']
outputs['Bypass Fraction'] = m.fs.splitter.split_fraction[0, 'bypass']

# Define the filename
output_filename = 'output/sweep_flowsheet_limited_%s_%s_%s.csv' % (kwargs_flowsheet['RO_type'], kwargs_flowsheet['RO_base'], kwargs_flowsheet['RO_level'])

global_results = parameter_sweep(m, sweep_params, outputs, output_filename, optimize_function=optimize)

print(global_results)
