from idaes.core.util import get_solver
from proteuslib.tools.parameter_sweep import UniformSample, NormalSample, parameter_sweep

from RO_with_energy_recovery import (build,
	set_operating_conditions,
	initialize_system,
	solve,
	optimize)

# Set up the solver
solver = get_solver(options={'nlp_scaling_method': 'user-scaling'})

# Build, set, and initialize the system (these steps will change depending on the underlying model)
m = build()
set_operating_conditions(m, water_recovery=0.5, over_pressure=0.3, solver=solver)
initialize_system(m, solver=solver)

# Simulate once outside the parameter sweep to ensure everything is appropriately initialized 
solve(m, solver=solver)

# Define the sampling type and ranges for three different variables
sweep_params = {}
sweep_params['A_comp'] = NormalSample(m.fs.RO.A_comp, 4.0e-12, 0.5e-12)
sweep_params['B_comp'] = NormalSample(m.fs.RO.B_comp, 3.5e-8, 0.5e-8)
sweep_params['Porosity'] = UniformSample(m.fs.RO.spacer_porosity, 0.95, 0.99)

# Run the parameter sweep study using num_samples randomly drawn from the above range
output_filename = 'output/monte_carlo_results.csv'
num_samples = 100

global_results = parameter_sweep(m, sweep_params, {'EC':m.fs.specific_energy_consumption, 'LCOW': m.fs.costing.LCOW},
	output_filename, optimize_function=optimize, optimize_kwargs={'solver':solver}, num_samples=num_samples)

# Print the results for a quick check (also written to monte_carlo_results.csv)
print(global_results)
