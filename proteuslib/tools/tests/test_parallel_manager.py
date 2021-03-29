import pytest
import numpy as np
import sys
import os

from sim_LSRRO_Nstage import build_model, simulate
from proteuslib.tools.parallel_manager import run_param_sweep

# ================================================================
# Build and simulate the model once as an initial condition
# ================================================================

if len(sys.argv) > 1:
    num_stages = int(sys.argv[1])
else:
    num_stages = 2

m = build_model(N=num_stages)
m = simulate(m, N=num_stages);

# ================================================================
# Set parameter sweep options
# sweep_params['Short/Pretty-print Name'] =
# (path.to.model.variable, lower_limit, upper_limit, num_samples)
# ================================================================

sweep_params = dict()

# sweep_params['fs.P1.control_volume.properties_out[0].pressure'] = [55e5, 65e5, 4]
# sweep_params['fs.P2.control_volume.properties_out[0].pressure'] = [65e5, 75e5, 4]

#sweep_params['m.fs.M1.feed.flow_mass_comp'] = [0.005, 0.155, 72]
#sweep_params['m.fs.M1.feed.flow_mass_comp'] = [0.025, 0.15, 72]
#sweep_params['m.fs.recovery'] = [0.3, 0.7, 72]

sweep_params['Membrane Cost'] = (m.fs.costing_param.mem_cost, 15, 60, 4)
sweep_params['Electricity Cost']  = (m.fs.costing_param.electricity_cost, 0.03, 0.15, 4)
# sweep_params['Water Permeability'] = (m.fs.Stage[:].A, 2e-12, 10e-12, 4)

# ================================================================
# Run the parameter sweep
# ================================================================

run_param_sweep(m, num_stages, sweep_params, output_dir='output')
