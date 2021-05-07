import pytest
import numpy as np
import sys
import os

from sim_LSRRO_Nstage import build_model, simulate, optimization, display_metrics
from proteuslib.tools.parameter_sweep import parameter_sweep

# ================================================================
# Build and simulate the model once as an initial condition
# ================================================================

if len(sys.argv) > 1:
    num_stages = int(sys.argv[1])
else:
    num_stages = 2

m = build_model(N=num_stages)
m = simulate(m, N=num_stages)

# ================================================================
# Set parameter sweep options
# sweep_params['Short/Pretty-print Name'] =
# (path.to.model.variable, lower_limit, upper_limit, num_samples)
# ================================================================

sweep_params = dict()
sweep_params['Membrane Cost'] = (m.fs.costing_param.mem_cost, 15, 60, 4)
sweep_params['Electricity Cost']  = (m.fs.costing_param.electricity_cost, 0.03, 0.15, 4)
# sweep_params['Water Permeability'] = (m.fs.Stage[:].A, 2e-12, 10e-12, 4)

# ================================================================
# Run the parameter sweep
# ================================================================

parameter_sweep(m, sweep_params, ['EC', 'LCOW'], objective='LCOW', output_dir='output_z',
    num_stages=num_stages, optimization=optimization, display_metrics=display_metrics)
