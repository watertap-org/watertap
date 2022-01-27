###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

import pytest
import os
import numpy as np
import pyomo.environ as pyo

from pyomo.environ import value

from watertap.tools.parameter_sweep import (_init_mpi,
                                               _build_combinations,
                                               _divide_combinations,
                                               _update_model_values,
                                               _aggregate_results,
                                               _interp_nan_values,
                                               _process_sweep_params,
                                               _write_output_to_h5,
                                               _read_output_h5,
                                               _create_local_output_skeleton,
                                               _create_global_output,
                                               parameter_sweep,
                                               LinearSample,
                                               UniformSample,
                                               NormalSample,
                                               SamplingType)
from watertap.tools.recursive_parameter_sweep import *
from watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (build,
    set_operating_conditions,
    initialize_system,
    solve,
    optimize)

# -----------------------------------------------------------------------------

@pytest.mark.component
def test_recursive_parameter_sweep():

    seed = None
    results_file = "monte_carlo_results.csv"
    read_from_file = False

    # Set up the solver
    solver = get_solver()

    # Build, set, and initialize the system (these steps will change depending on the underlying model)
    m = build()
    set_operating_conditions(m, water_recovery=0.5, over_pressure=0.3, solver=solver)
    initialize_system(m, solver=solver)

    # Simulate once outside the parameter sweep to ensure everything is appropriately initialized
    solve(m, solver=solver)

    # Define the sampling type and ranges for three different variables
    if read_from_file:
        sweep_params = get_sweep_params_from_yaml(m, 'mc_sweep_params.yaml')
    else:
        sweep_params = get_sweep_params(m, use_LHS=use_LHS)

    # Define the outputs to be saved
    outputs = {}
    outputs['EC'] = m.fs.specific_energy_consumption
    outputs['LCOW'] = m.fs.costing.LCOW

    # Run the parameter sweep study using num_samples randomly drawn from the above range
    num_samples = 2

    # Run the parameter sweep
    recursive_parameter_sweep(m, sweep_params, outputs, results_file=results_file,
        optimize_function=optimize, optimize_kwargs={'solver':solver}, req_num_samples=num_samples,
        seed=seed, reinitialize_before_sweep=False, reinitialize_function=initialize_system,
        reinitialize_kwargs={'solver':solver}
