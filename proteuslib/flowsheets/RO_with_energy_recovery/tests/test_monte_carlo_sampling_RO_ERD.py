###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################
import pytest
import numpy as np

from idaes.core.util import get_solver
from proteuslib.tools.parameter_sweep import UniformSample, NormalSample, parameter_sweep

from proteuslib.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (build,
    set_operating_conditions,
    initialize_system,
    solve,
    optimize)


@pytest.mark.component
def test_monte_carlo_sampling():

    # Define truth values
    truth_values = np.array([[2.84923065e-12, 2.95005437e-08, 9.80005773e-01, 2.36410287e+00, 4.44141629e-01],
                         [3.46351569e-12, 3.30797282e-08, 9.83385027e-01, 2.24418284e+00, 4.27201202e-01],
                         [3.61939655e-12, 3.06107079e-08, 9.79926626e-01, 2.21839505e+00, 4.23913360e-01],
                         [3.69412179e-12, 2.46992965e-08, 9.62620625e-01, 2.21024752e+00, 4.22844107e-01],
                         [3.73591412e-12, 3.33879140e-08, 9.77460037e-01, 2.18944444e+00, 4.20825636e-01],
                         [3.87531481e-12, 3.79140761e-08, 9.81571173e-01, 2.15808848e+00, 4.17338985e-01],
                         [4.15951955e-12, 3.52110687e-08, 9.61217760e-01, 2.11424952e+00, 4.11325437e-01],
                         [4.43270381e-12, 4.06688472e-08, 9.50731531e-01, 2.08852919e+00, 4.05242154e-01],
                         [4.81217268e-12, 4.23105397e-08, 9.77675085e-01, 2.11018817e+00, 4.00395308e-01],
                         [4.87240588e-12, 3.41378590e-08, 9.89554444e-01, 2.04159947e+00, 4.00711632e-01]])

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
    sweep_params['Spacer_porosity'] = UniformSample(m.fs.RO.spacer_porosity, 0.95, 0.99)

    # Run the parameter sweep study using num_samples randomly drawn from the above range
    num_samples = 10

    # Suppress the CSV-file output for this test
    output_filename = None 

    global_results = parameter_sweep(m, sweep_params, {'EC':m.fs.specific_energy_consumption, 'LCOW': m.fs.costing.LCOW},
        output_filename, optimize_function=optimize, optimize_kwargs={'solver':solver}, num_samples=num_samples, seed=1)

    # Loop through individual values for specificity
    for value, truth_value in zip(truth_values.flatten(), global_results.flatten()):
        assert value == pytest.approx(truth_value)

