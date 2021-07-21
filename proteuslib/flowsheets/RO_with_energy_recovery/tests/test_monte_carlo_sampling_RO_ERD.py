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

from proteuslib.flowsheets.RO_with_energy_recovery.monte_carlo_sampling_RO_ERD import run_parameter_sweep

@pytest.mark.component
def test_monte_carlo_sampling():

    # Define truth values
    truth_values = np.array([[2.84923065e-12, 2.95005437e-08, 9.80005773e-01, 2.73387066e+00, 4.58469101e-01],
                             [3.46351569e-12, 3.30797282e-08, 9.83385027e-01, 2.72795706e+00, 4.48196806e-01],
                             [3.61939655e-12, 3.06107079e-08, 9.79926626e-01, 2.72796370e+00, 4.46124238e-01],
                             [3.69412179e-12, 2.46992965e-08, 9.62620625e-01, 2.73246138e+00, 4.44857668e-01],
                             [3.73591412e-12, 3.33879140e-08, 9.77460037e-01, 2.72789671e+00, 4.44539702e-01],
                             [3.87531481e-12, 3.79140761e-08, 9.81571173e-01, 2.72602000e+00, 4.42998209e-01],
                             [4.15951955e-12, 3.52110687e-08, 9.61217760e-01, 2.72989529e+00, 4.39492394e-01],
                             [4.43270381e-12, 4.06688472e-08, 9.50731531e-01, 2.73127570e+00, 4.36625830e-01],
                             [4.81217268e-12, 4.23105397e-08, 9.77675085e-01, 2.72300670e+00, 4.34671785e-01],
                             [4.87240588e-12, 3.41378590e-08, 9.89554444e-01, 2.72040293e+00, 4.35064581e-01]])

    global_results = run_parameter_sweep(None)

    # Loop through individual values for specificity
    for value, truth_value in zip(truth_values.flatten(), global_results.flatten()):
        assert value == pytest.approx(truth_value)
