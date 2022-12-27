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

from watertap.examples.flowsheets.RO_with_energy_recovery.monte_carlo_sampling_RO_ERD import (
    run_parameter_sweep,
)


@pytest.mark.integration
def test_monte_carlo_sampling():

    # Define truth values (defined with rng_seed=1)
    truth_values = np.array(
        [
            [2.849231e-12, 2.950054e-08, 9.800058e-01, 2.754121e00, 4.949683e-01],
            [3.463516e-12, 3.307973e-08, 9.833850e-01, 2.747771e00, 4.854287e-01],
            [3.619397e-12, 3.061071e-08, 9.799266e-01, 2.748404e00, 4.835434e-01],
            [3.694122e-12, 2.469930e-08, 9.626206e-01, 2.755873e00, 4.821188e-01],
            [3.735914e-12, 3.338791e-08, 9.774600e-01, 2.748674e00, 4.818635e-01],
            [3.875315e-12, 3.791408e-08, 9.815712e-01, 2.746096e00, 4.803844e-01],
            [4.159520e-12, 3.521107e-08, 9.612178e-01, 2.753369e00, 4.766497e-01],
            [4.432704e-12, 4.066885e-08, 9.507315e-01, 2.756508e00, 4.734972e-01],
            [4.812173e-12, 4.231054e-08, 9.776751e-01, 2.743801e00, 4.725490e-01],
            [4.872406e-12, 3.413786e-08, 9.895544e-01, 2.739614e00, 4.737730e-01],
        ]
    )

    # Run the parameter sweep
    global_results = run_parameter_sweep(None, seed=1)

    # Compare individual values for specificity
    for value, truth_value in zip(global_results.flatten(), truth_values.flatten()):
        assert value == pytest.approx(truth_value)


@pytest.mark.integration
def test_monte_carlo_sampling_with_files():

    # Define truth values (defined with rng_seed=1)
    truth_values = np.array(
        [
            [2.849231e-12, 2.950054e-08, 9.800058e-01, 2.754121e00, 4.949683e-01],
            [3.463516e-12, 3.307973e-08, 9.833850e-01, 2.747771e00, 4.854287e-01],
            [3.619397e-12, 3.061071e-08, 9.799266e-01, 2.748404e00, 4.835434e-01],
            [3.694122e-12, 2.469930e-08, 9.626206e-01, 2.755873e00, 4.821188e-01],
            [3.735914e-12, 3.338791e-08, 9.774600e-01, 2.748674e00, 4.818635e-01],
            [3.875315e-12, 3.791408e-08, 9.815712e-01, 2.746096e00, 4.803844e-01],
            [4.159520e-12, 3.521107e-08, 9.612178e-01, 2.753369e00, 4.766497e-01],
            [4.432704e-12, 4.066885e-08, 9.507315e-01, 2.756508e00, 4.734972e-01],
            [4.812173e-12, 4.231054e-08, 9.776751e-01, 2.743801e00, 4.725490e-01],
            [4.872406e-12, 3.413786e-08, 9.895544e-01, 2.739614e00, 4.737730e-01],
        ]
    )

    # Define the file of a
    cwd = os.path.dirname(__file__)
    sweep_params_fpath = os.path.join(cwd, "..", "mc_sweep_params.yaml")
    default_config_fpath = os.path.join(cwd, "..", "default_configuration.yaml")

    print("sweep_params_fpath = ", sweep_params_fpath)
    print("default_config_fpath = ", default_config_fpath)

    # Run the parameter sweep
    global_results = run_parameter_sweep(
        seed=1,
        read_sweep_params_from_file=True,
        sweep_params_fname=sweep_params_fpath,
        read_model_defauls_from_file=True,
        defaults_fname=default_config_fpath,
    )

    # Compare individual values for specificity
    for value, truth_value in zip(global_results.flatten(), truth_values.flatten()):
        assert value == pytest.approx(truth_value)


@pytest.mark.integration
def test_lhs_sampling():

    # Define truth values (defined with rng_seed=1)
    truth_values = np.array(
        [
            [
                6.61414143e-13,
                7.40569553e-09,
                9.53202978e-01,
                2.90551565e00,
                7.13296099e-01,
            ],
            [
                1.06113138e-12,
                2.83866433e-08,
                9.74340177e-01,
                2.80854292e00,
                5.92574507e-01,
            ],
            [
                1.42905375e-12,
                2.99443415e-08,
                9.89512570e-01,
                2.77545192e00,
                5.50904323e-01,
            ],
            [
                1.83480893e-12,
                3.37424165e-08,
                9.78156219e-01,
                2.77045843e00,
                5.25076864e-01,
            ],
            [
                2.21768149e-12,
                1.79885975e-08,
                9.82679322e-01,
                2.76125308e00,
                5.12188143e-01,
            ],
            [
                2.54863544e-12,
                1.05788392e-08,
                9.69505557e-01,
                2.76404284e00,
                5.03126360e-01,
            ],
            [
                2.84418360e-12,
                2.29178372e-08,
                9.59253697e-01,
                2.76514112e00,
                4.94806891e-01,
            ],
            [
                3.29995997e-12,
                1.57480856e-08,
                9.64769290e-01,
                2.75818442e00,
                4.88041144e-01,
            ],
            [
                3.39788643e-12,
                2.33656477e-08,
                9.57873046e-01,
                2.76040149e00,
                4.85649494e-01,
            ],
            [
                3.85404230e-12,
                1.23239305e-08,
                9.73578427e-01,
                2.75058745e00,
                4.82138899e-01,
            ],
        ]
    )

    # Run the parameter sweep
    global_results = run_parameter_sweep(seed=1, use_LHS=True)

    # Compare individual values for specificity
    for value, truth_value in zip(global_results.flatten(), truth_values.flatten()):
        assert value == pytest.approx(truth_value)
