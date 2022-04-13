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
            [2.849231e-12, 2.950054e-08, 9.800058e-01, 2.735876e00, 4.587524e-01],
            [3.463516e-12, 3.307973e-08, 9.833850e-01, 2.729647e00, 4.484156e-01],
            [3.619397e-12, 3.061071e-08, 9.799266e-01, 2.729649e00, 4.463418e-01],
            [3.694122e-12, 2.469930e-08, 9.626206e-01, 2.734365e00, 4.451184e-01],
            [3.735914e-12, 3.338791e-08, 9.774600e-01, 2.729575e00, 4.447564e-01],
            [3.875315e-12, 3.791408e-08, 9.815712e-01, 2.727602e00, 4.431959e-01],
            [4.159520e-12, 3.521107e-08, 9.612178e-01, 2.731659e00, 4.397268e-01],
            [4.432704e-12, 4.066885e-08, 9.507315e-01, 2.733099e00, 4.368726e-01],
            [4.812173e-12, 4.231054e-08, 9.776751e-01, 2.724421e00, 4.348367e-01],
            [4.872406e-12, 3.413786e-08, 9.895544e-01, 2.721691e00, 4.352008e-01],
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
            [2.849231e-12, 2.950054e-08, 9.800058e-01, 2.735876e00, 4.587524e-01],
            [3.463516e-12, 3.307973e-08, 9.833850e-01, 2.729647e00, 4.484156e-01],
            [3.619397e-12, 3.061071e-08, 9.799266e-01, 2.729649e00, 4.463418e-01],
            [3.694122e-12, 2.469930e-08, 9.626206e-01, 2.734365e00, 4.451184e-01],
            [3.735914e-12, 3.338791e-08, 9.774600e-01, 2.729575e00, 4.447564e-01],
            [3.875315e-12, 3.791408e-08, 9.815712e-01, 2.727602e00, 4.431959e-01],
            [4.159520e-12, 3.521107e-08, 9.612178e-01, 2.731659e00, 4.397268e-01],
            [4.432704e-12, 4.066885e-08, 9.507315e-01, 2.733099e00, 4.368726e-01],
            [4.812173e-12, 4.231054e-08, 9.776751e-01, 2.724421e00, 4.348367e-01],
            [4.872406e-12, 3.413786e-08, 9.895544e-01, 2.721691e00, 4.352008e-01],
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
                2.88439747e00,
                6.79583725e-01,
            ],
            [
                1.06113138e-12,
                2.83866433e-08,
                9.74340177e-01,
                2.79162022e00,
                5.60514329e-01,
            ],
            [
                1.42905375e-12,
                2.99443415e-08,
                9.89512570e-01,
                2.75991573e00,
                5.17460908e-01,
            ],
            [
                1.83480893e-12,
                3.37424165e-08,
                9.78156219e-01,
                2.75310862e00,
                4.91161174e-01,
            ],
            [
                2.21768149e-12,
                1.79885975e-08,
                9.82679322e-01,
                2.74356449e00,
                4.76269356e-01,
            ],
            [
                2.54863544e-12,
                1.05788392e-08,
                9.69505557e-01,
                2.74396258e00,
                4.66585039e-01,
            ],
            [
                2.84418360e-12,
                2.29178372e-08,
                9.59253697e-01,
                2.74369493e00,
                4.58896025e-01,
            ],
            [
                3.29995997e-12,
                1.57480856e-08,
                9.64769290e-01,
                2.73699331e00,
                4.50942843e-01,
            ],
            [
                3.39788643e-12,
                2.33656477e-08,
                9.57873046e-01,
                2.73830993e00,
                4.49066094e-01,
            ],
            [
                3.85404230e-12,
                1.23239305e-08,
                9.73578427e-01,
                2.73028618e00,
                4.43937384e-01,
            ],
        ]
    )

    # Run the parameter sweep
    global_results = run_parameter_sweep(seed=1, use_LHS=True)

    # Compare individual values for specificity
    for value, truth_value in zip(global_results.flatten(), truth_values.flatten()):
        assert value == pytest.approx(truth_value)
