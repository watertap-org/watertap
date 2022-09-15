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
import ast

from pyomo.environ import value

from watertap.tools.parameter_sweep.sampling_types import *
from watertap.tools.parameter_sweep import (
    RecursiveParameterSweep,
    recursive_parameter_sweep,
)
from watertap.tools.parameter_sweep.parameter_sweep_writer import *
from watertap.tools.parameter_sweep.tests.test_parameter_sweep import (
    _read_output_h5,
    _get_rank0_path,
    _assert_dictionary_correctness,
)

import watertap.tools.MPI as MPI

# -----------------------------------------------------------------------------


@pytest.fixture
def model():

    """
    # Example Usage:
    # Set the number of trials
    nn = 50
    results = np.zeros((nn, 2))

    for k in range(nn):
        # Set a random value for the parameter a
        m.fs.a.set_value(np.random.fs.rand())

        # Attempt to solve the model (infeasible if m.fs.a > success_prob)
        solver.solve(m)

        # Store the value of a and the solution x
        a = pyo.value(m.fs.a)
        x = pyo.value(m.fs.x)

        if np.abs(pyo.value(m.fs.err)) < 1e-6:
            results[k, :] = [a, x]
        else:
            results[k, :] = [a, np.nan]
    """

    # Initialize a model and solver
    m = pyo.ConcreteModel()

    m.fs = pyo.Block()

    # Declare decision variable and param
    m.fs.x = pyo.Var()
    m.fs.a = pyo.Param(mutable=True)
    m.fs.success_prob = pyo.Param(initialize=0.5, mutable=True)

    # Define expressions and constraints:
    # Numbers must sum to success_prob and x must be positive
    m.fs.err = pyo.Expression(expr=m.fs.a + m.fs.x - m.fs.success_prob)
    m.fs.sum = pyo.Constraint(expr=m.fs.err == 0.0)
    m.fs.pos = pyo.Constraint(expr=m.fs.x >= 0.0)

    return m


@pytest.mark.component
def test_aggregate_filtered_input_arr():

    ps = RecursiveParameterSweep()

    input_dict = {
        "outputs": {
            "x_val": {
                "lower bound": None,
                "units": "None",
                "upper bound": None,
                "value": np.array(
                    [
                        0.31655848,
                        0.2763452,
                        0.26241279,
                        0.15511682,
                        0.1511865,
                        0.09723662,
                        0.05410589,
                        0.6797816,
                        0.62896394,
                        0.6128707,
                        0.23852064,
                        0.17110508,
                        0.13195544,
                    ]
                ),
            }
        },
        "solve_successful": [True] * 13,
        "sweep_params": {
            "fs.a": {
                "units": "None",
                "value": np.array(
                    [
                        0.38344152,
                        0.4236548,
                        0.43758721,
                        0.54488318,
                        0.5488135,
                        0.60276338,
                        0.64589411,
                        0.0202184,
                        0.07103606,
                        0.0871293,
                        0.46147936,
                        0.52889492,
                        0.56804456,
                    ]
                ),
            }
        },
    }

    truth_arr = np.array(
        [
            0.38344152,
            0.4236548,
            0.43758721,
            0.54488318,
            0.5488135,
            0.60276338,
            0.64589411,
            0.0202184,
            0.07103606,
            0.0871293,
            0.46147936,
            0.52889492,
            0.56804456,
        ]
    )

    req_num_samples = 10
    values_arr = ps._aggregate_filtered_input_arr(input_dict, req_num_samples)

    assert np.allclose(values_arr[:, 0], truth_arr[0:req_num_samples], equal_nan=True)


@pytest.mark.component
def test_recursive_parameter_sweep(model, tmp_path):

    comm = MPI.COMM_WORLD

    tmp_path = _get_rank0_path(comm, tmp_path)

    results_fname = os.path.join(tmp_path, "global_results")
    csv_results_file = str(results_fname) + ".csv"
    h5_results_file = str(results_fname) + ".h5"

    ps = RecursiveParameterSweep(
        csv_results_file_name=csv_results_file,
        h5_results_file_name=h5_results_file,
        debugging_data_dir=tmp_path,
    )

    m = model

    solver = pyo.SolverFactory("ipopt")

    sweep_params = {}
    sweep_params["a_val"] = UniformSample(m.fs.a, 0.0, 1.0)

    outputs = {}
    outputs["x_val"] = m.fs.x

    # Run the parameter sweep study using num_samples randomly drawn from the above range
    num_samples = 10
    seed = 0

    # Run the parameter sweep
    data = ps.parameter_sweep(
        m,
        sweep_params,
        outputs=outputs,
        req_num_samples=num_samples,
        seed=seed,
    )

    reference_save_data = np.array(
        [
            [0.38344152, 0.11655848],
            [0.4236548, 0.0763452],
            [0.43758721, 0.06241279],
            [0.0187898, 0.4812102],
            [0.0202184, 0.4797816],
            [0.06022547, 0.43977453],
            [0.07103606, 0.42896394],
            [0.0871293, 0.4128707],
            [0.10204481, 0.39795519],
            [0.11827443, 0.38172557],
        ]
    )

    assert np.shape(data) == (10, 2)
    assert np.allclose(reference_save_data, data, equal_nan=True)
    assert np.allclose(np.sum(data, axis=1), value(m.fs.success_prob))

    if ps.rank == 0:
        # Check that the global results file is created
        assert os.path.isfile(csv_results_file)

        # Check that all local output files have been created
        for k in range(ps.num_procs):
            assert os.path.isfile(os.path.join(tmp_path, f"local_results_{k:03}.h5"))
            assert os.path.isfile(os.path.join(tmp_path, f"local_results_{k:03}.csv"))

        csv_data = np.genfromtxt(csv_results_file, skip_header=1, delimiter=",")

        # Compare the last row of the imported data to truth
        assert np.allclose(data[-1, :], reference_save_data[-1, :], equal_nan=True)

        truth_dict = {
            "outputs": {
                "x_val": {
                    "lower bound": -1.7976931348623157e308,
                    "units": "None",
                    "upper bound": 1.7976931348623157e308,
                    "value": np.array(
                        [
                            0.11655848,
                            0.0763452,
                            0.06241279,
                            0.4812102,
                            0.4797816,
                            0.43977453,
                            0.42896394,
                            0.4128707,
                            0.39795519,
                            0.38172557,
                        ]
                    ),
                }
            },
            "solve_successful": [
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
            ],
            "sweep_params": {
                "fs.a": {
                    "units": "None",
                    "value": np.array(
                        [
                            0.38344152,
                            0.4236548,
                            0.43758721,
                            0.0187898,
                            0.0202184,
                            0.06022547,
                            0.07103606,
                            0.0871293,
                            0.10204481,
                            0.11827443,
                        ]
                    ),
                }
            },
        }

        read_dict = _read_output_h5(h5_results_file)
        _assert_dictionary_correctness(truth_dict, read_dict)

        assert np.allclose(
            data[:, -1], read_dict["outputs"]["x_val"]["value"][:num_samples]
        )

        # Check for the companion text file
        txt_fpath = os.path.join(tmp_path, "{0}.txt".format(h5_results_file))
        assert os.path.exists(txt_fpath)

        truth_txt_dict = {
            "outputs": ["x_val"],
            "sweep_params": ["fs.a"],
        }

        with open(txt_fpath, "r") as f:
            f_contents = f.read()
            read_txt_dict = ast.literal_eval(f_contents)
        assert read_txt_dict == truth_txt_dict


@pytest.mark.component
def test_recursive_parameter_sweep_function(model, tmp_path):
    comm = MPI.COMM_WORLD
    tmp_path = _get_rank0_path(comm, tmp_path)

    m = model

    solver = pyo.SolverFactory("ipopt")

    sweep_params = {}
    sweep_params["a_val"] = UniformSample(m.fs.a, 0.0, 1.0)

    outputs = {}
    outputs["x_val"] = m.fs.x

    # Run the parameter sweep study using num_samples randomly drawn from the above range
    num_samples = 10
    seed = 0

    results_fname = os.path.join(tmp_path, "global_results")
    csv_results_file = str(results_fname) + ".csv"
    h5_results_file = str(results_fname) + ".h5"

    # Run the parameter sweep
    data = recursive_parameter_sweep(
        m,
        sweep_params,
        outputs=outputs,
        csv_results_file_name=csv_results_file,
        h5_results_file_name=h5_results_file,
        req_num_samples=num_samples,
        debugging_data_dir=tmp_path,
        seed=seed,
    )

    reference_save_data = np.array(
        [
            [0.38344152, 0.11655848],
            [0.4236548, 0.0763452],
            [0.43758721, 0.06241279],
            [0.0187898, 0.4812102],
            [0.0202184, 0.4797816],
            [0.06022547, 0.43977453],
            [0.07103606, 0.42896394],
            [0.0871293, 0.4128707],
            [0.10204481, 0.39795519],
            [0.11827443, 0.38172557],
        ]
    )

    assert np.shape(data) == (10, 2)
    assert np.allclose(reference_save_data, data, equal_nan=True)
    assert np.allclose(np.sum(data, axis=1), value(m.fs.success_prob))

    if comm.rank == 0:
        # Check that the global results file is created
        assert os.path.isfile(csv_results_file)

        # Check that all local output files have been created
        for k in range(comm.size):
            assert os.path.isfile(os.path.join(tmp_path, f"local_results_{k:03}.h5"))
            assert os.path.isfile(os.path.join(tmp_path, f"local_results_{k:03}.csv"))

        csv_data = np.genfromtxt(csv_results_file, skip_header=1, delimiter=",")

        # Compare the last row of the imported data to truth
        assert np.allclose(data[-1, :], reference_save_data[-1, :], equal_nan=True)

        # Check for the h5 output
        truth_dict = {
            "outputs": {
                "x_val": {
                    "lower bound": -1.7976931348623157e308,
                    "units": "None",
                    "upper bound": 1.7976931348623157e308,
                    "value": np.array(
                        [
                            0.11655848,
                            0.0763452,
                            0.06241279,
                            0.4812102,
                            0.4797816,
                            0.43977453,
                            0.42896394,
                            0.4128707,
                            0.39795519,
                            0.38172557,
                        ]
                    ),
                }
            },
            "solve_successful": [
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
            ],
            "sweep_params": {
                "fs.a": {
                    "units": "None",
                    "value": np.array(
                        [
                            0.38344152,
                            0.4236548,
                            0.43758721,
                            0.0187898,
                            0.0202184,
                            0.06022547,
                            0.07103606,
                            0.0871293,
                            0.10204481,
                            0.11827443,
                        ]
                    ),
                }
            },
        }

        read_dict = _read_output_h5(h5_results_file)
        _assert_dictionary_correctness(truth_dict, read_dict)
        assert np.allclose(
            data[:, -1], read_dict["outputs"]["x_val"]["value"][:num_samples]
        )

        # Check for the companion text file
        txt_fpath = os.path.join(tmp_path, "{0}.txt".format(h5_results_file))
        assert os.path.exists(txt_fpath)

        truth_txt_dict = {
            "outputs": ["x_val"],
            "sweep_params": ["fs.a"],
        }

        with open(txt_fpath, "r") as f:
            f_contents = f.read()
            read_txt_dict = ast.literal_eval(f_contents)
        assert read_txt_dict == truth_txt_dict
