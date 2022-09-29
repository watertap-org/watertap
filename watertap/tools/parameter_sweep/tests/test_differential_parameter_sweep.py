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
from watertap.tools.parameter_sweep import DifferentialParameterSweep
from watertap.tools.parameter_sweep.parameter_sweep_writer import *
from watertap.tools.parameter_sweep.tests.test_parameter_sweep import (
    _read_output_h5,
    _get_rank0_path,
    _assert_dictionary_correctness,
    _assert_h5_csv_agreement,
    _optimization,
    _reinitialize,
)
import watertap.tools.MPI as MPI


@pytest.fixture
def model():
    m = pyo.ConcreteModel()
    m.fs = fs = pyo.Block()

    fs.input = pyo.Var(["a", "b"], within=pyo.UnitInterval, initialize=0.5)
    fs.output = pyo.Var(["c", "d"], within=pyo.UnitInterval, initialize=0.5)

    fs.slack = pyo.Var(["ab_slack", "cd_slack"], bounds=(0, 0), initialize=0.0)
    fs.slack_penalty = pyo.Param(default=1000.0, mutable=True, within=pyo.PositiveReals)

    fs.ab_constr = pyo.Constraint(
        expr=(fs.output["c"] + fs.slack["ab_slack"] == 2 * fs.input["a"])
    )
    fs.cd_constr = pyo.Constraint(
        expr=(fs.output["d"] + fs.slack["cd_slack"] == 3 * fs.input["b"])
    )

    fs.performance = pyo.Expression(expr=pyo.summation(fs.output))

    m.objective = pyo.Objective(
        expr=m.fs.performance - m.fs.slack_penalty * pyo.summation(m.fs.slack),
        sense=pyo.maximize,
    )
    return m


@pytest.mark.component
def test_create_differential_sweep_params_normal(model):

    m = model

    differential_sweep_specs = {
        "fs.a": {
            "diff_sample_type": NormalSample,
            "relative_std_dev": 0.01,
            "pyomo_object": m.fs.input["a"],
        },
        "fs.b": {
            "diff_sample_type": NormalSample,
            "relative_std_dev": 0.5,
            "pyomo_object": m.fs.input["b"],
        },
    }

    ps = DifferentialParameterSweep()
    local_values = np.array([0.0, 1.0, 2.0])

    diff_sweep_param_dict = ps._create_differential_sweep_params(
        local_values, differential_sweep_specs
    )

    expected_dict = {
        "fs.a": NormalSample(m.fs.input["a"], 0.0, 0.01),
        "fs.b": NormalSample(m.fs.input["b"], 1.0, 0.5),
    }

    for key, value in diff_sweep_param_dict.items():
        print(value.mean)
        assert value.mean == expected_dict[key].mean
        assert value.sd == expected_dict[key].sd


@pytest.mark.component
def test_create_differential_sweep_params_others(model):

    m = model

    differential_sweep_specs = {
        "fs.a": {
            "diff_mode": "sum",
            "diff_sample_type": GeomSample,
            "relative_lb": 0.01,
            "relative_ub": 10.0,
            "pyomo_object": m.fs.input["a"],
        },
        "fs.b": {
            "diff_mode": "product",
            "diff_sample_type": UniformSample,
            "relative_lb": 0.01,
            "relative_ub": 0.1,
            "pyomo_object": m.fs.input["b"],
        },
    }

    ps = DifferentialParameterSweep()
    local_values = np.array([0.1, 1.0, 2.0])

    diff_sweep_param_dict = ps._create_differential_sweep_params(
        local_values, differential_sweep_specs
    )

    expected_dict = {
        "fs.a": GeomSample(m.fs.input["a"], 0.099, 0.101, 1),
        "fs.b": UniformSample(m.fs.input["b"], 0.01, 0.1),
    }

    for key, value in diff_sweep_param_dict.items():
        assert value.lower_limit == expected_dict[key].lower_limit
        assert value.upper_limit == expected_dict[key].upper_limit


@pytest.mark.component
def test_differential_parameter_sweep(model, tmp_path):

    comm = MPI.COMM_WORLD
    tmp_path = _get_rank0_path(comm, tmp_path)

    results_fname = os.path.join(tmp_path, "global_results")
    csv_results_file_name = str(results_fname) + ".csv"
    h5_results_file_name = str(results_fname) + ".h5"

    m = model

    differential_sweep_specs = {
        "fs.a": {
            "diff_mode": "sum",
            "diff_sample_type": NormalSample,
            "relative_std_dev": 0.01,
            "pyomo_object": m.fs.input["a"],
        },
        "fs.b": {
            "diff_mode": "product",
            "diff_sample_type": UniformSample,
            "relative_lb": 0.01,
            "relative_ub": 0.01,
            "pyomo_object": m.fs.input["b"],
        },
    }

    ps = DifferentialParameterSweep(
        comm=comm,
        csv_results_file_name=csv_results_file_name,
        h5_results_file_name=h5_results_file_name,
        debugging_data_dir=tmp_path,
        interpolate_nan_outputs=True,
        optimize_function=_optimization,
        reinitialize_function=_reinitialize,
        reinitialize_kwargs={"slack_penalty": 10.0},
    )

    m = model
    m.fs.slack_penalty = 1000.0
    m.fs.slack.setub(0)

    A = m.fs.input["a"]
    B = m.fs.input["b"]
    sweep_params = {A.name: (A, 0.1, 0.9, 3), B.name: (B, 0.0, 0.5, 3)}

    # Call the parameter_sweep function
    global_results_dict, _ = ps.parameter_sweep(
        m,
        sweep_params,
        differential_sweep_specs,
        outputs=None,
        # optimize_function=_optimization,
        # reinitialize_function=_reinitialize,
        # reinitialize_kwargs={"slack_penalty": 10.0},
        seed=0,
    )

    import pprint

    print()
    pprint.pprint(global_results_dict)

    if ps.rank == 0:
        truth_dict = {
            "outputs": {
                "fs.output[c]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array(
                        [
                            0.2,
                            0.2,
                            0.2,
                            1.0,
                            1.0,
                            1.0,
                            1.0,
                            1.0,
                            1.0,
                            0.23528105,
                            0.23528105,
                            0.23528106,
                            1.0,
                            1.0,
                            1.0,
                            1.0,
                            1.0,
                            1.0,
                        ]
                    ),
                },
                "fs.output[d]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array(
                        [
                            0.0,
                            0.75,
                            1.0,
                            9.98996974e-09,
                            0.75,
                            1.0,
                            9.92236517e-09,
                            0.75,
                            1.0,
                            4.23516474e-22,
                            7.50000000e-03,
                            1.50000098e-02,
                            9.77884773e-09,
                            7.50000977e-03,
                            1.50000098e-02,
                            9.77884942e-09,
                            7.50000977e-03,
                            1.50000098e-02,
                        ]
                    ),
                },
                "fs.performance": {
                    "value": np.array(
                        [
                            0.2,
                            0.95,
                            1.2,
                            1.0,
                            1.75,
                            2.0,
                            1.0,
                            1.75,
                            2.0,
                            0.23528105,
                            0.24278105,
                            0.25028107,
                            1.0,
                            1.00750001,
                            1.01500001,
                            1.0,
                            1.00750001,
                            1.01500001,
                        ]
                    )
                },
                "fs.slack[ab_slack]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 0,
                    "value": np.array(
                        [
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.8,
                            0.8,
                            0.8,
                            0.0,
                            0.0,
                            0.0,
                            0.03528104,
                            0.03528104,
                            0.03528104,
                            0.83528104,
                            0.83528104,
                            0.83528104,
                        ]
                    ),
                },
                "fs.slack[cd_slack]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 0,
                    "value": np.array(
                        [
                            0.0,
                            0.0,
                            0.5,
                            0.0,
                            0.0,
                            0.5,
                            0.0,
                            0.0,
                            0.5,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                        ]
                    ),
                },
                "objective": {
                    "value": np.array(
                        [
                            0.2,
                            0.95,
                            -3.79999989,
                            1.0,
                            1.75,
                            -2.9999999,
                            -6.99999989,
                            -6.24999989,
                            -10.9999998,
                            0.23528105,
                            0.24278105,
                            0.25028107,
                            0.64718964,
                            0.65468964,
                            0.66218964,
                            -7.35281036,
                            -7.34531036,
                            -7.33781036,
                        ]
                    )
                },
            },
            "sweep_params": {
                "fs.input[a]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array(
                        [
                            0.1,
                            0.1,
                            0.1,
                            0.5,
                            0.5,
                            0.5,
                            0.9,
                            0.9,
                            0.9,
                            0.11764052,
                            0.11764052,
                            0.11764052,
                            0.51764052,
                            0.51764052,
                            0.51764052,
                            0.91764052,
                            0.91764052,
                            0.91764052,
                        ]
                    ),
                },
                "fs.input[b]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array(
                        [
                            0.0,
                            0.25,
                            0.5,
                            0.0,
                            0.25,
                            0.5,
                            0.0,
                            0.25,
                            0.5,
                            0.0,
                            0.0025,
                            0.005,
                            0.0,
                            0.0025,
                            0.005,
                            0.0,
                            0.0025,
                            0.005,
                        ]
                    ),
                },
            },
        }

        read_dict = _read_output_h5(h5_results_file_name)
        _assert_dictionary_correctness(truth_dict, read_dict)
        _assert_h5_csv_agreement(csv_results_file_name, read_dict)
