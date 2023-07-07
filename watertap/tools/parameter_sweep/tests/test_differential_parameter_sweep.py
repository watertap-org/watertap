#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pytest
import os
import copy
import numpy as np
import pyomo.environ as pyo

from watertap.tools.parameter_sweep.sampling_types import (
    NormalSample,
    GeomSample,
    UniformSample,
)
from watertap.tools.parameter_sweep import (
    DifferentialParameterSweep,
    differential_parameter_sweep,
)
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
def test_check_differential_sweep_key_validity(model):

    m = model

    A = m.fs.input["a"]
    B = m.fs.input["b"]
    sweep_params = {A.name: (A, 0.1, 0.9, 3), B.name: (B, 0.0, 0.5, 3)}

    differential_sweep_specs = {
        A.name: {
            "diff_mode": "sum",
            "diff_sample_type": NormalSample,
            "std_dev": 0.01,
            "pyomo_object": m.fs.input["a"],
        },
        B.name: {
            "diff_mode": "product",
            "diff_sample_type": UniformSample,
            "relative_lb": 0.01,
            "relative_ub": 0.01,
            "pyomo_object": m.fs.input["b"],
        },
    }

    ps = DifferentialParameterSweep(differential_sweep_specs=differential_sweep_specs)
    sweep_params, _ = ps._process_sweep_params(sweep_params)
    ps.outputs = None
    ps._check_differential_sweep_key_validity(sweep_params)

    assert ps.diff_spec_index == [0, 1]


@pytest.mark.component
def test_create_differential_sweep_params_normal(model):

    m = model

    differential_sweep_specs = {
        "fs.a": {
            "diff_sample_type": NormalSample,
            "std_dev": 0.01,
            "pyomo_object": m.fs.input["a"],
        },
        "fs.b": {
            "diff_sample_type": NormalSample,
            "std_dev": 0.5,
            "pyomo_object": m.fs.input["b"],
        },
    }

    ps = DifferentialParameterSweep(differential_sweep_specs=differential_sweep_specs)
    local_values = np.array([0.0, 1.0, 2.0])

    ps.diff_spec_index = [0, 1]
    diff_sweep_param_dict = ps._create_differential_sweep_params(local_values)

    expected_dict = {
        "fs.a": NormalSample(m.fs.input["a"], 0.0, 0.01, 1),
        "fs.b": NormalSample(m.fs.input["b"], 1.0, 0.5, 1),
    }

    for key, value in diff_sweep_param_dict.items():
        assert value.mean == expected_dict[key].mean
        assert value.sd == expected_dict[key].sd


@pytest.mark.component
def test_create_differential_sweep_params_sum_prod(model):

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

    ps = DifferentialParameterSweep(differential_sweep_specs=differential_sweep_specs)
    local_values = np.array([0.1, 1.0, 2.0])

    ps.diff_spec_index = [0, 1]
    diff_sweep_param_dict = ps._create_differential_sweep_params(local_values)

    expected_dict = {
        "fs.a": GeomSample(m.fs.input["a"], 0.099, 1.1, 1),
        "fs.b": UniformSample(m.fs.input["b"], 0.01, 0.1, 1),
    }

    for key, value in diff_sweep_param_dict.items():
        assert value.lower_limit == expected_dict[key].lower_limit
        assert value.upper_limit == expected_dict[key].upper_limit


@pytest.mark.component
def test_create_differential_sweep_params_percentile(model):

    m = model

    differential_sweep_specs = {
        "fs.b": {
            "diff_mode": "percentile",
            "diff_sample_type": UniformSample,
            "relative_lb": 0.01,
            "relative_ub": 0.1,
            "nominal_lb": 0.0,
            "nominal_ub": 1.0,
            "pyomo_object": m.fs.input["b"],
        },
    }

    ps = DifferentialParameterSweep(differential_sweep_specs=differential_sweep_specs)
    local_values = np.array([0.1, 1.0, 2.0])

    ps.diff_spec_index = [0, 1]
    diff_sweep_param_dict = ps._create_differential_sweep_params(local_values)

    expected_dict = {
        "fs.b": UniformSample(m.fs.input["b"], 0.11, 0.2, 1),
    }

    for key, value in diff_sweep_param_dict.items():
        assert value.lower_limit == expected_dict[key].lower_limit
        assert value.upper_limit == expected_dict[key].upper_limit


@pytest.mark.component
def test_bad_differential_sweep_specs(model, tmp_path):

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

    A = m.fs.input["a"]
    B = m.fs.input["b"]
    sweep_params = {A.name: (A, 0.1, 0.9, 3), B.name: (B, 0.0, 0.5, 3)}

    ps = DifferentialParameterSweep(differential_sweep_specs=differential_sweep_specs)
    with pytest.raises(ValueError):
        ps.parameter_sweep(
            m,
            sweep_params,
            outputs=None,
            seed=0,
        )


@pytest.mark.component
def test_differential_sweep_outputs(model):

    comm = MPI.COMM_WORLD

    m = model
    m.fs.slack_penalty = 1000.0
    m.fs.slack.setub(0)

    A = m.fs.input["a"]
    B = m.fs.input["b"]
    sweep_params = {A.name: (A, 0.1, 0.9, 3), B.name: (B, 0.0, 0.5, 3)}

    differential_sweep_specs = {
        B.name: {
            "diff_mode": "product",
            "diff_sample_type": UniformSample,
            "relative_lb": 0.01,
            "relative_ub": 0.01,
            "pyomo_object": m.fs.input["b"],
        },
    }

    outputs = {"fs.output[c]": m.fs.output["c"]}

    ps = DifferentialParameterSweep(
        comm=comm,
        optimize_function=_optimization,
        reinitialize_function=_reinitialize,
        reinitialize_kwargs={"slack_penalty": 10.0},
        differential_sweep_specs=differential_sweep_specs,
    )

    sweep_params, _ = ps._process_sweep_params(sweep_params)

    ps.outputs = outputs
    ps._define_differential_sweep_outputs(sweep_params)

    # Finally test for the keys
    expected_keys = ["fs.output[c]", "fs.input[a]"]
    assert expected_keys == list(ps.differential_outputs.keys())


@pytest.mark.component
def test_differential_parameter_sweep(model, tmp_path):

    comm = MPI.COMM_WORLD
    tmp_path = "./output/"  # _get_rank0_path(comm, tmp_path)

    results_fname = os.path.join(tmp_path, "global_results")
    csv_results_file_name = str(results_fname) + ".csv"
    h5_results_file_name = str(results_fname) + ".h5"

    m = model
    m.fs.slack_penalty = 1000.0
    m.fs.slack.setub(0)

    A = m.fs.input["a"]
    B = m.fs.input["b"]
    sweep_params = {A.name: (A, 0.1, 0.9, 3), B.name: (B, 0.0, 0.5, 3)}

    differential_sweep_specs = {
        A.name: {
            "diff_mode": "sum",
            "diff_sample_type": NormalSample,
            "std_dev": 0.01,
            "pyomo_object": m.fs.input["a"],
        },
        B.name: {
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
        differential_sweep_specs=differential_sweep_specs,
    )

    # Call the parameter_sweep function
    global_results_dict, _ = ps.parameter_sweep(
        m,
        sweep_params,
        outputs=None,
        seed=0,
    )

    if ps.rank == 0:

        truth_dict = {
            "differential_idx": np.array(
                [
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    0.0,
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                    5.0,
                    6.0,
                    7.0,
                    8.0,
                ]
            ),
            "nominal_idx": np.array(
                [
                    0.0,
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                    5.0,
                    6.0,
                    7.0,
                    8.0,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                ]
            ),
            "outputs": {
                "fs.input[a]": {
                    "full_name": "fs.input[a]",
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
                    "full_name": "fs.input[b]",
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
                "fs.output[c]": {
                    "full_name": "fs.output[c]",
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array(
                        [
                            0.2,
                            0.2,
                            0.20000001,
                            1.00000001,
                            1.00000001,
                            1.00000001,
                            1.00000001,
                            1.00000001,
                            1.00000001,
                            0.23528105,
                            0.23528105,
                            0.23528106,
                            1.00000001,
                            1.00000001,
                            1.00000001,
                            1.00000001,
                            1.00000001,
                            1.00000001,
                        ]
                    ),
                },
                "fs.output[d]": {
                    "full_name": "fs.output[d]",
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array(
                        [
                            0.00000000e00,
                            7.50000000e-01,
                            1.00000001e00,
                            9.98996974e-09,
                            7.50000010e-01,
                            1.00000001e00,
                            9.92236517e-09,
                            7.50000010e-01,
                            1.00000001e00,
                            1.09912079e-22,
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
                    "full_name": "fs.performance",
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
                            1.0075,
                            1.0150,
                            1.0,
                            1.0075,
                            1.015,
                        ]
                    ),
                },
                "fs.slack[ab_slack]": {
                    "full_name": "fs.slack[ab_slack]",
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 0,
                    "value": np.array(
                        [
                            0.0,
                            0.0,
                            -9.77208577e-09,
                            -9.54438112e-09,
                            -9.54432508e-09,
                            -9.54432509e-09,
                            0.8,
                            0.8,
                            0.8,
                            0.0,
                            0.0,
                            -9.77210838e-09,
                            3.52810371e-02,
                            3.52810371e-02,
                            3.52810371e-02,
                            8.35281037e-01,
                            8.35281037e-01,
                            8.35281037e-01,
                        ]
                    ),
                },
                "fs.slack[cd_slack]": {
                    "full_name": "fs.slack[cd_slack]",
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 0,
                    "value": np.array(
                        [
                            0.0,
                            0.0,
                            0.5,
                            -9.98996974e-09,
                            -9.77226571e-09,
                            0.5,
                            -9.92236517e-09,
                            -9.77226571e-09,
                            0.5,
                            0.0,
                            0.0,
                            -9.77036119e-09,
                            -9.77884773e-09,
                            -9.76848697e-09,
                            -9.77036119e-09,
                            -9.77884942e-09,
                            -9.76848697e-09,
                            -9.77036119e-09,
                        ]
                    ),
                },
                "fs.slack_penalty": {
                    "full_name": "fs.slack_penalty",
                    "units": "None",
                    "value": np.array(
                        [
                            1000.0,
                            1000.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            1000.0,
                            1000.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                        ]
                    ),
                },
                "objective": {
                    "full_name": "objective",
                    "value": np.array(
                        [
                            0.2,
                            0.95,
                            -3.79999979,
                            1.00000021,
                            1.75000021,
                            -2.99999979,
                            -6.99999978,
                            -6.24999979,
                            -10.99999979,
                            0.23528105,
                            0.24278105,
                            0.25028126,
                            0.64718975,
                            0.65468975,
                            0.66218975,
                            -7.35281025,
                            -7.34531025,
                            -7.33781025,
                        ]
                    ),
                },
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
                "fs.input[a]": {
                    "full_name": "fs.input[a]",
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
                    "full_name": "fs.input[b]",
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
        _assert_h5_csv_agreement(csv_results_file_name, read_dict)
        if ps.parallel_manager.number_of_processes() > 1:
            # Compare the sorted dictionary. We need to work with a sorted dictionary
            # because the differential parameter sweep produces a global dictionary
            # that is jumbled by the number of procs.
            if ps.parallel_manager.rank == 0:
                sorted_truth_dict = sort_output_dict(truth_dict)
                sorted_read_dict = sort_output_dict(read_dict)
                sorted_global_results_dict = sort_output_dict(global_results_dict)
                _assert_dictionary_correctness(sorted_truth_dict, sorted_read_dict)
                _assert_dictionary_correctness(
                    sorted_truth_dict, sorted_global_results_dict
                )
        else:
            _assert_dictionary_correctness(truth_dict, global_results_dict)
            _assert_dictionary_correctness(truth_dict, read_dict)


@pytest.mark.component
def test_differential_parameter_sweep_selective(model, tmp_path):

    comm = MPI.COMM_WORLD
    tmp_path = _get_rank0_path(comm, tmp_path)

    results_fname = os.path.join(tmp_path, "global_results")
    csv_results_file_name = str(results_fname) + ".csv"
    h5_results_file_name = str(results_fname) + ".h5"

    m = model
    m.fs.slack_penalty = 1000.0
    m.fs.slack.setub(0)

    A = m.fs.input["a"]
    B = m.fs.input["b"]
    sweep_params = {A.name: (A, 0.1, 0.9, 3), B.name: (B, 0.0, 0.5, 3)}

    differential_sweep_specs = {
        B.name: {
            "diff_mode": "product",
            "diff_sample_type": UniformSample,
            "relative_lb": 0.25,
            "relative_ub": 0.75,
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
        differential_sweep_specs=differential_sweep_specs,
        num_diff_samples=2,
    )

    # Call the parameter_sweep function
    global_results_dict, _ = ps.parameter_sweep(
        m,
        sweep_params,
        outputs=None,
        seed=0,
    )

    if ps.rank == 0:

        truth_dict = {
            "differential_idx": np.array(
                [
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    2.0,
                    2.0,
                    3.0,
                    3.0,
                    4.0,
                    4.0,
                    5.0,
                    5.0,
                    6.0,
                    6.0,
                    7.0,
                    7.0,
                    8.0,
                    8.0,
                ]
            ),
            "nominal_idx": np.array(
                [
                    0.0,
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                    5.0,
                    6.0,
                    7.0,
                    8.0,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                ]
            ),
            "outputs": {
                "fs.input[a]": {
                    "full_name": "fs.input[a]",
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
                            0.1,
                            0.1,
                            0.1,
                            0.1,
                            0.1,
                            0.1,
                            0.5,
                            0.5,
                            0.5,
                            0.5,
                            0.5,
                            0.5,
                            0.9,
                            0.9,
                            0.9,
                            0.9,
                            0.9,
                            0.9,
                        ]
                    ),
                },
                "fs.input[b]": {
                    "full_name": "fs.input[b]",
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
                            0.0,
                            0.13110169,
                            0.15189867,
                            0.26220338,
                            0.30379734,
                            0.0,
                            0.0,
                            0.13110169,
                            0.15189867,
                            0.26220338,
                            0.30379734,
                            0.0,
                            0.0,
                            0.13110169,
                            0.15189867,
                            0.26220338,
                            0.30379734,
                        ]
                    ),
                },
                "fs.output[c]": {
                    "full_name": "fs.output[c]",
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
                            0.2,
                            0.2,
                            0.2,
                            0.2,
                            0.2,
                            0.2,
                            1.0,
                            1.0,
                            1.0,
                            1.0,
                            1.0,
                            1.0,
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
                    "full_name": "fs.output[d]",
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array(
                        [
                            0.0,
                            0.75,
                            1.0,
                            9.89148858e-09,
                            7.50000010e-01,
                            1.00000001e00,
                            9.77798881e-09,
                            7.50000010e-01,
                            1.0,
                            0.0,
                            0.0,
                            3.93305064e-01,
                            4.55696012e-01,
                            7.86610138e-01,
                            9.11392035e-01,
                            9.88979157e-09,
                            9.88979157e-09,
                            3.93305074e-01,
                            4.55696022e-01,
                            7.86610138e-01,
                            9.11392035e-01,
                            9.88962016e-09,
                            9.88962016e-09,
                            3.93305074e-01,
                            4.55696022e-01,
                            7.86610138e-01,
                            9.11392035e-01,
                        ]
                    ),
                },
                "fs.performance": {
                    "full_name": "fs.performance",
                    "value": np.array(
                        [
                            0.2,
                            0.95,
                            1.20,
                            1.00,
                            1.75,
                            2.0,
                            1.00,
                            1.75,
                            2.0,
                            0.2,
                            0.2,
                            0.59330506,
                            0.65569601,
                            0.98661015,
                            1.11139204,
                            1.00000002,
                            1.00000002,
                            1.39330508,
                            1.45569603,
                            1.78661015,
                            1.91139204,
                            1.0,
                            1.0,
                            1.39330508,
                            1.45569603,
                            1.78661015,
                            1.91139204,
                        ]
                    ),
                },
                "fs.slack[ab_slack]": {
                    "full_name": "fs.slack[ab_slack]",
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 0,
                    "value": np.array(
                        [
                            0.0,
                            0.0,
                            -9.77208577e-09,
                            -9.54438112e-09,
                            -9.54432513e-09,
                            -9.54432506e-09,
                            0.8,
                            0.8,
                            0.8,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            -9.77208594e-09,
                            -9.77208576e-09,
                            -9.54438113e-09,
                            -9.54438113e-09,
                            -9.54432505e-09,
                            -9.54432507e-09,
                            -9.54432502e-09,
                            -9.54432513e-09,
                            0.8,
                            0.8,
                            0.8,
                            0.8,
                            0.8,
                            0.8,
                        ]
                    ),
                },
                "fs.slack[cd_slack]": {
                    "full_name": "fs.slack[cd_slack]",
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 0,
                    "value": np.array(
                        [
                            0.0,
                            0.0,
                            0.5,
                            -9.89148858e-09,
                            -9.77226571e-09,
                            0.5,
                            -9.77798881e-09,
                            -9.77226571e-09,
                            0.5,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            -9.77228650e-09,
                            -9.77247711e-09,
                            -9.88979157e-09,
                            -9.88979157e-09,
                            -9.77216575e-09,
                            -9.77218082e-09,
                            -9.77228649e-09,
                            -9.77247712e-09,
                            -9.88962016e-09,
                            -9.88962016e-09,
                            -9.77216575e-09,
                            -9.77218081e-09,
                            -9.77228649e-09,
                            -9.77247712e-09,
                        ]
                    ),
                },
                "fs.slack_penalty": {
                    "full_name": "fs.slack_penalty",
                    "units": "None",
                    "value": np.array(
                        [
                            1000.0,
                            1000.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            1000.0,
                            1000.0,
                            1000.0,
                            1000.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                        ]
                    ),
                },
                "objective": {
                    "full_name": "objective",
                    "value": np.array(
                        [
                            0.2,
                            0.95,
                            -3.79999979,
                            1.0,
                            1.75,
                            -3.0,
                            -7.0,
                            -6.24999979,
                            -10.99999979,
                            0.2,
                            0.2,
                            0.59330506,
                            0.65569601,
                            0.98661034,
                            1.11139224,
                            1.0,
                            1.0,
                            1.39330528,
                            1.45569622,
                            1.78661034,
                            1.91139224,
                            -6.99999978,
                            -6.99999978,
                            -6.60669472,
                            -6.54430377,
                            -6.21338966,
                            -6.08860776,
                        ]
                    ),
                },
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
                True,
                True,
                True,
                True,
                True,
                True,
                True,
            ],
            "sweep_params": {
                "fs.input[a]": {
                    "full_name": "fs.input[a]",
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
                            0.1,
                            0.1,
                            0.1,
                            0.1,
                            0.1,
                            0.1,
                            0.5,
                            0.5,
                            0.5,
                            0.5,
                            0.5,
                            0.5,
                            0.9,
                            0.9,
                            0.9,
                            0.9,
                            0.9,
                            0.9,
                        ]
                    ),
                },
                "fs.input[b]": {
                    "full_name": "fs.input[b]",
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
                            0.0,
                            0.13110169,
                            0.15189867,
                            0.26220338,
                            0.30379734,
                            0.0,
                            0.0,
                            0.13110169,
                            0.15189867,
                            0.26220338,
                            0.30379734,
                            0.0,
                            0.0,
                            0.13110169,
                            0.15189867,
                            0.26220338,
                            0.30379734,
                        ]
                    ),
                },
            },
        }

        read_dict = _read_output_h5(h5_results_file_name)
        _assert_h5_csv_agreement(csv_results_file_name, read_dict)
        if ps.parallel_manager.number_of_processes() > 1:
            # Compare the sorted dictionary. We need to work with a sorted dictionary
            # because the differential parameter sweep produces a global dictionary
            # that is jumbled by the number of procs.
            if ps.parallel_manager.rank == 0:
                sorted_truth_dict = sort_output_dict(truth_dict)
                sorted_read_dict = sort_output_dict(read_dict)
                sorted_global_results_dict = sort_output_dict(global_results_dict)
                _assert_dictionary_correctness(sorted_truth_dict, sorted_read_dict)
                _assert_dictionary_correctness(
                    sorted_truth_dict, sorted_global_results_dict
                )
        else:
            _assert_dictionary_correctness(truth_dict, global_results_dict)
            _assert_dictionary_correctness(truth_dict, read_dict)


@pytest.mark.component
def test_differential_parameter_sweep_function(model, tmp_path):

    comm = MPI.COMM_WORLD
    tmp_path = _get_rank0_path(comm, tmp_path)

    results_fname = os.path.join(tmp_path, "global_results")
    csv_results_file_name = str(results_fname) + ".csv"
    h5_results_file_name = str(results_fname) + ".h5"

    m = model
    m.fs.slack_penalty = 1000.0
    m.fs.slack.setub(0)

    A = m.fs.input["a"]
    B = m.fs.input["b"]
    sweep_params = {A.name: (A, 0.1, 0.9, 3), B.name: (B, 0.0, 0.5, 3)}

    differential_sweep_specs = {
        A.name: {
            "diff_mode": "sum",
            "diff_sample_type": NormalSample,
            "std_dev": 0.01,
            "pyomo_object": m.fs.input["a"],
        },
        B.name: {
            "diff_mode": "product",
            "diff_sample_type": UniformSample,
            "relative_lb": 0.01,
            "relative_ub": 0.01,
            "pyomo_object": m.fs.input["b"],
        },
    }

    # Call the parameter_sweep function
    global_results_dict, _ = differential_parameter_sweep(
        model,
        sweep_params,
        differential_sweep_specs,
        outputs=None,
        mpi_comm=comm,
        csv_results_file_name=csv_results_file_name,
        h5_results_file_name=h5_results_file_name,
        debugging_data_dir=tmp_path,
        interpolate_nan_outputs=True,
        optimize_function=_optimization,
        reinitialize_function=_reinitialize,
        reinitialize_kwargs={"slack_penalty": 10.0},
        seed=0,
    )

    if comm.rank == 0:

        truth_dict = {
            "differential_idx": np.array(
                [
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    0.0,
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                    5.0,
                    6.0,
                    7.0,
                    8.0,
                ]
            ),
            "nominal_idx": np.array(
                [
                    0.0,
                    1.0,
                    2.0,
                    3.0,
                    4.0,
                    5.0,
                    6.0,
                    7.0,
                    8.0,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                ]
            ),
            "outputs": {
                "fs.input[a]": {
                    "full_name": "fs.input[a]",
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
                    "full_name": "fs.input[b]",
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
                "fs.output[c]": {
                    "full_name": "fs.output[c]",
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
                    "full_name": "fs.output[d]",
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
                            7.50000010e-01,
                            1.0,
                            1.09912079e-22,
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
                    "full_name": "fs.performance",
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
                            1.0075,
                            1.0150,
                            1.0,
                            1.0075,
                            1.015,
                        ]
                    ),
                },
                "fs.slack[ab_slack]": {
                    "full_name": "fs.slack[ab_slack]",
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 0,
                    "value": np.array(
                        [
                            0.0,
                            0.0,
                            -9.77208577e-09,
                            -9.54438112e-09,
                            -9.54432508e-09,
                            -9.54432509e-09,
                            0.8,
                            0.8,
                            0.8,
                            0.0,
                            0.0,
                            -9.77210838e-09,
                            3.52810371e-02,
                            3.52810371e-02,
                            3.52810371e-02,
                            8.35281037e-01,
                            8.35281037e-01,
                            8.35281037e-01,
                        ]
                    ),
                },
                "fs.slack[cd_slack]": {
                    "full_name": "fs.slack[cd_slack]",
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 0,
                    "value": np.array(
                        [
                            0.0,
                            0.0,
                            0.5,
                            -9.98996974e-09,
                            -9.77226571e-09,
                            4.99999990e-01,
                            -9.92236517e-09,
                            -9.77226571e-09,
                            4.99999990e-01,
                            0.0,
                            0.0,
                            -9.77036119e-09,
                            -9.77884773e-09,
                            -9.76848697e-09,
                            -9.77036119e-09,
                            -9.77884942e-09,
                            -9.76848697e-09,
                            -9.77036119e-09,
                        ]
                    ),
                },
                "fs.slack_penalty": {
                    "full_name": "fs.slack_penalty",
                    "units": "None",
                    "value": np.array(
                        [
                            1000.0,
                            1000.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            1000.0,
                            1000.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                            10.0,
                        ]
                    ),
                },
                "objective": {
                    "full_name": "objective",
                    "value": np.array(
                        [
                            0.2,
                            0.95,
                            -3.8,
                            1.0,
                            1.75,
                            -3.0,
                            -6.99999978,
                            -6.24999979,
                            -10.99999979,
                            0.23528105,
                            0.24278105,
                            0.25028126,
                            0.64718975,
                            0.65468975,
                            0.66218975,
                            -7.35281025,
                            -7.34531025,
                            -7.33781025,
                        ]
                    ),
                },
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
                "fs.input[a]": {
                    "full_name": "fs.input[a]",
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
                    "full_name": "fs.input[b]",
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
        _assert_h5_csv_agreement(csv_results_file_name, read_dict)
        if comm.Get_size() > 1:
            # Compare the sorted dictionary. We need to work with a sorted dictionary
            # because the differential parameter sweep produces a global dictionary
            # that is jumbled by the number of procs.
            if comm.rank == 0:
                sorted_truth_dict = sort_output_dict(truth_dict)
                sorted_read_dict = sort_output_dict(read_dict)
                sorted_global_results_dict = sort_output_dict(global_results_dict)
                _assert_dictionary_correctness(sorted_truth_dict, sorted_read_dict)
                _assert_dictionary_correctness(
                    sorted_truth_dict, sorted_global_results_dict
                )
        else:
            _assert_dictionary_correctness(truth_dict, global_results_dict)
            _assert_dictionary_correctness(truth_dict, read_dict)


def sort_output_dict(input_dict):
    """Simple utility function to sort all values in ascending order"""

    sorted_dict = copy.deepcopy(input_dict)
    for key, item in input_dict.items():
        if key in ["sweep_params", "outputs"]:  # != "solve_successful":
            for subkey, subitem in item.items():
                sorted_dict[key][subkey]["value"] = np.sort(subitem["value"])
        elif key == "solve_successful":
            sorted_dict["solve_successful"] = np.sort(input_dict[key]).tolist()
        elif key in ["nominal_idx", "differential_idx"]:
            sorted_dict[key] = np.sort(item)

    return sorted_dict
