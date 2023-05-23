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
import numpy as np
import pyomo.environ as pyo

from pyomo.environ import value

from watertap.tools.parameter_sweep.sampling_types import *
from watertap.tools.parameter_sweep import ParameterSweep, parameter_sweep
from watertap.tools.parameter_sweep.parameter_sweep_writer import *

from watertap.tools.parameter_sweep.rayio_param_sweep_func import do_rayio_sweep
from watertap.tools.parameter_sweep.tests.test_parameter_sweep import (
    _read_output_h5,
    _assert_dictionary_correctness,
    _assert_h5_csv_agreement,
    _optimization,
)

# -----------------------------------------------------------------------------


def build(**kwargs):
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


def reinitialization(**kwargs):
    pass


@pytest.mark.component
def test_parameter_sweep(tmp_path):
    results_fname = os.path.join(tmp_path, "global_results")
    csv_results_file_name = str(results_fname) + ".csv"
    h5_results_file_name = str(results_fname) + ".h5"

    custom_do_param_sweep_kwargs = {
        "build_function": build,
        "build_kwargs": {},
        "num_workers": 10,
        "load_form_json": True,
    }
    ps = ParameterSweep(
        optimize_function=_optimization,
        csv_results_file_name=csv_results_file_name,
        h5_results_file_name=h5_results_file_name,
        debugging_data_dir=tmp_path,
        interpolate_nan_outputs=True,
        custom_do_param_sweep=do_rayio_sweep,
        custom_do_param_sweep_kwargs=custom_do_param_sweep_kwargs,
    )

    m = build()
    m.fs.slack_penalty = 1000.0
    m.fs.slack.setub(0)

    A = m.fs.input["a"]
    B = m.fs.input["b"]
    sweep_params = {A.name: (A, 0.1, 0.9, 3), B.name: (B, 0.0, 0.5, 3)}
    outputs = {
        "output_c": m.fs.output["c"],
        "output_d": m.fs.output["d"],
        "performance": m.fs.performance,
    }

    # Call the parameter_sweep function
    _ = ps.parameter_sweep(
        m,
        sweep_params,
        outputs=outputs,
    )

    # NOTE: rank 0 "owns" tmp_path, so it needs to be
    #       responsible for doing any output file checking
    #       tmp_path can be deleted as soon as this method
    #       returns
    if ps.rank == 0:
        # Check that the global results file is created
        assert os.path.isfile(csv_results_file_name)
        assert os.path.isfile(os.path.join(tmp_path, "interpolated_global_results.csv"))

        # Check that all local output files have been created
        for k in range(ps.num_procs):
            assert os.path.isfile(os.path.join(tmp_path, f"local_results_{k:03}.h5"))
            assert os.path.isfile(os.path.join(tmp_path, f"local_results_{k:03}.csv"))

        # Attempt to read in the data
        data = np.genfromtxt(csv_results_file_name, skip_header=1, delimiter=",")
        print(data[-1])

        print(data)
        # Compare the last row of the imported data to truth
        truth_data = [0.9, 0.5, 1, 1, 2]
        assert np.allclose(data[-1], truth_data, equal_nan=True)

    # Check for the h5 output

    truth_dict = {
        "outputs": {
            "output_c": {
                "lower bound": 0,
                "units": "None",
                "upper bound": 1,
                "value": np.array(
                    [0.2, 0.2, np.nan, 1.0, 1.0, np.nan, np.nan, np.nan, np.nan]
                ),
            },
            "output_d": {
                "lower bound": 0,
                "units": "None",
                "upper bound": 1,
                "value": np.array(
                    [
                        0.0,
                        0.75,
                        np.nan,
                        0.0,
                        0.75,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                    ]
                ),
            },
            "performance": {
                "value": np.array(
                    [
                        0.2,
                        0.95,
                        np.nan,
                        1.0,
                        1.75,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                    ]
                )
            },
        },
        "solve_successful": [
            True,
            True,
            False,
            True,
            True,
            False,
            False,
            False,
            False,
        ],
        "sweep_params": {
            "fs.input[a]": {
                "lower bound": 0,
                "units": "None",
                "upper bound": 1,
                "value": np.array([0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.9, 0.9, 0.9]),
            },
            "fs.input[b]": {
                "lower bound": 0,
                "units": "None",
                "upper bound": 1,
                "value": np.array([0.0, 0.25, 0.5, 0.0, 0.25, 0.5, 0.0, 0.25, 0.5]),
            },
        },
    }

    read_dict = _read_output_h5(h5_results_file_name)

    _assert_dictionary_correctness(truth_dict, read_dict)
    _assert_h5_csv_agreement(csv_results_file_name, read_dict)

    # Check if there is a text file created
    import ast

    truth_txt_dict = {
        "outputs": ["output_c", "output_d", "performance"],
        "sweep_params": ["fs.input[a]", "fs.input[b]"],
    }

    txt_fpath = h5_results_file_name + ".txt"
    assert os.path.exists(txt_fpath)
    f = open(txt_fpath)
    f_contents = f.read()
    read_txt_dict = ast.literal_eval(f_contents)
    assert read_txt_dict == truth_txt_dict
