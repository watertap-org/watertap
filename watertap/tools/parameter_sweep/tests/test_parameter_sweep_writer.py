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
import copy

from watertap.tools.parameter_sweep.parameter_sweep import *
from watertap.tools.parameter_sweep.parameter_sweep_writer import *
from watertap.tools.parameter_sweep.tests.test_parameter_sweep import (
    _get_rank0_path,
    _read_output_h5,
    _assert_dictionary_correctness,
)

# ------------------------------------------------------------------------------


class TestParallelWriterManager:
    @pytest.fixture(scope="class")
    def model(self):
        m = pyo.ConcreteModel()
        m.fs = fs = pyo.Block()

        fs.input = pyo.Var(["a", "b"], within=pyo.UnitInterval, initialize=0.5)
        fs.output = pyo.Var(["c", "d"], within=pyo.UnitInterval, initialize=0.5)

        fs.slack = pyo.Var(["ab_slack", "cd_slack"], bounds=(0, 0), initialize=0.0)
        fs.slack_penalty = pyo.Param(
            default=1000.0, mutable=True, within=pyo.PositiveReals
        )

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

    @pytest.mark.unit
    def test_interp_nan_values(self):
        ps = ParameterSweep()
        ps_writer = ParameterSweepWriter(
            ps.parallel_manager,
            csv_results_file_name=None,
            h5_results_file_name=None,
            debugging_data_dir=None,
            interpolate_nan_outputs=True,
        )

        global_values = np.array(
            [
                [0, 0, 0],
                [0, 0, 1],
                [0, 1, 0],
                [0, 1, 1],
                [1, 0, 0],
                [1, 0, 1],
                [1, 1, 0],
                [1, 1, 1],
                [0.5, 0.5, 0.5],
                [1, 1, 1],
            ]
        )

        global_results = np.array([0, 1, 2, 3, 4, 5, 6, 7, np.nan, np.nan])[
            np.newaxis
        ].T

        global_results_clean = ps_writer._interp_nan_values(
            global_values, global_results
        )

        assert np.shape(global_results_clean)[1] == np.shape(global_results)[1]
        assert np.shape(global_results_clean)[0] == np.shape(global_results)[0]

        assert (global_results_clean[8]) == pytest.approx(np.mean(global_results[0:8]))
        assert (global_results_clean[9]) == pytest.approx(global_results[7])

    @pytest.mark.unit
    def test_h5_read_write(self, tmp_path):
        ps = ParameterSweep()

        tmp_path = _get_rank0_path(None, tmp_path)
        h5_fname = "h5_test_{0}.h5".format(ps.parallel_manager.get_rank())

        ps_writer = ParameterSweepWriter(
            ps.parallel_manager,
            csv_results_file_name=None,
            h5_results_file_name=h5_fname,
            debugging_data_dir=None,
            interpolate_nan_outputs=True,
        )

        input_dict = {
            "outputs": {
                "fs.input[a]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.input[b]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.output[c]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.2, 0.2, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.output[d]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.75, 0.0, 0.0, 0.75, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.slack[ab_slack]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.slack[cd_slack]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
            },
            "solve_successful": [True] * 9,
            "sweep_params": {
                "fs.input[a]": {
                    "lower bound": None,
                    "units": "None",
                    "upper bound": None,
                    "value": np.array([0.1, 0.1, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.input[b]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.25, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0]),
                },
            },
        }

        reference_dict = copy.deepcopy(input_dict)
        reference_dict["sweep_params"]["fs.input[a]"]["lower bound"] = np.finfo("d").min
        reference_dict["sweep_params"]["fs.input[a]"]["upper bound"] = np.finfo("d").max

        ps_writer._write_output_to_h5(input_dict, os.path.join(tmp_path, h5_fname))

        read_dictionary = _read_output_h5(os.path.join(tmp_path, h5_fname))
        _assert_dictionary_correctness(reference_dict, read_dictionary)

    @pytest.mark.unit
    def test_h5_parents(self, tmp_path):
        ps = ParameterSweep()

        tmp_path = _get_rank0_path(None, tmp_path)
        h5_fname = "h5_test_{0}.h5".format(ps.parallel_manager.get_rank())
        h5_parent_groups = ["loop1", "loop2"]

        embedded_dict = {
            "outputs": {
                "fs.input[a]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.input[b]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.output[c]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.2, 0.2, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.output[d]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.75, 0.0, 0.0, 0.75, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.slack[ab_slack]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.slack[cd_slack]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                },
            },
            "solve_successful": [True] * 9,
            "sweep_params": {
                "fs.input[a]": {
                    "lower bound": None,
                    "units": "None",
                    "upper bound": None,
                    "value": np.array([0.1, 0.1, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0]),
                },
                "fs.input[b]": {
                    "lower bound": 0,
                    "units": "None",
                    "upper bound": 1,
                    "value": np.array([0.0, 0.25, 0.0, 0.0, 0.25, 0.0, 0.0, 0.0, 0.0]),
                },
            },
        }

        ps_writer_dict = dict()
        reference_dict = dict()
        for h5_parent_group in h5_parent_groups:
            ps_writer_dict[h5_parent_group] = ParameterSweepWriter(
                ps.parallel_manager,
                csv_results_file_name=None,
                h5_results_file_name=h5_fname,
                h5_parent_group_name=h5_parent_group,
                debugging_data_dir=None,
                interpolate_nan_outputs=True,
            )

            reference_dict[h5_parent_group] = embedded_dict
            reference_dict[h5_parent_group]["sweep_params"]["fs.input[a]"][
                "lower bound"
            ] = np.finfo("d").min
            reference_dict[h5_parent_group]["sweep_params"]["fs.input[a]"][
                "upper bound"
            ] = np.finfo("d").max

            ps_writer_dict[h5_parent_group]._write_output_to_h5(
                embedded_dict, os.path.join(tmp_path, h5_fname)
            )

        # Now read the dictionary
        read_dictionary = _read_loop_h5(
            os.path.join(tmp_path, h5_fname), h5_parent_groups
        )
        _assert_loop_correctness(reference_dict, read_dictionary)


def _read_loop_h5(filepath, parent_list):

    f = h5py.File(filepath, "r")

    output_dict = {}
    for parent in parent_list:
        output_dict[parent] = _read_output_h5(f[parent])

    f.close()

    return output_dict


def _assert_loop_correctness(truth_dict, test_dict):
    assert truth_dict.keys() == test_dict.keys()

    for key, item in truth_dict.items():
        _assert_dictionary_correctness(item, test_dict[key])
