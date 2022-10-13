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
import yaml

from pyomo.environ import value

from watertap.tools.parameter_sweep.parameter_sweep_reader import (
    ParameterSweepReader,
    get_sweep_params_from_yaml,
    set_defaults_from_yaml,
)

# Imports for conditional fails
from idaes.config import bin_directory as idaes_bin_directory


class TestInputParser:
    @pytest.fixture(scope="class")
    def model(self):
        m = pyo.ConcreteModel()
        m.fs = pyo.Block()
        m.fs.x = pyo.Var(initialize=-0.5)
        m.fs.a = pyo.Param(initialize=1.0, mutable=True)
        m.fs.b = pyo.Param(initialize=2.0, mutable=False)
        m.fs.sum = pyo.Expression(expr=m.fs.a + m.fs.x)
        m.fs.lim = pyo.Constraint(expr=m.fs.sum >= 0.0)
        m.fs.obj = pyo.Objective(expr=m.fs.x, sense=pyo.minimize)

        return m

    @pytest.fixture(scope="class")
    def get_simple_yaml_file(self):
        filename = "temp.yaml"
        dict_to_write = {"a_val": 1.1, "b_val": 2.2, "c_val": 3.3, "d_val": 4.4}

        with open(filename, "w") as fp:
            yaml.dump(dict_to_write, fp, sort_keys=False)

        return filename, dict_to_write

    @pytest.fixture(scope="class")
    def get_param_yaml_file(self):
        filename = "temp.yaml"
        dict_to_write = {
            "a_val": {
                "type": "LinearSample",
                "param": "fs.a",
                "lower_limit": 10,
                "upper_limit": 20,
                "num_samples": 5,
            },
            "b_val": {
                "type": "UniformSample",
                "param": "fs.a",
                "lower_limit": 10,
                "upper_limit": 20,
            },
            "c_val": {"type": "NormalSample", "param": "fs.a", "mean": 10, "std": 2},
            "d_val": {
                "type": "LatinHypercubeSample",
                "param": "fs.a",
                "lower_limit": 10,
                "upper_limit": 20,
            },
            "e_val": {
                "type": "GeomSample",
                "param": "fs.a",
                "lower_limit": 1,
                "upper_limit": 10,
                "num_samples": 3,
            },
            "f_val": {
                "type": "ReverseGeomSample",
                "param": "fs.a",
                "lower_limit": 1,
                "upper_limit": 10,
                "num_samples": 3,
            },
        }

        with open(filename, "w") as fp:
            yaml.dump(dict_to_write, fp, sort_keys=False)

        return filename, dict_to_write

    @pytest.fixture(scope="class")
    def get_param_yaml_file_error(self):
        filename = "temp.yaml"
        dict_to_write = {"a_val": {"param": "nonexistent_variable"}}

        with open(filename, "w") as fp:
            yaml.dump(dict_to_write, fp, sort_keys=False)

        return filename, dict_to_write

    @pytest.fixture(scope="class")
    def get_default_yaml_file(self, model):
        filename = "temp.yaml"

        m = model

        dict_to_write = {}

        for pyo_obj in m.component_data_objects((pyo.Var, pyo.Param), active=True):
            for k in range(10):
                dict_to_write[pyo_obj.name] = np.random.rand()

        with open(filename, "w") as fp:
            yaml.dump(dict_to_write, fp, sort_keys=False)

        return filename, dict_to_write

    @pytest.fixture(scope="class")
    def get_default_yaml_file_error(self, model):
        filename = "temp.yaml"

        dict_to_write = {"nonexistent_variable": 1.0}

        with open(filename, "w") as fp:
            yaml.dump(dict_to_write, fp, sort_keys=False)

        return filename, dict_to_write

    @pytest.mark.unit
    def test_yaml_to_dict(self, get_simple_yaml_file):
        filename, truth_values = get_simple_yaml_file

        psr = ParameterSweepReader()
        values = psr._yaml_to_dict(filename)

        for key, truth_value in truth_values.items():
            assert key in values
            assert values[key] == truth_value

    @pytest.mark.unit
    def test_yaml_to_dict_error(self):
        psr = ParameterSweepReader()
        with pytest.raises(Exception):
            values = psr._yaml_to_dict("nonexistent_file.yaml")

    @pytest.mark.unit
    def test_get_sweep_params_from_yaml(self, model, get_param_yaml_file):
        m = model
        filename, truth_values = get_param_yaml_file

        values = get_sweep_params_from_yaml(m, filename)

        for key, truth_value in truth_values.items():
            assert key in values

    @pytest.mark.unit
    def test_get_sweep_params_from_yaml_error(self, model, get_param_yaml_file_error):
        m = model
        filename, truth_values = get_param_yaml_file_error

        psr = ParameterSweepReader()
        with pytest.raises(ValueError):
            psr.get_sweep_params_from_yaml(m, filename)

    @pytest.mark.unit
    def test_set_defaults_from_yaml(self, model, get_default_yaml_file):
        m = model
        filename, truth_values = get_default_yaml_file

        set_defaults_from_yaml(m, filename)

        for key, truth_value in truth_values.items():
            component = m.find_component(key)
            if component is m.fs.b:  # don't set non-mutable params
                assert value(component) == 2.0
            else:
                assert value(component) == truth_value

    @pytest.mark.unit
    def test_set_defaults_from_yaml_error(self, model, get_default_yaml_file_error):
        m = model
        filename, truth_values = get_default_yaml_file_error

        psr = ParameterSweepReader()
        with pytest.raises(ValueError):
            psr.set_defaults_from_yaml(m, filename)

        os.remove(filename)
