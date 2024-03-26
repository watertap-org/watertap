#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
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
import pyomo.environ as pyo

from pyomo.environ import ConcreteModel

from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import (
    unscaled_constraints_generator,
    calculate_scaling_factors,
)
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

import watertap.property_models.NaCl_prop_pack as props

import idaes.logger as idaeslog

__author__ = "Marcus Holly"

_log = idaeslog.getLogger(__name__)

# Set up solver
solver = get_solver()


# Start test class
# TODO: Consider using dummy metadata rather than importing property package
class TestScaling:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = props.NaClParameterBlock()

        return m

    @pytest.mark.unit
    def test_constraint_scaling(self, model):
        m = model

        # Create the state block and pull in the default metadata
        m.fs.stream = m.fs.properties.build_state_block([0])
        metadata = m.fs.properties.get_metadata().properties

        # Check that constraints are unscaled
        for p in metadata.list_supported_properties():
            var_str = p.name
            if p.method is not None and m.fs.stream[0].is_property_constructed(var_str):
                if hasattr(self, "eq_" + var_str):
                    con = getattr(self, "eq_" + var_str)
                    unscaled_constraint_list = list(unscaled_constraints_generator(con))
                    assert len(unscaled_constraint_list) == 1

        # Scale constraints
        # #TODO: need to replace with transform_property_constraints
        for p in metadata.list_supported_properties():
            var_str = p.name
            if p.method is not None and m.fs.stream[0].is_property_constructed(var_str):
                var = getattr(self, var_str)
                if not isinstance(var, pyo.Var):
                    continue  # properties that are not vars do not have constraints
                # adding a conditional to check if a constraint exists for property; in the case when we only add and object reference, there would not be a constraint
                if hasattr(self, "eq_" + var_str):
                    con = getattr(self, "eq_" + var_str)
                    for ind, c in con.items():
                        sf = iscale.get_scaling_factor(
                            var[ind], default=1, warning=True
                        )
                        iscale.constraint_scaling_transform(c, sf)

        # Check that constraints are scaled
        for p in metadata.list_supported_properties():
            var_str = p.name
            if p.method is not None and m.fs.stream[0].is_property_constructed(var_str):
                if hasattr(self, "eq_" + var_str):
                    con = getattr(self, "eq_" + var_str)
                    unscaled_constraint_list = list(unscaled_constraints_generator(con))
                    assert len(unscaled_constraint_list) == 0

    @pytest.mark.unit
    def test_logger_warning(self, model, caplog):
        caplog.set_level(idaeslog.INFO, logger="watertap.core.util.scaling")

        m = model

        # Create the state block and pull in the default metadata
        metadata = m.fs.properties.get_metadata().properties

        # Scale constraints
        # TODO: need to replace with transform_property_constraints and find a way to force it to fail - use fail_flag or dummy data?
        for p in metadata.list_supported_properties():
            var_str = p.name
            if p.method is not None and m.fs.stream[0].is_property_constructed(var_str):
                var = getattr(self, var_str)
                if not isinstance(var, pyo.Var):
                    continue  # properties that are not vars do not have constraints
                # adding a conditional to check if a constraint exists for property; in the case when we only add and object reference, there would not be a constraint
                if hasattr(self, "eq_" + var_str):
                    con = getattr(self, "eq_" + var_str)
                    for ind, c in con.items():
                        sf = iscale.get_scaling_factor(
                            var[ind], default=1, warning=True
                        )
                        iscale.constraint_scaling_transform(c, sf)
                else:
                    msg = (
                        f"If there was a property constraint written for the variable, {var}, that constraint was not "
                        f"scaled. The transform_property_constraints tool expects constraints to have the following naming "
                        f"convention: 'eq_' + '{var_str}'. This suggests that the user may have defined a property in "
                        f"metadata but failed to follow the naming convention for its constraint. If there is no property "
                        f"constraint associated with the {var_str}, this warning can be ignored."
                    )
                    _log.warning(msg)

                assert (
                    "If there was a property constraint written for the variable, {var}, that constraint was not "
                    "scaled. The transform_property_constraints tool expects constraints to have the following naming "
                    "convention: 'eq_' + '{var_str}'. This suggests that the user may have defined a property in "
                    "metadata but failed to follow the naming convention for its constraint. If there is no property "
                    "constraint associated with the {var_str}, this warning can be ignored."
                    in caplog.text
                )

    @pytest.mark.unit
    def test_no_logger_warning_for_params_set_as_fixed_vars(self, model, caplog):
        caplog.set_level(idaeslog.INFO, logger="watertap.core.util.scaling")

        m = model

        calculate_scaling_factors(m.fs.stream[0])

        assert (
            "If there was a property constraint written for the variable"
            not in caplog.text
        )
