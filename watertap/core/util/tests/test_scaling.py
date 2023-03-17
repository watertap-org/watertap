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

from idaes.core import FlowsheetBlock

import watertap.property_models.NaCl_prop_pack as props

from pyomo.environ import ConcreteModel, Var, Constraint, Block, SolverFactory

from idaes.core.solvers import get_solver
from watertap.core.util.initialization import (
    check_dof,
    assert_degrees_of_freedom,
    assert_no_degrees_of_freedom,
    check_solve,
)
import idaes.logger as idaeslog

__author__ = "Marcus Holly"

_log = idaeslog.getLogger(__name__)

# Set up solver
solver = get_solver()

# add test for checking that a constraint is made for eq + var_str
# add test for checking that warning displays when name is wrong
# add test for checking that warning displays when the var has no property constraint


class TestScaling:
    @pytest.fixture(scope="class")
    def m(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = props.NaClParameterBlock()
        return m
