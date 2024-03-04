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
import pyomo.environ as pyo
import idaes.core.util.model_statistics as istat

from pyomo.util.check_units import assert_units_consistent
from idaes.core import UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.testing import initialization_tester
from watertap.costing import WaterTAPCosting
from watertap.unit_models.tests.test_electrolyzer import build

__author__ = "Hunter Barber"

solver = get_solver()


class TestElectrolyzerCosting:
    @pytest.fixture(scope="class")
    def build_costing(self):
        m = build()
        initialization_tester(m)
        solver.solve(m)

        # build costing model block
        m.fs.costing = WaterTAPCosting()
        m.fs.costing.base_currency = pyo.units.USD_2020

        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.unit.costing.initialize()

        return m

    @pytest.mark.unit
    def test_costing(self, build_costing):

        m = build_costing

        # testing gac costing block dof and initialization
        assert assert_units_consistent(m) is None
        assert istat.degrees_of_freedom(m) == 0
        m.fs.unit.costing.initialize()

        # solve
        results = solver.solve(m)

        # Check for optimal solution
        assert pyo.check_optimal_termination(results)

        # check solution values
        assert pytest.approx(2.0 * 17930, rel=1e-3) == pyo.value(
            m.fs.unit.costing.capital_cost
        )
        assert pytest.approx(82.50, rel=1e-3) == pyo.value(
            m.fs.costing.aggregate_flow_electricity
        )
        assert pytest.approx(50040, rel=1e-3) == pyo.value(
            m.fs.costing.aggregate_flow_costs["electricity"]
        )
