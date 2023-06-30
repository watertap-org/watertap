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
import idaes.core as idc
import idaes.core.util.model_statistics as istat

from pyomo.util.check_units import assert_units_consistent
from idaes.core.solvers import get_solver
from watertap.costing.watertap_costing_package import WaterTAPCosting
from watertap.costing.units.rectifier import cost_rectifier

solver = get_solver()


@pytest.mark.component
def test_rectifier_costing():

    # build generic models
    m = pyo.ConcreteModel()
    m.fs = idc.FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()
    # change base units for value assertions
    m.fs.costing.base_currency = pyo.units.USD_2021
    m.fs.unit = idc.UnitModelBlock()

    # build power variable required for cost_rectifier
    m.fs.unit.power = pyo.Var(
        initialize=100,
        domain=pyo.NonNegativeReals,
        units=pyo.units.kW,
    )
    m.fs.unit.power.fix()

    # cost rectifier
    m.fs.unit.costing = idc.UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=cost_rectifier,
    )

    # cost process for cost_flow to be applied
    m.fs.costing.cost_process()

    # check model and solve
    assert_units_consistent(m)
    assert istat.degrees_of_freedom(m) == 0
    results = solver.solve(m)
    assert pyo.check_optimal_termination(results)

    # check values
    assert pytest.approx(59320, rel=1e-3) == pyo.value(m.fs.unit.costing.capital_cost)
    assert pytest.approx(111.1, rel=1e-3) == pyo.value(m.fs.unit.costing.ac_power)
    assert pytest.approx(111.1, rel=1e-3) == pyo.value(
        m.fs.costing.aggregate_flow_electricity
    )
    assert pytest.approx(80040, rel=1e-3) == pyo.value(
        m.fs.costing.aggregate_flow_costs["electricity"]
    )
