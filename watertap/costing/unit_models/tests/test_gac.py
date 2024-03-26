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
import idaes.core.util.model_statistics as istat

from pyomo.util.check_units import assert_units_consistent
from idaes.core import (
    UnitModelCostingBlock,
)
from idaes.core.solvers import get_solver
from idaes.core.util.testing import initialization_tester
from watertap.costing import WaterTAPCosting
from watertap.unit_models.tests.test_gac import build_crittenden


__author__ = "Hunter Barber"

solver = get_solver()


class TestGACCosting:
    @pytest.fixture(scope="class")
    def build(self):
        m = build_crittenden()
        initialization_tester(m)
        solver.solve(m)

        return m

    @pytest.mark.component
    def test_robust_costing_pressure(self, build):
        m = build

        m.fs.costing = WaterTAPCosting()
        m.fs.costing.base_currency = pyo.units.USD_2020

        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()

        # testing gac costing block dof and initialization
        assert assert_units_consistent(m) is None
        assert istat.degrees_of_freedom(m) == 0
        m.fs.unit.costing.initialize()

        # solve
        results = solver.solve(m)

        # Check for optimal solution
        assert pyo.check_optimal_termination(results)

        cost = m.fs.unit.costing
        # Check for known cost solution of default twin alternating contactors
        assert pyo.value(m.fs.costing.gac_pressure.num_contactors_op) == 1
        assert pyo.value(m.fs.costing.gac_pressure.num_contactors_redundant) == 1
        assert pytest.approx(56900, rel=1e-3) == pyo.value(cost.contactor_cost)
        assert pytest.approx(4.359, rel=1e-3) == pyo.value(cost.adsorbent_unit_cost)
        assert pytest.approx(17450, rel=1e-3) == pyo.value(cost.adsorbent_cost)
        assert pytest.approx(81690, rel=1e-3) == pyo.value(cost.other_process_cost)
        assert pytest.approx(2.0 * 156000, rel=1e-3) == pyo.value(cost.capital_cost)
        assert pytest.approx(12680, rel=1e-3) == pyo.value(cost.gac_makeup_cost)
        assert pytest.approx(27660, rel=1e-3) == pyo.value(cost.gac_regen_cost)
        assert pytest.approx(0.01631, rel=1e-3) == pyo.value(cost.energy_consumption)
        assert pytest.approx(40370, rel=1e-3) == pyo.value(cost.fixed_operating_cost)

    @pytest.mark.component
    def test_robust_costing_gravity(self, build):
        mr_grav = build

        mr_grav.fs.costing = WaterTAPCosting()
        mr_grav.fs.costing.base_currency = pyo.units.USD_2020

        mr_grav.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=mr_grav.fs.costing,
            costing_method_arguments={"contactor_type": "gravity"},
        )
        mr_grav.fs.costing.cost_process()

        # testing gac costing block dof and initialization
        assert assert_units_consistent(mr_grav) is None
        assert istat.degrees_of_freedom(mr_grav) == 0
        mr_grav.fs.unit.costing.initialize()

        # solve
        results = solver.solve(mr_grav)

        # Check for optimal solution
        assert pyo.check_optimal_termination(results)

        cost = mr_grav.fs.unit.costing
        # Check for known cost solution of default twin alternating contactors
        assert pyo.value(mr_grav.fs.costing.gac_gravity.num_contactors_op) == 1
        assert pyo.value(mr_grav.fs.costing.gac_gravity.num_contactors_redundant) == 1
        assert pytest.approx(163200, rel=1e-3) == pyo.value(cost.contactor_cost)
        assert pytest.approx(4.359, rel=1e-3) == pyo.value(cost.adsorbent_unit_cost)
        assert pytest.approx(17450, rel=1e-3) == pyo.value(cost.adsorbent_cost)
        assert pytest.approx(159500, rel=1e-3) == pyo.value(cost.other_process_cost)
        assert pytest.approx(2.0 * 340200, rel=1e-3) == pyo.value(cost.capital_cost)
        assert pytest.approx(12680, rel=1e-3) == pyo.value(cost.gac_makeup_cost)
        assert pytest.approx(27660, rel=1e-3) == pyo.value(cost.gac_regen_cost)
        assert pytest.approx(2.476, rel=1e-3) == pyo.value(cost.energy_consumption)
        assert pytest.approx(40370, rel=1e-3) == pyo.value(cost.fixed_operating_cost)

    @pytest.mark.component
    def test_robust_costing_modular_contactors(self, build):
        m = build

        m.fs.costing = WaterTAPCosting()
        m.fs.costing.base_currency = pyo.units.USD_2020

        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()

        m.fs.costing.gac_pressure.num_contactors_op.fix(4)
        m.fs.costing.gac_pressure.num_contactors_redundant.fix(2)

        results = solver.solve(m)

        cost = m.fs.unit.costing
        # Check for known cost solution when changing volume scale of vessels in parallel
        assert pyo.value(m.fs.costing.gac_pressure.num_contactors_op) == 4
        assert pyo.value(m.fs.costing.gac_pressure.num_contactors_redundant) == 2
        assert pytest.approx(89040, rel=1e-3) == pyo.value(cost.contactor_cost)
        assert pytest.approx(69690, rel=1e-3) == pyo.value(cost.other_process_cost)
        assert pytest.approx(2.0 * 176200, rel=1e-3) == pyo.value(cost.capital_cost)

    @pytest.mark.component
    def test_robust_costing_max_gac_ref(self, build):
        m = build

        # scale flow up 10x
        unit_feed = m.fs.unit.process_flow.properties_in[0]
        unit_feed.flow_mol_phase_comp["Liq", "H2O"].fix(10 * 824.0736620370348)
        unit_feed.flow_mol_phase_comp["Liq", "TCE"].fix(10 * 5.644342973110135e-05)

        m.fs.costing = WaterTAPCosting()
        m.fs.costing.base_currency = pyo.units.USD_2020

        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        # not necessarily an optimum solution because poor scaling
        # but just checking the conditional
        results = solver.solve(m)

        # Check for bed_mass_gac_cost_ref to be overwritten
        # if bed_mass_gac is greater than bed_mass_gac_cost_max_ref
        assert pyo.value(m.fs.unit.bed_mass_gac) > pyo.value(
            m.fs.costing.gac_pressure.bed_mass_max_ref
        )
        assert pyo.value(m.fs.unit.costing.bed_mass_gac_ref) == (
            pytest.approx(pyo.value(m.fs.costing.gac_pressure.bed_mass_max_ref), 1e-5)
        )
