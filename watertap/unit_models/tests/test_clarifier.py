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
"""
Tests for costing clarifier.
"""
__author__ = "Chenyu Wang"

from watertap.unit_models.clarifier import Clarifier
from idaes.models.unit_models.separator import SplittingType
import pytest
from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    value,
    units,
)
from idaes.core import FlowsheetBlock
from watertap.property_models.activated_sludge.asm1_properties import ASM1ParameterBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import UnitModelCostingBlock
from idaes.core.util.testing import initialization_tester
from pyomo.util.check_units import assert_units_consistent
from watertap.costing import WaterTAPCosting
from watertap.costing.unit_models.clarifier import (
    cost_circular_clarifier,
    cost_rectangular_clarifier,
    cost_primary_clarifier,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestClarifierCosting:
    @pytest.fixture(scope="class")
    def clarifier_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = ASM1ParameterBlock()

        m.fs.unit = Clarifier(
            property_package=m.fs.properties,
            outlet_list=["underflow", "effluent"],
            split_basis=SplittingType.componentFlow,
        )

        m.fs.unit.inlet.temperature.fix(298.15 * units.K)
        m.fs.unit.inlet.pressure.fix(1 * units.atm)
        m.fs.unit.inlet.flow_vol.fix(18446 * units.m**3 / units.day)

        m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(27 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_S"].fix(58 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(92 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_S"].fix(363 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_BH"].fix(50 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_BA"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_P"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_O"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NO"].fix(0 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NH"].fix(23 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_ND"].fix(5 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_ND"].fix(16 * units.g / units.m**3)

        # Alkalinity was givien in mg/L based on C
        m.fs.unit.inlet.alkalinity[0].fix(7 * units.mol / units.m**3)

        # Unit option
        m.fs.unit.split_fraction[0, "effluent", "H2O"].fix(0.993)
        m.fs.unit.split_fraction[0, "effluent", "S_I"].fix(0.993)
        m.fs.unit.split_fraction[0, "effluent", "S_S"].fix(0.993)
        m.fs.unit.split_fraction[0, "effluent", "X_I"].fix(0.5192)
        m.fs.unit.split_fraction[0, "effluent", "X_S"].fix(0.5192)
        m.fs.unit.split_fraction[0, "effluent", "X_BH"].fix(0.5192)
        m.fs.unit.split_fraction[0, "effluent", "X_BA"].fix(0.5192)
        m.fs.unit.split_fraction[0, "effluent", "X_P"].fix(0.5192)
        m.fs.unit.split_fraction[0, "effluent", "S_O"].fix(0.993)
        m.fs.unit.split_fraction[0, "effluent", "S_NO"].fix(0.993)
        m.fs.unit.split_fraction[0, "effluent", "S_NH"].fix(0.993)
        m.fs.unit.split_fraction[0, "effluent", "S_ND"].fix(0.993)
        m.fs.unit.split_fraction[0, "effluent", "X_ND"].fix(0.5192)
        m.fs.unit.split_fraction[0, "effluent", "S_ALK"].fix(0.993)

        return m

    @pytest.mark.unit
    def test_dof(self, clarifier_frame):
        m = clarifier_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_units(self, clarifier_frame):
        assert_units_consistent(clarifier_frame)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, clarifier_frame):
        initialization_tester(clarifier_frame)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, clarifier_frame):
        m = clarifier_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_circular_clarifier_costing(self, clarifier_frame):
        m = clarifier_frame

        m.fs.costing = WaterTAPCosting()

        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing, costing_method=cost_circular_clarifier
        )
        m.fs.unit.surface_area.fix(1500 * units.m**2)
        m.fs.costing.cost_process()
        results = solver.solve(m)

        assert_optimal_termination(results)

        # Check solutions
        assert pytest.approx(1681573 * 2, rel=1e-5) == value(
            m.fs.unit.costing.capital_cost
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_rectangular_clarifier_costing(self, clarifier_frame):
        m = clarifier_frame

        m.fs.costing = WaterTAPCosting()

        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=cost_rectangular_clarifier,
        )
        m.fs.unit.surface_area.fix(1500 * units.m**2)
        m.fs.costing.cost_process()
        results = solver.solve(m)

        assert_optimal_termination(results)

        # Check solutions
        assert pytest.approx(2131584 * 2, rel=1e-5) == value(
            m.fs.unit.costing.capital_cost
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_primary_clarifier_costing(self, clarifier_frame):
        m = clarifier_frame

        m.fs.costing = WaterTAPCosting()

        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=cost_primary_clarifier,
        )
        m.fs.costing.cost_process()
        results = solver.solve(m)

        assert_optimal_termination(results)

        # Check solutions
        assert pytest.approx(1390570 * 2, rel=1e-5) == value(
            m.fs.unit.costing.capital_cost
        )
