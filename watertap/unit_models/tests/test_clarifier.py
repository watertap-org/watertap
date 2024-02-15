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

from pyomo.environ import (
    ConcreteModel,
    units,
)

from idaes.core import (
    FlowsheetBlock,
)

from idaes.core.solvers import get_solver

from watertap.unit_models.tests.unit_test_harness import UnitTestHarness
import idaes.core.util.scaling as iscale

from idaes.core import UnitModelCostingBlock
from watertap.costing import WaterTAPCosting

from watertap.unit_models.clarifier import Clarifier
from idaes.models.unit_models.separator import SplittingType

from watertap.property_models.activated_sludge.asm1_properties import ASM1ParameterBlock
from watertap.costing.unit_models.clarifier import (
    cost_circular_clarifier,
    cost_rectangular_clarifier,
    cost_primary_clarifier,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
def build():
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
    m.fs.unit.inlet.conc_mass_comp[0, "X_BA"].fix(1e-3 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_P"].fix(1e-3 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_O"].fix(1e-3 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_NO"].fix(1e-3 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_NH"].fix(23 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "S_ND"].fix(5 * units.g / units.m**3)
    m.fs.unit.inlet.conc_mass_comp[0, "X_ND"].fix(16 * units.g / units.m**3)

    # Alkalinity was given in mg/L based on C
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

    # Set scaling factors for badly scaled variables
    iscale.set_scaling_factor(
        m.fs.unit.underflow_state[0.0].conc_mass_comp["X_BA"], 1e3
    )
    iscale.set_scaling_factor(m.fs.unit.underflow_state[0.0].conc_mass_comp["X_P"], 1e3)
    iscale.set_scaling_factor(m.fs.unit.underflow_state[0.0].conc_mass_comp["S_O"], 1e3)
    iscale.set_scaling_factor(
        m.fs.unit.underflow_state[0.0].conc_mass_comp["S_NO"], 1e3
    )
    iscale.set_scaling_factor(m.fs.unit.effluent_state[0.0].conc_mass_comp["X_BA"], 1e7)
    iscale.set_scaling_factor(m.fs.unit.effluent_state[0.0].conc_mass_comp["X_P"], 1e7)
    iscale.set_scaling_factor(m.fs.unit.effluent_state[0.0].conc_mass_comp["S_O"], 1e3)
    iscale.set_scaling_factor(m.fs.unit.effluent_state[0.0].conc_mass_comp["S_NO"], 1e3)

    return m


class TestCircularCosting(UnitTestHarness):
    def configure(self):
        m = build()

        # Add unit model costing
        m.fs.costing = WaterTAPCosting()

        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing, costing_method=cost_circular_clarifier
        )
        m.fs.unit.surface_area.fix(1500 * units.m**2)

        iscale.set_scaling_factor(m.fs.unit.costing.capital_cost, 1e-6)

        m.fs.costing.cost_process()

        self.unit_model_block = m.fs.unit

        self.unit_solutions[m.fs.unit.costing.capital_cost] = 1681573 * 2


class TestRectangularCosting(UnitTestHarness):
    def configure(self):
        m = build()

        # Add unit model costing
        m.fs.costing = WaterTAPCosting()

        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=cost_rectangular_clarifier,
        )
        m.fs.unit.surface_area.fix(1500 * units.m**2)

        iscale.set_scaling_factor(m.fs.unit.costing.capital_cost, 1e-6)

        m.fs.costing.cost_process()

        self.unit_model_block = m.fs.unit

        self.unit_solutions[m.fs.unit.costing.capital_cost] = 2131584 * 2


class TestPrimaryClarifierCosting(UnitTestHarness):
    def configure(self):
        m = build()

        # Add unit model costing
        m.fs.costing = WaterTAPCosting()

        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing, costing_method=cost_primary_clarifier
        )
        m.fs.unit.surface_area.fix(1500 * units.m**2)

        iscale.set_scaling_factor(m.fs.unit.costing.capital_cost, 1e-6)

        m.fs.costing.cost_process()

        self.unit_model_block = m.fs.unit

        self.unit_solutions[m.fs.unit.costing.capital_cost] = 1390570 * 2
