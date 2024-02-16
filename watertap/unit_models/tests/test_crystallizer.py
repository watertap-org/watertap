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
from watertap.costing import WaterTAPCosting, CrystallizerCostType

from watertap.unit_models.crystallizer import Crystallization
import watertap.property_models.cryst_prop_pack as props

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = Crystallization(property_package=m.fs.properties)

    # fully specify system
    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.2126
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    feed_pressure = 101325
    feed_temperature = 273.15 + 20
    eps = 1e-6
    crystallizer_temperature = 273.15 + 55
    crystallizer_yield = 0.40

    # Fully define feed
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(eps)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(eps)
    m.fs.unit.inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)

    # Define operating conditions
    m.fs.unit.temperature_operating.fix(crystallizer_temperature)
    m.fs.unit.crystallization_yield["NaCl"].fix(crystallizer_yield)

    # Fix growth rate, crystal length and Sounders brown constant to default values
    m.fs.unit.crystal_growth_rate.fix()
    m.fs.unit.souders_brown_constant.fix()
    m.fs.unit.crystal_median_length.fix()

    iscale.set_scaling_factor(m.fs.unit.pressure_operating, 1e-3)
    iscale.set_scaling_factor(m.fs.unit.properties_out[0].pressure, 1e-5)
    iscale.set_scaling_factor(m.fs.unit.properties_solids[0].pressure, 1e-5)
    iscale.set_scaling_factor(m.fs.unit.properties_vapor[0].pressure, 1e-3)
    iscale.set_scaling_factor(m.fs.unit.properties_out[0].flow_vol_phase["Liq"], 1e3)
    iscale.set_scaling_factor(m.fs.unit.properties_out[0].flow_vol_phase["Vap"], 1e8)
    iscale.set_scaling_factor(m.fs.unit.properties_out[0].pressure_sat, 1e-3)
    iscale.set_scaling_factor(m.fs.unit.properties_solids[0].flow_vol_phase["Vap"], 1e8)
    iscale.set_scaling_factor(m.fs.unit.properties_solids[0].flow_vol_phase["Sol"], 1e8)
    iscale.set_scaling_factor(
        m.fs.unit.properties_vapor[0].dens_mass_solvent["Vap"], 1e1
    )
    iscale.set_scaling_factor(m.fs.unit.properties_out[0].flow_vol_phase["Sol"], 1e12)
    iscale.set_scaling_factor(
        m.fs.unit.properties_solids[0].flow_vol_phase["Liq"], 1e11
    )
    iscale.set_scaling_factor(m.fs.unit.properties_vapor[0].flow_vol_phase["Liq"], 1e12)
    iscale.set_scaling_factor(m.fs.unit.properties_vapor[0].flow_vol_phase["Sol"], 1e12)

    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestCrystallizer(UnitTestHarness):
    def configure(self):
        m = build()

        # Add unit model costing
        m.fs.costing = WaterTAPCosting()

        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"cost_type": CrystallizerCostType.mass_basis},
        )

        iscale.set_scaling_factor(m.fs.unit.costing.capital_cost, 1e-5)

        m.fs.costing.cost_process()

        m.fs.costing.initialize()

        self.unit_model_block = m.fs.unit

        self.unit_solutions[
            m.fs.unit.solids.flow_mass_phase_comp[0, "Sol", "NaCl"]
        ] = 0.08504096
        self.unit_solutions[
            m.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
        ] = 0.12756
        self.unit_solutions[
            m.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.34570767
        self.unit_solutions[
            m.fs.unit.properties_out[0].mass_frac_phase_comp["Liq", "NaCl"]
        ] = 0.26953035
        self.unit_solutions[m.fs.unit.pressure_operating] = 11991.82426
        self.unit_solutions[m.fs.unit.work_mechanical[0]] = 1127.222711
        self.unit_solutions[m.fs.unit.diameter_crystallizer] = 1.204904
        self.unit_solutions[m.fs.unit.volume_suspension] = 1.6194469
        self.unit_solutions[m.fs.unit.t_res] = 1.022821
        self.unit_solutions[m.fs.unit.costing.capital_cost] = 600166.82793
        self.unit_solutions[
            m.fs.costing.aggregate_flow_costs["electricity"]
        ] = 8596.50939
        self.unit_solutions[m.fs.costing.aggregate_flow_costs["steam"]] = 2720557.6394
        self.unit_solutions[m.fs.costing.aggregate_flow_costs["NaCl"]] = 0
