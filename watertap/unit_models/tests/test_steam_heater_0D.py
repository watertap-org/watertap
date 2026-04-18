#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

from pyomo.environ import (
    ConcreteModel,
)

from idaes.core import FlowsheetBlock
from watertap.costing import WaterTAPCosting
import watertap.property_models.water_prop_pack as props_w
import watertap.property_models.seawater_prop_pack as props_sw
from watertap.unit_models.steam_heater_0D import SteamHeater0D
from watertap.core.solvers import get_solver
from watertap.unit_models.mvc.components.lmtd_chen_callback import (
    delta_temperature_chen_callback,
)
from idaes.models.unit_models.heat_exchanger import (
    HeatExchangerFlowPattern,
)
import idaes.core.util.scaling as iscale
import pytest
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

from idaes.core import UnitModelCostingBlock

solver = get_solver()


def build(estimate_cooling_water=False, test_with_steam_cost=False):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.steam_properties = props_w.WaterParameterBlock()
    m.fs.cold_side_properties = props_sw.SeawaterParameterBlock()
    m.fs.unit = SteamHeater0D(
        hot_side_name="hot",
        cold_side_name="cold",
        hot={
            "property_package": m.fs.steam_properties,
        },
        cold={
            "property_package": m.fs.cold_side_properties,
        },
        delta_temperature_callback=delta_temperature_chen_callback,
        flow_pattern=HeatExchangerFlowPattern.countercurrent,
        estimate_cooling_water=estimate_cooling_water,
    )

    if test_with_steam_cost:
        m.fs.costing = WaterTAPCosting()
        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"cost_steam_flow": True},
        )
        m.fs.costing.cost_process()
    m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0.5)
    m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(0)
    m.fs.unit.hot_side_inlet.temperature.fix(273.15 + 140)
    m.fs.unit.hot_side_inlet.pressure[0].fix(201325)
    m.fs.unit.cold_side_inlet.pressure.fix(101325)
    m.fs.unit.cold_side_inlet.temperature.fix(273.15 + 25)
    total_flow_mass = 10
    TDS_mass_frac = 0.035
    m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        TDS_mass_frac * total_flow_mass
    )
    m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        (1 - TDS_mass_frac) * total_flow_mass
    )
    m.fs.unit.area.fix(10)
    m.fs.unit.overall_heat_transfer_coefficient.fix(2e3)
    iscale.calculate_scaling_factors(m.fs.unit)
    return m


@pytest.mark.requires_idaes_solver
class TestSteamHeater0D(UnitTestHarness):
    def configure(self):
        m = build()

        self.unit_solutions[
            m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"]
        ] = 0.5
        self.unit_solutions[
            m.fs.unit.hot_side_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.5
        self.unit_solutions[m.fs.unit.hot_side_outlet.temperature[0]] = 341.356646567839
        self.unit_solutions[m.fs.unit.cold_side_outlet.temperature[0]] = (
            329.521210682437
        )

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"]
                + m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "TDS"],
                "out": m.fs.unit.hot_side_outlet.flow_mass_phase_comp[0, "Vap", "H2O"]
                + m.fs.unit.hot_side_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_side_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.cold_side_outlet.flow_mass_phase_comp[0, "Liq", "H2O"],
            },
        }
        return m


@pytest.mark.requires_idaes_solver
class TestSteamHeater0D_with_costing(UnitTestHarness):
    def configure(self):
        m = build(test_with_steam_cost=True)
        m.fs.costing.initialize()
        self.unit_solutions[
            m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"]
        ] = 0.5
        self.unit_solutions[
            m.fs.unit.hot_side_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.5
        self.unit_solutions[m.fs.unit.hot_side_outlet.temperature[0]] = 341.356646567839
        self.unit_solutions[m.fs.unit.cold_side_outlet.temperature[0]] = (
            329.521210682437
        )
        # just checks these are properly constructed
        self.unit_solutions[m.fs.costing.steam_cost] = 0.004
        self.unit_solutions[m.fs.costing.aggregate_flow_steam] = 0.299398766101876
        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"]
                + m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "TDS"],
                "out": m.fs.unit.hot_side_outlet.flow_mass_phase_comp[0, "Vap", "H2O"]
                + m.fs.unit.hot_side_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_side_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.cold_side_outlet.flow_mass_phase_comp[0, "Liq", "H2O"],
            },
        }
        return m


@pytest.mark.requires_idaes_solver
class TestCondenserwithEstimation(UnitTestHarness):
    def configure(self):
        m = build(estimate_cooling_water=True)
        outlet_temperature = 340
        m.fs.unit.cold_side_outlet.temperature[0].fix(outlet_temperature)
        m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0.5)
        total_flow_mass = 12
        TDS_mass_frac = 0.035
        m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].unfix()
        m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].unfix()
        m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
            TDS_mass_frac * total_flow_mass
        )
        m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].set_value(
            (1 - TDS_mass_frac) * total_flow_mass
        )

        self.unit_solutions[m.fs.unit.hot_side_outlet.temperature[0]] = (
            347.6197595677517
        )
        self.unit_solutions[m.fs.unit.cold_side_outlet.temperature[0]] = 340
        self.unit_solutions[m.fs.unit.area] = 10
        self.unit_solutions[
            m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.26018868152521274
        self.unit_solutions[
            m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 7.173773647766579
        self.unit_solutions[
            m.fs.unit.cold_side.properties_in[0].mass_frac_phase_comp["Liq", "TDS"]
        ] = 0.035

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"]
                + m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "TDS"],
                "out": m.fs.unit.hot_side_outlet.flow_mass_phase_comp[0, "Vap", "H2O"]
                + m.fs.unit.hot_side_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_side_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.cold_side_outlet.flow_mass_phase_comp[0, "Liq", "H2O"],
            },
        }

        return m
