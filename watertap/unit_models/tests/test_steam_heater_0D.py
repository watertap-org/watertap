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

from pyomo.environ import (
    ConcreteModel,
)

from idaes.core import FlowsheetBlock
import watertap.property_models.water_prop_pack as props_w
import watertap.property_models.seawater_prop_pack as props_sw
from watertap.unit_models.steam_heater_0D import SteamHeater0D, Mode
from watertap.core.solvers import get_solver
from watertap.unit_models.mvc.components.lmtd_chen_callback import (
    delta_temperature_chen_callback,
)
from idaes.models.unit_models.heat_exchanger import (
    HeatExchangerFlowPattern,
)
import idaes.core.util.scaling as iscale
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness


solver = get_solver()


def build(mode, estimate_cooling_water=False):
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
        mode=mode,
        estimate_cooling_water=estimate_cooling_water,
    )

    m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].set_value(1)
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


class TestSteamHeater0D(UnitTestHarness):
    def configure(self):
        m = build(Mode.HEATER)

        self.unit_solutions[
            m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"]
        ] = 0.6913572208080987
        self.unit_solutions[m.fs.unit.hot_side_outlet.temperature[0]] = (
            381.3881358453978
        )
        self.unit_solutions[m.fs.unit.cold_side_outlet.temperature[0]] = (
            338.631650369243
        )
        return m


class TestCondenserNoEstimation(UnitTestHarness):
    def configure(self):
        m = build(Mode.CONDENSER, estimate_cooling_water=False)
        m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0.5)
        m.fs.unit.area.unfix()

        self.unit_solutions[m.fs.unit.hot_side_outlet.temperature[0]] = (
            375.3825125826731
        )
        self.unit_solutions[m.fs.unit.cold_side_outlet.temperature[0]] = (
            327.67799780959757
        )
        self.unit_solutions[m.fs.unit.area] = 7.089052938081363
        self.unit_solutions[
            m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"]
        ] = 0.5
        return m


class TestCondenserwithEstimation(UnitTestHarness):
    def configure(self):
        m = build(Mode.CONDENSER, estimate_cooling_water=True)
        m.fs.unit.area.unfix()

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
            373.69526509227114
        )
        self.unit_solutions[m.fs.unit.cold_side_outlet.temperature[0]] = 340
        self.unit_solutions[m.fs.unit.area] = 7.774904144204368
        self.unit_solutions[
            m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.24841897743647612
        self.unit_solutions[
            m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 6.849266092188223
        self.unit_solutions[
            m.fs.unit.cold_side.properties_in[0].mass_frac_phase_comp["Liq", "TDS"]
        ] = 0.035

        return m
