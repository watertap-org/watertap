#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
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

from pyomo.environ import ConcreteModel

from idaes.core import FlowsheetBlock
from idaes.models.unit_models.heat_exchanger import HeatExchangerFlowPattern
import idaes.core.util.scaling as iscale

import watertap.property_models.water_prop_pack as props_w
import watertap.property_models.seawater_prop_pack as props_sw
from watertap.unit_models.steam_heater_0D import SteamHeater0D
from watertap.unit_models.mvc.components.lmtd_chen_callback import (
    delta_temperature_chen_callback,
)
from watertap.core.solvers import get_solver
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness


solver = get_solver()


def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.steam_properties = props_w.WaterParameterBlock()
    m.fs.cold_side_properties = props_sw.SeawaterParameterBlock()

    m.fs.unit = SteamHeater0D(
        hot_side_name="hot",
        cold_side_name="cold",
        hot={"property_package": m.fs.steam_properties},
        cold={"property_package": m.fs.cold_side_properties},
        delta_temperature_callback=delta_temperature_chen_callback,
        flow_pattern=HeatExchangerFlowPattern.countercurrent,
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


@pytest.mark.requires_idaes_solver
class TestSteamHeater0D(UnitTestHarness):
    def configure(self):
        m = build()

        self.unit_solutions[
            m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"]
        ] = 0.7474967732076081
        self.unit_solutions[m.fs.unit.hot_side_outlet.temperature[0]] = (
            393.5672256519623
        )
        self.unit_solutions[m.fs.unit.cold_side_outlet.temperature[0]] = (
            340.9522230977077
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
class TestCondenserNoEstimation(UnitTestHarness):
    def configure(self):
        m = build()
        m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0.5)
        m.fs.unit.area.unfix()

        self.unit_solutions[m.fs.unit.hot_side_outlet.temperature[0]] = (
            393.5672256519623
        )
        self.unit_solutions[m.fs.unit.cold_side_outlet.temperature[0]] = (
            326.68592002256395
        )
        self.unit_solutions[m.fs.unit.area] = 6.129593261211377
        self.unit_solutions[
            m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"]
        ] = 0.5

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
        m = build()
        m.fs.unit.area.unfix()
        m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0.5)

        outlet_temperature = 340
        m.fs.unit.cold_side_outlet.temperature[0].fix(outlet_temperature)

        TDS_mass_frac = 0.035
        total_flow_mass_guess = 10
        m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].unfix()
        m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].unfix()
        m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].set_value(
            (1 - TDS_mass_frac) * total_flow_mass_guess
        )
        m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].set_value(
            TDS_mass_frac * total_flow_mass_guess
        )
        m.fs.unit.cold_side.properties_in[0].mass_frac_phase_comp["Liq", "TDS"].fix(
            TDS_mass_frac
        )

        self.unit_solutions[
            m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"]
        ] = 0.5
        self.unit_solutions[m.fs.unit.cold_side_outlet.temperature[0]] = 340
        self.unit_solutions[m.fs.unit.hot_side_outlet.temperature[0]] = (
            393.5672256519623
        )
        self.unit_solutions[
            m.fs.unit.cold_side.properties_in[0].mass_frac_phase_comp["Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.area] = 6.6472767942455455
        self.unit_solutions[
            m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.23938598366941452
        self.unit_solutions[
            m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 6.600213549742428

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
