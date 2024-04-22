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
from pyomo.environ import (
    ConcreteModel,
)

from idaes.core import FlowsheetBlock
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
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness


# Set up solver
solver = get_solver()


def build():
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
    )

    m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].set_value(6)
    m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(0)
    m.fs.unit.hot_side_inlet.temperature.fix(273.15 + 140)
    m.fs.unit.hot_side_inlet.pressure[0].fix(201325)
    m.fs.unit.cold_side_inlet.pressure.fix(101325)
    m.fs.unit.cold_side_inlet.temperature.fix(273.15 + 70)
    m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0.035 * 6)
    m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(6)
    m.fs.unit.area.fix(50)
    m.fs.unit.overall_heat_transfer_coefficient.fix(2e3)

    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestSteamHeater0D(UnitTestHarness):
    def configure(self):
        m = build()

        self.unit_solutions[
            m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"]
        ] = 3.044918783049
        self.unit_solutions[m.fs.unit.hot_side_outlet.temperature[0]] = 349.976261897
        return m
