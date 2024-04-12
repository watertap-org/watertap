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
    value,
)

from idaes.core import FlowsheetBlock
import watertap.property_models.water_prop_pack as props_w
from watertap.unit_models.steam_heater import FWHCondensing0D
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver
from watertap.unit_models.mvc.components.lmtd_chen_callback import (
    delta_temperature_chen_callback,
)
from idaes.models.unit_models.heat_exchanger import (
    HeatExchangerFlowPattern,
)

# Set up solver
solver = get_solver()


@pytest.mark.unit
def test_unit_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.default_property_package = props_w.WaterParameterBlock()

    m.fs.unit = FWHCondensing0D(
        hot_side_name="hot",
        cold_side_name="cold",
        hot={
            "property_package": m.fs.default_property_package,
            "has_pressure_change": False,
        },
        cold={
            "property_package": m.fs.default_property_package,
            "has_pressure_change": False,
        },
        delta_temperature_callback=delta_temperature_chen_callback,
        flow_pattern=HeatExchangerFlowPattern.countercurrent,
    )

    m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].set_value(6)
    m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(0)
    m.fs.unit.hot_side_inlet.temperature.fix(273.15 + 140)
    m.fs.unit.cold_side_inlet.pressure.fix(101325)
    m.fs.unit.cold_side_inlet.temperature.fix(273.15 + 70)
    m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    m.fs.unit.cold_side_outlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    m.fs.unit.cold_side_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(6)
    m.fs.unit.area.fix(50)
    m.fs.unit.overall_heat_transfer_coefficient.fix(2e3)

    assert degrees_of_freedom(m) == 0

    m.fs.unit.initialize()

    # assert degrees_of_freedom(m) == 0

    assert pytest.approx(0.687041189050, rel=1e-5) == value(
        m.fs.unit.hot_side_inlet.flow_mass_phase_comp[0, "Vap", "H2O"]
    )
