###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################
import sys
import pytest
from io import StringIO

from pyomo.environ import ConcreteModel, assert_optimal_termination
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.unit_models.mvc.components import Evaporator
import watertap.property_models.seawater_prop_pack as props_sw
import watertap.property_models.water_prop_pack as props_w

solver = get_solver()

# -----------------------------------------------------------------------------
@pytest.mark.requires_idaes_solver
@pytest.mark.component
def test_evaporator():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_feed = props_sw.SeawaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.evaporator = Evaporator(
        property_package_feed=m.fs.properties_feed,
        property_package_vapor=m.fs.properties_vapor,
    )

    # scaling
    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.properties_feed.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )
    m.fs.properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Vap", "H2O")
    )
    m.fs.properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    iscale.set_scaling_factor(m.fs.evaporator.area, 1e-3)
    iscale.set_scaling_factor(m.fs.evaporator.U, 1e-3)
    iscale.set_scaling_factor(m.fs.evaporator.delta_temperature_in, 1e-1)
    iscale.set_scaling_factor(m.fs.evaporator.delta_temperature_out, 1e-1)
    iscale.set_scaling_factor(m.fs.evaporator.lmtd, 1e-1)
    iscale.calculate_scaling_factors(m)

    # state variables
    # Feed inlet
    m.fs.evaporator.inlet_feed.flow_mass_phase_comp[0, "Liq", "H2O"].fix(10)
    m.fs.evaporator.inlet_feed.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0.05)
    m.fs.evaporator.inlet_feed.temperature[0].fix(273.15 + 50.52)  # K
    m.fs.evaporator.inlet_feed.pressure[0].fix(1e5)  # Pa

    # Evaporator specifications
    m.fs.evaporator.outlet_brine.temperature[0].fix(273.15 + 60)
    m.fs.evaporator.U.fix(1e3)  # W/K-m^2
    m.fs.evaporator.area.fix(100)  # m^2
    m.fs.evaporator.delta_temperature_out.fix(5)
    m.fs.evaporator.delta_temperature_in.fix(30)

    # check build
    assert_units_consistent(m)
    assert degrees_of_freedom(m) == 0

    # initialize
    m.fs.evaporator.initialize_build()

    # solve
    solver = get_solver()
    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)

    # check values, TODO: make a report for the evaporator
    vapor_blk = m.fs.evaporator.properties_vapor[0]
    assert vapor_blk.flow_mass_phase_comp["Vap", "H2O"].value == pytest.approx(
        0.4171, rel=1e-3
    )
    assert m.fs.evaporator.lmtd.value == pytest.approx(13.79, rel=1e-3)
    assert m.fs.evaporator.heat_transfer.value == pytest.approx(1.379e6, rel=1e-3)
