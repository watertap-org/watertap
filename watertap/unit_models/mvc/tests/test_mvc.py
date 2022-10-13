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

from pyomo.environ import (
    ConcreteModel,
    TransformationFactory,
    assert_optimal_termination,
)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.util.initialization import propagate_state

# Import components
from watertap.unit_models.mvc.components import Evaporator
from watertap.unit_models.mvc.components import Compressor
from watertap.unit_models.mvc.components import Condenser

# Import property packages
import watertap.property_models.seawater_prop_pack as props_sw
import watertap.property_models.water_prop_pack as props_w


def build(m):
    # Properties
    m.fs.properties_feed = props_sw.SeawaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()

    # Evaporator
    m.fs.evaporator = Evaporator(
        property_package_feed=m.fs.properties_feed,
        property_package_vapor=m.fs.properties_vapor,
    )
    # Compressor
    m.fs.compressor = Compressor(property_package=m.fs.properties_vapor)

    # Condenser
    m.fs.condenser = Condenser(property_package=m.fs.properties_vapor)

    # Connections
    m.fs.s01 = Arc(
        source=m.fs.evaporator.outlet_vapor, destination=m.fs.compressor.inlet
    )
    m.fs.s02 = Arc(source=m.fs.compressor.outlet, destination=m.fs.condenser.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)
    m.fs.evaporator.connect_to_condenser(m.fs.condenser)


def scale(m):
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

    # Evaporator
    iscale.set_scaling_factor(m.fs.evaporator.area, 1e-3)
    iscale.set_scaling_factor(m.fs.evaporator.U, 1e-3)
    iscale.set_scaling_factor(m.fs.evaporator.delta_temperature_in, 1e-1)
    iscale.set_scaling_factor(m.fs.evaporator.delta_temperature_out, 1e-1)
    iscale.set_scaling_factor(m.fs.evaporator.lmtd, 1e-1)

    # Compressor
    iscale.set_scaling_factor(m.fs.compressor.control_volume.work, 1e-6)

    # Condenser
    iscale.set_scaling_factor(m.fs.condenser.control_volume.heat, 1e-6)

    iscale.calculate_scaling_factors(m)


def specify(m):
    # Feed inlet
    m.fs.evaporator.inlet_feed.flow_mass_phase_comp[0, "Liq", "H2O"].fix(10)
    m.fs.evaporator.inlet_feed.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0.05)
    m.fs.evaporator.inlet_feed.temperature[0].fix(273.15 + 50.52)  # K
    m.fs.evaporator.inlet_feed.pressure[0].fix(1e5)  # Pa

    # Evaporator
    # m.fs.evaporator.outlet_vapor.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0.5)
    m.fs.evaporator.outlet_brine.temperature[0].fix(273.15 + 60)
    m.fs.evaporator.U.fix(1e3)  # W/K-m^2
    m.fs.evaporator.area.fix(400)  # m^2

    # Compressor
    m.fs.compressor.pressure_ratio = 2
    m.fs.compressor.control_volume.work.fix(5.8521e05)
    m.fs.compressor.efficiency.fix(0.8)


def initialize(m, solver=None):
    if solver is None:
        solver = get_solver()
    optarg = solver.options

    # initialize evaporator
    m.fs.evaporator.initialize_build(
        delta_temperature_in=30, delta_temperature_out=5, outlvl=idaeslog.INFO_HIGH
    )

    # initialize compressor
    propagate_state(m.fs.s01)
    m.fs.compressor.initialize(outlvl=idaeslog.INFO_HIGH)

    # initialize condenser
    propagate_state(m.fs.s02)
    m.fs.condenser.initialize_build(heat=-m.fs.evaporator.heat_transfer.value)


@pytest.mark.requires_idaes_solver
@pytest.mark.component
def test_mvc():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    build(m)
    assert_units_consistent(m)
    scale(m)
    specify(m)
    assert degrees_of_freedom(m) == 0

    solver = get_solver()
    initialize(m, solver=solver)

    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)

    m.fs.compressor.report()
    m.fs.condenser.report()
    m.fs.evaporator.display()
    brine_blk = m.fs.evaporator.properties_brine[0]
    # evaporator values
    assert brine_blk.pressure.value == pytest.approx(1.9849e4, rel=1e-3)
    assert m.fs.evaporator.lmtd.value == pytest.approx(30.44, rel=1e-3)
    assert m.fs.evaporator.heat_transfer.value == pytest.approx(1.2176e7, rel=1e-3)

    # compressor values
    compressed_blk = m.fs.compressor.control_volume.properties_out[0]
    assert m.fs.compressor.control_volume.work[0].value == pytest.approx(
        5.8521e5, rel=1e-3
    )
    assert compressed_blk.pressure.value == pytest.approx(3.9720e4, rel=1e-3)
    assert compressed_blk.temperature.value == pytest.approx(407.70, rel=1e-3)

    # condenser values
    condensed_blk = m.fs.condenser.control_volume.properties_out[0]
    assert m.fs.condenser.control_volume.heat[0].value == pytest.approx(
        -1.2176e7, rel=1e-3
    )
    assert condensed_blk.temperature.value == pytest.approx(342.20, rel=1e-3)
