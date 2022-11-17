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
import pytest
from io import StringIO

from pyomo.environ import ConcreteModel, assert_optimal_termination
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale

from watertap.unit_models.mvc.components import Condenser
import watertap.property_models.water_prop_pack as props

solver = get_solver()

# -----------------------------------------------------------------------------
@pytest.mark.component
def test_complete_condense():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.WaterParameterBlock()
    m.fs.unit = Condenser(property_package=m.fs.properties)

    # scaling
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Vap", "H2O"))
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    iscale.set_scaling_factor(m.fs.unit.control_volume.heat, 1e-6)
    iscale.calculate_scaling_factors(m)

    # state variables
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(1)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(1e-8)
    m.fs.unit.inlet.temperature[0].fix(400)  # K
    m.fs.unit.inlet.pressure[0].fix(0.5e5)  # Pa

    m.fs.unit.outlet.temperature[0].fix(340)  # K

    # solving
    assert_units_consistent(m)
    assert degrees_of_freedom(m) == 0

    m.fs.unit.initialize_build()

    solver = get_solver()
    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)

    # Check values
    distillate_blk = m.fs.unit.control_volume.properties_out[0]
    assert distillate_blk.flow_mass_phase_comp["Liq", "H2O"].value == pytest.approx(
        1.0000, rel=1e-3
    )
    assert distillate_blk.temperature.value == pytest.approx(340, rel=1e-3)
    assert distillate_blk.pressure.value == pytest.approx(50000, rel=1e-3)
    assert m.fs.unit.control_volume.heat[0].value == pytest.approx(-2.4358e6, rel=1e-3)
