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

from pyomo.environ import ConcreteModel, assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.unit_models.mvc.components import Compressor
import watertap.property_models.water_prop_pack as props

solver = get_solver()

# -----------------------------------------------------------------------------
@pytest.mark.component
def test_compressor():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.WaterParameterBlock()
    m.fs.compressor = Compressor(property_package=m.fs.properties)
    # m.fs.compressor.control_volume.display()

    # scaling
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Vap", "H2O"))
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    iscale.set_scaling_factor(m.fs.compressor.control_volume.work, 1e-6)
    iscale.calculate_scaling_factors(m)

    # state variables
    m.fs.compressor.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(1)
    m.fs.compressor.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(1e-8)
    m.fs.compressor.inlet.temperature[0].fix(350)  # K
    m.fs.compressor.inlet.pressure[0].fix(0.5e5)  # Pa

    # specifications
    m.fs.compressor.pressure_ratio.fix(2)
    m.fs.compressor.efficiency.fix(0.8)

    # check build
    assert_units_consistent(m)
    assert degrees_of_freedom(m) == 0

    # initialize
    m.fs.compressor.initialize()

    # solve
    solver = get_solver()
    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)

    assert pytest.approx(1.1534e5, rel=1e-4) == value(
        m.fs.compressor.control_volume.work[0]
    )
    assert pytest.approx(1e-08, rel=1e-4) == value(
        m.fs.compressor.outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
    )
    assert pytest.approx(1.0, rel=1e-4) == value(
        m.fs.compressor.outlet.flow_mass_phase_comp[0, "Vap", "H2O"]
    )
    assert pytest.approx(429.57, rel=1e-4) == value(
        m.fs.compressor.outlet.temperature[0]
    )
    assert pytest.approx(1.0e5, rel=1e-4) == value(m.fs.compressor.outlet.pressure[0])

    m.fs.compressor.report()
    perf_dict = m.fs.compressor._get_performance_contents()
    assert perf_dict == {
        "vars": {
            "Pressure ratio": m.fs.compressor.pressure_ratio,
            "Efficiency": m.fs.compressor.efficiency,
            "Work": m.fs.compressor.control_volume.work[0],
        }
    }
