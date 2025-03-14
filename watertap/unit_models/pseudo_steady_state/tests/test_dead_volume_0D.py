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
    check_optimal_termination,
)

from idaes.core import (
    FlowsheetBlock,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
)

from pyomo.util.check_units import assert_units_consistent
import watertap.property_models.NaCl_prop_pack as props
from idaes.core.util.testing import initialization_tester
from watertap.core.solvers import get_solver
from watertap.unit_models.pseudo_steady_state.dead_volume_0D import DeadVolume0D

from idaes.core.util.scaling import (
    calculate_scaling_factors,
)
from idaes.core.util.model_diagnostics import DiagnosticsToolbox

__author__ = "Alexander V. Dudchenko"


@pytest.fixture
def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = DeadVolume0D(property_package=m.fs.properties)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(1)
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(0.01)
    m.fs.unit.inlet.temperature.fix(293)
    m.fs.unit.inlet.pressure.fix(5e5)
    m.fs.unit.volume[0, "Liq"].fix(0.0010074962222595425)

    m.fs.unit.delta_state.volume[0, "Liq"].fix(0.0010074962222595425)
    m.fs.unit.accumulation_time.fix(1)
    m.fs.unit.delta_state.dens_mass_phase[0, "Liq"].fix(1002.4851485147615)
    # m.fs.unit.delta_state.mass_frac_phase_comp[0, "Liq", "H2O"].fix(1)
    m.fs.unit.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"].fix(
        0.009900990099009898
    )
    # m.fs.unit.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"].fix(0.01)

    for key, obj in m.fs.unit.inlet.flow_mass_phase_comp.items():
        val = value(obj)
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1 / val, index=("Liq", key[-1])
        )

    calculate_scaling_factors(m)
    return m


@pytest.mark.component
def test_build(build):
    m = build
    # m.fs.unit.display()
    assert degrees_of_freedom(m) == 0
    assert assert_units_consistent(m) is None


@pytest.mark.component
def test_solve(build):
    m = build
    # m.fs.unit.display()
    assert degrees_of_freedom(m) == 0

    initialization_tester(m)

    # test that dead volume did not change outflow comp
    assert (
        pytest.approx(1, rel=1e-5)
        == m.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "H2O"].value
    )
    assert (
        pytest.approx(0.01, rel=1e-5)
        == m.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "NaCl"].value
    )

    dead_volume_water_mass = 0.9999999999999194
    dead_volume_salt_mass = 0.009999999999999183
    # Test now, but double our salt concentration
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(0.02)

    solver = get_solver()
    result = solver.solve(m, tee=True)
    assert check_optimal_termination(result)
    assert (
        pytest.approx(m.fs.unit.inlet.pressure[0].value, rel=1e-5)
        == m.fs.unit.outlet.pressure[0].value
    )

    new_fraction = (dead_volume_salt_mass + 0.02) / (
        1 + dead_volume_water_mass + dead_volume_salt_mass + 0.02
    )
    assert (
        pytest.approx(new_fraction, rel=1e-5)
        == m.fs.unit.dead_volume.properties_out[0.0]
        .mass_frac_phase_comp["Liq", "NaCl"]
        .value
    )
    m.fs.unit.display()
