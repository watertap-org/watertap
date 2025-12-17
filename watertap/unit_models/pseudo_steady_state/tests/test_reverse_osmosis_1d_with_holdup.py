#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
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
from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.reverse_osmosis_0D import (
    ConcentrationPolarizationType,
    MassTransferCoefficient,
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
from watertap.unit_models.pseudo_steady_state.reverse_osmosis_1D_with_holdup import (
    ReverseOsmosis1DwithHoldUp,
)

from idaes.core.util.scaling import (
    calculate_scaling_factors,
)
from idaes.core.util.model_diagnostics import DiagnosticsToolbox


import idaes.core.util.scaling as iscale

__author__ = "Alexander V. Dudchenko"
# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()
# -----------------------------------------------------------------------------


# @pytest.fixture
def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = ReverseOsmosis1DwithHoldUp(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        finite_elements=2,
        module_type="spiral_wound",
        has_full_reporting=True,
    )

    # fully specify system
    feed_flow_mass = 1000 / 3600
    feed_mass_frac_NaCl = 0.034283
    feed_pressure = 70e5

    feed_temperature = 273.15 + 25
    A = 4.2e-12
    B = 3.5e-8
    pressure_atmospheric = 1e5
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )

    m.fs.unit.inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
    m.fs.unit.feed_side.N_Re[0, 0].fix(400)
    m.fs.unit.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(0.5)
    m.fs.unit.feed_side.spacer_porosity.fix(0.97)
    m.fs.unit.feed_side.channel_height.fix(0.001)

    ## accumulation terms
    m.fs.unit.feed_side.volume.fix(0.01)
    m.fs.unit.feed_side.accumulation_time.fix(1)
    for dv in m.fs.unit.feed_side.delta_state.node_dens_mass_phase:
        print(dv)
        m.fs.unit.feed_side.delta_state.node_dens_mass_phase[dv].fix(1000)

    for mfp in m.fs.unit.feed_side.delta_state.node_mass_phase_comp:

        if "NaCl" in mfp:
            print(mfp)
            m.fs.unit.feed_side.delta_state.node_mass_frac_phase_comp[mfp].fix(
                feed_mass_frac_NaCl * 0.5
            )

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )

    iscale.calculate_scaling_factors(m.fs.unit)
    return m


# @pytest.mark.component
def test_build(build):
    m = build
    assert_units_consistent(m)
    assert degrees_of_freedom(m) == 0


if __name__ == "__main__":
    m = build()
    # assert False
    m.fs.unit.initialize()
    m.fs.unit.feed_side.display()
    m.fs.unit.feed_side.material_flow_linking_constraints.pprint()
    m.fs.unit.feed_side.material_balances.pprint()
    # m.fs.unit.feed_side.material_holdup_calculation.pprint()
    m.fs.unit.feed_side.material_flow_dx_disc_eq.pprint()
    # initialization_tester(m)
    for prop in m.fs.unit.feed_side.properties:
        print(prop)
        for index, obj in m.fs.unit.feed_side.properties[
            prop
        ].mass_frac_phase_comp.items():
            print(f"{index}: {value(obj)}")
