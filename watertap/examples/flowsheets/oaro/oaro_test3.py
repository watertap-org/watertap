#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
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
import idaes.logger as idaeslog
from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    Constraint,
    assert_optimal_termination,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.network import Port
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    FlowDirection,
)
from watertap.unit_models.osmotically_assisted_reverse_osmosis_0D import (
    OsmoticallyAssistedReverseOsmosis0D,
)
import watertap.property_models.NaCl_prop_pack as props

from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)

from watertap.core import (
    MembraneChannel0DBlock,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
    FrictionFactor,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


def main():
    # set up solver
    # solver = get_solver()

    # build, set, and initialize
    m = build()
    # set_operating_conditions(m, number_of_stages=number_of_stages)
    # initialize_system(m, number_of_stages, solver=solver)

    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)

    # print("\n***---Optimization results---***")
    # display_system(m)
    # display_design(m)
    # if erd_type == ERDtype.pump_as_turbine:
    #     display_state(m)
    # else:
    #     pass

    return m


def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = OsmoticallyAssistedReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        pressure_change_type=PressureChangeType.calculated,
    )

    # fully specify system
    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.035
    feed_pressure = 50e5
    feed_temperature = 273.15 + 25
    membrane_area = 50
    A = 4.2e-12
    B = 1.3e-8
    pressure_atmospheric = 101325

    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)

    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        0.7719873688757197
    )
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        0.004970483059191647
    )
    m.fs.unit.permeate_inlet.pressure[0].fix(405467.2130804832)
    m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)

    m.fs.unit.area.fix(membrane_area)

    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)

    m.fs.unit.structural_parameter.fix(300e-6)

    m.fs.unit.permeate_side.channel_height.fix(0.001)
    m.fs.unit.permeate_side.spacer_porosity.fix(0.75)
    m.fs.unit.feed_side.channel_height.fix(0.002)
    m.fs.unit.feed_side.spacer_porosity.fix(0.75)
    m.fs.unit.feed_side.velocity[0, 0].fix(0.1)

    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    calculate_scaling_factors(m)

    m.fs.unit.initialize()

    return m


def display_state(m):
    print("--------state---------")

    def print_state(s, b):
        feed_flow_mass = sum(
            m.fs.feed.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "NaCl"]
        )
        flow_mass = sum(
            b.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "NaCl"]
        )
        normalized_flow_mass = flow_mass / feed_flow_mass * 100
        mass_frac_ppm = b.flow_mass_phase_comp[0, "Liq", "NaCl"].value / flow_mass * 1e3
        pressure_bar = b.pressure[0].value / 1e5
        print(
            s.ljust(20)
            + ": %.0f, %.3f g/L, %.1f bar"
            % (normalized_flow_mass, mass_frac_ppm, pressure_bar)
        )

        print_state(f"OARO feed inlet", m.fs.unit.feed_inlet)
        print_state(f"OARO permeate inlet", m.fs.unit.permeate_inlet)
        print_state(f"OARO feed outlet", m.fs.unit.feed_outlet)
        print_state(f"OARO permeate outlet", m.fs.unit.permeate_outlet)


if __name__ == "__main__":
    m = main()
