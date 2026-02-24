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
    assert_optimal_termination,
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


@pytest.fixture
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
        finite_elements=10,
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
    m.fs.unit.feed_side.volume.fix(0.001)
    m.fs.unit.feed_side.accumulation_time.fix(5)
    for dv in m.fs.unit.feed_side.delta_state.node_dens_mass_phase:
        m.fs.unit.feed_side.delta_state.node_dens_mass_phase[dv].fix(1000)

    for mfp in m.fs.unit.feed_side.delta_state.node_mass_phase_comp:

        if "NaCl" in mfp:
            m.fs.unit.feed_side.delta_state.node_mass_frac_phase_comp[mfp].fix(
                feed_mass_frac_NaCl  # * 0.9
            )

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )

    iscale.calculate_scaling_factors(m.fs.unit)
    return m


@pytest.mark.component
def test_build(build):
    m = build
    assert_units_consistent(m)
    assert degrees_of_freedom(m) == 0


@pytest.mark.component
def test_no_holdup(build):
    m = build

    m.fs.unit.feed_side.volume.fix(0.0)
    m.fs.unit.initialize()
    results = solver.solve(m)
    assert check_optimal_termination(results)

    # check that holdup terms are zero
    for idxx in m.fs.unit.feed_side.accumulation_mass_transfer_term:
        assert pytest.approx(0, abs=1e-8) == value(
            m.fs.unit.feed_side.accumulation_mass_transfer_term[idxx]
        )


@pytest.mark.component
def test_holdup(build):
    m = build
    m.fs.unit.feed_side.volume.fix(0.01)
    m.fs.unit.initialize()
    for dv in m.fs.unit.feed_side.delta_state.node_dens_mass_phase:
        dense = m.fs.unit.feed_side.properties[0, 0].dens_mass_phase["Liq"].value
        m.fs.unit.feed_side.delta_state.node_dens_mass_phase[dv].fix(dense)

    for mfp in m.fs.unit.feed_side.delta_state.node_mass_frac_phase_comp:

        if "NaCl" in mfp:
            frac_nacl = (
                m.fs.unit.feed_side.properties[0, 0]
                .mass_frac_phase_comp["Liq", "NaCl"]
                .value
            )
            m.fs.unit.feed_side.delta_state.node_mass_frac_phase_comp[mfp].fix(
                frac_nacl  # * 0.9
            )

    results = solver.solve(m)
    assert check_optimal_termination(results)

    def get_manual_node_mass_fraction(m, inflow_node=0, outflow_node=1):
        delta_state = m.fs.unit.feed_side.delta_state
        inflow_mass = m.fs.unit.feed_side.properties[
            0, inflow_node
        ].flow_mass_phase_comp  # [ "Liq", "NaCl"]
        outflow_mass = m.fs.unit.feed_side.properties[
            0, outflow_node
        ].flow_mass_phase_comp  # [ "Liq", "NaCl"]

        mem_transfer = m.fs.unit.feed_side.mass_transfer_term
        new_node_fraction = (
            delta_state.node_mass_phase_comp[0.0, outflow_node, "Liq", "NaCl"]
            + (inflow_mass["Liq", "NaCl"] - outflow_mass["Liq", "NaCl"])
            * m.fs.unit.feed_side.accumulation_time[0]
        ) / (
            delta_state.node_mass_phase_comp[0.0, outflow_node, "Liq", "H2O"]
            + delta_state.node_mass_phase_comp[0.0, outflow_node, "Liq", "NaCl"]
            + (inflow_mass["Liq", "NaCl"] - outflow_mass["Liq", "NaCl"])
            * m.fs.unit.feed_side.accumulation_time[0]
            + (inflow_mass["Liq", "H2O"] - outflow_mass["Liq", "H2O"])
            * m.fs.unit.feed_side.accumulation_time[0]
            + mem_transfer[0.0, outflow_node, "Liq", "NaCl"]
            * m.fs.unit.feed_side.length
            / m.fs.unit.nfe
            * m.fs.unit.feed_side.accumulation_time[0]
            + mem_transfer[0.0, outflow_node, "Liq", "H2O"]
            * m.fs.unit.feed_side.length
            / m.fs.unit.nfe
            * m.fs.unit.feed_side.accumulation_time[0]
        )
        return value(new_node_fraction)

    prior_lg = 0
    for lg in m.fs.unit.feed_side.length_domain:
        print(f"Length domain {lg}:")
        if lg != 0:
            new_node_fraction = get_manual_node_mass_fraction(
                m, inflow_node=prior_lg, outflow_node=lg
            )
            actual_fraction = m.fs.unit.feed_side.properties[
                0, lg
            ].mass_frac_phase_comp[
                "Liq", "NaCl"
            ]  # [ "Liq", "NaCl"]

            print(lg, new_node_fraction, value(actual_fraction))
            # verify manually calculate mass fraction is close to actual. There is a small
            # difference as we manually calculate the dx flow mass, rather then use euler method.
            assert pytest.approx(new_node_fraction, rel=1e-3) == value(actual_fraction)
        prior_lg = lg
    assert value(m.fs.unit.feed_side.node_volume) == 0.001

    ## mass balance test!

    mass_in = m.fs.unit.feed_side.properties[0, 0].flow_mass_phase_comp
    mass_retentate = m.fs.unit.feed_side.properties[0, 1].flow_mass_phase_comp
    mass_permeate = m.fs.unit.mixed_permeate[0].flow_mass_phase_comp

    mass_start = m.fs.unit.feed_side.delta_state.node_mass_phase_comp

    mass_end = m.fs.unit.feed_side.node_mass_phase_comp

    for comp in ["H2O", "NaCl"]:
        total_mass_in = value(
            mass_in["Liq", comp] * m.fs.unit.feed_side.accumulation_time[0]
        )
        total_mass_out = value(
            mass_retentate["Liq", comp] * m.fs.unit.feed_side.accumulation_time[0]
        ) + value(mass_permeate["Liq", comp] * m.fs.unit.feed_side.accumulation_time[0])
        total_mass_start = sum(
            value(mass_start[0, idx, "Liq", comp])
            for idx in m.fs.unit.feed_side.difference_elements
        )
        total_mass_end = sum(
            value(mass_end[0, idx, "Liq", comp])
            for idx in m.fs.unit.feed_side.difference_elements
        )
        print(
            comp,
            total_mass_in,
            total_mass_out,
            total_mass_start,
            total_mass_end,
            total_mass_start - total_mass_end,
            total_mass_out - total_mass_in,
        )
        # verify that mass in - mass out == change of mass in holdup volume (e.g. conservation of mass)
        assert pytest.approx(
            total_mass_start - total_mass_end, rel=1e-5
        ) == pytest.approx(
            total_mass_out - total_mass_in,
            rel=1e-5,
        )


if __name__ == "__main__":
    m = build()
    # assert False
    # m.fs.unit.feed_side.volume.fix(0.1)
    m.fs.unit.initialize()
    for dv in m.fs.unit.feed_side.delta_state.node_dens_mass_phase:
        dense = m.fs.unit.feed_side.properties[0, 0].dens_mass_phase["Liq"].value
        m.fs.unit.feed_side.delta_state.node_dens_mass_phase[dv].fix(dense)

    for mfp in m.fs.unit.feed_side.delta_state.node_mass_frac_phase_comp:

        if "NaCl" in mfp:
            frac_nacl = (
                m.fs.unit.feed_side.properties[0, 0]
                .mass_frac_phase_comp["Liq", "NaCl"]
                .value
            )
            m.fs.unit.feed_side.delta_state.node_mass_frac_phase_comp[mfp].fix(
                frac_nacl  # * 0.9
            )

    results = solver.solve(m)
    assert_optimal_termination(results)
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
    print(value(m.fs.unit.feed_side.node_volume))
    print(value(m.fs.unit.feed_side.volume))

    # manually compute expected change in salt simply based on flow
    def get_manual_node_mass_fraction(m, inflow_node=0, outflow_node=1):
        delta_state = m.fs.unit.feed_side.delta_state
        inflow_mass = m.fs.unit.feed_side.properties[
            0, inflow_node
        ].flow_mass_phase_comp  # [ "Liq", "NaCl"]
        outflow_mass = m.fs.unit.feed_side.properties[
            0, outflow_node
        ].flow_mass_phase_comp  # [ "Liq", "NaCl"]

        mem_transfer = m.fs.unit.feed_side.mass_transfer_term
        new_node_fraction = (
            delta_state.node_mass_phase_comp[0.0, outflow_node, "Liq", "NaCl"]
            + (inflow_mass["Liq", "NaCl"] - outflow_mass["Liq", "NaCl"])
            * m.fs.unit.feed_side.accumulation_time[0]
        ) / (
            delta_state.node_mass_phase_comp[0.0, outflow_node, "Liq", "H2O"]
            + delta_state.node_mass_phase_comp[0.0, outflow_node, "Liq", "NaCl"]
            + (inflow_mass["Liq", "NaCl"] - outflow_mass["Liq", "NaCl"])
            * m.fs.unit.feed_side.accumulation_time[0]
            + (inflow_mass["Liq", "H2O"] - outflow_mass["Liq", "H2O"])
            * m.fs.unit.feed_side.accumulation_time[0]
            + mem_transfer[0.0, outflow_node, "Liq", "NaCl"]
            * m.fs.unit.feed_side.length
            / m.fs.unit.nfe
            * m.fs.unit.feed_side.accumulation_time[0]
            + mem_transfer[0.0, outflow_node, "Liq", "H2O"]
            * m.fs.unit.feed_side.length
            / m.fs.unit.nfe
            * m.fs.unit.feed_side.accumulation_time[0]
        )
        return value(new_node_fraction)

    prior_lg = 0
    for lg in m.fs.unit.feed_side.length_domain:
        print(f"Length domain {lg}:")
        if lg != 0:
            new_node_fraction = get_manual_node_mass_fraction(
                m, inflow_node=prior_lg, outflow_node=lg
            )
            actual_fraction = m.fs.unit.feed_side.properties[
                0, lg
            ].mass_frac_phase_comp[
                "Liq", "NaCl"
            ]  # [ "Liq", "NaCl"]

            print(lg, new_node_fraction, value(actual_fraction))
        prior_lg = lg
