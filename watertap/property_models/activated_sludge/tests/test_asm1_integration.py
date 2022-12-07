#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Tests for ASM1 property and reaction packages.

Verified againts results from:

[1] J. Alex, L. Benedetti, J. Copp, K.V. Gernaey, U. Jeppsson, I. Nopens, M.N. Pons,
J.P. Steyer and P. Vanrolleghem, "Benchmark Simulation Model no. 1 (BSM1)", 2018
"""

# Some more inforation about this module
__author__ = "Andrew Lee"

import pytest

import pyomo.environ as pyo
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.models.unit_models import CSTR
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state

from watertap.property_models.activated_sludge.asm1_properties import ASM1ParameterBlock
from watertap.property_models.activated_sludge.asm1_reactions import (
    ASM1ReactionParameterBlock,
)


@pytest.mark.integration
def test_ASM1_reactor():
    m = pyo.ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = ASM1ParameterBlock()
    m.fs.rxn_props = ASM1ReactionParameterBlock(property_package=m.fs.props)

    m.fs.R1 = CSTR(property_package=m.fs.props, reaction_package=m.fs.rxn_props)
    m.fs.R2 = CSTR(property_package=m.fs.props, reaction_package=m.fs.rxn_props)
    m.fs.R3 = CSTR(property_package=m.fs.props, reaction_package=m.fs.rxn_props)
    m.fs.R4 = CSTR(property_package=m.fs.props, reaction_package=m.fs.rxn_props)
    m.fs.R5 = CSTR(property_package=m.fs.props, reaction_package=m.fs.rxn_props)

    # Link units
    m.fs.stream1 = Arc(source=m.fs.R1.outlet, destination=m.fs.R2.inlet)
    m.fs.stream2 = Arc(source=m.fs.R2.outlet, destination=m.fs.R3.inlet)
    m.fs.stream3 = Arc(source=m.fs.R3.outlet, destination=m.fs.R4.inlet)
    m.fs.stream4 = Arc(source=m.fs.R4.outlet, destination=m.fs.R5.inlet)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    assert_units_consistent(m)

    # Feed conditions based on manual mass balance of inlet and recycle streams
    m.fs.R1.inlet.flow_vol.fix(92230 * pyo.units.m**3 / pyo.units.day)
    m.fs.R1.inlet.temperature.fix(298.15 * pyo.units.K)
    m.fs.R1.inlet.pressure.fix(1 * pyo.units.atm)
    m.fs.R1.inlet.conc_mass_comp[0, "S_I"].fix(30 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "S_S"].fix(14.6112 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "X_I"].fix(1149 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "X_S"].fix(89.324 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "X_BH"].fix(
        2542.03 * pyo.units.g / pyo.units.m**3
    )
    m.fs.R1.inlet.conc_mass_comp[0, "X_BA"].fix(148.6 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "X_P"].fix(448 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "S_O"].fix(0.3928 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "S_NO"].fix(8.32 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "S_NH"].fix(7.696 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "S_ND"].fix(1.9404 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.conc_mass_comp[0, "X_ND"].fix(5.616 * pyo.units.g / pyo.units.m**3)
    m.fs.R1.inlet.alkalinity.fix(4.704 * pyo.units.mol / pyo.units.m**3)

    m.fs.R1.volume.fix(1000 * pyo.units.m**3)
    m.fs.R2.volume.fix(1000 * pyo.units.m**3)
    m.fs.R3.volume.fix(1333 * pyo.units.m**3)
    m.fs.R4.volume.fix(1333 * pyo.units.m**3)
    m.fs.R5.volume.fix(1333 * pyo.units.m**3)

    assert degrees_of_freedom(m) == 0

    # Initialize flowsheet
    m.fs.R1.initialize()
    propagate_state(m.fs.stream1)
    m.fs.R2.initialize()
    propagate_state(m.fs.stream2)
    m.fs.R3.initialize()
    propagate_state(m.fs.stream3)
    m.fs.R4.initialize()
    propagate_state(m.fs.stream4)
    m.fs.R5.initialize()

    # For aerobic reactors, need to fix the oxygen concentration in outlet
    # To do this, we also need to deactivate the constraint linking O2 from
    # the previous unit
    # Doing this before initialization will cause issues with DoF however
    m.fs.R3.outlet.conc_mass_comp[0, "S_O"].fix(1.72 * pyo.units.g / pyo.units.m**3)
    m.fs.stream2.expanded_block.conc_mass_comp_equality[0, "S_O"].deactivate()

    m.fs.R4.outlet.conc_mass_comp[0, "S_O"].fix(2.43 * pyo.units.g / pyo.units.m**3)
    m.fs.stream3.expanded_block.conc_mass_comp_equality[0, "S_O"].deactivate()
    m.fs.R5.outlet.conc_mass_comp[0, "S_O"].fix(0.491 * pyo.units.g / pyo.units.m**3)
    m.fs.stream4.expanded_block.conc_mass_comp_equality[0, "S_O"].deactivate()

    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert pyo.check_optimal_termination(results)

    # Verify results against reference
    # First reactor (anoxic)
    assert pyo.value(m.fs.R1.outlet.flow_vol[0]) == pytest.approx(1.0675, rel=1e-4)
    assert pyo.value(m.fs.R1.outlet.temperature[0]) == pytest.approx(298.15, rel=1e-4)
    assert pyo.value(m.fs.R1.outlet.pressure[0]) == pytest.approx(101325, rel=1e-4)
    assert pyo.value(m.fs.R1.outlet.conc_mass_comp[0, "S_I"]) == pytest.approx(
        30e-3, rel=1e-5
    )
    assert pyo.value(m.fs.R1.outlet.conc_mass_comp[0, "S_S"]) == pytest.approx(
        2.81e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R1.outlet.conc_mass_comp[0, "X_I"]) == pytest.approx(
        1149e-3, rel=1e-3
    )
    assert pyo.value(m.fs.R1.outlet.conc_mass_comp[0, "X_S"]) == pytest.approx(
        82.1e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R1.outlet.conc_mass_comp[0, "X_BH"]) == pytest.approx(
        2552e-3, rel=1e-3
    )
    assert pyo.value(m.fs.R1.outlet.conc_mass_comp[0, "X_BA"]) == pytest.approx(
        149e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R1.outlet.conc_mass_comp[0, "X_P"]) == pytest.approx(
        449e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R1.outlet.conc_mass_comp[0, "S_O"]) == pytest.approx(
        4.3e-6, rel=1e-2
    )
    assert pyo.value(m.fs.R1.outlet.conc_mass_comp[0, "S_NO"]) == pytest.approx(
        5.36e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R1.outlet.conc_mass_comp[0, "S_NH"]) == pytest.approx(
        7.92e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R1.outlet.conc_mass_comp[0, "S_ND"]) == pytest.approx(
        1.22e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R1.outlet.conc_mass_comp[0, "X_ND"]) == pytest.approx(
        5.29e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R1.outlet.alkalinity[0]) == pytest.approx(4.93e-3, rel=1e-2)

    # Second reactor (anoixic)
    assert pyo.value(m.fs.R2.outlet.flow_vol[0]) == pytest.approx(1.0675, rel=1e-4)
    assert pyo.value(m.fs.R2.outlet.temperature[0]) == pytest.approx(298.15, rel=1e-4)
    assert pyo.value(m.fs.R2.outlet.pressure[0]) == pytest.approx(101325, rel=1e-4)
    assert pyo.value(m.fs.R2.outlet.conc_mass_comp[0, "S_I"]) == pytest.approx(
        30e-3, rel=1e-5
    )
    assert pyo.value(m.fs.R2.outlet.conc_mass_comp[0, "S_S"]) == pytest.approx(
        1.46e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R2.outlet.conc_mass_comp[0, "X_I"]) == pytest.approx(
        1149e-3, rel=1e-3
    )
    assert pyo.value(m.fs.R2.outlet.conc_mass_comp[0, "X_S"]) == pytest.approx(
        76.4e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R2.outlet.conc_mass_comp[0, "X_BH"]) == pytest.approx(
        2553e-3, rel=1e-3
    )
    assert pyo.value(m.fs.R2.outlet.conc_mass_comp[0, "X_BA"]) == pytest.approx(
        148e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R2.outlet.conc_mass_comp[0, "X_P"]) == pytest.approx(
        449e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R2.outlet.conc_mass_comp[0, "S_O"]) == pytest.approx(
        6.31e-8, rel=1e-2
    )
    assert pyo.value(m.fs.R2.outlet.conc_mass_comp[0, "S_NO"]) == pytest.approx(
        3.65e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R2.outlet.conc_mass_comp[0, "S_NH"]) == pytest.approx(
        8.34e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R2.outlet.conc_mass_comp[0, "S_ND"]) == pytest.approx(
        0.882e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R2.outlet.conc_mass_comp[0, "X_ND"]) == pytest.approx(
        5.03e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R2.outlet.alkalinity[0]) == pytest.approx(5.08e-3, rel=1e-2)

    # Third reactor (aerobic)
    assert pyo.value(m.fs.R3.outlet.flow_vol[0]) == pytest.approx(1.0675, rel=1e-4)
    assert pyo.value(m.fs.R3.outlet.temperature[0]) == pytest.approx(298.15, rel=1e-4)
    assert pyo.value(m.fs.R3.outlet.pressure[0]) == pytest.approx(101325, rel=1e-4)
    assert pyo.value(m.fs.R3.outlet.conc_mass_comp[0, "S_I"]) == pytest.approx(
        30e-3, rel=1e-5
    )
    assert pyo.value(m.fs.R3.outlet.conc_mass_comp[0, "S_S"]) == pytest.approx(
        1.15e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R3.outlet.conc_mass_comp[0, "X_I"]) == pytest.approx(
        1149e-3, rel=1e-3
    )
    assert pyo.value(m.fs.R3.outlet.conc_mass_comp[0, "X_S"]) == pytest.approx(
        64.9e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R3.outlet.conc_mass_comp[0, "X_BH"]) == pytest.approx(
        2557e-3, rel=1e-3
    )
    assert pyo.value(m.fs.R3.outlet.conc_mass_comp[0, "X_BA"]) == pytest.approx(
        149e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R3.outlet.conc_mass_comp[0, "X_P"]) == pytest.approx(
        450e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R3.outlet.conc_mass_comp[0, "S_O"]) == pytest.approx(
        1.72e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R3.outlet.conc_mass_comp[0, "S_NO"]) == pytest.approx(
        6.54e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R3.outlet.conc_mass_comp[0, "S_NH"]) == pytest.approx(
        5.55e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R3.outlet.conc_mass_comp[0, "S_ND"]) == pytest.approx(
        0.829e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R3.outlet.conc_mass_comp[0, "X_ND"]) == pytest.approx(
        4.39e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R3.outlet.alkalinity[0]) == pytest.approx(4.67e-3, rel=1e-2)

    # Fourth reactor (aerobic)
    assert pyo.value(m.fs.R4.outlet.flow_vol[0]) == pytest.approx(1.0675, rel=1e-4)
    assert pyo.value(m.fs.R4.outlet.temperature[0]) == pytest.approx(298.15, rel=1e-4)
    assert pyo.value(m.fs.R4.outlet.pressure[0]) == pytest.approx(101325, rel=1e-4)
    assert pyo.value(m.fs.R4.outlet.conc_mass_comp[0, "S_I"]) == pytest.approx(
        30e-3, rel=1e-5
    )
    assert pyo.value(m.fs.R4.outlet.conc_mass_comp[0, "S_S"]) == pytest.approx(
        0.995e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R4.outlet.conc_mass_comp[0, "X_I"]) == pytest.approx(
        1149e-3, rel=1e-3
    )
    assert pyo.value(m.fs.R4.outlet.conc_mass_comp[0, "X_S"]) == pytest.approx(
        55.7e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R4.outlet.conc_mass_comp[0, "X_BH"]) == pytest.approx(
        2559e-3, rel=1e-3
    )
    assert pyo.value(m.fs.R4.outlet.conc_mass_comp[0, "X_BA"]) == pytest.approx(
        150e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R4.outlet.conc_mass_comp[0, "X_P"]) == pytest.approx(
        451e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R4.outlet.conc_mass_comp[0, "S_O"]) == pytest.approx(
        2.43e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R4.outlet.conc_mass_comp[0, "S_NO"]) == pytest.approx(
        9.30e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R4.outlet.conc_mass_comp[0, "S_NH"]) == pytest.approx(
        2.97e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R4.outlet.conc_mass_comp[0, "S_ND"]) == pytest.approx(
        0.767e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R4.outlet.conc_mass_comp[0, "X_ND"]) == pytest.approx(
        3.88e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R4.outlet.alkalinity[0]) == pytest.approx(4.29e-3, rel=1e-2)

    # Fifth reactor (aerobic)
    assert pyo.value(m.fs.R5.outlet.flow_vol[0]) == pytest.approx(1.0675, rel=1e-4)
    assert pyo.value(m.fs.R5.outlet.temperature[0]) == pytest.approx(298.15, rel=1e-4)
    assert pyo.value(m.fs.R5.outlet.pressure[0]) == pytest.approx(101325, rel=1e-4)
    assert pyo.value(m.fs.R5.outlet.conc_mass_comp[0, "S_I"]) == pytest.approx(
        30e-3, rel=1e-5
    )
    assert pyo.value(m.fs.R5.outlet.conc_mass_comp[0, "S_S"]) == pytest.approx(
        0.889e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R5.outlet.conc_mass_comp[0, "X_I"]) == pytest.approx(
        1149e-3, rel=1e-3
    )
    assert pyo.value(m.fs.R5.outlet.conc_mass_comp[0, "X_S"]) == pytest.approx(
        49.3e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R5.outlet.conc_mass_comp[0, "X_BH"]) == pytest.approx(
        2559e-3, rel=1e-3
    )
    assert pyo.value(m.fs.R5.outlet.conc_mass_comp[0, "X_BA"]) == pytest.approx(
        150e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R5.outlet.conc_mass_comp[0, "X_P"]) == pytest.approx(
        452e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R5.outlet.conc_mass_comp[0, "S_O"]) == pytest.approx(
        0.491e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R5.outlet.conc_mass_comp[0, "S_NO"]) == pytest.approx(
        10.4e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R5.outlet.conc_mass_comp[0, "S_NH"]) == pytest.approx(
        1.73e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R5.outlet.conc_mass_comp[0, "S_ND"]) == pytest.approx(
        0.688e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R5.outlet.conc_mass_comp[0, "X_ND"]) == pytest.approx(
        3.53e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R5.outlet.alkalinity[0]) == pytest.approx(4.13e-3, rel=1e-2)
