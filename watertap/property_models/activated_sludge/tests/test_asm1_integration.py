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
Tests for ASM1 proeprty and reaction packages.

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

    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.props = ASM1ParameterBlock()
    m.fs.rxn_props = ASM1ReactionParameterBlock(
        default={"property_package": m.fs.props}
    )

    m.fs.R1 = CSTR(
        default={
            "property_package": m.fs.props,
            "reaction_package": m.fs.rxn_props,
        }
    )
    m.fs.R2 = CSTR(
        default={
            "property_package": m.fs.props,
            "reaction_package": m.fs.rxn_props,
        }
    )

    # Link units
    m.fs.stream1 = Arc(source=m.fs.R1.outlet, destination=m.fs.R2.inlet)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    assert_units_consistent(m)

    # Feed conditions based on manual mass balance of inlet and recycle streams
    m.fs.R1.inlet.flow_vol.fix(92230 * pyo.units.m**3 / pyo.units.day)
    m.fs.R1.inlet.temperature.fix(298.15 * pyo.units.K)
    m.fs.R1.inlet.pressure.fix(1 * pyo.units.atm)
    m.fs.R1.inlet.conc_mass_comp[0, "inert_soluble"].fix(
        30 * pyo.units.g / pyo.units.m**3
    )
    m.fs.R1.inlet.conc_mass_comp[0, "substrate"].fix(
        14.6112 * pyo.units.g / pyo.units.m**3
    )
    m.fs.R1.inlet.conc_mass_comp[0, "inert_particulate"].fix(
        1149 * pyo.units.g / pyo.units.m**3
    )
    m.fs.R1.inlet.conc_mass_comp[0, "biodegradable"].fix(
        89.324 * pyo.units.g / pyo.units.m**3
    )
    m.fs.R1.inlet.conc_mass_comp[0, "heterotrophic"].fix(
        2542.03 * pyo.units.g / pyo.units.m**3
    )
    m.fs.R1.inlet.conc_mass_comp[0, "autotrophic"].fix(
        148.6 * pyo.units.g / pyo.units.m**3
    )
    m.fs.R1.inlet.conc_mass_comp[0, "decay_particulate"].fix(
        448 * pyo.units.g / pyo.units.m**3
    )
    m.fs.R1.inlet.conc_mass_comp[0, "oxygen"].fix(
        0.3928 * pyo.units.g / pyo.units.m**3
    )
    m.fs.R1.inlet.conc_mass_comp[0, "nitrates"].fix(
        8.32 * pyo.units.g / pyo.units.m**3
    )
    m.fs.R1.inlet.conc_mass_comp[0, "ammonium"].fix(
        7.696 * pyo.units.g / pyo.units.m**3
    )
    m.fs.R1.inlet.conc_mass_comp[0, "nitrogen_soluble"].fix(
        1.9404 * pyo.units.g / pyo.units.m**3
    )
    m.fs.R1.inlet.conc_mass_comp[0, "nitrogen_particulate"].fix(
        5.616 * pyo.units.g / pyo.units.m**3
    )
    m.fs.R1.inlet.conc_mass_comp[0, "alkalinity"].fix(
        4.704 * pyo.units.g / pyo.units.m**3
    )

    m.fs.R1.volume.fix(1000 * pyo.units.m**3)
    m.fs.R2.volume.fix(1000 * pyo.units.m**3)

    assert degrees_of_freedom(m) == 0

    # Initialize flowsheet
    m.fs.R1.initialize()
    propagate_state(m.fs.stream1)
    m.fs.R2.initialize()

    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert pyo.check_optimal_termination(results)

    # Verify results against reference
    assert pyo.value(m.fs.R1.outlet.flow_vol[0]) == pytest.approx(1.0675, rel=1e-4)
    assert pyo.value(m.fs.R1.outlet.temperature[0]) == pytest.approx(298.15, rel=1e-4)
    assert pyo.value(m.fs.R1.outlet.pressure[0]) == pytest.approx(101325, rel=1e-4)
    assert pyo.value(
        m.fs.R1.outlet.conc_mass_comp[0, "inert_soluble"]
    ) == pytest.approx(30e-3, rel=1e-5)
    assert pyo.value(m.fs.R1.outlet.conc_mass_comp[0, "substrate"]) == pytest.approx(
        2.81e-3, rel=1e-2
    )
    assert pyo.value(
        m.fs.R1.outlet.conc_mass_comp[0, "inert_particulate"]
    ) == pytest.approx(1149e-3, rel=1e-3)
    assert pyo.value(
        m.fs.R1.outlet.conc_mass_comp[0, "biodegradable"]
    ) == pytest.approx(82.1e-3, rel=1e-2)
    assert pyo.value(
        m.fs.R1.outlet.conc_mass_comp[0, "heterotrophic"]
    ) == pytest.approx(2552e-3, rel=1e-3)
    assert pyo.value(m.fs.R1.outlet.conc_mass_comp[0, "autotrophic"]) == pytest.approx(
        149e-3, rel=1e-2
    )
    assert pyo.value(
        m.fs.R1.outlet.conc_mass_comp[0, "decay_particulate"]
    ) == pytest.approx(449e-3, rel=1e-2)
    assert pyo.value(m.fs.R1.outlet.conc_mass_comp[0, "oxygen"]) == pytest.approx(
        4.3e-6, rel=1e-2
    )
    assert pyo.value(m.fs.R1.outlet.conc_mass_comp[0, "nitrates"]) == pytest.approx(
        5.36e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R1.outlet.conc_mass_comp[0, "ammonium"]) == pytest.approx(
        7.92e-3, rel=1e-2
    )
    assert pyo.value(
        m.fs.R1.outlet.conc_mass_comp[0, "nitrogen_soluble"]
    ) == pytest.approx(1.22e-3, rel=1e-2)
    assert pyo.value(
        m.fs.R1.outlet.conc_mass_comp[0, "nitrogen_particulate"]
    ) == pytest.approx(5.29e-3, rel=1e-2)
    assert pyo.value(m.fs.R1.outlet.conc_mass_comp[0, "alkalinity"]) == pytest.approx(
        4.93e-3, rel=1e-2
    )

    assert pyo.value(m.fs.R2.outlet.flow_vol[0]) == pytest.approx(1.0675, rel=1e-4)
    assert pyo.value(m.fs.R2.outlet.temperature[0]) == pytest.approx(298.15, rel=1e-4)
    assert pyo.value(m.fs.R2.outlet.pressure[0]) == pytest.approx(101325, rel=1e-4)
    assert pyo.value(
        m.fs.R2.outlet.conc_mass_comp[0, "inert_soluble"]
    ) == pytest.approx(30e-3, rel=1e-5)
    assert pyo.value(m.fs.R2.outlet.conc_mass_comp[0, "substrate"]) == pytest.approx(
        1.46e-3, rel=1e-2
    )
    assert pyo.value(
        m.fs.R2.outlet.conc_mass_comp[0, "inert_particulate"]
    ) == pytest.approx(1149e-3, rel=1e-3)
    assert pyo.value(
        m.fs.R2.outlet.conc_mass_comp[0, "biodegradable"]
    ) == pytest.approx(76.4e-3, rel=1e-2)
    assert pyo.value(
        m.fs.R2.outlet.conc_mass_comp[0, "heterotrophic"]
    ) == pytest.approx(2553e-3, rel=1e-3)
    assert pyo.value(m.fs.R2.outlet.conc_mass_comp[0, "autotrophic"]) == pytest.approx(
        148e-3, rel=1e-2
    )
    assert pyo.value(
        m.fs.R2.outlet.conc_mass_comp[0, "decay_particulate"]
    ) == pytest.approx(449e-3, rel=1e-2)
    assert pyo.value(m.fs.R2.outlet.conc_mass_comp[0, "oxygen"]) == pytest.approx(
        6.31e-8, rel=1e-2
    )
    assert pyo.value(m.fs.R2.outlet.conc_mass_comp[0, "nitrates"]) == pytest.approx(
        3.65e-3, rel=1e-2
    )
    assert pyo.value(m.fs.R2.outlet.conc_mass_comp[0, "ammonium"]) == pytest.approx(
        8.34e-3, rel=1e-2
    )
    assert pyo.value(
        m.fs.R2.outlet.conc_mass_comp[0, "nitrogen_soluble"]
    ) == pytest.approx(0.882e-3, rel=1e-2)
    assert pyo.value(
        m.fs.R2.outlet.conc_mass_comp[0, "nitrogen_particulate"]
    ) == pytest.approx(5.03e-3, rel=1e-2)
    assert pyo.value(m.fs.R2.outlet.conc_mass_comp[0, "alkalinity"]) == pytest.approx(
        5.08e-3, rel=1e-2
    )
