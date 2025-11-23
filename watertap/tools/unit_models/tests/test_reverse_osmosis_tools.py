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

from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Feed

from watertap.property_models.seawater_prop_pack import (
    SeawaterParameterBlock,
    SeawaterStateBlockData,
)
from watertap.tools.unit_models import calculate_operating_pressure
from watertap.unit_models import ReverseOsmosis0D
from watertap.core.solvers import get_solver

solver = get_solver()


def build_seawater_prop_model():
    """
    Create feed model using seawater property package
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.feed.properties[0].pressure_osm_phase
    m.fs.feed.properties[0].conc_mass_phase_comp
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(0.965)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(0.035)
    m.fs.feed.properties[0].temperature.fix(273 + 25)
    m.fs.feed.properties[0].pressure.fix(101325)

    # Initialize and solve for the initial conditions
    m.fs.feed.initialize()
    results = solver.solve(m)
    print(f"Solve termination {results.solver.termination_condition}")
    # Unfix TDS mass flow
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].unfix()

    return m


m = build_seawater_prop_model()
osm = calculate_operating_pressure(
    m.fs.feed.properties[0], water_recovery_mass=0.5, solver=solver
)
print(f"Calculated osmotic pressure: {osm:.2f} Pa")
