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
from pyomo.environ import (
    ConcreteModel,
    value,
    Constraint,
    Expression,
    Objective,
    Param,
    TransformationFactory,
    units,
    assert_optimal_termination,
)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Mixer, Separator, Product, Feed
from idaes.core import UnitModelCostingBlock
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

import watertap.property_models.seawater_prop_pack as properties


def main():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = properties.SeawaterParameterBlock()
    m.fs.state_block = m.fs.properties.build_state_block([0], default={})

    # display the state block, it only has the state variables
    m.fs.state_block[0].display()

    # set the state variables
    m.fs.state_block[0].flow_mass_phase_comp["Liq", "H2O"].fix(0.965)
    m.fs.state_block[0].flow_mass_phase_comp["Liq", "TDS"].fix(0.035)
    m.fs.state_block[0].temperature.fix(273 + 25)
    m.fs.state_block[0].pressure.fix(101325)
    m.fs.state_block[0].display()

    # supported properties are located in the metadata
    metadata = m.fs.properties.get_metadata().properties
    for variable_name in metadata:
        print(variable_name)

    # ***include a picture of the properties table*** from
    # https://watertap.readthedocs.io/en/latest/technical_reference/property_models/seawater.html

    # properties are created on demand for the state block
    m.fs.state_block[0].mass_frac_phase_comp
    m.fs.state_block[0].display()

    # variable and constraint are created for the property, but not yet solved
    solver = get_solver()
    solver.solve(m.fs.state_block[0])
    m.fs.state_block[0].display()

    # *** this could be hidden *** it's also might not be necessary
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )
    iscale.calculate_scaling_factors(m)

    # determine the osmotic pressure
    m.fs.state_block[0].pressure_osm
    solver.solve(m.fs.state_block[0])
    m.fs.state_block[0].display()

    # print value and change units
    print(
        "Osmotic pressure: %.1f bar"
        % value(units.convert(m.fs.state_block[0].pressure_osm, to_units=units.bar))
    )

    # equation oriented modeling can solve using different variables
    m.fs.state_block[0].flow_vol_phase["Liq"].fix(1 * (units.m**3 / units.hr))
    m.fs.state_block[0].mass_frac_phase_comp["Liq", "TDS"].fix(0.05)

    m.fs.state_block[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    m.fs.state_block[0].flow_mass_phase_comp["Liq", "TDS"].unfix()

    solver.solve(m.fs.state_block[0])
    m.fs.state_block[0].display()
    print(
        "Osmotic pressure: %.1f bar"
        % value(units.convert(m.fs.state_block[0].pressure_osm, to_units=units.bar))
    )

    # another example
    m.fs.state_block[0].conc_mass_phase_comp["Liq", "TDS"].fix(100 * units.g / units.L)
    m.fs.state_block[0].mass_frac_phase_comp["Liq", "TDS"].unfix()
    solver.solve(m.fs.state_block[0])
    print(
        "Osmotic pressure: %.1f bar"
        % value(units.convert(m.fs.state_block[0].pressure_osm, to_units=units.bar))
    )

    # even can specify the property and solve for the state
    m.fs.state_block[0].pressure_osm.fix(65 * units.bar)
    m.fs.state_block[0].conc_mass_phase_comp["Liq", "TDS"].unfix()

    solver.solve(m.fs.state_block[0])
    m.fs.state_block[0].display()

    # make graphs through simulating multiple points
    m.fs.state_block[0].pressure_osm.unfix()
    conc_mass = []
    pressure_osm = []
    for i in range(0, 251):
        m.fs.state_block[0].conc_mass_phase_comp["Liq", "TDS"].fix(i)
        solver.solve(m.fs.state_block[0])
        conc_mass.append(value(m.fs.state_block[0].conc_mass_phase_comp["Liq", "TDS"]))
        pressure_osm.append(
            value(units.convert(m.fs.state_block[0].pressure_osm, to_units=units.bar))
        )

    return conc_mass, pressure_osm


if __name__ == "__main__":
    c, p = main()
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 8))
    plt.rc("font", size=18)
    plt.plot(c, p, "red", lw=2.0)
    plt.xlabel("Concentration (g/L)")
    plt.ylabel("Osmotic Pressure (bar)")
    plt.xlim(0, 250)
    plt.ylim(0, 250)
