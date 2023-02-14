from pyomo.environ import ConcreteModel, assert_optimal_termination, units
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom

import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

import watertap.examples.SOP.SOP_prop_pack as props
from watertap.examples.SOP.selective_oil_permeation import SelectiveOilPermeation


def main():
    # create model, flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    # attach property package
    m.fs.properties = props.SopParameterBlock()
    # build the unit model
    m.fs.SOP = SelectiveOilPermeation(
        property_package=m.fs.properties, has_pressure_change=True
    )

    # scale model
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "oil"))
    iscale.calculate_scaling_factors(m.fs)

    # print DOF before specifying
    print("DOF before specifying:", degrees_of_freedom(m.fs))

    # specify model
    # feed
    m.fs.SOP.feed_side.properties_in[0].temperature.fix(298.15)
    m.fs.SOP.feed_side.properties_in[0].pressure.fix(2 * units.bar)
    m.fs.SOP.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
    m.fs.SOP.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "oil"].fix(0.2)
    # unit model
    m.fs.SOP.mass_transfer_oil.fix(0.15)
    m.fs.SOP.area.fix(10)
    m.fs.SOP.properties_permeate[0].pressure.fix(1 * units.bar)
    m.fs.SOP.deltaP.fix(-0.5 * units.bar)

    print("DOF after specifying:", degrees_of_freedom(m.fs))

    # solve model
    assert_units_consistent(m)  # check that units are consistent
    assert (
        degrees_of_freedom(m) == 0
    )  # check that the degrees of freedom are what we expect

    m.fs.SOP.initialize()
    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)

    # display metrics
    m.fs.SOP.report()
    return m


if __name__ == "__main__":
    main()
