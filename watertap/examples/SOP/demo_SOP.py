from pyomo.environ import ConcreteModel, assert_optimal_termination, units, value
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom

import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

import watertap.examples.SOP.SOP_prop_pack as props
from watertap.examples.SOP.selective_oil_permeation import SelectiveOilPermeation


def detect_badly_scaled_vars(m):
    print()
    badly_scaled_var_list = list(
        iscale.badly_scaled_var_generator(m, large=1e2, small=1e-2, zero=1e-8)
    )
    if len(badly_scaled_var_list) == 0:
        print("No badly scaled variables detected")
    else:
        print("Badly scaled variables:")
        for (var, sv) in badly_scaled_var_list:
            print(f"    Name: {var.name}, value = {value(var)}, scaled value = {sv}")


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
    print("\n--- display---")
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    print("\n - display scaling factors -")
    m.fs.SOP.scaling_factor.display()
    detect_badly_scaled_vars(m)

    # display metrics
    m.fs.SOP.report()
    return m


if __name__ == "__main__":
    main()
