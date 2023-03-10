from pyomo.environ import ConcreteModel, assert_optimal_termination, value
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

import watertap.examples.SOP.SOP_prop_pack as props


def detect_badly_scaled_vars(m):
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
    # build a state block, must specify a time which by convention for steady state models is just 0
    m.fs.stream = m.fs.properties.build_state_block([0])

    # display the state block, it only has the state variables and they are all unfixed
    print("\n---first display---")
    m.fs.stream[0].display()

    # attempt to access properties so that they are built
    for p in m.fs.properties.get_metadata().properties.list_supported_properties():
        getattr(m.fs.stream[0], p.name)
    print("\n---second display---")
    m.fs.stream[0].display()

    # fix state variables
    m.fs.stream[0].temperature.fix(298.15)
    m.fs.stream[0].pressure.fix(101325)
    m.fs.stream[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
    m.fs.stream[0].flow_mass_phase_comp["Liq", "oil"].fix(1e-2)

    # the user should provide the scale for the flow rate, so that our tools can ensure the model is well scaled
    iscale.set_scaling_factor(m.fs.stream[0].flow_mass_phase_comp["Liq", "H2O"], 1)
    iscale.set_scaling_factor(m.fs.stream[0].flow_mass_phase_comp["Liq", "oil"], 100)
    iscale.calculate_scaling_factors(m.fs)  # this utility scales the model

    # solving
    assert_units_consistent(m)  # check that units are consistent
    assert (
        degrees_of_freedom(m) == 0
    )  # check that the degrees of freedom are what we expect

    solver = get_solver()
    print("\n- solve model -")
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    print("\n - display scaling factors -")
    m.fs.stream[0].scaling_factor.display()
    print()
    detect_badly_scaled_vars(m)

    # display results
    print("\n---third display---")
    m.fs.stream[0].display()

    # try unfixing some variables, and fix some different ones, and then re-solve
    m.fs.stream[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    m.fs.stream[0].flow_mass_phase_comp["Liq", "oil"].unfix()
    m.fs.stream[0].flow_vol_phase_comp["Liq", "H2O"].fix(1e-4)
    m.fs.stream[0].flow_vol_phase_comp["Liq", "oil"].fix(1e-4)

    # Reset default scaling factors and recalculate all scaling factors
    iscale.set_scaling_factor(m.fs.stream[0].flow_mass_phase_comp["Liq", "H2O"], 10)
    iscale.set_scaling_factor(m.fs.stream[0].flow_mass_phase_comp["Liq", "oil"], 10)
    iscale.calculate_scaling_factors(m.fs)

    assert_units_consistent(m)  # check that units are consistent
    assert (
        degrees_of_freedom(m) == 0
    )  # check that the degrees of freedom are what we expect
    print("\n- solve model -")
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    print("\n - display scaling factors -")
    m.fs.stream[0].scaling_factor.display()
    print()
    detect_badly_scaled_vars(m)

    # display resolved results
    print("\n---fourth display---")
    m.fs.stream[0].display()

    return m


if __name__ == "__main__":
    main()
