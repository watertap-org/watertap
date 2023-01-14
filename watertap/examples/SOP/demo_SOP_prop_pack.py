from pyomo.environ import ConcreteModel, assert_optimal_termination
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

import watertap.examples.SOP.SOP_prop_pack as props


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
    m.fs.stream[0].visc_d_phase_comp
    m.fs.stream[0].mass_frac_phase_comp
    m.fs.stream[0].dens_mass_phase
    m.fs.stream[0].flow_vol_phase
    m.fs.stream[0].conc_mass_phase_comp
    print("\n---second display---")
    m.fs.stream[0].display()

    # fix state variables
    m.fs.stream[0].temperature.fix(298.15)
    m.fs.stream[0].pressure.fix(101325)
    m.fs.stream[0].flow_mass_phase_comp["Liq", "H2O"].fix(0.05)
    m.fs.stream[0].flow_mass_phase_comp["Liq", "oil"].fix(0.05)

    # the user should provide the scale for the flow rate, so that our tools can ensure the model is well scaled
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "oil")
    )
    iscale.calculate_scaling_factors(m.fs)  # this utility scales the model

    # solving
    assert_units_consistent(m)  # check that units are consistent
    assert (
        degrees_of_freedom(m) == 0
    )  # check that the degrees of freedom are what we expect

    solver = get_solver()
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)

    # display results
    print("\n---third display---")
    m.fs.stream[0].display()

    # try unfixing some variables, and fix some different ones, and then resolve
    m.fs.stream[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    m.fs.stream[0].flow_mass_phase_comp["Liq", "oil"].unfix()
    m.fs.stream[0].flow_vol_phase["Liq"].fix(1e-4)
    m.fs.stream[0].mass_frac_phase_comp["Liq", "oil"].fix(0.5)
    assert_units_consistent(m)  # check that units are consistent
    assert (
        degrees_of_freedom(m) == 0
    )  # check that the degrees of freedom are what we expect
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)

    # display resolved results
    print("\n---fourth display---")
    m.fs.stream[0].display()

    return m


if __name__ == "__main__":
    main()
