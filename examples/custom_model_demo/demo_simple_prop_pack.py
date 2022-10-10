from pyomo.environ import ConcreteModel, assert_optimal_termination
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

import watertap.examples.custom_model_demo.simple_prop_pack as props


def main():
    # create model, flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    # attach property package
    m.fs.properties = props.PropParameterBlock()
    # build a state block, must specify a time which by convention for steady state models is just 0
    m.fs.stream = m.fs.properties.build_state_block([0])

    # display the state block, it only has the state variables and they are all unfixed
    print("\n---first display---")
    m.fs.stream[0].display()

    # attempt to access properties so that they are built
    m.fs.stream[0].mass_frac_phase_comp
    # after touching the property, the state block automatically builds it,
    # note the mass_frac_phase_comp variable and the constraint to calculate it
    print("\n---second display---")
    m.fs.stream[0].display()

    # touch another property
    m.fs.stream[0].conc_mass_phase_comp
    # after touching this property, the state block automatically builds it AND any other properties that are necessary,
    # note that now there is the conc_mass_phase_comp and dens_mass_phase variable and associated constraints
    print("\n---third display---")
    m.fs.stream[0].display()

    # touch another property
    m.fs.stream[0].flow_vol_phase

    # now that we have a state block, we can fix the state variables and solve for the properties
    m.fs.stream[0].temperature.fix(273.15 + 25)
    m.fs.stream[0].pressure.fix(101325)
    m.fs.stream[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
    m.fs.stream[0].flow_mass_phase_comp["Liq", "NaCl"].fix(0.035)
    m.fs.stream[0].flow_mass_phase_comp["Liq", "TSS"].fix(120e-6)

    # the user should provide the scale for the flow rate, so that our tools can ensure the model is well scaled
    # generally scaling factors should be such that if it is multiplied by the variable it will range between 0.01 and 100
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e4, index=("Liq", "TSS")
    )
    iscale.calculate_scaling_factors(m.fs)  # this utility scales the model

    # solving
    assert_units_consistent(m)  # check that units are consistent
    assert (
        degrees_of_freedom(m) == 0
    )  # check that the degrees of freedom are what we expect

    solver = get_solver()

    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)

    # display results
    print("\n---fourth display---")
    m.fs.stream[0].display()
    # note that the properties are solved, and the body of the constraints are small (residual)

    # equation oriented modeling has several advantages, one of them is that we can unfix variables and fix others
    # instead of setting the mass flow rates, we can set the volumetric flow rate and mass fractions
    m.fs.stream[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    m.fs.stream[0].flow_mass_phase_comp["Liq", "NaCl"].unfix()
    m.fs.stream[0].flow_mass_phase_comp["Liq", "TSS"].unfix()

    m.fs.stream[0].flow_vol_phase["Liq"].fix(1.5e-3)
    m.fs.stream[0].mass_frac_phase_comp["Liq", "NaCl"].fix(0.05)
    m.fs.stream[0].mass_frac_phase_comp["Liq", "TSS"].fix(80e-6)

    # resolve
    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)

    print("\n---fifth display---")
    m.fs.stream[0].display()

    return m


if __name__ == "__main__":
    main()
