from pyomo.environ import ConcreteModel, assert_optimal_termination
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

import watertap.examples.SOP.SOP_prop_pack as props

# TODO need to add vol_frac_phase_comp in here?
# TODO include flow_mass_phase in here?


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
    metadata = m.fs.properties.get_metadata().properties
    for v_name in metadata:
        getattr(m.fs.stream[0], v_name)
    print("\n---second display---")
    m.fs.stream[0].display()

    # fix state variables
    m.fs.stream[0].temperature.fix(298.15)
    m.fs.stream[0].pressure.fix(101325)
    m.fs.stream[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
    m.fs.stream[0].flow_mass_phase_comp["Liq", "oil"].fix(1e-2)

    # the user should provide the scale for the flow rate, so that our tools can ensure the model is well scaled
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
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
    m.fs.stream[0].flow_vol_phase_comp["Liq", "H2O"].fix(1e-4)
    m.fs.stream[0].flow_vol_phase_comp["Liq", "oil"].fix(1e-4)
    assert_units_consistent(m)  # check that units are consistent
    assert (
        degrees_of_freedom(m) == 0
    )  # check that the degrees of freedom are what we expect
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)

    badly_scaled_var_list = list(
        iscale.badly_scaled_var_generator(m, large=1e2, small=1e-2, zero=1e-8)
    )
    if len(badly_scaled_var_list) != 0:
        lst = []
        for (var, val) in badly_scaled_var_list:
            lst.append((var.name, val))
        print("badly scaled variables:", lst)
        assert False

    m.fs.stream[0].scaling_factor.display()
    # display resolved results
    print("\n---fourth display---")
    m.fs.stream[0].display()

    return m


if __name__ == "__main__":
    main()
