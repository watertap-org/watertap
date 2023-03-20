from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    units,
    value,
    Objective,
    SolverFactory,
)
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DegeneracyHunter

import watertap.property_models.SOP_prop_pack as props
from watertap.unit_models.selective_oil_permeation import SelectiveOilPermeation


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
    # build the unit model
    m.fs.SOP = SelectiveOilPermeation(property_package=m.fs.properties)

    # print DOF before specifying
    print("DOF before specifying:", degrees_of_freedom(m.fs))

    # specify model
    m.fs.SOP.feed_side.properties_in[0].temperature.fix(298.15)  # temp in K
    m.fs.SOP.feed_side.properties_in[0].pressure.fix(1.48 * units.bar)
    # The below feed data corresponds roughly to 3.8 L/min flow with 500 ppm oil
    m.fs.SOP.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        6.3e-2
    )  # H2O flow in kg/s
    m.fs.SOP.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "oil"].fix(
        3.2e-5
    )  # oil flow in kg/s
    m.fs.SOP.feed_side.properties_out[0].pressure.fix(1.20 * units.bar)
    m.fs.SOP.area.fix(1.4)  # membrane area in m^2
    m.fs.SOP.properties_permeate[0].pressure.fix(1 * units.bar)
    print("DOF after specifying:", degrees_of_freedom(m.fs))
    print()

    # scale model
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e5, index=("Liq", "oil")
    )
    iscale.calculate_scaling_factors(m)
    print(
        "\n------- display badly scaled variables - after initialize, before solve -------"
    )
    detect_badly_scaled_vars(m)
    # solve model
    assert_units_consistent(m)  # check that units are consistent

    assert (
        degrees_of_freedom(m) == 0
    )  # check that the degrees of freedom are what we expect
    solver = get_solver()

    # Use of Degeneracy Hunter for troubleshooting model.
    m.fs.dummy_objective = Objective(expr=0)
    solver.options["max_iter"] = 0
    solver.solve(m, tee=True)
    dh = DegeneracyHunter(m, solver=SolverFactory("cbc"))
    dh.check_residuals(tol=0.1)

    print("------- initialize model -------")
    solver.options["max_iter"] = 2000
    m.fs.SOP.initialize()

    print("\n------- solve model -------")
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)

    print("\n------- display scaling factors -------")
    m.fs.SOP.scaling_factor.display()

    print("\n------- display badly scaled variables - after solve -------")
    detect_badly_scaled_vars(m)  # TODO do I still want this here?

    # display metrics
    print("\n------- display report -------")
    m.fs.SOP.report()

    return m


if __name__ == "__main__":
    m = main()
