from pyomo.environ import ConcreteModel, assert_optimal_termination
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom

import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

import watertap.examples.custom_model_demo.simple_prop_pack as props
from watertap.examples.custom_model_demo.simple_filter import Filtration


def main():
    # create model, flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    # attach property package
    m.fs.properties = props.PropParameterBlock()
    # build the unit model
    m.fs.filter = Filtration(property_package=m.fs.properties)

    # display model
    # note that there are the recovery and removal fraction variables are on m.fs.filter
    # any variable that starts with a _ can be ignored, they are references
    # note that are three separate state blocks on the model (properties_in, properties_out, properties_waste)
    print("first display")
    m.fs.filter.display()

    # now specify the model
    # note there are 8 degrees of freedom
    print("DOF before specifying:", degrees_of_freedom(m.fs))

    # specify the feed
    m.fs.filter.properties_in[0].pressure.fix(2e5)
    m.fs.filter.properties_in[0].temperature.fix(273.15 + 25)
    m.fs.filter.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(1)
    m.fs.filter.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].fix(0.035)
    m.fs.filter.properties_in[0].flow_mass_phase_comp["Liq", "TSS"].fix(120e-6)
    # an alternative to setting the state variables at the state block is to use the port like below
    # note that the time domain 0, is now accessed with the other indices, this is the case for ports
    # m.fs.filter.inlet.pressure[0].fix(2e5)
    # m.fs.filter.inlet.temperature[0].fix(273.15 + 25)
    # m.fs.filter.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(1)
    # m.fs.filter.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].fix(0.035)
    # m.fs.filter.inlet.flow_mass_phase_comp[0, 'Liq', 'TSS'].fix(120e-6)

    # specify the recovery or removal
    m.fs.filter.removal_fraction_mass_phase_comp["Liq", "TSS"].fix(0.9)
    m.fs.filter.recovery_mass_phase_comp["Liq", "H2O"].fix(0.97)
    m.fs.filter.recovery_mass_phase_comp["Liq", "NaCl"].fix(0.97)

    # Currently the outlet pressure of the waste is unused (i.e. not used in any constraint) so it isn't counted in the
    # degrees of freedom, but if we connected the waste to another unit model then the pressure would be used.
    # So really the unit model has 9 DOF and I'm fixing the last one here.
    m.fs.filter.properties_waste[0].pressure.fix(101325)
    print("DOF after specifying:", degrees_of_freedom(m.fs))

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

    print("second display")
    m.fs.filter.display()

    return m


if __name__ == "__main__":
    main()
