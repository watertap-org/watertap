from pyomo.environ import ConcreteModel, assert_optimal_termination
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

    # display model
    print(
        "--------------------------------- first display ---------------------------------"
    )
    m.fs.SOP.display()

    # print DOF before specifying
    print("DOF before specifying:", degrees_of_freedom(m.fs))

    # now specify the model
    # fix vairables
    m.fs.SOP.feed_side.properties_in[0].temperature.fix(273.15 + 25)

    # display model
    print(
        "--------------------------------- second display ---------------------------------"
    )
    m.fs.SOP.display()

    # print DOF after specifying
    print("DOF after specifying:", degrees_of_freedom(m.fs))
    # TODO

    return m


if __name__ == "__main__":
    main()
