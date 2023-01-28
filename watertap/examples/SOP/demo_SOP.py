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
    m.fs.SOP = SelectiveOilPermeation(property_package=m.fs.properties)

    # display model
    print("first display")
    m.fs.SOP.display()

    # now specify the model
    print("DOF before specifying:", degrees_of_freedom(m.fs))

    # TODO

    return m


if __name__ == "__main__":
    main()
