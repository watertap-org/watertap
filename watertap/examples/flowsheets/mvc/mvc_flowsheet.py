from pyomo.environ import (ConcreteModel, value, Param, Var, Constraint, Expression, Objective, TransformationFactory,
                           Block, NonNegativeReals, RangeSet, check_optimal_termination,
                           units as pyunits)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
# Import components
from watertap.examples.flowsheets.mvc.components.evaporator import Evaporator
from watertap.examples.flowsheets.mvc.components.compressor import Compressor
# Import property packages
import watertap.property_models.seawater_prop_pack as props_sw
import watertap.property_models.water_prop_pack as props_w

def build(m):
    # Properties
    m.fs.properties_feed = props_sw.SeawaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()

    # Evaporator
    m.fs.evaporator = Evaporator(default={"property_package_feed": m.fs.properties_feed,
                                          "property_package_vapor": m.fs.properties_vapor,})
    # Compressor
    m.fs.compressor = Compressor(default={"property_package": m.fs.properties_vapor})

    return m

def solve_flowsheet():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    build(m)
    TransformationFactory("network.expand_arcs").apply_to(m)


if __name__ == "__main__":
    m = solve_flowsheet()
