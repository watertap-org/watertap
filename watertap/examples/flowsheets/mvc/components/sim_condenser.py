from pyomo.environ import (ConcreteModel,
                           SolverFactory,
                           assert_optimal_termination)
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.examples.flowsheets.mvc.components.complete_condenser import Condenser
import watertap.property_models.water_prop_pack as props

def main():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={'dynamic': False})
    m.fs.properties = props.WaterParameterBlock()
    m.fs.unit = Condenser(default={"property_package": m.fs.properties})

    # scaling
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Vap', 'H2O'))
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
    iscale.calculate_scaling_factors(m)

    # state variables
    m.fs.unit.inlet.flow_mass_phase_comp[0, 'Vap', 'H2O'].fix(1)
    m.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(1e-8)
    m.fs.unit.inlet.temperature[0].fix(350)  # K
    m.fs.unit.inlet.pressure[0].fix(0.5e5)  # Pa

    # solving
    assert_units_consistent(m)
    m.fs.unit.display()
    print(degrees_of_freedom(m))
    assert(degrees_of_freedom(m) == 0)

    m.fs.unit.initialize()

    assert False
    solver = get_solver()
    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)

    m.fs.compressor.report()

if __name__ == "__main__":
    main()