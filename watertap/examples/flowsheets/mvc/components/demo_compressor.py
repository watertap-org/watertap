from pyomo.environ import (ConcreteModel,
                           SolverFactory,
                           assert_optimal_termination)
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.examples.flowsheets.mvc.components.compressor import Compressor
import watertap.property_models.water_prop_pack as props

def main():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={'dynamic': False})
    m.fs.properties = props.WaterParameterBlock()
    m.fs.compressor = Compressor(default={"property_package": m.fs.properties})
    #m.fs.compressor.control_volume.display()

    # scaling
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Vap', 'H2O'))
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
    iscale.set_scaling_factor(m.fs.compressor.control_volume.work, 1e-6)
    iscale.calculate_scaling_factors(m)

    # state variables
    m.fs.compressor.inlet.flow_mass_phase_comp[0,'Vap','H2O'].fix(1)
    m.fs.compressor.inlet.flow_mass_phase_comp[0,'Liq','H2O'].fix(1e-8)
    m.fs.compressor.inlet.temperature[0].fix(350) # K
    m.fs.compressor.inlet.pressure[0].fix(0.5e5) # Pa

    # specifications
    m.fs.compressor.pressure_ratio.fix(2)
    m.fs.compressor.efficiency.fix(0.8)


    # m.fs.compressor.eq_compressor_work.pprint()
    # m.fs.compressor.control_volume.enthalpy_balances.pprint()
    # assert False

    # solving
    assert_units_consistent(m)
    assert(degrees_of_freedom(m) == 0)

    m.fs.compressor.initialize(outlvl=idaeslog.INFO_HIGH)
    solver = get_solver()
    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)

    m.fs.compressor.report()

    return m

if __name__ == "__main__":
    m = main()