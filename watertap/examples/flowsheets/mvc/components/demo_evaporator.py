from pyomo.environ import (ConcreteModel,
                           SolverFactory,
                           assert_optimal_termination)
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.examples.flowsheets.mvc.components.evaporator import Evaporator
import watertap.property_models.seawater_prop_pack as props_sw
import watertap.property_models.water_prop_pack as props_w

def main():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={'dynamic': False})
    m.fs.properties_feed = props_sw.SeawaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.evaporator = Evaporator(default={"property_package_feed": m.fs.properties_feed,
                                          "property_package_vapor": m.fs.properties_vapor,})

    # scaling
    m.fs.properties_feed.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
    m.fs.properties_feed.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'TDS'))
    m.fs.properties_vapor.set_default_scaling('flow_mass_phase_comp', 1, index=('Vap', 'H2O'))
    m.fs.properties_vapor.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
    iscale.set_scaling_factor(m.fs.evaporator.area, 1e-3)
    iscale.set_scaling_factor(m.fs.evaporator.U, 1e-3)
    iscale.set_scaling_factor(m.fs.evaporator.delta_temperature_in, 1e-1)
    iscale.set_scaling_factor(m.fs.evaporator.delta_temperature_out, 1e-1)
    iscale.set_scaling_factor(m.fs.evaporator.lmtd, 1e-1)
    #iscale.set_scaling_factor(m.fs.evaporator.heat_transfer, 1e-6)
    iscale.calculate_scaling_factors(m)

    #assert False

    # state variables
    # Feed inlet
    m.fs.evaporator.inlet_feed.flow_mass_phase_comp[0,'Liq','H2O'].fix(10)
    m.fs.evaporator.inlet_feed.flow_mass_phase_comp[0,'Liq','TDS'].fix(0.05)
    m.fs.evaporator.inlet_feed.temperature[0].fix(273.15+50.52) # K
    m.fs.evaporator.inlet_feed.pressure[0].fix(1e5) # Pa

    # Condenser inlet
    m.fs.evaporator.inlet_condenser.flow_mass_phase_comp[0, 'Vap', 'H2O'].fix(0.5)
    m.fs.evaporator.inlet_condenser.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(1e-8)
    m.fs.evaporator.inlet_condenser.temperature[0].fix(400)  # K
    m.fs.evaporator.inlet_condenser.pressure[0].fix(0.5e5)  # Pa

    # Evaporator/condenser specifications
    # m.fs.evaporator.outlet_brine.flow_mass_phase_comp[0,'Liq','TDS'].fix(0.1)
    #m.fs.evaporator.outlet_vapor.flow_mass_phase_comp[0,'Vap','H2O'].fix(0.5)
    m.fs.evaporator.outlet_brine.temperature[0].fix(273.15+60)
    #m.fs.evaporator.outlet_vapor.pressure[0].fix(0.3e5)
    m.fs.evaporator.U.fix(1e3) # W/K-m^2
    m.fs.evaporator.area.fix(100) # m^2
    #m.fs.evaporator.outlet_condenser.temperature[0].fix(340)
    #m.fs.evaporator.heat_transfer.fix(24e6)

    # solving
    # assert_units_consistent(m)
    # print(degrees_of_freedom(m))
    # assert(degrees_of_freedom(m) == 0)
    m.fs.evaporator.initialize(outlvl=idaeslog.INFO_HIGH)
    #assert False
    solver = get_solver()
    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)

    return m

if __name__ == "__main__":
    m = main()