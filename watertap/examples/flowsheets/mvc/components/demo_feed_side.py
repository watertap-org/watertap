from pyomo.environ import (ConcreteModel,
                           SolverFactory,
                           assert_optimal_termination)
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.examples.flowsheets.mvc.components.feed_side import Feed_side
import watertap.property_models.seawater_prop_pack as props_sw
import watertap.property_models.water_prop_pack as props_w

def main():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={'dynamic': False})
    m.fs.properties_feed = props_sw.SeawaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.evaporator = Feed_side(default={"property_package_feed": m.fs.properties_feed,
                                          "property_package_vapor": m.fs.properties_vapor,})

    # scaling
    m.fs.properties_feed.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
    m.fs.properties_feed.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'TDS'))
    m.fs.properties_vapor.set_default_scaling('flow_mass_phase_comp', 1, index=('Vap', 'H2O'))
    m.fs.properties_vapor.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
    #iscale.set_scaling_factor(m.fs.evaporator.heat_transfer, 1e-6) # found automatically now based on enthalpy flows
    iscale.calculate_scaling_factors(m)

    # state variables
    m.fs.evaporator.inlet_feed.flow_mass_phase_comp[0,'Liq','H2O'].fix(1)
    m.fs.evaporator.inlet_feed.flow_mass_phase_comp[0,'Liq','TDS'].fix(0.05)
    m.fs.evaporator.inlet_feed.temperature[0].fix(310) # K
    m.fs.evaporator.inlet_feed.pressure[0].fix(1e5) # Pa

    # specifications
    m.fs.evaporator.outlet_brine.temperature[0].fix(333)
    m.fs.evaporator.outlet_vapor.flow_mass_phase_comp[0,'Vap','H2O'].fix(0.5)
    #m.fs.evaporator.heat_transfer.fix(24e6)

    # solving
    assert_units_consistent(m)
    assert(degrees_of_freedom(m) == 0)
    # m.fs.compressor.initialize(outlvl=idaeslog.INFO_HIGH)
    m.fs.evaporator.initialize()
    solver = get_solver()
    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)
    # bad_scale_var_gen = iscale.badly_scaled_var_generator(m)
    # for (var, val) in bad_scale_var_gen:
    #     print(var.name, val)

    return m

if __name__ == "__main__":
    m = main()