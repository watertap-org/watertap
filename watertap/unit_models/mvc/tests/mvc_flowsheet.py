from pyomo.environ import (ConcreteModel, TransformationFactory,
                           assert_optimal_termination)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.util.initialization import propagate_state

# Import components
from watertap.unit_models.mvc.components.evaporator import Evaporator
from watertap.unit_models.mvc.components.compressor import Compressor
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

    # Connections
    m.fs.s01 = Arc(source=m.fs.evaporator.outlet_vapor, destination=m.fs.compressor.inlet)
    m.fs.s02 = Arc(source=m.fs.compressor.outlet, destination=m.fs.evaporator.inlet_condenser)

    TransformationFactory("network.expand_arcs").apply_to(m)


def scale(m):
    m.fs.properties_feed.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
    m.fs.properties_feed.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'TDS'))
    m.fs.properties_vapor.set_default_scaling('flow_mass_phase_comp', 1, index=('Vap', 'H2O'))
    m.fs.properties_vapor.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
    # Evaporator
    iscale.set_scaling_factor(m.fs.evaporator.area, 1e-3)
    iscale.set_scaling_factor(m.fs.evaporator.U, 1e-3)
    iscale.set_scaling_factor(m.fs.evaporator.delta_temperature_in, 1e-1)
    iscale.set_scaling_factor(m.fs.evaporator.delta_temperature_out, 1e-1)
    iscale.set_scaling_factor(m.fs.evaporator.lmtd, 1e-1)
    # iscale.set_scaling_factor(m.fs.evaporator.heat_transfer, 1e-6)
    # Compressor
    iscale.set_scaling_factor(m.fs.compressor.control_volume.work, 1e-6)

    iscale.calculate_scaling_factors(m)


def specify(m):
    # state variables
    # Feed inlet
    m.fs.evaporator.inlet_feed.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(1)
    m.fs.evaporator.inlet_feed.flow_mass_phase_comp[0, 'Liq', 'TDS'].fix(0.05)
    m.fs.evaporator.inlet_feed.temperature[0].fix(273.15 + 50.52)  # K
    m.fs.evaporator.inlet_feed.pressure[0].fix(1e5)  # Pa

    m.fs.evaporator.outlet_vapor.flow_mass_phase_comp[0,'Vap','H2O'].fix(0.5)
    # m.fs.evaporator.outlet_brine.temperature[0].fix(273.15 + 60)
    # m.fs.evaporator.outlet_vapor.pressure[0].fix(0.3e5)
    m.fs.evaporator.U.fix(1e3)  # W/K-m^2
    m.fs.evaporator.area.fix(100)  # m^2

    # specifications
    m.fs.compressor.pressure_ratio.fix(2)
    m.fs.compressor.efficiency.fix(0.8)

def initialize(m):
    m.fs.evaporator.inlet_condenser.flow_mass_phase_comp[0, 'Vap', 'H2O'].fix(0.5)
    m.fs.evaporator.inlet_condenser.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(1e-8)
    m.fs.evaporator.inlet_condenser.temperature[0].fix(400)  # K
    m.fs.evaporator.inlet_condenser.pressure[0].fix(0.5e5)  # Pa

    m.fs.evaporator.initialize(outlvl=idaeslog.INFO_HIGH)

    m.fs.evaporator.inlet_condenser.flow_mass_phase_comp[0, 'Vap', 'H2O'].unfix()
    m.fs.evaporator.inlet_condenser.flow_mass_phase_comp[0, 'Liq', 'H2O'].unfix()
    m.fs.evaporator.inlet_condenser.temperature[0].unfix()  # K
    m.fs.evaporator.inlet_condenser.pressure[0].unfix()  # Pa


    propagate_state(m.fs.s01)
    m.fs.compressor.control_volume.properties_in[0].display()

    m.fs.compressor.initialize(outlvl=idaeslog.INFO_HIGH)


def solve_flowsheet():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    build(m)
    scale(m)
    specify(m)
    initialize(m)
    print(degrees_of_freedom(m))

    solver = get_solver()
    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)
    m.display()

if __name__ == "__main__":
    m = solve_flowsheet()
