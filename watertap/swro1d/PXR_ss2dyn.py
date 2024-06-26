from pyomo.environ import (
    ConcreteModel,
    Var,
    Expression,
    Constraint,
    TerminationCondition,
    SolverFactory,
    SolverStatus,
    value,
    units as pyunits,
)
from pyomo.network import Port
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from watertap.unit_models.pressure_exchanger import (
    PressureExchanger,
    PressureExchangeType,
)
from watertap.core.solvers import get_solver
import watertap.property_models.seawater_prop_pack as props


m = ConcreteModel()
m.fs = FlowsheetBlock(
    dynamic=True,
    time_set=[0, 3],
    time_units=pyunits.s
)
# m.fs = FlowsheetBlock(dynamic=False)
m.fs.properties = props.SeawaterParameterBlock()
m.fs.unit = PressureExchanger(dynamic=True, has_holdup=True, property_package=m.fs.properties)

# Specify inlet conditions
temperature = 25 + 273.15
flow_vol = 1e-3
lowP_mass_frac_TDS = 0.035
lowP_pressure = 101325
highP_mass_frac_TDS = 0.07
highP_pressure = 50e5

m.fs.unit.feed_side.properties_in[0].flow_vol_phase["Liq"].fix(flow_vol)
m.fs.unit.feed_side.properties_in[0].mass_frac_phase_comp["Liq", "TDS"].fix(
    lowP_mass_frac_TDS
)

m.fs.unit.feed_side.properties_in[0].pressure.fix(lowP_pressure)
m.fs.unit.feed_side.properties_in[0].temperature.fix(temperature)

m.fs.unit.brine_side.properties_in[0].flow_vol_phase["Liq"].fix(flow_vol)
m.fs.unit.brine_side.properties_in[0].mass_frac_phase_comp["Liq", "TDS"].fix(
    highP_mass_frac_TDS
)
m.fs.unit.brine_side.properties_in[0].pressure.fix(highP_pressure)
m.fs.unit.brine_side.properties_in[0].temperature.fix(temperature)

# solve inlet conditions and only fix state variables (i.e. unfix flow_vol and mass_frac_phase)

solver = get_solver()

results = solver.solve(m.fs.unit.brine_side.properties_in[0])
assert results.solver.termination_condition == TerminationCondition.optimal
print(value(m.fs.unit.brine_side.properties_in[0].flow_mass_phase_comp["Liq", "TDS"]))
m.fs.unit.brine_side.properties_in[0].flow_mass_phase_comp["Liq", "TDS"].fix()
m.fs.unit.brine_side.properties_in[0].flow_vol_phase["Liq"].unfix()
m.fs.unit.brine_side.properties_in[0].mass_frac_phase_comp["Liq", "TDS"].unfix()

results = solver.solve(m.fs.unit.feed_side.properties_in[0])
assert results.solver.termination_condition == TerminationCondition.optimal
print(value(m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "TDS"]))
m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix()
m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "TDS"].fix()
m.fs.unit.feed_side.properties_in[0].flow_vol_phase["Liq"].unfix()
m.fs.unit.feed_side.properties_in[0].mass_frac_phase_comp["Liq", "TDS"].unfix()

opt = SolverFactory("ipopt")
solve_status = opt.solve(m, tee=True)

m.fs.unit.report()

# Specify unit
# efficiency = 0.95
# m.fs.unit.efficiency_pressure_exchanger.fix(efficiency)