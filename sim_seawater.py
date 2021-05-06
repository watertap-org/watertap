from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition, value, Constraint
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.environ import units as pyunits
from idaes.core.util.constants import Constants
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent, check_units_equivalent
import idaes.core.util.scaling as iscale

import proteuslib.property_models.seawater_prop_pack as props

print('--------------------------------------------------------------------------------------------------------------!')
m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.properties = props.SeawaterParameterBlock()
m.fs.stream = m.fs.properties.build_state_block([0], default={})

# specify conditions
mass_flow = 1e11
x_TDS = 0.035

m.fs.stream[0].flow_mass_phase_comp['Liq', 'TDS'].fix(x_TDS * mass_flow)
m.fs.stream[0].flow_mass_phase_comp['Liq', 'H2O'].fix((1 - x_TDS) * mass_flow)
m.fs.stream[0].temperature.fix(273.15 + 25)
m.fs.stream[0].pressure.fix(101325)

m.fs.stream[0].mass_frac_phase_comp
m.fs.stream[0].dens_mass_phase
m.fs.stream[0].dens_mass_solvent
m.fs.stream[0].flow_vol_phase     #### Problem: when mass flow too high, solver fails when activating this touch variable--> works after setting scaling factors for mass flow as reciprocal of fixed vals
m.fs.stream[0].conc_mass_phase_comp
m.fs.stream[0].flow_mol_phase_comp #### Problem: when mass flow too high, solver fails when activating this touch variable
m.fs.stream[0].mole_frac_phase_comp #### Problem: when mass flow too high, solver fails when activating this touch variable
m.fs.stream[0].molality_comp
m.fs.stream[0].visc_d_phase
m.fs.stream[0].osm_coeff
m.fs.stream[0].pressure_osm
m.fs.stream[0].enth_mass_phase
m.fs.stream[0].enth_flow
m.fs.stream[0].pressure_sat
m.fs.stream[0].cp_phase
m.fs.stream[0].therm_cond_phase
m.fs.stream[0].dh_vap

# scaling
m.fs.properties.set_default_scaling('flow_mass_phase_comp', (1/((1-x_TDS) * mass_flow)), index=('Liq', 'H2O'))
m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1/((x_TDS) * mass_flow), index=('Liq', 'TDS'))
iscale.calculate_scaling_factors(m.fs)

m.fs.stream[0].scaling_factor.display()
bad_var_list = iscale.badly_scaled_var_generator(m)
type(bad_var_list)
print("\n\n\nBad variables List----------------------------------------------------")
[[print(v), print(k)] for v, k in bad_var_list]
# [print(v) for v in bad_var_list]

print("Bad variables List End----------------------------------------------------\n\n\n")

# solving
assert (degrees_of_freedom(m) == 0)

solver = SolverFactory('ipopt')
solver.options = {'nlp_scaling_method': 'user-scaling'}
results = solver.solve(m, tee=True)
assert results.solver.termination_condition == TerminationCondition.optimal

# display results
m.fs.stream[0].display()
