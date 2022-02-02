from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale

import cryst_prop_pack as props


from io import StringIO
from pyomo.util.infeasible import (log_active_constraints, log_close_to_bounds,
                                   log_infeasible_bounds,
                                   log_infeasible_constraints)
from pyomo.common.log import LoggingIntercept
import logging


# create model, flowsheet
m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
# attach property package
m.fs.properties = props.NaClParameterBlock(default=
                                           {"heat_of_crystallization_model": props.HeatOfCrystallizationModel.constant})
# build a state block, must specify a time which by convention for steady state models is just 0
m.fs.stream = m.fs.properties.build_state_block([0], default={'defined_state': True})

# display the state block, it only has the state variables and they are all unfixed
print('\n---first display---')
m.fs.stream[0].display()

# Attempt to access properties so that they are built
m.fs.stream[0].mass_frac_phase_comp
m.fs.stream[0].solubility_mass_phase_comp
m.fs.stream[0].dens_mass_solvent
m.fs.stream[0].dens_mass_solute
m.fs.stream[0].solubility_mass_frac_phase_comp
m.fs.stream[0].dens_mass_phase
m.fs.stream[0].dh_vap_solvent
m.fs.stream[0].cp_solvent
m.fs.stream[0].cp_solute
m.fs.stream[0].cp_phase
m.fs.stream[0].flow_vol_phase
m.fs.stream[0].flow_vol
m.fs.stream[0].pressure_sat
m.fs.stream[0].temperature_sat_solvent
m.fs.stream[0].conc_mass_phase_comp
m.fs.stream[0].enth_mass_solvent
m.fs.stream[0].enth_mass_solute
m.fs.stream[0].enth_mass_phase
m.fs.stream[0].dh_crystallization
m.fs.stream[0].flow_mol_phase_comp
m.fs.stream[0].mole_frac_phase_comp

# ========
# Case 1a: 
# ========
# Fix the state variables and solve for the properties
m.fs.stream[0].temperature.fix(273.15 + 25)
m.fs.stream[0].pressure.fix(5e5)
m.fs.stream[0].flow_mass_phase_comp['Liq', 'H2O'].fix(1)
m.fs.stream[0].flow_mass_phase_comp['Liq', 'NaCl'].fix(0.27)
m.fs.stream[0].flow_mass_phase_comp['Sol', 'NaCl'].fix(0)
m.fs.stream[0].flow_mass_phase_comp['Vap', 'H2O'].fix(0)

m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))
iscale.calculate_scaling_factors(m.fs) 
assert_units_consistent(m)  
assert(degrees_of_freedom(m) == 0) 
solver = SolverFactory('ipopt')
solver.options = {'tol': 1e-8, 'nlp_scaling_method': 'user-scaling'}
results = solver.solve(m, tee=True)
assert results.solver.termination_condition == TerminationCondition.optimal
m.fs.stream[0].display()
m.fs.stream[0].flow_vol.display()
m.fs.stream[0].enth_flow.display()


# ========
# Case 1b: 
# ========
print('\n==== New solve ====\n')
# # Initialize/re-solve with other properties fixed
m.fs.stream[0].flow_mass_phase_comp['Liq', 'H2O'].unfix()
m.fs.stream[0].flow_mass_phase_comp['Liq', 'NaCl'].unfix()
m.fs.stream[0].flow_mass_phase_comp['Sol', 'NaCl'].unfix()
m.fs.stream[0].flow_mass_phase_comp['Vap', 'H2O'].unfix()
m.fs.stream[0].flow_vol_phase['Liq'].fix(2e-2)
m.fs.stream[0].flow_vol_phase['Sol'].fix(0)
m.fs.stream[0].flow_vol_phase['Vap'].fix(0)
m.fs.stream[0].mass_frac_phase_comp['Liq', 'NaCl'].fix(0.05)
assert(degrees_of_freedom(m) == 0) 
results = solver.solve(m, tee=True)
assert results.solver.termination_condition == TerminationCondition.optimal
print('\n==== New solve ====\n')
m.fs.stream[0].display()