from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale

from idaes.core.util.scaling import badly_scaled_var_generator

import coagulation_prop_pack as props

# create model, flowsheet
m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
# attach property package
m.fs.properties = props.CoagulationParameterBlock()
# build a state block, must specify a time which by convention for steady state models is just 0
m.fs.stream = m.fs.properties.build_state_block([0], default={})

# display the state block, it only has the state variables and they are all unfixed
print('\n---first display---')
m.fs.stream[0].display()

# attempt to access properties so that they are built
m.fs.stream[0].mass_frac_phase_comp
# after touching the property, the state block automatically builds it,
# note the mass_frac_phase_comp variable and the constraint to calculate it
print('\n---second display---')
m.fs.stream[0].display()

# touch another property
m.fs.stream[0].conc_mass_phase_comp
# after touching this property, the state block automatically builds it AND any other properties that are necessary,
# note that now there is the conc_mass_phase_comp and dens_mass_phase variable and associated constraints
print('\n---third display---')
m.fs.stream[0].display()

# touch a property to construct
m.fs.stream[0].flow_vol_phase
m.fs.stream[0].flow_mol_phase_comp
m.fs.stream[0].mol_frac_phase_comp
m.fs.stream[0].conc_mol_phase_comp

# now that we have a state block, we can fix the state variables and solve for the properties
m.fs.stream[0].temperature.fix(273.15 + 25)
m.fs.stream[0].pressure.fix(101325)
m.fs.stream[0].flow_mass_phase_comp['Liq', 'H2O'].fix(10)
m.fs.stream[0].flow_mass_phase_comp['Liq', 'TSS'].fix(0.1)
m.fs.stream[0].flow_mass_phase_comp['Liq', 'TDS'].fix(0.1)
m.fs.stream[0].flow_mass_phase_comp['Sol', 'Sludge'].fix(0.001)

# the user should provide the scale for the flow rate, so that our tools can ensure the model is well scaled
# generally scaling factors should be such that if it is multiplied by the variable it will range between 0.01 and 100
m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e-1, index=('Liq', 'H2O'))
m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e1, index=('Liq', 'TSS'))
m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e1, index=('Liq', 'TDS'))
m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e3, index=('Sol', 'Sludge'))
iscale.calculate_scaling_factors(m.fs)  # this utility scales the model

# Looking for poor scaling
badly_scaled_var_list = list(badly_scaled_var_generator(m, large=1e2, small=1e-2))
if len(badly_scaled_var_list) != 0:
    lst = []
    for (var, val) in badly_scaled_var_list:
        lst.append((var.name, val))
        print(var.name, var.value)
    print("The following variable(s) are poorly scaled: {lst}".format(lst=lst))

# solving
assert_units_consistent(m)  # check that units are consistent
assert(degrees_of_freedom(m) == 0)  # check that the degrees of freedom are what we expect

solver = SolverFactory('ipopt')
solver.options = {'tol': 1e-8, 'nlp_scaling_method': 'user-scaling'}

results = solver.solve(m, tee=True)
assert results.solver.termination_condition == TerminationCondition.optimal

# Looking for poor scaling
badly_scaled_var_list = list(badly_scaled_var_generator(m, large=1e2, small=1e-2))
if len(badly_scaled_var_list) != 0:
    lst = []
    for (var, val) in badly_scaled_var_list:
        lst.append((var.name, val))
        print(var.name, var.value)
    print("The following variable(s) are poorly scaled: {lst}".format(lst=lst))


# display results
print('\n---fourth display---')
m.fs.stream[0].display()
# note that the properties are solved, and the body of the constraints are small (residual)

# equation oriented modeling has several advantages, one of them is that we can unfix variables and fix others
# instead of setting the mass flow rates, we can set the volumetric flow rate and mass fractions
m.fs.stream[0].flow_mass_phase_comp['Liq', 'H2O'].unfix()
m.fs.stream[0].flow_mass_phase_comp['Liq', 'TSS'].unfix()
m.fs.stream[0].flow_mass_phase_comp['Liq', 'TDS'].unfix()

m.fs.stream[0].flow_vol_phase['Liq'].fix(1.5e-3)
m.fs.stream[0].mass_frac_phase_comp['Liq', 'TSS'].fix(80e-4)
m.fs.stream[0].mass_frac_phase_comp['Liq', 'TDS'].fix(80e-4)

# resolve
results = solver.solve(m, tee=True)
assert results.solver.termination_condition == TerminationCondition.optimal

print('\n---fifth display---')
m.fs.stream[0].display()

# Display contraint with scaling factors
#m.fs.stream[0].eq_flow_mol_phase_comp.pprint()

# Looking for poor scaling
badly_scaled_var_list = list(badly_scaled_var_generator(m, large=1e2, small=1e-2))
if len(badly_scaled_var_list) != 0:
    lst = []
    for (var, val) in badly_scaled_var_list:
        lst.append((var.name, val))
        print(var.name, var.value)
    print("The following variable(s) are poorly scaled: {lst}".format(lst=lst))
