from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition, \
    value, Constraint, Var, Objective, Expression, assert_optimal_termination
from pyomo.environ import units as pyunits
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent, check_units_equivalent
import pyomo.util.infeasible as infeas
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom, number_variables,number_unfixed_variables, number_total_constraints,report_statistics
import idaes.core.util.model_statistics as stats
from idaes.core.util.constants import Constants
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

#from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock
from ion_DSPMDE_prop_pack_inwork import DSPMDEParameterBlock # testing with the updated property package

from electrodialysis_0d import Electrodialysis0D

from idaes.core.util import get_solver
import pandas as pd
solver = get_solver()

# create model, flowsheet
m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
# create dict to define ions (the prop pack of Adam requires this)
ion_dict = {
    "solute_list": ["Na_+", "Cl_-"],
    "diffusivity_data": {("Liq", "Na_+"): 1.33e-9,
                         ("Liq", "Cl_-"): 2.03e-9
                         },
    "mw_data": {"H2O": 18e-3,
                "Na_+": 23e-3,
                "Cl_-": 35.5e-3},
    "electrical_mobility_data": {"Na_+":5.19e-8, "Cl_-":7.92e-8},
    "stokes_radius_data": {"Na_+": 0.184e-9,
                           "Cl_-": 0.121e-9},
    "charge": {"Na_+": 1,
               "Cl_-": -1
               },
}
# attach prop pack to flowsheet
m.fs.properties = DSPMDEParameterBlock(default=ion_dict)
m.fs.stream = m.fs.properties.build_state_block([0], default={})
#print('\n---Before Specifying---')
#m.fs.stream[0].display()

m.fs.stream[0].temperature.fix(298.15)
m.fs.stream[0].pressure.fix(101325)
m.fs.stream[0].flow_mol_phase_comp['Liq', 'H2O'].fix(0.013)
m.fs.stream[0].flow_mol_phase_comp['Liq', 'Na_+'].fix(2.46e-5)
m.fs.stream[0].flow_mol_phase_comp['Liq', 'Cl_-'].fix(2.46e-5)
m.fs.stream[0].conc_mol_phase_comp
m.fs.stream[0].molality_comp
m.fs.stream.initialize()
print('\n---Before Solving---')
m.fs.stream[0].display()

m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e1, index=('Liq', 'H2O'))
m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e4, index=('Liq', 'Na_+'))
m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e4, index=('Liq', 'Cl_-'))

print('\n report model statistics before solving')
#print(report_statistics(m.fs))
assert_units_consistent(m)  # check that units are consistent
assert(degrees_of_freedom(m) == 0)  # check that the degrees of freedom are what we expect

solver = SolverFactory('ipopt')
solver.options = {'tol': 1e-8, 'nlp_scaling_method': 'user-scaling'}
iscale.calculate_scaling_factors(m.fs)
print("REPORT BADLY SCALED VARS & CONSTRAINS")
badly_scaled_var_values = {
        var.name: val
        for (var, val) in iscale.badly_scaled_var_generator(m, large=100, small=0.01, zero = 1e-10)
        }
for j, k in badly_scaled_var_values.items():
    print(j, ':', k)
results = solver.solve(m, tee=False)
assert results.solver.termination_condition == TerminationCondition.optimal

#print('\n report model statistics after solving')
#print(report_statistics(m.fs))

# display results
print('\n---After Solving---')
m.fs.stream[0].display()
# note that the properties are solved, and the body of the constraints are small (residual)

