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

from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock
#from ion_DSPMDE_prop_pack_inwork import DSPMDEParameterBlock # testing with the updated property package

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

ion_tran_num = {('cem','Na_+'): 1, ('cem', 'Cl_-'):0, ('aem','Na_+'): 0, ('aem', 'Cl_-'):1}
water_trans_number = {'cem': 5, 'aem': 5}
#water_permeability = 2.16e-14

# attach prop pack to flowsheet
m.fs.properties = DSPMDEParameterBlock(default=ion_dict)

# build the unit model, pass prop pack to the model
m.fs.unit = Electrodialysis0D(default = {"property_package": m.fs.properties})
m.fs.unit.diluate_channel.display()


print('----------------------------------------------')
#print('DOF before specifying:', degrees_of_freedom(m.fs))
#print('number of var beofre specifying:', number_variables(m.fs))
#print('number of constr beofre specifying:', number_total_constraints(m.fs))
print('report model statistics before specifying',report_statistics(m.fs))

assert_units_consistent(m)


# specify the feed for each inlet stream
m.fs.unit.inlet_diluate.pressure.fix(101325)
m.fs.unit.inlet_diluate.temperature.fix(298.15)
#m.fs.unit.outlet_diluate.temperature.fix(298.15)
m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, 'Liq', 'H2O'].fix(0.013)
m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, 'Liq', 'Na_+'].fix(2.46e-5)
m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, 'Liq', 'Cl_-'].fix(2.46e-5)

#m.fs.unit.inlet_diluate.conc_mol_phase_comp[0, 'Liq', 'Na_+'].fix(102.67)
#m.fs.unit.inlet_diluate.conc_mol_phase_comp[0, 'Liq', 'Cl_-'].fix(102.67)


m.fs.unit.inlet_concentrate.pressure.fix(101325)
m.fs.unit.inlet_concentrate.temperature.fix(298.15)
#m.fs.unit.outlet_concentrate.temperature.fix(298.15)
m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'H2O'].fix(0.013)
m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'Na_+'].fix(2.46e-5)
m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'Cl_-'].fix(2.46e-5)

#m.fs.unit.inlet_concentrate.conc_mol_phase_comp[0, 'Liq', 'Na_+'].fix(102.67)
#m.fs.unit.inlet_concentrate.conc_mol_phase_comp[0, 'Liq', 'Cl_-'].fix(102.67)



m.fs.unit.water_trans_number_membrane.fix(1.9)
m.fs.unit.water_permeability_membrane['cem'].fix(2.16e-14)
m.fs.unit.water_permeability_membrane['aem'].fix(1.75e-14)
m.fs.unit.voltage.fix(0.4)
#m.fs.unit.current.fix(30)
m.fs.unit.current_utilization.fix(1)
m.fs.unit.spacer_thickness.fix(1.5e-4)
m.fs.unit.membrane_surface_resistence['cem'].fix(1.89e-4)
m.fs.unit.membrane_surface_resistence['aem'].fix(1.77e-4)
#m.fs.unit.slt_eq_conductivity.fix(0.001)
m.fs.unit.cell_width.fix(0.1)
m.fs.unit.cell_length.fix(0.43)
#m.fs.unit.T = 298.15
m.fs.unit.membrane_thickness['aem'].fix(1.3e-4)
m.fs.unit.membrane_thickness['cem'].fix(1.3e-4)
m.fs.unit.ion_diffusivity_membrane.fix(7e-9)
m.fs.unit.ion_trans_number_membrane['cem','Na_+'].fix(1)
m.fs.unit.ion_trans_number_membrane['aem','Na_+'].fix(0)
m.fs.unit.ion_trans_number_membrane['cem','Cl_-'].fix(0)
m.fs.unit.ion_trans_number_membrane['aem','Cl_-'].fix(1)


#m.fs.unit.R.fix()
#m.fs.unit.vHcoef.fix()
#m.fs.unit.osmotic_coef.fix()
#m.fs.unit.faraday_const.fix()
#m.fs.unit.water_density.fix()
#m.fs.unit.water_MW.fix()


print('----------------------------------------------')
#print('DOF after specifying:', degrees_of_freedom(m.fs)) #should be zero after specifying everything 
#print('number of var after specifyig', number_variables(m.fs))
#print('number of unfixed var after specifying', number_unfixed_variables(m.fs))
#print('number of constr after specifying', number_total_constraints(m.fs))
print('report model statistics after specifying',report_statistics(m.fs))

#m.fs.unit.pprint()
if degrees_of_freedom(m.fs) != 0:
    print("error")
    exit()


# set scaling factors for state vars and call the 'calculate_scaling_factors' function
m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e1, index=('Liq', 'H2O'))
m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e4, index=('Liq', 'Na_+'))
m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e4, index=('Liq', 'Cl_-'))


# NOTE: We have to skip this step for now due to an error in Adams' Prop Pack
iscale.calculate_scaling_factors(m.fs)

print("REPORT UNSCALED VARS & CONSTRAINS")
for i in iscale.unscaled_constraints_generator(m):
    print('unscaled constaints:',  i.name)

print("REPORT BADLY SCALED VARS & CONSTRAINS")
badly_scaled_var_values = {
        var.name: val
        for (var, val) in iscale.badly_scaled_var_generator(m, large=100, small=0.01, zero = 1e-10)
        }

for j, k in badly_scaled_var_values.items():
    print(j, ':', k)
# Intialize the model
m.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)

# Solve the model
results = solver.solve(m, tee=True)
assert_optimal_termination(results)
# Display results 'cleanly'
# NOTE: This doesn't work because the 'report' function is
#       expecting our ports to be named 'inlet' and 'outlet',
#       respectively
#m.fs.unit.report()

# Display full set of model info on the unit
#m.fs.unit.pprint()
#print('Scaling Inspect')
#iscale.badly_scaled_var_generator(m.fs,large=10000.0, small=0.001, zero=1e-10, descend_into=True, include_fixed=False)




# Display the material balance constraints
#m.fs.unit.diluate_channel.material_balances.pprint()
#m.fs.unit.concentrate_channel.material_balances.pprint()
#print(iscale.extreme_jacobian_entries(m))
# Display the mass transfer terms 
print('===check property information===')
m.fs.unit.diluate_channel.properties_in[0].flow_mol_phase_comp['Liq','H2O'].pprint()
m.fs.unit.diluate_channel.properties_in[0].flow_mol_phase_comp['Liq','Na_+'].pprint()
m.fs.unit.diluate_channel.properties_in[0].flow_mol_phase_comp['Liq','Cl_-'].pprint()
m.fs.unit.diluate_channel.properties_in[0].conc_mol_phase_comp['Liq','H2O'].pprint()
m.fs.unit.diluate_channel.properties_in[0].conc_mol_phase_comp['Liq','Na_+'].pprint()
m.fs.unit.diluate_channel.properties_in[0].conc_mol_phase_comp['Liq','Cl_-'].pprint()

m.fs.unit.diluate_channel.properties_out[0].conc_mol_phase_comp['Liq','H2O'].pprint()
m.fs.unit.diluate_channel.properties_out[0].conc_mol_phase_comp['Liq','Na_+'].pprint()
m.fs.unit.diluate_channel.properties_out[0].conc_mol_phase_comp['Liq','Cl_-'].pprint()





m.fs.unit.diluate_channel.mass_transfer_term.pprint()
m.fs.unit.concentrate_channel.mass_transfer_term.pprint()

print('----------------------------------------------')
m.fs.unit.elec_migration_flux_in.pprint()
#m.fs.unit.elec_migration_flux_in.pprint()
m.fs.unit.elec_migration_flux_out.pprint()
#m.fs.unit.concentrate_channel.elec_migration_flux_out.pprint()
print('----------------------------------------------')
m.fs.unit.nonelec_flux_in.pprint()
#m.fs.unit.concentrate_channel.nonelec_flux_in.pprint()
m.fs.unit.nonelec_flux_out.pprint()
#m.fs.unit.concentrate_channel.nonelec_flux_out.pprint()
print('----------------------------------------------')
m.fs.unit.report()

print(Constants.gas_constant)
