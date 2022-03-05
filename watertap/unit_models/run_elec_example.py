from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition, \
    value, Constraint, Var, Objective, Expression
from pyomo.environ import units as pyunits
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent, check_units_equivalent
import pyomo.util.infeasible as infeas
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom, number_variables
import idaes.core.util.model_statistics as stats
from idaes.core.util.constants import Constants
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock
from electrodialysis_0D import Electrodialysis0D

from idaes.core.util import get_solver

solver = get_solver()

# create model, flowsheet
m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
# create dict to define ions (the prop pack of Adam requires this)
ion_dict = {
    "solute_list": ["Na_+", "Cl_-"],
    "diffusivity_data": {("Liq", "Na_+"): 1.33e-9,
                         ("Liq", "Cl_-"): 2.03e-9},
    "mw_data": {"H2O": 18e-3,
                "Na_+": 23e-3,
                "Cl_-": 35e-3},
    "stokes_radius_data": {"Na_+": 0.184e-9,
                           "Cl_-": 0.121e-9},
    "charge": {"Na_+": 1,
               "Cl_-": -1},
}
# attach prop pack to flowsheet
m.fs.properties = DSPMDEParameterBlock(default=ion_dict)

# build the unit model, pass prop pack to the model
m.fs.unit = Electrodialysis0D(default={
    "property_package": m.fs.properties })


print('----------------------------------------------')
print('DOF before specifying:', degrees_of_freedom(m.fs))

assert_units_consistent(m)


# specify the feed for each inlet stream
m.fs.unit.inlet_dilute.pressure.fix(101325)
m.fs.unit.inlet_dilute.temperature.fix(298.15)
m.fs.unit.inlet_dilute.flow_mol_phase_comp[0, 'Liq', 'H2O'].fix(1)
m.fs.unit.inlet_dilute.flow_mol_phase_comp[0, 'Liq', 'Na_+'].fix(0.01)
m.fs.unit.inlet_dilute.flow_mol_phase_comp[0, 'Liq', 'Cl_-'].fix(0.01)

m.fs.unit.inlet_concentrate.pressure.fix(101325)
m.fs.unit.inlet_concentrate.temperature.fix(298.15)
m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'H2O'].fix(1)
m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'Na_+'].fix(0.01)
m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'Cl_-'].fix(0.01)


print('----------------------------------------------')
print('DOF after specifying:', degrees_of_freedom(m.fs))

# set scaling factors for state vars and call the 'calculate_scaling_factors' function
m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1, index=('Liq', 'H2O'))
m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e2, index=('Liq', 'Na_+'))
m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e2, index=('Liq', 'Cl_-'))

# NOTE: We have to skip this step for now due to an error in Adams' Prop Pack
#iscale.calculate_scaling_factors(m.fs)

# Intialize the model
m.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)

# Solve the model
results = solver.solve(m, tee=True)

# Display results 'cleanly'
# NOTE: This doesn't work because the 'report' function is
#       expecting our ports to be named 'inlet' and 'outlet',
#       respectively
#m.fs.unit.report()

# Display full set of model info on the unit
#m.fs.unit.pprint()

# Display the material balance constraints
m.fs.unit.dilute_side.material_balances.pprint()
m.fs.unit.concentrate_side.material_balances.pprint()

# Display the mass transfer terms
m.fs.unit.dilute_side.mass_transfer_term.pprint()
m.fs.unit.concentrate_side.mass_transfer_term.pprint()
