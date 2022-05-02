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
from electrodialysis_1D import Electrodialysis1D

from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock,
)

from idaes.core import AqueousPhase
from idaes.core.components import Solvent, Solute, Cation, Anion
from idaes.core.phases import PhaseType as PT

from idaes.generic_models.properties.core.pure.ConstantProperties import Constant
from idaes.generic_models.properties.core.state_definitions import FpcTP
from idaes.generic_models.properties.core.eos.ideal import Ideal

from idaes.core.util import get_solver

solver = get_solver()

def build_model():
    # create model, flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    # create dict to define ions (the prop pack of Adam requires this)
    ion_dict = {
        "solute_list": ["Na_+", "Cl_-", "NaCl"],
        "diffusivity_data": {("Liq", "Na_+"): 1.33e-9,
                             ("Liq", "Cl_-"): 2.03e-9,
                             ("Liq", "NaCl"): 1.70e-9},
        "mw_data": {"H2O": 18e-3,
                    "Na_+": 23e-3,
                    "Cl_-": 35e-3,
                    "NaCl": 58e-3},
        "stokes_radius_data": {"Na_+": 0.184e-9,
                               "Cl_-": 0.121e-9,
                               "NaCl": 0.305e-9},
        "charge": {"Na_+": 1,
                   "Cl_-": -1,
                   "NaCl": 0},
    }
    # attach prop pack to flowsheet
    m.fs.properties = DSPMDEParameterBlock(default=ion_dict)

    # build the unit model, pass prop pack to the model
    m.fs.unit = Electrodialysis1D(default={
        "property_package": m.fs.properties })


    print('----------------------------------------------')
    print('DOF before specifying:', degrees_of_freedom(m.fs))

    assert_units_consistent(m)

    return m

def build_model_generic():
    # create model, flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    # Configuration dictionary for generic
    thermo_config = {
        "components": {
            "H2O": {
                "type": Solvent,
                # Define the methods used to calculate the following properties
                "dens_mol_liq_comp": Constant,
                "enth_mol_liq_comp": Constant,
                "cp_mol_liq_comp": Constant,
                "entr_mol_liq_comp": Constant,
                # Parameter data is always associated with the methods defined above
                "parameter_data": {
                    "mw": (18.0153, pyunits.g / pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.K / pyunits.mol),
                },
                # End parameter_data
            },
            "Na_+": {
                "type": Cation,
                "charge": 1,
                # Define the methods used to calculate the following properties
                "dens_mol_liq_comp": Constant,
                "enth_mol_liq_comp": Constant,
                "cp_mol_liq_comp": Constant,
                "entr_mol_liq_comp": Constant,
                # Parameter data is always associated with the methods defined above
                "parameter_data": {
                    "mw": (22.989769, pyunits.g / pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.K / pyunits.mol),
                },
                # End parameter_data
            },
            "Cl_-": {
                "type": Anion,
                "charge": -1,
                # Define the methods used to calculate the following properties
                "dens_mol_liq_comp": Constant,
                "enth_mol_liq_comp": Constant,
                "cp_mol_liq_comp": Constant,
                "entr_mol_liq_comp": Constant,
                # Parameter data is always associated with the methods defined above
                "parameter_data": {
                    "mw": (35.453, pyunits.g / pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.K / pyunits.mol),
                },
                # End parameter_data
            },
            "NaCl": {
                "type": Solute,
                "valid_phase_types": PT.aqueousPhase,
                # Define the methods used to calculate the following properties
                "dens_mol_liq_comp": Constant,
                "enth_mol_liq_comp": Constant,
                "cp_mol_liq_comp": Constant,
                "entr_mol_liq_comp": Constant,
                # Parameter data is always associated with the methods defined above
                "parameter_data": {
                    "mw": (58.442, pyunits.g / pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (75.312, pyunits.J / pyunits.mol / pyunits.K),
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.K / pyunits.mol),
                },
                # End parameter_data
            },
        },
        # End Component list
        "phases": {
            "Liq": {"type": AqueousPhase, "equation_of_state": Ideal},
        },
        "state_definition": FpcTP,
        "state_bounds": {
            "temperature": (273.15, 300, 650),
            "pressure": (5e4, 1e5, 1e6),
        },
        "pressure_ref": 1e5,
        "temperature_ref": 300,
        "base_units": {
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        }
    }
    # End thermo_config definition

    # attach prop pack to flowsheet
    m.fs.properties = GenericParameterBlock(default=thermo_config)

    # build the unit model, pass prop pack to the model
    m.fs.unit = Electrodialysis1D(default={
        "property_package": m.fs.properties })


    print('----------------------------------------------')
    print('DOF before specifying:', degrees_of_freedom(m.fs))

    assert_units_consistent(m)

    return m

def fix_inlets_and_vars(m):
    # specify the feed for each inlet stream
    m.fs.unit.inlet_dilute.pressure.fix(101325)
    m.fs.unit.inlet_dilute.temperature.fix(298.15)
    m.fs.unit.inlet_dilute.flow_mol_phase_comp[0, 'Liq', 'H2O'].fix(2)
    m.fs.unit.inlet_dilute.flow_mol_phase_comp[0, 'Liq', 'Na_+'].fix(0.2)
    m.fs.unit.inlet_dilute.flow_mol_phase_comp[0, 'Liq', 'Cl_-'].fix(0.2)
    m.fs.unit.inlet_dilute.flow_mol_phase_comp[0, 'Liq', 'NaCl'].fix(0.002)

    m.fs.unit.inlet_concentrate.pressure.fix(101325)
    m.fs.unit.inlet_concentrate.temperature.fix(298.15)
    m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'H2O'].fix(1)
    m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'Na_+'].fix(0.15)
    m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'Cl_-'].fix(0.15)
    m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'NaCl'].fix(0.0015)

    m.fs.unit.cell_length.fix(0.5)
    m.fs.unit.cell_width.fix(0.1)
    m.fs.unit.membrane_thickness['cem'].fix()
    m.fs.unit.membrane_thickness['aem'].fix()

    m.fs.unit.ion_diffusivity_membrane.fix()
    m.fs.unit.water_permeability_membrane.fix()

    print('----------------------------------------------')
    print('DOF after specifying:', degrees_of_freedom(m.fs))

    if degrees_of_freedom(m.fs) != 0:
        print("ERROR! Problem not square")
        exit(-1)

def scale_model(m):
    # set scaling factors for state vars and call the 'calculate_scaling_factors' function

    # Special function to attempt to fix scaling issues from prop pack
    #   by manually applying the inlet condition as the initial guess for the states.
    #m.fs.unit.propogate_initial_state()

    # # TODO: Figure out if this is NOT the proper way to  set scaling factors in 1D unit model
    m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1, index=('Liq', 'H2O'))
    m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e1, index=('Liq', 'Na_+'))
    m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e1, index=('Liq', 'Cl_-'))
    m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e3, index=('Liq', 'NaCl'))

    # NOTE: We have to skip this step for now due to an error in Adams' Prop Pack

    # # TODO: Adam's scaling may need to be revisted

    # # TODO: Figure out why 'flow_mass_phase_comp' is being constructed when it is
    #           not being called for... this is because the way in which the prop pack
    #           is written is that the 'conc_mol_phase_comp' is a function of the
    #           'flow_mass_phase_comp' (even though it doesn't need to be).
    iscale.calculate_scaling_factors(m)

    unscaled_constraint_list = list(iscale.unscaled_constraints_generator(m))
    print("List of unscaled constraints")
    print("----------------------------")
    for j in unscaled_constraint_list:
        print(j)
    print()

    # check that all variables have scaling factors
    unscaled_var_list = list(iscale.unscaled_variables_generator(m))
    print("List of unscaled variables")
    print("--------------------------")
    for j in unscaled_var_list:
        print(j)
    print()

    # check if any variables are badly scaled
    badly_scaled_var_values = {
        var.name: val
        for (var, val) in iscale.badly_scaled_var_generator(
            m, large=1e3, small=1e-3
        )
    }
    print("List of poorly scaled variables")
    print("-------------------------------")
    for j in badly_scaled_var_values:
        print(str(j) + "\t" + str(badly_scaled_var_values[j]))
    print()

def initialize_model(m):
    # Intialize the model
    m.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)

def run_model(m):
    # Check scaling one last time
    unscaled_constraint_list = list(iscale.unscaled_constraints_generator(m))
    print("List of unscaled constraints")
    print("----------------------------")
    for j in unscaled_constraint_list:
        print(j)
    print()

    # check that all variables have scaling factors
    unscaled_var_list = list(iscale.unscaled_variables_generator(m))
    print("List of unscaled variables")
    print("--------------------------")
    for j in unscaled_var_list:
        print(j)
    print()

    # check if any variables are badly scaled
    badly_scaled_var_values = {
        var.name: val
        for (var, val) in iscale.badly_scaled_var_generator(
            m, large=1e3, small=1e-3
        )
    }
    print("List of poorly scaled variables")
    print("-------------------------------")
    for j in badly_scaled_var_values:
        print(str(j) + "\t" + str(badly_scaled_var_values[j]))
    print()

    # Solve the model
    results = solver.solve(m, tee=True)

def display_results(m):
    # Display results 'cleanly'
    m.fs.unit.report()

    # Display full set of model info on the unit
    #m.fs.unit.pprint()

def view_model_constraints(m):
    # Display the material balance constraints
    m.fs.unit.dilute_side.material_balances.pprint()
    m.fs.unit.dilute_side.material_flow_dx_disc_eq.pprint()
    m.fs.unit.dilute_side.mass_transfer_term.pprint()
    m.fs.unit.eq_mass_transfer_term_dilute.pprint()

    print()

    # Display the mass transfer terms
    m.fs.unit.concentrate_side.material_balances.pprint()
    m.fs.unit.concentrate_side.material_flow_dx_disc_eq.pprint()
    m.fs.unit.concentrate_side.mass_transfer_term.pprint()
    m.fs.unit.eq_mass_transfer_term_concentrate.pprint()

    m.fs.unit.cell_length.pprint()

    print()

    m.fs.unit.nonelec_flux.pprint()
    m.fs.unit.eq_nonelec_flux.pprint()

def view_model_control_volumes(m):
    # Display the full control volume equation set
    m.fs.unit.dilute_side.pprint()
    print()
    m.fs.unit.concentrate_side.pprint()

def view_model_properties(m):
    # Display the full control volume equation set
    m.fs.unit.dilute_side.properties.pprint()
    print()
    m.fs.unit.concentrate_side.properties.pprint()

## Run for testing purposes ##
if __name__ == "__main__":

   #m = build_model()
   m = build_model_generic()

   fix_inlets_and_vars(m)
   scale_model(m)

   # ONLY use these line to solve the model (requires proper solvers)
   initialize_model(m)
   run_model(m)
   display_results(m)

   m.fs.unit.dilute_side.scaling_factor.pprint()
   m.fs.unit.dilute_side.constraint_transformed_scaling_factor.pprint()

   print(iscale.get_scaling_factor(m.fs.unit.dilute_side._flow_terms[0.0,0.0,"Liq","NaCl"]))

   print(value(m.fs.unit.eq_nonelec_flux[0,0.1,"Liq","H2O"]))

   #m.fs.unit.dilute_side.properties[0,1].conc_mol_phase_comp.pprint()
   #m.fs.unit.concentrate_side.properties[0,1].conc_mol_phase_comp.pprint()

   #view_model_constraints(m)
   #view_model_control_volumes(m)
   #view_model_properties(m)

   # Evaluate the osm_pressure from generic prop pack
   #print(value(m.fs.unit.concentrate_side.properties[0,1].pressure_osm_phase['Liq']))
