###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################

""" Simple Softening process with Lime"""

# Importing the object for units from pyomo
from pyomo.environ import units as pyunits

# Imports from idaes core
from idaes.core import AqueousPhase
from idaes.core.components import Solvent, Solute
from idaes.core.phases import PhaseType as PT

# Imports from idaes generic models
import idaes.generic_models.properties.core.pure.Perrys as Perrys
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ideal import Ideal

# Importing the enum for concentration unit basis used in the 'get_concentration_term' function
from idaes.generic_models.properties.core.generic.generic_reaction import ConcentrationForm

# Import the object/function for heat of reaction
from idaes.generic_models.properties.core.reactions.dh_rxn import constant_dh_rxn

# Import safe log power law equation
from idaes.generic_models.properties.core.reactions.equilibrium_forms import log_power_law_equil

# Import specific pyomo objects
from pyomo.environ import (ConcreteModel,
                           Constraint,
                           SolverStatus,
                           TerminationCondition,
                           value,
                           Suffix)

from idaes.core.util import scaling as iscale

# Import pyomo methods to check the system units
from pyomo.util.check_units import assert_units_consistent


from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof
from idaes.core.util import get_solver

# Import the idaes objects for Generic Properties and Reactions
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)
from idaes.generic_models.properties.core.generic.generic_reaction import (
        GenericReactionParameterBlock)

# Import the idaes object for the StoichiometricReactor unit model
from idaes.generic_models.unit_models.stoichiometric_reactor import \
    StoichiometricReactor

# Import the core idaes objects for Flowsheets and types of balances
from idaes.core import FlowsheetBlock

# Import log10 function from pyomo
from pyomo.environ import log10

import idaes.logger as idaeslog

# Grab the scaling utilities
from proteuslib.flowsheets.full_treatment_train.electrolyte_scaling_utils import (
    approximate_chemical_state_args,
    calculate_chemical_scaling_factors)

__author__ = "Srikanth Allu"

# Configuration dictionary
simple_softening_thermo_config = {
    "components": {
        'H2O': {"type": Solvent,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (18.0153, pyunits.g/pyunits.mol),
                    "pressure_crit": (220.64E5, pyunits.Pa),
                    "temperature_crit": (647, pyunits.K),
                    # Comes from Perry's Handbook:  p. 2-98
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.459, pyunits.kmol*pyunits.m**-3),
                        '2': (0.30542, pyunits.dimensionless),
                        '3': (647.13, pyunits.K),
                        '4': (0.081, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-285.830, pyunits.kJ/pyunits.mol),
                    "enth_mol_form_vap_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    # Comes from Perry's Handbook:  p. 2-174
                    "cp_mol_liq_comp_coeff": {
                        '1': (2.7637E5, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (-2.0901E3, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (8.125, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (-1.4116E-2, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (9.3701E-6, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "cp_mol_ig_comp_coeff": {
                        'A': (30.09200, pyunits.J/pyunits.mol/pyunits.K),
                        'B': (6.832514, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**-1),
                        'C': (6.793435, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**-2),
                        'D': (-2.534480, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**-3),
                        'E': (0.082139, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**2),
                        'F': (-250.8810, pyunits.kJ/pyunits.mol),
                        'G': (223.3967, pyunits.J/pyunits.mol/pyunits.K),
                        'H': (0, pyunits.kJ/pyunits.mol)},
                    "entr_mol_form_liq_comp_ref": (69.95, pyunits.J/pyunits.K/pyunits.mol),
                    "pressure_sat_comp_coeff": {
                        'A': (4.6543, None),  # [1], temperature range 255.9 K - 373 K
                        'B': (1435.264, pyunits.K),
                        'C': (-64.848, pyunits.K)}
                                },
                    # End parameter_data
                    },
        'Ca(OH)2': {  "type": Solute,  "valid_phase_types": PT.aqueousPhase,
                    # Define the methods used to calculate the following properties
                    "dens_mol_liq_comp": Constant,
                    "enth_mol_liq_comp": Constant,
                    "cp_mol_liq_comp": Constant,
                    "entr_mol_liq_comp": Constant,
                    # Parameter data is always associated with the methods defined above
                    "parameter_data": {
                        "mw": (74.093, pyunits.g/pyunits.mol),
                        "dens_mol_liq_comp_coeff": (55, pyunits.kmol*pyunits.m**-3),
                        "enth_mol_form_liq_comp_ref": (-945.53, pyunits.kJ/pyunits.mol),
                        "cp_mol_liq_comp_coeff": (167039, pyunits.J/pyunits.kmol/pyunits.K),
                        "entr_mol_form_liq_comp_ref": (100, pyunits.J/pyunits.K/pyunits.mol)
                            },
                    # End parameter_data
                    },
        'CaCO3': {  "type": Solute,  "valid_phase_types": PT.aqueousPhase,
                    # Define the methods used to calculate the following properties
                    "dens_mol_liq_comp": Constant,
                    "enth_mol_liq_comp": Constant,
                    "cp_mol_liq_comp": Constant,
                    "entr_mol_liq_comp": Constant,
                    # Parameter data is always associated with the methods defined above
                    "parameter_data": {
                        "dens_mol_liq_comp_coeff": (55, pyunits.kmol*pyunits.m**-3),
                        "enth_mol_form_liq_comp_ref": (-945.53, pyunits.kJ/pyunits.mol),
                        "cp_mol_liq_comp_coeff": (167039, pyunits.J/pyunits.kmol/pyunits.K),
                        "entr_mol_form_liq_comp_ref": (100, pyunits.J/pyunits.K/pyunits.mol)
                            },
                    # End parameter_data
                    },
        'Ca(HCO3)2': {  "type": Solute,  "valid_phase_types": PT.aqueousPhase,
                    # Define the methods used to calculate the following properties
                    "dens_mol_liq_comp": Constant,
                    "enth_mol_liq_comp": Constant,
                    "cp_mol_liq_comp": Constant,
                    "entr_mol_liq_comp": Constant,
                    # Parameter data is always associated with the methods defined above
                    "parameter_data": {
                        "dens_mol_liq_comp_coeff": (55, pyunits.kmol*pyunits.m**-3),
                        "enth_mol_form_liq_comp_ref": (-945.53, pyunits.kJ/pyunits.mol),
                        "cp_mol_liq_comp_coeff": (167039, pyunits.J/pyunits.kmol/pyunits.K),
                        "entr_mol_form_liq_comp_ref": (100, pyunits.J/pyunits.K/pyunits.mol)
                            },
                    # End parameter_data
                    },
        'Mg(OH)2': {  "type": Solute,  "valid_phase_types": PT.aqueousPhase,
                    # Define the methods used to calculate the following properties
                    "dens_mol_liq_comp": Constant,
                    "enth_mol_liq_comp": Constant,
                    "cp_mol_liq_comp": Constant,
                    "entr_mol_liq_comp": Constant,
                    # Parameter data is always associated with the methods defined above
                    "parameter_data": {
                        "mw": (74.093, pyunits.g/pyunits.mol),
                        "dens_mol_liq_comp_coeff": (55, pyunits.kmol*pyunits.m**-3),
                        "enth_mol_form_liq_comp_ref": (-945.53, pyunits.kJ/pyunits.mol),
                        "cp_mol_liq_comp_coeff": (167039, pyunits.J/pyunits.kmol/pyunits.K),
                        "entr_mol_form_liq_comp_ref": (100, pyunits.J/pyunits.K/pyunits.mol)
                            },
                    # End parameter_data
                    },
        'Mg(HCO3)2': {  "type": Solute,  "valid_phase_types": PT.aqueousPhase,
                    # Define the methods used to calculate the following properties
                    "dens_mol_liq_comp": Constant,
                    "enth_mol_liq_comp": Constant,
                    "cp_mol_liq_comp": Constant,
                    "entr_mol_liq_comp": Constant,
                    # Parameter data is always associated with the methods defined above
                    "parameter_data": {
                        "dens_mol_liq_comp_coeff": (55, pyunits.kmol*pyunits.m**-3),
                        "enth_mol_form_liq_comp_ref": (-945.53, pyunits.kJ/pyunits.mol),
                        "cp_mol_liq_comp_coeff": (167039, pyunits.J/pyunits.kmol/pyunits.K),
                        "entr_mol_form_liq_comp_ref": (100, pyunits.J/pyunits.K/pyunits.mol)
                            },
                    # End parameter_data
                    },
              },
              # End Component list
        "phases":  {'Liq': {"type": AqueousPhase,
                            "equation_of_state": Ideal},
                    },

        "state_definition": FTPx,
        "state_bounds": {"flow_mol": (0, 50, 100),
                         "temperature": (273.15, 300, 650),
                         "pressure": (5e4, 1e5, 1e6)
                     },
        "pressure_ref": 1e5,
        "temperature_ref": 300,
        "base_units": {"time": pyunits.s,
                       "length": pyunits.m,
                       "mass": pyunits.kg,
                       "amount": pyunits.mol,
                       "temperature": pyunits.K},

    }
    # End simple_softening_thermo_config definition

# This config is REQUIRED to use StoichiometricReactor even if we have no equilibrium reactions
simple_softening_reaction_config = {
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},
    "equilibrium_reactions": {
        "H2O_Kw": {
                "stoichiometry": {("Liq", "H2O"): -1,
                                 ("Liq", "H_+"): 1,
                                 ("Liq", "OH_-"): 1},
               "heat_of_reaction": constant_dh_rxn,
               "equilibrium_constant": van_t_hoff,
               "equilibrium_form": log_power_law_equil,
               "concentration_form": ConcentrationForm.moleFraction,
               "parameter_data": {
                   "dh_rxn_ref": (55.830, pyunits.J/pyunits.mol),
                   "k_eq_ref": (10**-14/55.2/55.2, pyunits.dimensionless),
                   "T_eq_ref": (298, pyunits.K),

                   # By default, reaction orders follow stoichiometry
                   #    manually set reaction order here to override
                   "reaction_order": {("Liq", "H2O"): 0,
                                    ("Liq", "H_+"): 1,
                                    ("Liq", "OH_-"): 1}
                    }
                    # End parameter_data
               },
        "HOCl_Ka": {
                "stoichiometry": {("Liq", "HOCl"): -1,
                                 ("Liq", "H_+"): 1,
                                 ("Liq", "OCl_-"): 1},
               "heat_of_reaction": constant_dh_rxn,
               "equilibrium_constant": van_t_hoff,
               "equilibrium_form": log_power_law_equil,
               "concentration_form": ConcentrationForm.moleFraction,
               "parameter_data": {
                   "dh_rxn_ref": (13.8, pyunits.J/pyunits.mol),
                   "k_eq_ref": (10**-7.6422/55.2, pyunits.dimensionless),
                   "T_eq_ref": (298, pyunits.K),

                   # By default, reaction orders follow stoichiometry
                   #    manually set reaction order here to override
                   "reaction_order": {("Liq", "HOCl"): -1,
                                    ("Liq", "H_+"): 1,
                                    ("Liq", "OCl_-"): 1}
                    }
                    # End parameter_data
               }
               # End R2
         }
         # End equilibrium_reactions
    }
    # End reaction_config definition

# Get default solver for testing
solver = get_solver()

def build_simple_softening_unit(model,
                                mg_per_L_CaOH2_added = 2,
                                inlet_water_density_kg_per_L = 1,
                                inlet_temperature_K = 298,
                                inlet_pressure_Pa = 101325,
                                inlet_flow_mol_per_s = 10):
    model.fs.softening_thermo_params = GenericParameterBlock(default=simple_softening_thermo_config)
    model.fs.softening_rxn_params = GenericReactionParameterBlock(
            default={"property_package": model.fs.softening_thermo_params, **softening_reaction_config})
    model.fs.softening_unit = StoichiometricReactor(default={
                "property_package": model.fs.softening_thermo_params,
                "reaction_package": model.fs.softening_rxn_params,
                "has_heat_transfer": False,
                "has_heat_of_reaction": False,
                "energy_balance_type": EnergyBalanceType.none,
                "has_pressure_change": False})

    model.fs.unit.inlet.mole_frac_comp[0, "Mg(HCO3)2"].fix( 0.00003 )
    model.fs.unit.inlet.mole_frac_comp[0, "Ca(HCO3)2"].fix( 0.00003 )
    model.fs.unit.inlet.mole_frac_comp[0, "Mg(OH)2"].fix( 0. )
    model.fs.unit.inlet.mole_frac_comp[0, "CaCO3"].fix( 0. )

    total_molar_density = inlet_water_density_kg_per_L/18*1000 #mol/L

    model.fs.simple_softening_unit.inlet.mole_frac_comp[0, "OCl_-"].fix( total_chlorine_inlet/total_molar_density )
    model.fs.simple_softening_unit.inlet.mole_frac_comp[0, "Na_+"].fix( total_chlorine_inlet/total_molar_density )

    model.fs.unit.outlet.mole_frac_comp[0, "Ca(HCO3)2"].fix( 0.000015 )
    model.fs.unit.outlet.mole_frac_comp[0, "Mg(HCO3)2"].fix( 0.000015 )
    model.fs.unit.outlet.mole_frac_comp[0, "Ca(OH)2"].fix( 0.0000003 )

    # Perform a summation of all non-H2O molefractions to find the H2O molefraction
    sum = 0
    for i in model.fs.simple_softening_unit.inlet.mole_frac_comp:
        # NOTE: i will be a tuple with format (time, component)
        if i[1] != "H2O":
            sum += value(model.fs.simple_softening_unit.inlet.mole_frac_comp[i[0], i[1]])

    model.fs.simple_softening_unit.inlet.mole_frac_comp[0, "H2O"].fix( 1-sum )

    model.fs.simple_softening_unit.inlet.pressure.fix(inlet_pressure_Pa)
    model.fs.simple_softening_unit.inlet.temperature.fix(inlet_temperature_K)
    model.fs.simple_softening_unit.inlet.flow_mol.fix(inlet_flow_mol_per_s)

def initialize_chlorination_example(unit, state_args, user_scaling=True, debug_out=False):
    solver.options['bound_push'] = 1e-10
    solver.options['mu_init'] = 1e-6

    if user_scaling == True:
        solver.options["nlp_scaling_method"] = "user-scaling"

    unit.inlet.mole_frac_comp[0, "Ca(OH)2"].fix()

    if debug_out == True:
        unit.initialize(state_args=state_args, optarg=solver.options, outlvl=idaeslog.DEBUG)
    else:
        unit.initialize(state_args=state_args, optarg=solver.options)

    unit.inlet.mole_frac_comp[0, "Ca(OH)2"].unfix()

    iscale.constraint_autoscale_large_jac(unit)

def display_results_of_chlorination(softening_unit):
    print()
    print("=========== Softening Results ============")
    print("Outlet Temperature:       \t" + str(softening_unit.outlet.temperature[0].value))
    print("Outlet Pressure:          \t" + str(softening_unit.outlet.pressure[0].value))
    print("Outlet FlowMole:          \t" + str(softening_unit.outlet.flow_mol[0].value))
    print()
    total_molar_density = \
        value(softening_unit.control_volume.properties_out[0.0].dens_mol_phase['Liq'])/1000
    pH = -value(log10(softening_unit.outlet.mole_frac_comp[0, "H_+"]*total_molar_density))
    print("pH at Outlet:             \t" + str(pH))
    total_salt = value(softening_unit.outlet.mole_frac_comp[0, "Na_+"])*total_molar_density*23
    total_salt += value(softening_unit.outlet.mole_frac_comp[0, "Cl_-"])*total_molar_density*35.4
    psu = total_salt/(total_molar_density*18)*1000
    print("Salinity (PSU):           \t" + str(psu))
    print("NaOCl Dosing Rate (mg/s): \t" + str(softening_unit.dosing_rate.value))
    print("Free Chlorine (mg/L):     \t" + str(softening_unit.free_chlorine.value))
    print("\tDistribution:")
    hocl = (value(softening_unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq","HOCl"])/1000)/(chlorination_unit.free_chlorine.value/70900)
    print("\t % HOCl: \t" + str(hocl*100))
    ocl = (value(softening_unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq","OCl_-"])/1000)/(chlorination_unit.free_chlorine.value/70900)
    print("\t % OCl-: \t" + str(ocl*100))
    print("-------------------------------------------")
    print()

def run_chlorination_example():
    model = ConcreteModel()
    model.fs = FlowsheetBlock(default={"dynamic": False})

    build_simple_softening_unit(model, mg_per_L_CaOH2_added = 2)
    state_args, stoich_extents = approximate_chemical_state_args(model.fs.simple_softening_unit,
                                model.fs.simple_softening_rxn_params, simple_softening_reaction_config)

    calculate_chemical_scaling_factors(model.fs.simple_softening_unit,
                                model.fs.simple_softening_thermo_params,
                                model.fs.simple_softening_rxn_params, state_args)

    initialize_chlorination_example(model.fs.simple_softening_unit, state_args)

    solve_with_user_scaling(model, tee=True, bound_push=1e-10, mu_init=1e-6)

    display_results_of_chlorination(model.fs.simple_softening_unit)

    return model

# This example will try to find the dosing rate needed to yield 2mg/L free chlorine
def run_softening_constrained_outlet_example():
    model = ConcreteModel()
    model.fs = FlowsheetBlock(default={"dynamic": False})

    # Give bad initial guess for the dosing rate
    build_simple_softening_unit(model, mg_per_L_CaOH2_added = 1)
    state_args, stoich_extents = approximate_chemical_state_args(model.fs.simple_softening_unit,
                                model.fs.simple_softening_rxn_params, simple_softening_reaction_config)

    calculate_chemical_scaling_factors(model.fs.simple_softening_unit,
                                model.fs.simple_softening_thermo_params,
                                model.fs.simple_softening_rxn_params, state_args)

    initialize_chlorination_example(model.fs.simple_softening_unit, state_args)

    #Redefine the constraints and fixed vars here
    model.fs.simple_softening_unit.dosing_rate.unfix()
    model.fs.simple_softening_unit.free_chlorine.fix(2)

    solve_with_user_scaling(model, tee=True, bound_push=1e-10, mu_init=1e-6)

    display_results_of_chlorination(model.fs.softening_unit)

    return model

if __name__ == "__main__":
    model = run_chlorination_example()
    model = run_chlorination_constrained_outlet_example()
