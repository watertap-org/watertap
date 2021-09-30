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

""" Simple Softening process with addition of Lime"""

# Importing the object for units from pyomo
from pyomo.environ import units as pyunits

# Imports from idaes core
from idaes.core import AqueousPhase
from idaes.core.components import Solvent, Solute
from idaes.core.phases import PhaseType as PT

# Imports from idaes generic models
import idaes.generic_models.properties.core.pure.Perrys as Perrys
from idaes.generic_models.properties.core.pure.ConstantProperties import Constant
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ideal import Ideal
from idaes.generic_models.properties.core.reactions.rate_constant import arrhenius
from idaes.generic_models.properties.core.reactions.rate_forms import power_law_rate

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
from idaes.core import (AqueousPhase,
                        FlowsheetBlock,
                        EnergyBalanceType)

# Import log10 function from pyomo
from pyomo.environ import log10

import idaes.logger as idaeslog

# Grab the scaling utilities
from proteuslib.flowsheets.full_treatment_train.electrolyte_scaling_utils import (
    approximate_chemical_state_args,
    calculate_chemical_scaling_factors)

__author__ = "Srikanth Allu"

# Configuration dictionary
softening_thermo_config = {
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
    # End softening_thermo_config definition

# This config is REQUIRED to use StoichiometricReactor even if we have no equilibrium reactions
softening_reaction_config = {
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},
    "equilibrium_reactions": {
        "dummy": {
                "stoichiometry": {},
                "equilibrium_form": log_power_law_equil,
               }
         },
         # End equilibrium_reactions
    "rate_reactions": {
        "R1": {"stoichiometry": {("Liq", "Ca(HCO3)2"): -1,
                                 ("Liq", "Ca(OH)2"): -1,
                                 ("Liq", "CaCO3"): 2,
                                 ("Liq", "H2O"): 2},
               "heat_of_reaction": constant_dh_rxn,
               "rate_constant" : arrhenius,
               "rate_form" : power_law_rate,
               "concentration_form" : ConcentrationForm.moleFraction,
               "parameter_data": {
                   "arrhenius_const" : (1, pyunits.mol/pyunits.m**3/pyunits.s),
                   "energy_activation" : (0, pyunits.J/pyunits.mol),
                   "dh_rxn_ref": (0, pyunits.J/pyunits.mol)
              }
         },
        "R2": {"stoichiometry": {("Liq", "Mg(HCO3)2"): -1,
                                 ("Liq", "Ca(OH)2"): -2,
                                 ("Liq", "CaCO3"): 2,
                                 ("Liq", "Mg(OH)2"): 1,
                                 ("Liq", "H2O"): 2},
               "heat_of_reaction": constant_dh_rxn,
               "rate_constant" : arrhenius,
               "rate_form" : power_law_rate,
               "concentration_form" : ConcentrationForm.moleFraction,
               "parameter_data": {
                   "arrhenius_const" : (1, pyunits.mol/pyunits.m**3/pyunits.s),
                   "energy_activation" : (0, pyunits.J/pyunits.mol),
                   "dh_rxn_ref": (0, pyunits.J/pyunits.mol)
              }
         }
    }
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
    model.fs.softening_thermo_params = GenericParameterBlock(default=softening_thermo_config)
    model.fs.softening_rxn_params = GenericReactionParameterBlock(
            default={"property_package": model.fs.softening_thermo_params, **softening_reaction_config})
    model.fs.softening_unit = StoichiometricReactor(default={
                "property_package": model.fs.softening_thermo_params,
                "reaction_package": model.fs.softening_rxn_params,
                "has_heat_transfer": False,
                "has_heat_of_reaction": False,
                "energy_balance_type": EnergyBalanceType.none,
                "has_pressure_change": False})

    model.fs.softening_unit.inlet.mole_frac_comp[0, "Mg(HCO3)2"].fix( 0.00003 )
    model.fs.softening_unit.inlet.mole_frac_comp[0, "Ca(HCO3)2"].fix( 0.00003 )
    model.fs.softening_unit.inlet.mole_frac_comp[0, "Mg(OH)2"].fix( 0. )
    model.fs.softening_unit.inlet.mole_frac_comp[0, "CaCO3"].fix( 0. )
    model.fs.softening_unit.inlet.mole_frac_comp[0, "Ca(OH)2"].fix( 0.) # temperory fix for initializing

    total_molar_density = inlet_water_density_kg_per_L/18*1000 #mol/L

    model.fs.softening_unit.outlet.mole_frac_comp[0, "Ca(HCO3)2"].fix( 0.000015 )
    model.fs.softening_unit.outlet.mole_frac_comp[0, "Mg(HCO3)2"].fix( 0.000015 )
    model.fs.softening_unit.outlet.mole_frac_comp[0, "Ca(OH)2"].fix( 0.0000003 )

    # Perform a summation of all non-H2O molefractions to find the H2O molefraction
    sum = 0
    for i in model.fs.softening_unit.inlet.mole_frac_comp:
        # NOTE: i will be a tuple with format (time, component)
        if i[1] != "H2O":
            sum += value(model.fs.softening_unit.inlet.mole_frac_comp[i[0], i[1]])

    model.fs.softening_unit.inlet.mole_frac_comp[0, "H2O"].fix( 1-sum )
    model.fs.softening_unit.inlet.mole_frac_comp[0, "Ca(OH)2"].unfix() # unfix for initializing
    model.fs.softening_unit.inlet.mole_frac_comp.pprint()

    #model.fs.softening_unit.inlet.mole_frac_comp[0, "H2O"].fix( 0.99991 )

    model.fs.softening_unit.inlet.pressure.fix(inlet_pressure_Pa)
    model.fs.softening_unit.inlet.temperature.fix(inlet_temperature_K)
    model.fs.softening_unit.inlet.flow_mol.fix(inlet_flow_mol_per_s)

    model.fs.softening_unit.outlet.temperature.fix(inlet_temperature_K)
   
    check_dof(model)

def initialize_softening_example(unit, state_args, user_scaling=True, debug_out=True):
    check_dof(unit)
    solver.options['bound_push'] = 1e-10
    solver.options['mu_init'] = 1e-6

    if user_scaling == True:
        solver.options["nlp_scaling_method"] = "user-scaling"

    #unit.inlet.mole_frac_comp[0, "Ca(OH)2"].fix()

    if debug_out == True:
        solve_with_user_scaling(unit, tee=True, bound_push=1e-10, mu_init=1e-6)
    #    unit.initialize(state_args=state_args, optarg=solver.options, outlvl=idaeslog.DEBUG)
    else:
        solve_with_user_scaling(unit, tee=False, bound_push=1e-10, mu_init=1e-6)
    #    unit.initialize(state_args=state_args, optarg=solver.options)

    #unit.inlet.mole_frac_comp[0, "Ca(OH)2"].unfix()

    iscale.constraint_autoscale_large_jac(unit)
    check_dof(unit)

def display_results_of_softening(softening_unit):
    print()
    print("=========== Softening Results ============")
    print("Outlet Temperature:       \t" + str(softening_unit.outlet.temperature[0].value))
    print("Outlet Pressure:          \t" + str(softening_unit.outlet.pressure[0].value))
    print("Outlet FlowMole:          \t" + str(softening_unit.outlet.flow_mol[0].value))
    print()
    softening_unit.inlet.mole_frac_comp.pprint()
    softening_unit.outlet.mole_frac_comp.pprint()
    total_molar_density = \
        value(softening_unit.control_volume.properties_out[0.0].dens_mol_phase['Liq'])/1000

    total_hardness1 = 50000*2* softening_unit.outlet.mole_frac_comp[0, "Ca(HCO3)2"].value*total_molar_density
    total_hardness2 = 50000*2* softening_unit.outlet.mole_frac_comp[0, "Mg(HCO3)2"].value*total_molar_density
    print("hardness at Outlet:             \t" + str(total_hardness1+total_hardness2))
    print("-------------------------------------------")
    print()

def run_softening_example():
    model = ConcreteModel()
    model.fs = FlowsheetBlock(default={"dynamic": False})

    build_simple_softening_unit(model, mg_per_L_CaOH2_added = 2)
    state_args, stoich_extents = approximate_chemical_state_args(model.fs.softening_unit,
                                model.fs.softening_rxn_params, softening_reaction_config)

    calculate_chemical_scaling_factors(model.fs.softening_unit,
                                model.fs.softening_thermo_params,
                                model.fs.softening_rxn_params, state_args)

    initialize_softening_example(model.fs.softening_unit, state_args)

    solve_with_user_scaling(model, tee=True, bound_push=1e-10, mu_init=1e-6)

    display_results_of_softening(model.fs.softening_unit)

    return model

def run_softening_constrained_outlet_example():
    model = ConcreteModel()
    model.fs = FlowsheetBlock(default={"dynamic": False})

    build_simple_softening_unit(model, mg_per_L_CaOH2_added = 1)
    state_args, stoich_extents = approximate_chemical_state_args(model.fs.softening_unit,
                                model.fs.softening_rxn_params, softening_reaction_config)

    calculate_chemical_scaling_factors(model.fs.softening_unit,
                                model.fs.softening_thermo_params,
                                model.fs.softening_rxn_params, state_args)

    initialize_softening_example(model.fs.softening_unit, state_args)

    solve_with_user_scaling(model, tee=True, bound_push=1e-10, mu_init=1e-6)

    display_results_of_softening(model.fs.softening_unit)

    return model

if __name__ == "__main__":
    model = run_softening_example()
    #model = run_softening_constrained_outlet_example()
