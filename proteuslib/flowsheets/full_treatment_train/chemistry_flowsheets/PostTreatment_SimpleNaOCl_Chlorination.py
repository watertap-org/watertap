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

""" Flowsheet built as Separator and Chlorination process """

# Importing the object for units from pyomo
from pyomo.environ import units as pyunits

# Imports from idaes core
from idaes.core import AqueousPhase
from idaes.core.components import Solvent, Solute, Cation, Anion
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

# Import built-in van't Hoff function
from idaes.generic_models.properties.core.reactions.equilibrium_constant import van_t_hoff

# Import specific pyomo objects
from pyomo.environ import (ConcreteModel,
                           SolverStatus,
                           TerminationCondition,
                           value,
                           Suffix)

from idaes.core.util import scaling as iscale
from idaes.core.util.initialization import fix_state_vars, revert_state_vars

# Import pyomo methods to check the system units
from pyomo.util.check_units import assert_units_consistent


from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof
from idaes.core.util import get_solver

# Import the idaes objects for Generic Properties and Reactions
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)
from idaes.generic_models.properties.core.generic.generic_reaction import (
        GenericReactionParameterBlock)

# Import the idaes object for the EquilibriumReactor unit model
from idaes.generic_models.unit_models.equilibrium_reactor import EquilibriumReactor

# Import the core idaes objects for Flowsheets and types of balances
from idaes.core import FlowsheetBlock

# Import log10 function from pyomo
from pyomo.environ import log10

import idaes.logger as idaeslog

__author__ = "Austin Ladshaw"

# Configuration dictionary
simple_naocl_thermo_config = {
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
        'H_+': {"type": Cation, "charge": 1,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (1.00784, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.459, pyunits.kmol*pyunits.m**-3),
                        '2': (0.30542, pyunits.dimensionless),
                        '3': (647.13, pyunits.K),
                        '4': (0.081, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-230.000, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (2.7637E5, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (-2.0901E3, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (8.125, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (-1.4116E-2, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (9.3701E-6, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (-10.75, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'OH_-': {"type": Anion,
                "charge": -1,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (17.008, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.459, pyunits.kmol*pyunits.m**-3),
                        '2': (0.30542, pyunits.dimensionless),
                        '3': (647.13, pyunits.K),
                        '4': (0.081, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-230.000, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (2.7637E5, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (-2.0901E3, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (8.125, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (-1.4116E-2, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (9.3701E-6, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (-10.75, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'HOCl': {"type": Solute, "valid_phase_types": PT.aqueousPhase,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (52.46, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (4.985, pyunits.kmol*pyunits.m**-3),
                        '2': (0.36, pyunits.dimensionless),
                        '3': (1464.06, pyunits.K),
                        '4': (0.739, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-120.9, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (83993.8, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (142, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'OCl_-': {"type": Anion, "charge": -1,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (51.46, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (4.985, pyunits.kmol*pyunits.m**-3),
                        '2': (0.36, pyunits.dimensionless),
                        '3': (1464.06, pyunits.K),
                        '4': (0.739, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-107.1, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (83993.8, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (42, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    }
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

        # Inherent reactions
        "inherent_reactions": {
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
                   # End R1
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
                   # End R4
             }
             # End equilibrium_reactions
    }
    # End simple_naocl_thermo_config definition

# This config is REQUIRED to use EquilibriumReactor even if we have no equilibrium reactions
reaction_config = {
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
         }
         # End equilibrium_reactions
    }
    # End reaction_config definition

# Get default solver for testing
solver = get_solver()

def build_simple_naocl_chlorination_unit(model,
                                mg_per_L_NaOCl_added = 2,
                                inlet_water_density_kg_per_L = 1,
                                inlet_temperature_K = 298,
                                inlet_pressure_Pa = 101325,
                                inlet_flow_mol_per_s = 10):
    model.fs.simple_naocl_thermo_params = GenericParameterBlock(default=simple_naocl_thermo_config)
    model.fs.simple_naocl_rxn_params = GenericReactionParameterBlock(
            default={"property_package": model.fs.simple_naocl_thermo_params, **reaction_config})
    model.fs.simple_naocl_unit = EquilibriumReactor(default={
            "property_package": model.fs.simple_naocl_thermo_params,
            "reaction_package": model.fs.simple_naocl_rxn_params,
            "has_rate_reactions": False,
            "has_equilibrium_reactions": False,
            "has_heat_transfer": False,
            "has_heat_of_reaction": False,
            "has_pressure_change": False})

    model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "H_+"].fix( 0. )
    model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "OH_-"].fix( 0. )
    model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "HOCl"].fix( 0. )

    total_molar_density = inlet_water_density_kg_per_L/18*1000 #mol/L

    # Free Chlorine (mg-Cl2/L) = total_chlorine_inlet (mol/L) * 70,900
    #       Assumes chlorine is added as NaOCl (we don't track Na, because it doesn't matter)
    free_chlorine_added = mg_per_L_NaOCl_added #mg/L as Cl2
    total_chlorine_inlet = free_chlorine_added/70900 # mol/L
    total_molar_density+=total_chlorine_inlet

    model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "OCl_-"].fix( total_chlorine_inlet/total_molar_density )

    # Perform a summation of all non-H2O molefractions to find the H2O molefraction
    sum = 0
    for i in model.fs.simple_naocl_unit.inlet.mole_frac_comp:
        # NOTE: i will be a tuple with format (time, component)
        if i[1] != "H2O":
            sum += value(model.fs.simple_naocl_unit.inlet.mole_frac_comp[i[0], i[1]])

    model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "H2O"].fix( 1-sum )

    model.fs.simple_naocl_unit.inlet.pressure.fix(inlet_pressure_Pa)
    model.fs.simple_naocl_unit.inlet.temperature.fix(inlet_temperature_K)
    model.fs.simple_naocl_unit.inlet.flow_mol.fix(inlet_flow_mol_per_s)

# # TODO: Move this to utils?
def approximate_chemical_state_args(unit, rxn_params, reaction_config, contains_stoich_reactions=False):
    state_args = {}
    stoich_extents = {}

    # Set bulk values to their inlets
    state_args['pressure'] = unit.control_volume.properties_in[0.0].pressure.value
    state_args['temperature'] = unit.control_volume.properties_in[0.0].temperature.value
    state_args['flow_mol'] = unit.control_volume.properties_in[0.0].flow_mol.value

    # Set species based on inlets (and outlets for stoich reaction)
    state_args['mole_frac_comp'] = {}
    min = 1e-6
    for i in unit.control_volume.properties_in[0.0].mole_frac_comp:
        # Set state args to inlets on first pass
        if unit.inlet.mole_frac_comp[0, i].value > min:
            state_args['mole_frac_comp'][i] = unit.inlet.mole_frac_comp[0, i].value
        else:
            state_args['mole_frac_comp'][i] = min

    # Iterate through outlet mole fractions and note the fixed variables
    fixed = {}
    for i, species in unit.outlet.mole_frac_comp:
        if unit.outlet.mole_frac_comp[i, species].is_fixed():
            fixed[species] = True
            state_args['mole_frac_comp'][species] = unit.outlet.mole_frac_comp[0, species].value

    # Checking stoich reactions
    was_OH_changed = False
    was_H_changed = False
    if contains_stoich_reactions == True:
        for rid in rxn_params.rate_reaction_idx:
            #First loop establishes reaction extent
            extent = 0
            for phase, species in reaction_config["rate_reactions"][rid]["stoichiometry"]:
                # If a species here has its outlet fixed, then the difference between
                #   that species outlet and inlet values should serve as the basis for
                #   setting the values of the other species used in that reaction
                if species in fixed:
                    extent = unit.inlet.mole_frac_comp[0, species].value \
                            - unit.outlet.mole_frac_comp[0, species].value
                stoich_extents[rid] = extent

            # Loop again to set values based on extent
            for phase, species in reaction_config["rate_reactions"][rid]["stoichiometry"]:
                state_args['mole_frac_comp'][species] = unit.inlet.mole_frac_comp[0, species].value \
                        + extent*reaction_config["rate_reactions"][rid]["stoichiometry"][phase, species]
                if species == "H_+" and extent != 0.0:
                    was_H_changed = True
                if species == "OH_-" and extent != 0.0:
                    was_OH_changed = True

    # Lastly, we need for correct OH and/or H if they are changed by a stoich reaction
    if was_H_changed == False and was_OH_changed == False:
        state_args['mole_frac_comp']['H_+'] = 10**-7/55.6
        state_args['mole_frac_comp']['OH_-'] = 10**-7/55.6
    elif was_H_changed == False and was_OH_changed == True:
        state_args['mole_frac_comp']['H_+'] = 10**-14/(state_args['mole_frac_comp']['OH_-']*55.6)/55.6
    elif was_H_changed == True and was_OH_changed == False:
        state_args['mole_frac_comp']['OH_-'] = 10**-14/(state_args['mole_frac_comp']['H_+']*55.6)/55.6
    else:
        # Uncertain what to do at this point
        pass

    return state_args, stoich_extents

# # TODO: Move this to utils?
def calculate_chemical_scaling_factors_for_inherent_log_reactions(unit, thermo_params):
    try:
        # Iterate through the reactions to set appropriate eps values
        factor = 1e-4
        for rid in thermo_params.inherent_reaction_idx:
            scale = value(unit.control_volume.properties_out[0.0].k_eq[rid].expr)
            # Want to set eps in some fashion similar to this
            if scale < 1e-16:
                thermo_params.component("reaction_"+rid).eps.value = scale*factor
            else:
                thermo_params.component("reaction_"+rid).eps.value = 1e-16*factor

        for i in unit.control_volume.inherent_reaction_extent_index:
            scale = value(unit.control_volume.properties_out[0.0].k_eq[i[1]].expr)
            iscale.set_scaling_factor(unit.control_volume.inherent_reaction_extent[0.0,i[1]], 10/scale)
            iscale.constraint_scaling_transform(unit.control_volume.properties_out[0.0].
                    inherent_equilibrium_constraint[i[1]], 0.1)
    except:
        pass

# # TODO: Move this to utils?
def calculate_chemical_scaling_factors_for_equilibrium_log_reactions(unit, rxn_params):
    try:
        # Equilibrium reactions have eps in the 'simple_naocl_rxn_params'
        factor = 1e-4
        for rid in rxn_params.equilibrium_reaction_idx:
            if rid != "dummy":
                scale = value(unit.control_volume.reactions[0.0].k_eq[rid].expr)
                # Want to set eps in some fashion similar to this
                if scale < 1e-16:
                    rxn_params.component("reaction_"+rid).eps.value = scale*factor
                else:
                    rxn_params.component("reaction_"+rid).eps.value = 1e-16*factor

        for i in unit.control_volume.equilibrium_reaction_extent_index:
            if i[1] != "dummy":
                scale = value(unit.control_volume.reactions[0.0].k_eq[i[1]].expr)
                iscale.set_scaling_factor(unit.control_volume.equilibrium_reaction_extent[0.0,i[1]], 10/scale)
                iscale.constraint_scaling_transform(unit.control_volume.reactions[0.0].
                        equilibrium_constraint[i[1]], 0.1)
    except:
        pass

# # TODO: Move this to utils?
def calculate_chemical_scaling_factors_for_material_balances(unit):
    # Next, try adding scaling for species
    min = 1e-6
    for i in unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
        # i[0] = phase, i[1] = species
        if unit.inlet.mole_frac_comp[0, i[1]].value > min:
            scale = unit.inlet.mole_frac_comp[0, i[1]].value
        else:
            scale = min
        iscale.set_scaling_factor(unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]], 10/scale)
        iscale.set_scaling_factor(unit.control_volume.properties_out[0.0].mole_frac_phase_comp[i], 10/scale)
        iscale.set_scaling_factor(unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i], 10/scale)
        iscale.constraint_scaling_transform(
            unit.control_volume.properties_out[0.0].component_flow_balances[i[1]], 10/scale)
        iscale.constraint_scaling_transform(unit.control_volume.material_balances[0.0,i[1]], 10/scale)

# # TODO: Move this to utils?
def calculate_chemical_scaling_factors(unit, thermo_params, rxn_params, state_args, output_jac=False):
    calculate_chemical_scaling_factors_for_inherent_log_reactions(unit, thermo_params)
    calculate_chemical_scaling_factors_for_equilibrium_log_reactions(unit, rxn_params)
    calculate_chemical_scaling_factors_for_material_balances(unit)
    iscale.calculate_scaling_factors(unit)

    flags = fix_state_vars(unit.control_volume.properties_out, state_args)
    revert_state_vars(unit.control_volume.properties_out, flags)

    iscale.constraint_autoscale_large_jac(unit)

    if output_jac == True:
        jac, nlp = iscale.get_jacobian(unit, scaled=True)
        print("Extreme Jacobian entries:")
        for i in iscale.extreme_jacobian_entries(jac=jac, nlp=nlp, large=100):
            print(f"    {i[0]:.2e}, [{i[1]}, {i[2]}]")
        print("Unscaled constraints:")
        for c in iscale.unscaled_constraints_generator(unit):
            print(f"    {c}")
        print("Scaled constraints by factor:")
        for c, s in iscale.constraints_with_scale_factor_generator(unit):
            print(f"    {c}, {s}")
        print("Badly scaled variables:")
        for v, sv in iscale.badly_scaled_var_generator(unit, large=1e2, small=1e-2, zero=1e-12):
            print(f"    {v} -- {sv} -- {iscale.get_scaling_factor(v)}")
        print(f"Jacobian Condition Number: {iscale.jacobian_cond(jac=jac):.2e}")

def initialize_chlorination_example(unit, state_args, user_scaling=True, debug_out=False):
    solver.options['bound_push'] = 1e-10
    solver.options['mu_init'] = 1e-6

    if user_scaling == True:
        solver.options["nlp_scaling_method"] = "user-scaling"

    if debug_out == True:
        unit.initialize(state_args=state_args, optarg=solver.options, outlvl=idaeslog.DEBUG)
    else:
        unit.initialize(state_args=state_args, optarg=solver.options)

def display_results_of_chlorination(chlorination_unit):
    print()
    print("=========== Chlorination Results ============")
    print("Outlet Temperature:  \t" + str(chlorination_unit.outlet.temperature[0].value))
    print("Outlet Pressure:     \t" + str(chlorination_unit.outlet.pressure[0].value))
    print("Outlet FlowMole:     \t" + str(chlorination_unit.outlet.flow_mol[0].value))
    print()
    total_molar_density = \
        value(chlorination_unit.control_volume.properties_out[0.0].dens_mol_phase['Liq'])/1000
    pH = -value(log10(chlorination_unit.outlet.mole_frac_comp[0, "H_+"]*total_molar_density))
    print("pH at Outlet:        \t" + str(pH))
    hypo_remaining = value(chlorination_unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq","HOCl"])/1000
    hypo_remaining += value(chlorination_unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq","OCl_-"])/1000
    hypo_remaining = hypo_remaining*70900
    print("Free Chlorine (mg/L):\t" + str(hypo_remaining))
    print("\tDistribution:")
    hocl = (value(chlorination_unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq","HOCl"])/1000)/(hypo_remaining/70900)
    print("\t % HOCl: \t" + str(hocl*100))
    ocl = (value(chlorination_unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq","OCl_-"])/1000)/(hypo_remaining/70900)
    print("\t % OCl-: \t" + str(ocl*100))
    print("-------------------------------------------")
    print()

def run_chlorination_example():
    model = ConcreteModel()
    model.fs = FlowsheetBlock(default={"dynamic": False})

    build_simple_naocl_chlorination_unit(model, mg_per_L_NaOCl_added = 2)
    state_args, stoich_extents = approximate_chemical_state_args(model.fs.simple_naocl_unit,
                                model.fs.simple_naocl_rxn_params, reaction_config)

    calculate_chemical_scaling_factors(model.fs.simple_naocl_unit,
                                model.fs.simple_naocl_thermo_params,
                                model.fs.simple_naocl_rxn_params, state_args)

    initialize_chlorination_example(model.fs.simple_naocl_unit, state_args)

    display_results_of_chlorination(model.fs.simple_naocl_unit)


    return model


if __name__ == "__main__":
    model = run_chlorination_example()
