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

"""
    This test is to establish that the core chemistry packages in IDAES
    can solve a complex chlorination problem wherein the hypochlorite
    species get scavenged by NH3 in solution. NH3 binds with the chlorine
    in solution to dilute the effectiveness of the chlorination process,
    thus, it is important to have a model that can account for these
    impacts.

    Additionally, this test will establish that we can mix-and-match
    both inherent and equilibrium expressions in the system. This may
    be an important factor from a modeling perspective as not all
    unit processes will have all inherent reactions, thus, we may desire
    to add/remove reactions as equilibrium statements at will.

    Inherent Reactions:
        H2O <---> H + OH
        NH4 <---> H + NH3
        HOCl <---> H + OCl

    True NH3/HOCl Reactions (not included):
        NH3 + HOCl <---> NH2Cl + H2O
        NH2Cl + HOCl <---> NHCl2 + H2O
        NHCl2 + HOCl <---> NCl3 + H2O

    NOTE: Although there are several cascading reactions for NH3 to scavenge
        HOCl, we only care about the fact that every 1 NH3 can take upto
        3 HOCl, so here I am lumping all ammonia chloride formation reactions
        into a single reaction. This vastly improves the convergence behavior
        of the full system of equations.

    Effective NH3/HOCl Reaction (included):
        NH3 + 3 HOCl <---> NCl3 + 3 H2O
"""
# Importing testing libraries
import pytest

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

# Import k-value functions
from idaes.generic_models.properties.core.reactions.equilibrium_constant import (
    gibbs_energy,
    van_t_hoff)

# Import built-in van't Hoff function
from idaes.generic_models.properties.core.reactions.equilibrium_constant import van_t_hoff

# Import specific pyomo objects
from pyomo.environ import (ConcreteModel,
                           SolverStatus,
                           TerminationCondition,
                           value,
                           Suffix)

from idaes.core.util import scaling as iscale

# Import pyomo methods to check the system units
from pyomo.util.check_units import assert_units_consistent

# Import idaes methods to check the model during construction
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              fixed_variables_set,
                                              activated_constraints_set,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables)

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
thermo_config = {
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
        'NH3': {"type": Solute,
                "dens_mol_liq_comp": Perrys,
                "enth_mol_liq_comp": Perrys,
                "cp_mol_liq_comp": Perrys,
                "entr_mol_liq_comp": Perrys,
                "parameter_data": {
                    "mw": (17.031, pyunits.g/pyunits.mol),
                    "pressure_crit": (113E5, pyunits.Pa),
                    "temperature_crit": (405.4, pyunits.K),
                    "dens_mol_liq_comp_coeff": {
                        '1': (3.5383, pyunits.kmol*pyunits.m**-3),
                        '2': (0.25443, pyunits.dimensionless),
                        '3': (405.65, pyunits.K),
                        '4': (0.2888, pyunits.dimensionless)},
                    "cp_mol_ig_comp_coeff": {
                       'A': (19.99563, pyunits.J/pyunits.mol/pyunits.K),
                       'B': (49.77119, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**-1),
                       'C': (-15.37599, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**-2),
                       'D': (1.921168, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**-3),
                       'E': (0.189174, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**2),
                       'F': (-53.30667, pyunits.kJ/pyunits.mol),
                       'G': (203.8591, pyunits.J/pyunits.mol/pyunits.K),
                       'H': (-45.89806, pyunits.kJ/pyunits.mol)},
                    "cp_mol_liq_comp_coeff": {
                        '1': (71128, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "enth_mol_form_liq_comp_ref": (-80.29, pyunits.kJ/pyunits.mol),
                    "enth_mol_form_vap_comp_ref": (-45.9, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (111, pyunits.J/pyunits.K/pyunits.mol),
                    "entr_mol_form_vap_comp_ref": (192.77, pyunits.J/pyunits.mol),
                    "pressure_sat_comp_coeff": {
                        'A': (4.86886, None),
                        'B': (1113.928, pyunits.K),
                        'C': (-10.409, pyunits.K)}
                                }
                        # End parameter_data
                    },
        'NH4_+': {"type": Cation, "charge": 1,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (18.039, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (3.5383, pyunits.kmol*pyunits.m**-3),
                        '2': (0.25443, pyunits.dimensionless),
                        '3': (405.65, pyunits.K),
                        '4': (0.2888, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-132.5, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (71128, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (113.4, pyunits.J/pyunits.K/pyunits.mol)
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
                    },
        'NH2Cl': {"type": Solute, "valid_phase_types": PT.aqueousPhase,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (51.48, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (4.519, pyunits.kmol*pyunits.m**-3),
                        '2': (0.444, pyunits.dimensionless),
                        '3': (988.9, pyunits.K),
                        '4': (0.37, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (66.9, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (71128, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'NHCl2': {"type": Solute, "valid_phase_types": PT.aqueousPhase,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (85.92, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (4.519, pyunits.kmol*pyunits.m**-3),
                        '2': (0.444, pyunits.dimensionless),
                        '3': (988.9, pyunits.K),
                        '4': (0.37, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (176.6, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (71128, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'NCl3': {"type": Solute, "valid_phase_types": PT.aqueousPhase,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (120.365, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (4.519, pyunits.kmol*pyunits.m**-3),
                        '2': (0.444, pyunits.dimensionless),
                        '3': (988.9, pyunits.K),
                        '4': (0.37, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (299.2, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (71128, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
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
            "NH4_Ka": {
                    "stoichiometry": {("Liq", "NH4_+"): -1,
                                     ("Liq", "H_+"): 1,
                                     ("Liq", "NH3"): 1},
                   "heat_of_reaction": constant_dh_rxn,
                   "equilibrium_constant": van_t_hoff,
                   "equilibrium_form": log_power_law_equil,
                   "concentration_form": ConcentrationForm.moleFraction,
                   "parameter_data": {
                       "dh_rxn_ref": (52.21, pyunits.J/pyunits.mol),
                       "k_eq_ref": (10**-9.2767/55.2, pyunits.dimensionless),
                       "T_eq_ref": (298, pyunits.K),

                       # By default, reaction orders follow stoichiometry
                       #    manually set reaction order here to override
                       "reaction_order": {("Liq", "NH4_+"): -1,
                                        ("Liq", "H_+"): 1,
                                        ("Liq", "NH3"): 1}
                        }
                        # End parameter_data
                   },
                   # End R2
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
    # End thermo_config definition

# NOTE: These reactions are (usually) complete reactions, thus, it may be
#       better to model them as "stoichiometric" reactions for better
#       convergence behavior of the non-linear system
# Define the reaction_config for NH3/HOCl reaction
reaction_config = {
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},
    "equilibrium_reactions": {
        "NCl3_K": {
                "stoichiometry": {("Liq", "NH3"): 1,
                                 ("Liq", "HOCl"): 3,
                                 ("Liq", "NCl3"): -1,
                                 ("Liq", "H2O"): -3},
               "heat_of_reaction": constant_dh_rxn,
               "equilibrium_constant": van_t_hoff,
               "equilibrium_form": log_power_law_equil,
               "concentration_form": ConcentrationForm.moleFraction,
               "parameter_data": {
                   "dh_rxn_ref": (0, pyunits.J/pyunits.mol),
                   "k_eq_ref": (10**-11.4/55.2, pyunits.dimensionless),
                   "T_eq_ref": (298, pyunits.K),

                   # By default, reaction orders follow stoichiometry
                   #    manually set reaction order here to override
                   "reaction_order": {("Liq", "NH3"): 1,
                                    ("Liq", "HOCl"): 1,
                                    ("Liq", "NCl3"): -1,
                                    ("Liq", "H2O"): 0}
                    }
                    # End parameter_data
               },
               # End R1
         }
         # End equilibrium_reactions
    }
    # End reaction_config definition

# Get default solver for testing
solver = get_solver()

# Start test class
class TestChlorination():
    @pytest.fixture(scope="class")
    def chlorination_obj(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(default={"dynamic": False})
        model.fs.thermo_params = GenericParameterBlock(default=thermo_config)
        model.fs.rxn_params = GenericReactionParameterBlock(
                default={"property_package": model.fs.thermo_params, **reaction_config})
        model.fs.unit = EquilibriumReactor(default={
                "property_package": model.fs.thermo_params,
                "reaction_package": model.fs.rxn_params,
                "has_rate_reactions": False,
                "has_equilibrium_reactions": True,
                "has_heat_transfer": False,
                "has_heat_of_reaction": False,
                "has_pressure_change": False})

        model.fs.unit.inlet.mole_frac_comp[0, "H_+"].fix( 0. )
        model.fs.unit.inlet.mole_frac_comp[0, "OH_-"].fix( 0. )
        model.fs.unit.inlet.mole_frac_comp[0, "OCl_-"].fix( 0. )
        model.fs.unit.inlet.mole_frac_comp[0, "NH4_+"].fix( 0. )

        model.fs.unit.inlet.mole_frac_comp[0, "NH2Cl"].fix( 0. )
        model.fs.unit.inlet.mole_frac_comp[0, "NHCl2"].fix( 0. )
        model.fs.unit.inlet.mole_frac_comp[0, "NCl3"].fix( 0. )

        waste_stream_ammonia = 1 #mg/L
        total_molar_density = 54.8 #mol/L
        total_ammonia_inlet = waste_stream_ammonia/17000 # mol/L
        total_molar_density+=total_ammonia_inlet

        # Free Chlorine (mg-Cl2/L) = total_chlorine_inlet (mol/L) * 70,900
        free_chlorine_added = 15 #mg/L as Cl2
        total_chlorine_inlet = free_chlorine_added/70900 # mol/L
        total_molar_density+=total_chlorine_inlet

        model.fs.unit.inlet.mole_frac_comp[0, "NH3"].fix( total_ammonia_inlet/total_molar_density )
        model.fs.unit.inlet.mole_frac_comp[0, "HOCl"].fix( total_chlorine_inlet/total_molar_density )

        # Perform a summation of all non-H2O molefractions to find the H2O molefraction
        sum = 0
        for i in model.fs.unit.inlet.mole_frac_comp:
            # NOTE: i will be a tuple with format (time, component)
            if i[1] != "H2O":
                sum += value(model.fs.unit.inlet.mole_frac_comp[i[0], i[1]])

        model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix( 1-sum )

        model.fs.unit.inlet.pressure.fix(101325.0)
        model.fs.unit.inlet.temperature.fix(300.)
        model.fs.unit.inlet.flow_mol.fix(10)

        return model

    @pytest.mark.unit
    def test_build_model(self, chlorination_obj):
        model = chlorination_obj

        assert hasattr(model.fs.thermo_params, 'component_list')
        assert len(model.fs.thermo_params.component_list) == 10
        assert 'H2O' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.H2O, Solvent)
        assert 'H_+' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('H_+'), Cation)
        assert 'OH_-' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('OH_-'), Anion)

        assert 'OCl_-' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('OCl_-'), Anion)

        assert 'NH4_+' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('NH4_+'), Cation)

        assert 'NH2Cl' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('NH2Cl'), Solute)

        assert 'NHCl2' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('NHCl2'), Solute)

        assert 'NCl3' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.NCl3, Solute)

        assert 'NH3' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.NH3, Solute)

        assert 'HOCl' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.HOCl, Solute)

        assert hasattr(model.fs.thermo_params, 'phase_list')
        assert len(model.fs.thermo_params.phase_list) == 1
        assert isinstance(model.fs.thermo_params.Liq, AqueousPhase)

    @pytest.mark.component
    def test_units(self, chlorination_obj):
        model = chlorination_obj
        assert_units_consistent(model)

    @pytest.mark.unit
    def test_dof(self, chlorination_obj):
        model = chlorination_obj
        assert (degrees_of_freedom(model) == 0)

    @pytest.mark.unit
    def test_stats(self, chlorination_obj):
        model = chlorination_obj
        assert (number_variables(model) == 264)
        assert (number_total_constraints(model) == 85)
        assert (number_unused_variables(model) == 43)

    @pytest.mark.unit
    def test_custom_log_power_law_eps_options(self, chlorination_obj):
        model = chlorination_obj

        assert hasattr(model.fs.thermo_params, 'reaction_H2O_Kw')
        assert hasattr(model.fs.thermo_params.reaction_H2O_Kw, 'eps')

        assert hasattr(model.fs.rxn_params, 'reaction_NCl3_K')
        assert hasattr(model.fs.rxn_params.reaction_NCl3_K, 'eps')

    @pytest.mark.component
    def test_scaling(self, chlorination_obj):
        model = chlorination_obj

        # Iterate through the reactions to set appropriate eps values
        factor = 1e-4
        for rid in model.fs.thermo_params.inherent_reaction_idx:
            scale = value(model.fs.unit.control_volume.properties_out[0.0].k_eq[rid].expr)
            # Want to set eps in some fashion similar to this
            if scale < 1e-16:
                model.fs.thermo_params.component("reaction_"+rid).eps.value = scale*factor
            else:
                model.fs.thermo_params.component("reaction_"+rid).eps.value = 1e-16*factor

        for i in model.fs.unit.control_volume.inherent_reaction_extent_index:
            scale = value(model.fs.unit.control_volume.properties_out[0.0].k_eq[i[1]].expr)
            iscale.set_scaling_factor(model.fs.unit.control_volume.inherent_reaction_extent[0.0,i[1]], 10/scale)
            iscale.constraint_scaling_transform(model.fs.unit.control_volume.properties_out[0.0].
                    inherent_equilibrium_constraint[i[1]], 0.1)

        # Equilibrium reactions have eps in the 'rxn_params'

        factor = 1e-4
        for rid in model.fs.rxn_params.equilibrium_reaction_idx:
            scale = value(model.fs.unit.control_volume.reactions[0.0].k_eq[rid].expr)
            # Want to set eps in some fashion similar to this
            if scale < 1e-16:
                model.fs.rxn_params.component("reaction_"+rid).eps.value = scale*factor
            else:
                model.fs.rxn_params.component("reaction_"+rid).eps.value = 1e-16*factor

        for i in model.fs.unit.control_volume.equilibrium_reaction_extent_index:
            scale = value(model.fs.unit.control_volume.reactions[0.0].k_eq[i[1]].expr)
            iscale.set_scaling_factor(model.fs.unit.control_volume.equilibrium_reaction_extent[0.0,i[1]], 10/scale)
            iscale.constraint_scaling_transform(model.fs.unit.control_volume.reactions[0.0].
                    equilibrium_constraint[i[1]], 0.1)

        # Next, try adding scaling for species
        min = 1e-6
        for i in model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
            # i[0] = phase, i[1] = species
            if model.fs.unit.inlet.mole_frac_comp[0, i[1]].value > min:
                scale = model.fs.unit.inlet.mole_frac_comp[0, i[1]].value
            else:
                scale = min
            iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]], 10/scale)
            iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[i], 10/scale)
            iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i], 10/scale)
            iscale.constraint_scaling_transform(
                model.fs.unit.control_volume.properties_out[0.0].component_flow_balances[i[1]], 10/scale)
            iscale.constraint_scaling_transform(model.fs.unit.control_volume.material_balances[0.0,i[1]], 10/scale)

        iscale.calculate_scaling_factors(model.fs.unit)

        assert isinstance(model.fs.unit.control_volume.scaling_factor, Suffix)

        assert isinstance(model.fs.unit.control_volume.properties_out[0.0].scaling_factor, Suffix)

        assert isinstance(model.fs.unit.control_volume.properties_in[0.0].scaling_factor, Suffix)

        assert isinstance(model.fs.unit.control_volume.reactions[0.0].scaling_factor, Suffix)

    @pytest.mark.component
    def test_initialize(self, chlorination_obj):
        model = chlorination_obj

        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        solver.options['bound_push'] = 1e-10
        solver.options['mu_init'] = 1e-6
        model.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)

        fin_fixed_vars = fixed_variables_set(model)
        fin_act_consts = activated_constraints_set(model)

        assert degrees_of_freedom(model) == 0

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

    @pytest.mark.component
    def test_solve(self, chlorination_obj):
        model = chlorination_obj
        solver.options['max_iter'] = 200
        results = solver.solve(model, tee=True)
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution(self, chlorination_obj):
        model = chlorination_obj

        assert pytest.approx(300.002, rel=1e-5) == value(model.fs.unit.outlet.temperature[0])
        assert pytest.approx(10.0000, rel=1e-5) == value(model.fs.unit.outlet.flow_mol[0])
        assert pytest.approx(101325., rel=1e-5) == value(model.fs.unit.outlet.pressure[0])

        total_molar_density = \
            value(model.fs.unit.control_volume.properties_out[0.0].dens_mol_phase['Liq'])/1000
        assert pytest.approx(55.2044, rel=1e-5) == total_molar_density

        pH = -value(log10(model.fs.unit.outlet.mole_frac_comp[0, "H_+"]*total_molar_density))
        pOH = -value(log10(model.fs.unit.outlet.mole_frac_comp[0, "OH_-"]*total_molar_density))
        assert pytest.approx(6.0522334, rel=1e-4) == pH
        assert pytest.approx(7.94783425, rel=1e-4) == pOH

        hypo_remaining = value(model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq","HOCl"])/1000
        hypo_remaining += value(model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq","OCl_-"])/1000
        combined_chlorine = 0
        combined_chlorine += value(model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq","NH2Cl"])/1000
        combined_chlorine += 2*value(model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq","NHCl2"])/1000
        combined_chlorine += 3*value(model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq","NCl3"])/1000

        hypo_remaining = hypo_remaining*70900
        combined_chlorine = combined_chlorine*70900
        assert pytest.approx(2.50902829, rel=1e-3) == hypo_remaining
        assert pytest.approx(12.607332, rel=1e-3) == combined_chlorine
