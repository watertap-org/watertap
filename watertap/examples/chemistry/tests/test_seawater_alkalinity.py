###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

"""
    This test is to establish that the core chemistry packages in IDAES
    can solve a more complex aqueous speciation problem meant to mimic
    the salinity, pH, and alkalinity of real seawater. Note that not all
    possible reactions in seawater are being considered for this test,
    only those most relevant to salinity and alkalinity calculations.

    Reactions:
        H2O <---> H + OH
        H2CO3 <---> H + HCO3
        HCO3 <---> H + CO3
        NaHCO3 <---> Na + HCO3
    Other species:
        Cl
"""
# Importing testing libraries
import pytest

# Importing the object for units from pyomo
from pyomo.environ import units as pyunits

# Imports from idaes core
from idaes.core import AqueousPhase
from idaes.core.base.components import Solvent, Solute, Cation, Anion
from idaes.core.base.phases import PhaseType as PT

# Imports from idaes generic models
import idaes.models.properties.modular_properties.pure.Perrys as Perrys
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ideal import Ideal

# Importing the enum for concentration unit basis used in the 'get_concentration_term' function
from idaes.models.properties.modular_properties.base.generic_reaction import (
    ConcentrationForm,
)

# Import the object/function for heat of reaction
from idaes.models.properties.modular_properties.reactions.dh_rxn import constant_dh_rxn

# Import safe log power law equation
from idaes.models.properties.modular_properties.reactions.equilibrium_forms import (
    log_power_law_equil,
)

# Import k-value functions
from idaes.models.properties.modular_properties.reactions.equilibrium_constant import (
    gibbs_energy,
    van_t_hoff,
)

# Import specific pyomo objects
from pyomo.environ import (
    ConcreteModel,
    SolverStatus,
    TerminationCondition,
    value,
    Suffix,
)

from idaes.core.util import scaling as iscale

# Import pyomo methods to check the system units
from pyomo.util.check_units import assert_units_consistent

# Import idaes methods to check the model during construction
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    fixed_variables_set,
    activated_constraints_set,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)

# Import the idaes objects for Generic Properties and Reactions
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
)

# Import the idaes object for the EquilibriumReactor unit model
from idaes.models.unit_models.equilibrium_reactor import EquilibriumReactor

# Import the core idaes objects for Flowsheets and types of balances
from idaes.core import FlowsheetBlock

# Import log10 function from pyomo
from pyomo.environ import log10

__author__ = "Austin Ladshaw"

# Configuration dictionary
thermo_config = {
    "components": {
        "H2O": {
            "type": Solvent,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (18.0153, pyunits.g / pyunits.mol),
                "pressure_crit": (220.64e5, pyunits.Pa),
                "temperature_crit": (647, pyunits.K),
                # Comes from Perry's Handbook:  p. 2-98
                "dens_mol_liq_comp_coeff": {
                    "1": (5.459, pyunits.kmol * pyunits.m**-3),
                    "2": (0.30542, pyunits.dimensionless),
                    "3": (647.13, pyunits.K),
                    "4": (0.081, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-285.830, pyunits.kJ / pyunits.mol),
                "enth_mol_form_vap_comp_ref": (0, pyunits.kJ / pyunits.mol),
                # Comes from Perry's Handbook:  p. 2-174
                "cp_mol_liq_comp_coeff": {
                    "1": (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (8.125, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "cp_mol_ig_comp_coeff": {
                    "A": (30.09200, pyunits.J / pyunits.mol / pyunits.K),
                    "B": (
                        6.832514,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**-1,
                    ),
                    "C": (
                        6.793435,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**-2,
                    ),
                    "D": (
                        -2.534480,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**-3,
                    ),
                    "E": (
                        0.082139,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**2,
                    ),
                    "F": (-250.8810, pyunits.kJ / pyunits.mol),
                    "G": (223.3967, pyunits.J / pyunits.mol / pyunits.K),
                    "H": (0, pyunits.kJ / pyunits.mol),
                },
                "entr_mol_form_liq_comp_ref": (
                    69.95,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
                "pressure_sat_comp_coeff": {
                    "A": (4.6543, None),  # [1], temperature range 255.9 K - 373 K
                    "B": (1435.264, pyunits.K),
                    "C": (-64.848, pyunits.K),
                },
            },
            # End parameter_data
        },
        "H_+": {
            "type": Cation,
            "charge": 1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (1.00784, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": {
                    "1": (5.459, pyunits.kmol * pyunits.m**-3),
                    "2": (0.30542, pyunits.dimensionless),
                    "3": (647.13, pyunits.K),
                    "4": (0.081, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-230.000, pyunits.kJ / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (8.125, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "entr_mol_form_liq_comp_ref": (
                    -10.75,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "OH_-": {
            "type": Anion,
            "charge": -1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (17.008, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": {
                    "1": (5.459, pyunits.kmol * pyunits.m**-3),
                    "2": (0.30542, pyunits.dimensionless),
                    "3": (647.13, pyunits.K),
                    "4": (0.081, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-230.000, pyunits.kJ / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (8.125, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "entr_mol_form_liq_comp_ref": (
                    -10.75,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "Na_+": {
            "type": Cation,
            "charge": 1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (22.989769, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": {
                    "1": (5.252, pyunits.kmol * pyunits.m**-3),
                    "2": (0.347, pyunits.dimensionless),
                    "3": (1595.8, pyunits.K),
                    "4": (0.6598, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-240.1, pyunits.J / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (167039, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (0, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (0, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "entr_mol_form_liq_comp_ref": (59, pyunits.J / pyunits.K / pyunits.mol),
            },
            # End parameter_data
        },
        "Cl_-": {
            "type": Anion,
            "charge": -1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (35.453, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": {
                    "1": (4.985, pyunits.kmol * pyunits.m**-3),
                    "2": (0.36, pyunits.dimensionless),
                    "3": (1464.06, pyunits.K),
                    "4": (0.739, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-167.2, pyunits.kJ / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (83993.8, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (0, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (0, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "entr_mol_form_liq_comp_ref": (
                    56.5,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "H2CO3": {
            "type": Solute,
            "valid_phase_types": PT.aqueousPhase,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (62.03, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": {
                    "1": (5.4495, pyunits.kmol * pyunits.m**-3),
                    "2": (0.427, pyunits.dimensionless),
                    "3": (429.69, pyunits.K),
                    "4": (0.259, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-699.7, pyunits.kJ / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (135749.9, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (0, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (0, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "entr_mol_form_liq_comp_ref": (
                    187,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "HCO3_-": {
            "type": Anion,
            "charge": -1,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (61.0168, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": {
                    "1": (5.4495, pyunits.kmol * pyunits.m**-3),
                    "2": (0.427, pyunits.dimensionless),
                    "3": (429.69, pyunits.K),
                    "4": (0.259, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-692, pyunits.kJ / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (135749.9, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (0, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (0, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "entr_mol_form_liq_comp_ref": (
                    91.2,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "CO3_2-": {
            "type": Anion,
            "charge": -2,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (60.01, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": {
                    "1": (5.4495, pyunits.kmol * pyunits.m**-3),
                    "2": (0.427, pyunits.dimensionless),
                    "3": (429.69, pyunits.K),
                    "4": (0.259, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-677.1, pyunits.J / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (135749.9, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (0, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (0, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "entr_mol_form_liq_comp_ref": (
                    -56.9,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
        "NaHCO3": {
            "type": Solute,
            "valid_phase_types": PT.aqueousPhase,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (84.007, pyunits.g / pyunits.mol),
                "dens_mol_liq_comp_coeff": {
                    "1": (5.252, pyunits.kmol * pyunits.m**-3),
                    "2": (0.347, pyunits.dimensionless),
                    "3": (1595.8, pyunits.K),
                    "4": (0.6598, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-945.53, pyunits.kJ / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "1": (167039, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (0, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (0, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "entr_mol_form_liq_comp_ref": (
                    100,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
            },
            # End parameter_data
        },
    },
    # End Component list
    "phases": {
        "Liq": {"type": AqueousPhase, "equation_of_state": Ideal},
    },
    "state_definition": FTPx,
    "state_bounds": {
        "flow_mol": (0, 50, 100),
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
    },
    # Inherent reactions
    "inherent_reactions": {
        "H2O_Kw": {
            "stoichiometry": {
                ("Liq", "H2O"): -1,
                ("Liq", "H_+"): 1,
                ("Liq", "OH_-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (55.830, pyunits.J / pyunits.mol),
                "k_eq_ref": (10**-14 / 55.2 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "H2O"): 0,
                    ("Liq", "H_+"): 1,
                    ("Liq", "OH_-"): 1,
                },
            }
            # End parameter_data
        },
        # End R1
        "H2CO3_Ka1": {
            "stoichiometry": {
                ("Liq", "H2CO3"): -1,
                ("Liq", "H_+"): 1,
                ("Liq", "HCO3_-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (7.7, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (10**-6.35 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "H2CO3"): -1,
                    ("Liq", "H_+"): 1,
                    ("Liq", "HCO3_-"): 1,
                },
            }
            # End parameter_data
        },
        # End R2
        "H2CO3_Ka2": {
            "stoichiometry": {
                ("Liq", "HCO3_-"): -1,
                ("Liq", "H_+"): 1,
                ("Liq", "CO3_2-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (14.9, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (10**-10.33 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "HCO3_-"): -1,
                    ("Liq", "H_+"): 1,
                    ("Liq", "CO3_2-"): 1,
                },
            }
            # End parameter_data
        },
        # End R3
        "NaHCO3_K1": {
            "stoichiometry": {
                ("Liq", "NaHCO3"): -1,
                ("Liq", "Na_+"): 1,
                ("Liq", "HCO3_-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (13.43, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (10**0.27 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "NaHCO3"): -1,
                    ("Liq", "Na_+"): 1,
                    ("Liq", "HCO3_-"): 1,
                },
            }
            # End parameter_data
        }
        # End R4
    }
    # End equilibrium_reactions
}
# End thermo_config definition

# Define the reaction_config for water dissociation
reaction_config = {
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "equilibrium_reactions": {
        "H2O_Kw": {
            "stoichiometry": {
                ("Liq", "H2O"): -1,
                ("Liq", "H_+"): 1,
                ("Liq", "OH_-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (55.830, pyunits.J / pyunits.mol),
                "k_eq_ref": (10**-14 / 55.2 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "H2O"): 0,
                    ("Liq", "H_+"): 1,
                    ("Liq", "OH_-"): 1,
                },
            }
            # End parameter_data
        },
        # End R1
        "H2CO3_Ka1": {
            "stoichiometry": {
                ("Liq", "H2CO3"): -1,
                ("Liq", "H_+"): 1,
                ("Liq", "HCO3_-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (7.7, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (10**-6.35 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "H2CO3"): -1,
                    ("Liq", "H_+"): 1,
                    ("Liq", "HCO3_-"): 1,
                },
            }
            # End parameter_data
        },
        # End R2
        "H2CO3_Ka2": {
            "stoichiometry": {
                ("Liq", "HCO3_-"): -1,
                ("Liq", "H_+"): 1,
                ("Liq", "CO3_2-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (14.9, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (10**-10.33 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "HCO3_-"): -1,
                    ("Liq", "H_+"): 1,
                    ("Liq", "CO3_2-"): 1,
                },
            }
            # End parameter_data
        },
        # End R3
        "NaHCO3_K1": {
            "stoichiometry": {
                ("Liq", "NaHCO3"): -1,
                ("Liq", "Na_+"): 1,
                ("Liq", "HCO3_-"): 1,
            },
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (13.43, pyunits.kJ / pyunits.mol),
                "k_eq_ref": (10**0.27 / 55.2, pyunits.dimensionless),
                "T_eq_ref": (298, pyunits.K),
                # By default, reaction orders follow stoichiometry
                #    manually set reaction order here to override
                "reaction_order": {
                    ("Liq", "NaHCO3"): -1,
                    ("Liq", "Na_+"): 1,
                    ("Liq", "HCO3_-"): 1,
                },
            }
            # End parameter_data
        }
        # End R4
    }
    # End equilibrium_reactions
}
# End reaction_config definition

# Modify og configs to use either inherent or equilibrium reactions
thermo_only_config = {}
thermo_only_config.update(thermo_config)
del thermo_only_config["inherent_reactions"]

# Get default solver for testing
solver = get_solver()

# Start test class
class TestSeawaterAlkalinity:
    @pytest.fixture(scope="class")
    def inherent_reactions_config(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(dynamic=False)
        model.fs.thermo_params = GenericParameterBlock(**thermo_config)
        model.fs.rxn_params = GenericReactionParameterBlock(
            property_package=model.fs.thermo_params, **reaction_config
        )
        model.fs.unit = EquilibriumReactor(
            property_package=model.fs.thermo_params,
            reaction_package=model.fs.rxn_params,
            has_rate_reactions=False,
            has_equilibrium_reactions=False,
            has_heat_transfer=False,
            has_heat_of_reaction=False,
            has_pressure_change=False,
        )

        model.fs.unit.inlet.mole_frac_comp[0, "H_+"].fix(0.0)
        model.fs.unit.inlet.mole_frac_comp[0, "OH_-"].fix(0.0)
        model.fs.unit.inlet.mole_frac_comp[0, "HCO3_-"].fix(0.0)
        model.fs.unit.inlet.mole_frac_comp[0, "CO3_2-"].fix(0.0)

        total_nacl_inlet = 0.55  # mol/L
        total_carbonate_inlet = 0.00206  # mol/L
        frac_CO3_to_NaHCO3 = 1

        model.fs.unit.inlet.mole_frac_comp[0, "Na_+"].fix(total_nacl_inlet / 54.8)
        model.fs.unit.inlet.mole_frac_comp[0, "Cl_-"].fix(total_nacl_inlet / 54.8)

        model.fs.unit.inlet.mole_frac_comp[0, "NaHCO3"].fix(
            (total_carbonate_inlet * frac_CO3_to_NaHCO3) / 54.8
        )
        model.fs.unit.inlet.mole_frac_comp[0, "H2CO3"].fix(
            (total_carbonate_inlet * (1 - frac_CO3_to_NaHCO3)) / 54.8
        )

        # Perform a summation of all non-H2O molefractions to find the H2O molefraction
        sum = 0
        for i in model.fs.unit.inlet.mole_frac_comp:
            # NOTE: i will be a tuple with format (time, component)
            if i[1] != "H2O":
                sum += value(model.fs.unit.inlet.mole_frac_comp[i[0], i[1]])

        model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix(1 - sum)

        model.fs.unit.inlet.pressure.fix(101325.0)
        model.fs.unit.inlet.temperature.fix(298.0)
        model.fs.unit.inlet.flow_mol.fix(10)

        return model

    @pytest.fixture(scope="class")
    def equilibrium_reactions_config(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(dynamic=False)
        model.fs.thermo_params = GenericParameterBlock(**thermo_only_config)
        model.fs.rxn_params = GenericReactionParameterBlock(
            property_package=model.fs.thermo_params, **reaction_config
        )
        model.fs.unit = EquilibriumReactor(
            property_package=model.fs.thermo_params,
            reaction_package=model.fs.rxn_params,
            has_rate_reactions=False,
            has_equilibrium_reactions=True,
            has_heat_transfer=False,
            has_heat_of_reaction=False,
            has_pressure_change=False,
        )

        model.fs.unit.inlet.mole_frac_comp[0, "H_+"].fix(0.0)
        model.fs.unit.inlet.mole_frac_comp[0, "OH_-"].fix(0.0)
        model.fs.unit.inlet.mole_frac_comp[0, "HCO3_-"].fix(0.0)
        model.fs.unit.inlet.mole_frac_comp[0, "CO3_2-"].fix(0.0)

        total_nacl_inlet = 0.55  # mol/L
        total_carbonate_inlet = 0.00206  # mol/L
        frac_CO3_to_NaHCO3 = 1

        model.fs.unit.inlet.mole_frac_comp[0, "Na_+"].fix(total_nacl_inlet / 54.8)
        model.fs.unit.inlet.mole_frac_comp[0, "Cl_-"].fix(total_nacl_inlet / 54.8)

        model.fs.unit.inlet.mole_frac_comp[0, "NaHCO3"].fix(
            (total_carbonate_inlet * frac_CO3_to_NaHCO3) / 54.8
        )
        model.fs.unit.inlet.mole_frac_comp[0, "H2CO3"].fix(
            (total_carbonate_inlet * (1 - frac_CO3_to_NaHCO3)) / 54.8
        )

        # Perform a summation of all non-H2O molefractions to find the H2O molefraction
        sum = 0
        for i in model.fs.unit.inlet.mole_frac_comp:
            # NOTE: i will be a tuple with format (time, component)
            if i[1] != "H2O":
                sum += value(model.fs.unit.inlet.mole_frac_comp[i[0], i[1]])

        model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix(1 - sum)

        model.fs.unit.inlet.pressure.fix(101325.0)
        model.fs.unit.inlet.temperature.fix(298.0)
        model.fs.unit.inlet.flow_mol.fix(10)

        return model

    @pytest.mark.unit
    def test_build_model_inherent(self, inherent_reactions_config):
        model = inherent_reactions_config

        assert hasattr(model.fs.thermo_params, "component_list")
        assert len(model.fs.thermo_params.component_list) == 9
        assert "H2O" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.H2O, Solvent)
        assert "H_+" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("H_+"), Cation)
        assert "OH_-" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("OH_-"), Anion)
        assert "Na_+" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("Na_+"), Cation)
        assert "Cl_-" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("Cl_-"), Anion)
        assert "HCO3_-" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("HCO3_-"), Anion)
        assert "CO3_2-" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("CO3_2-"), Anion)
        assert "H2CO3" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.H2CO3, Solute)
        assert "NaHCO3" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.NaHCO3, Solute)

        assert hasattr(model.fs.thermo_params, "phase_list")
        assert len(model.fs.thermo_params.phase_list) == 1
        assert isinstance(model.fs.thermo_params.Liq, AqueousPhase)

    @pytest.mark.unit
    def test_build_model_equilibrium(self, equilibrium_reactions_config):
        model = equilibrium_reactions_config

        assert hasattr(model.fs.thermo_params, "component_list")
        assert len(model.fs.thermo_params.component_list) == 9
        assert "H2O" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.H2O, Solvent)
        assert "H_+" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("H_+"), Cation)
        assert "OH_-" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("OH_-"), Anion)
        assert "Na_+" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("Na_+"), Cation)
        assert "Cl_-" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("Cl_-"), Anion)
        assert "HCO3_-" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("HCO3_-"), Anion)
        assert "CO3_2-" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("CO3_2-"), Anion)
        assert "H2CO3" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.H2CO3, Solute)
        assert "NaHCO3" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.NaHCO3, Solute)

        assert hasattr(model.fs.thermo_params, "phase_list")
        assert len(model.fs.thermo_params.phase_list) == 1
        assert isinstance(model.fs.thermo_params.Liq, AqueousPhase)

    @pytest.mark.unit
    def test_units_inherent(self, inherent_reactions_config):
        model = inherent_reactions_config
        assert_units_consistent(model)

    @pytest.mark.unit
    def test_units_equilibrium(self, equilibrium_reactions_config):
        model = equilibrium_reactions_config
        assert_units_consistent(model)

    @pytest.mark.unit
    def test_dof_inherent(self, inherent_reactions_config):
        model = inherent_reactions_config
        assert degrees_of_freedom(model) == 0

    @pytest.mark.unit
    def test_dof_equilibrium(self, equilibrium_reactions_config):
        model = equilibrium_reactions_config
        assert degrees_of_freedom(model) == 0

    @pytest.mark.unit
    def test_stats_inherent(self, inherent_reactions_config):
        model = inherent_reactions_config
        assert number_variables(model) == 281
        assert number_total_constraints(model) == 72

    @pytest.mark.unit
    def test_stats_equilibrium(self, equilibrium_reactions_config):
        model = equilibrium_reactions_config
        assert number_variables(model) == 233
        assert number_total_constraints(model) == 72

    @pytest.mark.component
    def test_scaling_inherent(self, inherent_reactions_config):
        model = inherent_reactions_config

        for i in model.fs.unit.control_volume.inherent_reaction_extent_index:
            scale = value(
                model.fs.unit.control_volume.properties_out[0.0].k_eq[i[1]].expr
            )
            iscale.set_scaling_factor(
                model.fs.unit.control_volume.inherent_reaction_extent[0.0, i[1]],
                10 / scale,
            )
            iscale.constraint_scaling_transform(
                model.fs.unit.control_volume.properties_out[
                    0.0
                ].inherent_equilibrium_constraint[i[1]],
                0.1,
            )

        # Next, try adding scaling for species
        min = 1e-3
        for i in model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
            # i[0] = phase, i[1] = species
            if model.fs.unit.inlet.mole_frac_comp[0, i[1]].value > min:
                scale = model.fs.unit.inlet.mole_frac_comp[0, i[1]].value
            else:
                scale = min
            iscale.set_scaling_factor(
                model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]],
                10 / scale,
            )
            iscale.set_scaling_factor(
                model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
                    i
                ],
                10 / scale,
            )
            iscale.set_scaling_factor(
                model.fs.unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i],
                10 / scale,
            )
            iscale.constraint_scaling_transform(
                model.fs.unit.control_volume.properties_out[
                    0.0
                ].component_flow_balances[i[1]],
                10 / scale,
            )
            iscale.constraint_scaling_transform(
                model.fs.unit.control_volume.material_balances[0.0, i[1]], 10 / scale
            )

        iscale.calculate_scaling_factors(model.fs.unit)

        assert isinstance(model.fs.unit.control_volume.scaling_factor, Suffix)

        assert isinstance(
            model.fs.unit.control_volume.properties_out[0.0].scaling_factor, Suffix
        )

        assert isinstance(
            model.fs.unit.control_volume.properties_in[0.0].scaling_factor, Suffix
        )

    @pytest.mark.component
    def test_scaling_equilibrium(self, equilibrium_reactions_config):
        model = equilibrium_reactions_config

        for i in model.fs.unit.control_volume.equilibrium_reaction_extent_index:
            scale = value(model.fs.unit.control_volume.reactions[0.0].k_eq[i[1]].expr)
            iscale.set_scaling_factor(
                model.fs.unit.control_volume.equilibrium_reaction_extent[0.0, i[1]],
                10 / scale,
            )
            iscale.constraint_scaling_transform(
                model.fs.unit.control_volume.reactions[0.0].equilibrium_constraint[
                    i[1]
                ],
                0.1,
            )

        # Next, try adding scaling for species
        min = 1e-3
        for i in model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
            # i[0] = phase, i[1] = species
            if model.fs.unit.inlet.mole_frac_comp[0, i[1]].value > min:
                scale = model.fs.unit.inlet.mole_frac_comp[0, i[1]].value
            else:
                scale = min
            iscale.set_scaling_factor(
                model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]],
                10 / scale,
            )
            iscale.set_scaling_factor(
                model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
                    i
                ],
                10 / scale,
            )
            iscale.set_scaling_factor(
                model.fs.unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i],
                10 / scale,
            )
            iscale.constraint_scaling_transform(
                model.fs.unit.control_volume.properties_out[
                    0.0
                ].component_flow_balances[i[1]],
                10 / scale,
            )
            iscale.constraint_scaling_transform(
                model.fs.unit.control_volume.material_balances[0.0, i[1]], 10 / scale
            )

        iscale.calculate_scaling_factors(model.fs.unit)

        assert isinstance(model.fs.unit.control_volume.scaling_factor, Suffix)

        assert isinstance(
            model.fs.unit.control_volume.properties_out[0.0].scaling_factor, Suffix
        )

        assert isinstance(
            model.fs.unit.control_volume.properties_in[0.0].scaling_factor, Suffix
        )

        # When using equilibrium reactions, there are another set of scaling factors calculated
        assert isinstance(
            model.fs.unit.control_volume.reactions[0.0].scaling_factor, Suffix
        )

    @pytest.mark.component
    def test_initialize_inherent(self, inherent_reactions_config):
        model = inherent_reactions_config

        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.fs.unit.initialize(optarg=solver.options)

        fin_fixed_vars = fixed_variables_set(model)
        fin_act_consts = activated_constraints_set(model)

        assert degrees_of_freedom(model) == 0

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

    @pytest.mark.component
    def test_initialize_equilibrium(self, equilibrium_reactions_config):
        model = equilibrium_reactions_config

        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.fs.unit.initialize(optarg=solver.options)

        fin_fixed_vars = fixed_variables_set(model)
        fin_act_consts = activated_constraints_set(model)

        assert degrees_of_freedom(model) == 0

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

    @pytest.mark.component
    def test_solve_inherent(self, inherent_reactions_config):
        model = inherent_reactions_config
        solver.options["max_iter"] = 20
        results = solver.solve(model)
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solve_equilibrium(self, equilibrium_reactions_config):
        model = equilibrium_reactions_config
        solver.options["max_iter"] = 20
        results = solver.solve(model)
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution_inherent(self, inherent_reactions_config):
        model = inherent_reactions_config

        assert pytest.approx(297.9, rel=1e-5) == value(
            model.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(10.00029, rel=1e-5) == value(
            model.fs.unit.outlet.flow_mol[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            model.fs.unit.outlet.pressure[0]
        )

        total_molar_density = (
            value(
                model.fs.unit.control_volume.properties_out[0.0].dens_mol_phase["Liq"]
            )
            / 1000
        )
        assert pytest.approx(54.612040964248116, rel=1e-5) == total_molar_density

        total_salt = (
            value(model.fs.unit.outlet.mole_frac_comp[0, "Na_+"])
            * total_molar_density
            * 23
        )
        total_salt += (
            value(model.fs.unit.outlet.mole_frac_comp[0, "Cl_-"])
            * total_molar_density
            * 35.44
        )
        total_salt += (
            value(model.fs.unit.outlet.mole_frac_comp[0, "NaHCO3"])
            * total_molar_density
            * 84
        )
        psu = total_salt / (total_molar_density * 18) * 1000
        assert pytest.approx(32.66105, rel=1e-5) == psu

        total_carbonate = (
            value(model.fs.unit.outlet.mole_frac_comp[0, "NaHCO3"])
            * total_molar_density
        )
        total_carbonate += (
            value(model.fs.unit.outlet.mole_frac_comp[0, "H2CO3"]) * total_molar_density
        )
        total_carbonate += (
            value(model.fs.unit.outlet.mole_frac_comp[0, "HCO3_-"])
            * total_molar_density
        )
        total_carbonate += (
            value(model.fs.unit.outlet.mole_frac_comp[0, "CO3_2-"])
            * total_molar_density
        )
        assert pytest.approx(0.0020528746175810923, rel=1e-5) == total_carbonate

        carbonate_alk = (
            value(model.fs.unit.outlet.mole_frac_comp[0, "HCO3_-"])
            * total_molar_density
        )
        carbonate_alk += (
            2
            * value(model.fs.unit.outlet.mole_frac_comp[0, "CO3_2-"])
            * total_molar_density
        )
        carbonate_alk += (
            value(model.fs.unit.outlet.mole_frac_comp[0, "OH_-"]) * total_molar_density
        )
        carbonate_alk -= (
            value(model.fs.unit.outlet.mole_frac_comp[0, "H_+"]) * total_molar_density
        )
        carbonate_alk = carbonate_alk * 50000
        assert pytest.approx(79.389737, rel=1e-4) == carbonate_alk

        pH = -value(
            log10(model.fs.unit.outlet.mole_frac_comp[0, "H_+"] * total_molar_density)
        )
        pOH = -value(
            log10(model.fs.unit.outlet.mole_frac_comp[0, "OH_-"] * total_molar_density)
        )
        assert pytest.approx(8.31763834, rel=1e-4) == pH
        assert pytest.approx(5.6916662, rel=1e-4) == pOH

    @pytest.mark.component
    def test_solution_equilibrium(self, equilibrium_reactions_config):
        model = equilibrium_reactions_config

        assert pytest.approx(297.9, rel=1e-5) == value(
            model.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(10.0002, rel=1e-5) == value(
            model.fs.unit.outlet.flow_mol[0]
        )
        assert pytest.approx(101325, rel=1e-5) == value(
            model.fs.unit.outlet.pressure[0]
        )

        total_molar_density = (
            value(
                model.fs.unit.control_volume.properties_out[0.0].dens_mol_phase["Liq"]
            )
            / 1000
        )
        assert pytest.approx(54.612040964248116, rel=1e-5) == total_molar_density

        total_salt = (
            value(model.fs.unit.outlet.mole_frac_comp[0, "Na_+"])
            * total_molar_density
            * 23
        )
        total_salt += (
            value(model.fs.unit.outlet.mole_frac_comp[0, "Cl_-"])
            * total_molar_density
            * 35.44
        )
        total_salt += (
            value(model.fs.unit.outlet.mole_frac_comp[0, "NaHCO3"])
            * total_molar_density
            * 84
        )
        psu = total_salt / (total_molar_density * 18) * 1000
        assert pytest.approx(32.66105, rel=1e-5) == psu

        total_carbonate = (
            value(model.fs.unit.outlet.mole_frac_comp[0, "NaHCO3"])
            * total_molar_density
        )
        total_carbonate += (
            value(model.fs.unit.outlet.mole_frac_comp[0, "H2CO3"]) * total_molar_density
        )
        total_carbonate += (
            value(model.fs.unit.outlet.mole_frac_comp[0, "HCO3_-"])
            * total_molar_density
        )
        total_carbonate += (
            value(model.fs.unit.outlet.mole_frac_comp[0, "CO3_2-"])
            * total_molar_density
        )
        assert pytest.approx(0.0020528746175810923, rel=1e-5) == total_carbonate

        carbonate_alk = (
            value(model.fs.unit.outlet.mole_frac_comp[0, "HCO3_-"])
            * total_molar_density
        )
        carbonate_alk += (
            2
            * value(model.fs.unit.outlet.mole_frac_comp[0, "CO3_2-"])
            * total_molar_density
        )
        carbonate_alk += (
            value(model.fs.unit.outlet.mole_frac_comp[0, "OH_-"]) * total_molar_density
        )
        carbonate_alk -= (
            value(model.fs.unit.outlet.mole_frac_comp[0, "H_+"]) * total_molar_density
        )
        carbonate_alk = carbonate_alk * 50000
        assert pytest.approx(79.389737, rel=1e-4) == carbonate_alk

        pH = -value(
            log10(model.fs.unit.outlet.mole_frac_comp[0, "H_+"] * total_molar_density)
        )
        pOH = -value(
            log10(model.fs.unit.outlet.mole_frac_comp[0, "OH_-"] * total_molar_density)
        )
        assert pytest.approx(8.31763834, rel=1e-4) == pH
        assert pytest.approx(5.6916662, rel=1e-4) == pOH
