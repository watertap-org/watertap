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
    This test is to establish: (i) that IDAES can solve dilute system associated
    with remineralization processes, (ii) that IDAES can correctly assign ions
    to a set of apparent species, (iii) the solution of the system with the IDAES
    assigned apparent species set is the same as the prior solve used to establish
    that set, and (iv) that IDAES can mix-and-match inherent reactions with kinetic
    reactions to emulate more realistic solution chemistry. Solutions are checked
    with approximations to solution pH, alkalinity, and hardness after a typical
    remineralization process.

    NOTE: Scaling for the kinetics is not yet perfected. You will still see AMPL
            errors in the solve for kinetics due to the use of an 'unsafe' power
            law rate function. These errors occur mostly when values for molefractions
            become zero, which cannot be evaluated in a power law function.

    Inherent Reactions:
        H2O <---> H + OH
        H2CO3 <---> H + HCO3
        HCO3 <---> H + CO3
        H2O + CO2 <--> H2CO3

    Additional Apparent species:
        NaHCO3
        Ca(OH)2
        Ca(HCO3)2
        NaOH
        CaCO3

    Kinetic Reactions:
        NaHCO3 --> Na + HCO3
        Ca(OH)2 --> Ca + 2 OH
"""
# Importing testing libraries
import pytest

# Importing pyomo objects
from pyomo.environ import units as pyunits
from pyomo.environ import log10
from pyomo.util.check_units import assert_units_consistent

# Imports from idaes core
from idaes.core import AqueousPhase, VaporPhase, FlowsheetBlock, EnergyBalanceType
from idaes.core.components import Solvent, Solute, Cation, Anion, Apparent
from idaes.core.phases import PhaseType as PT

# Imports from idaes generic models
from idaes.generic_models.properties.core.pure import Perrys, NIST
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ideal import Ideal
from idaes.generic_models.properties.core.pure.ConstantProperties import Constant
from idaes.generic_models.properties.core.generic.generic_property import StateIndex
from idaes.generic_models.properties.core.phase_equil import SmoothVLE
from idaes.generic_models.properties.core.phase_equil.bubble_dew import IdealBubbleDew
from idaes.generic_models.properties.core.phase_equil.forms import fugacity

# Importing the generic model information and objects
from idaes.generic_models.properties.core.generic.generic_reaction import ConcentrationForm
from idaes.generic_models.properties.core.reactions.dh_rxn import constant_dh_rxn
from idaes.generic_models.properties.core.reactions.equilibrium_forms import log_power_law_equil
from idaes.generic_models.properties.core.reactions.rate_forms import power_law_rate
from idaes.generic_models.properties.core.reactions.equilibrium_constant import gibbs_energy, van_t_hoff
from idaes.generic_models.properties.core.reactions.rate_constant import arrhenius
from idaes.generic_models.properties.core.generic.generic_property import GenericParameterBlock
from idaes.generic_models.properties.core.generic.generic_reaction import GenericReactionParameterBlock
from idaes.generic_models.unit_models.equilibrium_reactor import EquilibriumReactor
from idaes.generic_models.unit_models.cstr import CSTR

# Import specific pyomo objects
from pyomo.environ import (ConcreteModel,
                           SolverStatus,
                           TerminationCondition,
                           value,
                           Suffix)

# Import idaes methods to check the model during construction
from idaes.core.util import scaling as iscale
from idaes.core.util.scaling import badly_scaled_var_generator
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom

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
              "enth_mol_ig_comp": NIST,
              "pressure_sat_comp": NIST,
              "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
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
                    },
        'CO2': {"type": Solute,
                "dens_mol_liq_comp": Perrys,
                "enth_mol_liq_comp": Perrys,
                "cp_mol_liq_comp": Perrys,
                "entr_mol_liq_comp": Perrys,
                "enth_mol_ig_comp": NIST,
                "pressure_sat_comp": NIST,
                "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                "parameter_data": {
                    "mw": (44.0095, pyunits.g/pyunits.mol),
                    "pressure_crit": (73.825E5, pyunits.Pa),
                    "temperature_crit": (304.23, pyunits.K),
                    "dens_mol_liq_comp_coeff": {
                        '1': (0.000789, pyunits.kmol * pyunits.m ** -3),
                        '2': (0.000956, pyunits.dimensionless),
                        '3': (500.78, pyunits.K),
                        '4': (0.94599, pyunits.dimensionless)},
                    "cp_mol_ig_comp_coeff": {
                       'A': (24.99735, pyunits.J/pyunits.mol/pyunits.K),
                       'B': (55.18696, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**-1),
                       'C': (-33.69137, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**-2),
                       'D': (7.948387, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**-3),
                       'E': (-0.136638, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**2),
                       'F': (-403.6075, pyunits.kJ/pyunits.mol),
                       'G': (228.2431, pyunits.J/pyunits.mol/pyunits.K),
                       'H': (0, pyunits.kJ/pyunits.mol)},
                    "cp_mol_liq_comp_coeff": {
                        '1': (-8.3043E6, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (1.0437E5, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (4.333E2, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (6.0052E-1, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "enth_mol_form_liq_comp_ref": (-285.83, pyunits.kJ/pyunits.mol),
                    "enth_mol_form_vap_comp_ref": (-393.52, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol),
                    "entr_mol_form_vap_comp_ref": (213.6, pyunits.J/pyunits.mol),
                    "pressure_sat_comp_coeff": {
                        'A': (6.81228, None),
                        'B': (1301.679, pyunits.K),
                        'C': (-3.494, pyunits.K)}
                                }
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
        'Na_+': {"type": Cation, "charge": 1,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (22.989769, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.252, pyunits.kmol*pyunits.m**-3),
                        '2': (0.347, pyunits.dimensionless),
                        '3': (1595.8, pyunits.K),
                        '4': (0.6598, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-240.1, pyunits.J/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (167039, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (59, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'Ca_2+': {"type": Cation, "charge": 2,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (40.078, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (13.5, pyunits.kmol*pyunits.m**-3),
                        '2': (1, pyunits.dimensionless),
                        '3': (1, pyunits.K),
                        '4': (1, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-542.83, pyunits.J/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (2.7637E5, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (-2.0901E3, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (8.125, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (-1.4116E-2, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (9.3701E-6, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (-53, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'H2CO3': {"type": Solute, "valid_phase_types": PT.aqueousPhase,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (62.03, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.4495, pyunits.kmol*pyunits.m**-3),
                        '2': (0.427, pyunits.dimensionless),
                        '3': (429.69, pyunits.K),
                        '4': (0.259, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-699.7, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (135749.9, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (187, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'HCO3_-': {"type": Anion, "charge": -1,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (61.0168, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.4495, pyunits.kmol*pyunits.m**-3),
                        '2': (0.427, pyunits.dimensionless),
                        '3': (429.69, pyunits.K),
                        '4': (0.259, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-692, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (135749.9, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (91.2, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'CO3_2-': {"type": Anion, "charge": -2,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (60.01, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.4495, pyunits.kmol*pyunits.m**-3),
                        '2': (0.427, pyunits.dimensionless),
                        '3': (429.69, pyunits.K),
                        '4': (0.259, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-677.1, pyunits.J/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (135749.9, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (-56.9, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'NaHCO3': {  "type": Apparent,  "valid_phase_types": PT.aqueousPhase,
                    "dissociation_species": {"Na_+":1, "HCO3_-":1},
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
        'Ca(OH)2': {  "type": Apparent,  "valid_phase_types": PT.aqueousPhase,
                    "dissociation_species": {"Ca_2+":1, "OH_-":2},
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
        'NaOH': {  "type": Apparent,  "valid_phase_types": PT.aqueousPhase,
                    "dissociation_species": {"Na_+":1, "OH_-":1},
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
        'CaCO3': {  "type": Apparent,  "valid_phase_types": PT.aqueousPhase,
                    "dissociation_species": {"Ca_2+":1, "CO3_2-":1},
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
        'Ca(HCO3)2': {  "type": Apparent,  "valid_phase_types": PT.aqueousPhase,
                    "dissociation_species": {"Ca_2+":1, "HCO3_-":2},
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
                            "equation_of_state": Ideal}},

        "state_definition": FTPx,
        "state_bounds": {"flow_mol": (0, 50, 100),
                         "temperature": (273.15, 300, 500),
                         "pressure": (5e4, 1e5, 1e6)
                     },
        "state_components": StateIndex.true,

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
            "H2CO3_Ka1": {
                    "stoichiometry": {("Liq", "H2CO3"): -1,
                                     ("Liq", "H_+"): 1,
                                     ("Liq", "HCO3_-"): 1},
                   "heat_of_reaction": constant_dh_rxn,
                   "equilibrium_constant": van_t_hoff,
                   "equilibrium_form": log_power_law_equil,
                   "concentration_form": ConcentrationForm.moleFraction,
                   "parameter_data": {
                       "dh_rxn_ref": (7.7, pyunits.kJ/pyunits.mol),
                       "k_eq_ref": (10**-6.35/55.2, pyunits.dimensionless),
                       "T_eq_ref": (298, pyunits.K),

                       # By default, reaction orders follow stoichiometry
                       #    manually set reaction order here to override
                       "reaction_order": {("Liq", "H2CO3"): -1,
                                        ("Liq", "H_+"): 1,
                                        ("Liq", "HCO3_-"): 1}
                        }
                        # End parameter_data
                   },
                   # End R2
            "H2CO3_Ka2": {
                    "stoichiometry": {("Liq", "HCO3_-"): -1,
                                     ("Liq", "H_+"): 1,
                                     ("Liq", "CO3_2-"): 1},
                   "heat_of_reaction": constant_dh_rxn,
                   "equilibrium_constant": van_t_hoff,
                   "equilibrium_form": log_power_law_equil,
                   "concentration_form": ConcentrationForm.moleFraction,
                   "parameter_data": {
                       "dh_rxn_ref": (14.9, pyunits.kJ/pyunits.mol),
                       "k_eq_ref": (10**-10.33/55.2, pyunits.dimensionless),
                       "T_eq_ref": (298, pyunits.K),

                       # By default, reaction orders follow stoichiometry
                       #    manually set reaction order here to override
                       "reaction_order": {("Liq", "HCO3_-"): -1,
                                        ("Liq", "H_+"): 1,
                                        ("Liq", "CO3_2-"): 1}
                        }
                        # End parameter_data
                   },
                   # End R3
            "CO2_to_H2CO3": {
                    "stoichiometry": {("Liq", "H2O"): -1,
                                     ("Liq", "CO2"): -1,
                                     ("Liq", "H2CO3"): 1},
                   "heat_of_reaction": constant_dh_rxn,
                   "equilibrium_constant": van_t_hoff,
                   "equilibrium_form": log_power_law_equil,
                   "concentration_form": ConcentrationForm.moleFraction,
                   "parameter_data": {
                       "dh_rxn_ref": (0, pyunits.kJ / pyunits.mol),
                       "k_eq_ref": (1.7*10**-3, pyunits.dimensionless),
                       "T_eq_ref": (298, pyunits.K),

                       # By default, reaction orders follow stoichiometry
                       #    manually set reaction order here to override
                       "reaction_order": {("Liq", "H2CO3"): 1,
                                        ("Liq", "CO2"): -1,
                                        ("Liq", "H2O"): 0}
                        }
                        # End parameter_data
                   },
             }
             # End equilibrium_reactions
    }
    # End thermo_config definition

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

# Start test class
class TestRemineralization():
    @pytest.fixture(scope="class")
    def remineralization_appr_equ(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(default={"dynamic": False})
        model.fs.thermo_params = GenericParameterBlock(default=thermo_config)
        model.fs.rxn_params = GenericReactionParameterBlock(
                default={"property_package": model.fs.thermo_params, **reaction_config })
        model.fs.unit = EquilibriumReactor(default={
                "property_package": model.fs.thermo_params,
                "reaction_package": model.fs.rxn_params,
                "has_rate_reactions": False,
                "has_equilibrium_reactions": False,
                "has_heat_transfer": False,
                "has_heat_of_reaction": False,
                "has_pressure_change": False,
                "energy_balance_type": EnergyBalanceType.none
                })

        model.fs.unit.inlet.mole_frac_comp[0, "H_+"].fix( 0. )
        model.fs.unit.inlet.mole_frac_comp[0, "CO3_2-"].fix( 0. )
        model.fs.unit.inlet.mole_frac_comp[0, "H2CO3"].fix( 0. )

        # Add in our species as if they were Ca(OH)2 and NaHCO3
        #       These are typical addatives for remineralization
        Ca_OH_2 = 6e-4  #mol/L
        NaHCO3 = 0.00206  #mol/L
        total = 55.6 + Ca_OH_2 + NaHCO3      #mol/L

        model.fs.unit.inlet.mole_frac_comp[0, "Na_+"].fix( NaHCO3/total )
        model.fs.unit.inlet.mole_frac_comp[0, "Ca_2+"].fix( Ca_OH_2/total )
        model.fs.unit.inlet.mole_frac_comp[0, "OH_-"].fix( 2*Ca_OH_2/total )
        model.fs.unit.inlet.mole_frac_comp[0, "HCO3_-"].fix( NaHCO3/total )

        # Add in a typical CO2 concentration for equilibrium reactor
        model.fs.unit.inlet.mole_frac_comp[0, "CO2"].fix( 0.0005 )

        # Perform a summation of all non-H2O molefractions to find the H2O molefraction
        sum = 0
        for i in model.fs.unit.inlet.mole_frac_comp:
            if i[1] != "H2O":
                sum += value(model.fs.unit.inlet.mole_frac_comp[i[0], i[1]])

        model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix( 1-sum )

        model.fs.unit.inlet.pressure.fix(101325.0)
        model.fs.unit.inlet.temperature.fix(298.)
        model.fs.unit.outlet.temperature.fix(298.)
        model.fs.unit.inlet.flow_mol.fix(10)

        return model

    @pytest.mark.unit
    def test_build_appr_equ(self, remineralization_appr_equ):
        model = remineralization_appr_equ

        assert hasattr(model.fs.thermo_params, 'component_list')
        assert len(model.fs.thermo_params.component_list) == 14
        assert 'H2O' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.H2O, Solvent)
        assert 'H_+' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('H_+'), Cation)
        assert 'OH_-' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('OH_-'), Anion)

        assert 'CO2' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('CO2'), Solute)

        assert 'Na_+' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('Na_+'), Cation)

        assert 'Ca_2+' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('Ca_2+'), Cation)

        assert 'H2CO3' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('H2CO3'), Solute)

        assert 'HCO3_-' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('HCO3_-'), Anion)

        assert 'CO3_2-' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('CO3_2-'), Anion)

        assert 'NaHCO3' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('NaHCO3'), Apparent)

        assert 'Ca(OH)2' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('Ca(OH)2'), Apparent)

        assert 'NaOH' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('NaOH'), Apparent)

        assert 'CaCO3' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('CaCO3'), Apparent)

        assert 'Ca(HCO3)2' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('Ca(HCO3)2'), Apparent)

        assert hasattr(model.fs.thermo_params, 'phase_list')
        assert len(model.fs.thermo_params.phase_list) == 1
        assert isinstance(model.fs.thermo_params.Liq, AqueousPhase)

    @pytest.mark.unit
    def test_units_appr_equ(self, remineralization_appr_equ):
        model = remineralization_appr_equ
        assert_units_consistent(model)

    @pytest.mark.unit
    def test_dof_appr_equ(self, remineralization_appr_equ):
        model = remineralization_appr_equ
        assert (degrees_of_freedom(model) == 0)

    @pytest.mark.component
    def test_scaling_appr_equ(self, remineralization_appr_equ):
        model = remineralization_appr_equ

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

        assert hasattr(model.fs.unit.control_volume, 'scaling_factor')
        assert isinstance(model.fs.unit.control_volume.scaling_factor, Suffix)

        assert hasattr(model.fs.unit.control_volume.properties_out[0.0], 'scaling_factor')
        assert isinstance(model.fs.unit.control_volume.properties_out[0.0].scaling_factor, Suffix)

        assert hasattr(model.fs.unit.control_volume.properties_in[0.0], 'scaling_factor')
        assert isinstance(model.fs.unit.control_volume.properties_in[0.0].scaling_factor, Suffix)

    @pytest.mark.component
    def test_initialize_solver_appr_equ(self, remineralization_appr_equ):
        model = remineralization_appr_equ
        solver.options['bound_push'] = 1e-10
        solver.options['mu_init'] = 1e-6
        model.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_solve_appr_equ(self, remineralization_appr_equ):
        model = remineralization_appr_equ
        solver.options['bound_push'] = 1e-10
        solver.options['mu_init'] = 1e-6
        results = solver.solve(model, tee=True)
        print(results.solver.termination_condition)
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution_appr_equ(self, remineralization_appr_equ):
        model = remineralization_appr_equ

        total_molar_density = value(model.fs.unit.control_volume.properties_out[0.0].dens_mol_phase['Liq'])/1000
        pH = -value(log10(model.fs.unit.outlet.mole_frac_comp[0, "H_+"]*total_molar_density))
        pOH = -value(log10(model.fs.unit.outlet.mole_frac_comp[0, "OH_-"]*total_molar_density))

        assert pytest.approx(8.2018656, rel=1e-4) == pH
        assert pytest.approx(5.7987456, rel=1e-4) == pOH

        # Calculate total hardness
        TH = 2*value(model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp[('Liq', 'Ca_2+')])/1000
        TH = TH*50000
        assert pytest.approx(59.524889, rel=1e-5) == TH

        # Calculating carbonate alkalinity to determine the split of total hardness
        CarbAlk = 2*value(model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp[('Liq', 'CO3_2-')])/1000
        CarbAlk += value(model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp[('Liq', 'HCO3_-')])/1000
        CarbAlk = 50000*CarbAlk
        assert pytest.approx(161.6301239, rel=1e-5) == CarbAlk

        # Non-Carbonate Hardness only exists if there is excess hardness above alkalinity
        if TH <= CarbAlk:
            NCH = 0
        else:
            NCH = TH - CarbAlk
        CH = TH - NCH
        assert pytest.approx(TH, rel=1e-5) == CH

    @pytest.mark.component
    def test_validation_appr_equ(self, remineralization_appr_equ):
        model = remineralization_appr_equ

        # Check the apparent species for valid distribution of ions
        #   i.e., if you use these as input values for apparent species to another
        #           process, then the model will result in same pH, hardness, and
        #           alkalinity, assuming same reaction sets apply
        nahco3 = value(model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp_apparent['Liq','NaHCO3'])
        assert pytest.approx( 3.6490406403736244e-05, rel=1e-3) == nahco3

        h2co3 = value(model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp_apparent['Liq','H2CO3'])
        assert pytest.approx( 8.127381781520279e-07, rel=1e-3) == h2co3

        caoh = value(model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp_apparent['Liq','Ca(OH)2'])
        assert pytest.approx( 5.303703532167982e-09, rel=1e-3) == caoh

        naoh = value(model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp_apparent['Liq','NaOH'])
        assert pytest.approx( 1.8209382127110066e-08, rel=1e-3) == naoh

        caco3 = value(model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp_apparent['Liq','CaCO3'])
        assert pytest.approx( 1.581439142347598e-07, rel=1e-3) == caco3

        cahco3 = value(model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp_apparent['Liq','Ca(HCO3)2'])
        assert pytest.approx( 1.0628273709826089e-05, rel=1e-3) == cahco3


# Configuration dictionary
thermo_config_cstr = {
    "components": {
        'H2O': {"type": Solvent,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              "enth_mol_ig_comp": NIST,
              "pressure_sat_comp": NIST,
              "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
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
                    },
        'CO2': {"type": Solute,
                "dens_mol_liq_comp": Perrys,
                "enth_mol_liq_comp": Perrys,
                "cp_mol_liq_comp": Perrys,
                "entr_mol_liq_comp": Perrys,
                "enth_mol_ig_comp": NIST,
                "pressure_sat_comp": NIST,
                "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                "parameter_data": {
                    "mw": (44.0095, pyunits.g/pyunits.mol),
                    "pressure_crit": (73.825E5, pyunits.Pa),
                    "temperature_crit": (304.23, pyunits.K),
                    "dens_mol_liq_comp_coeff": {
                        '1': (0.000789, pyunits.kmol * pyunits.m ** -3),
                        '2': (0.000956, pyunits.dimensionless),
                        '3': (500.78, pyunits.K),
                        '4': (0.94599, pyunits.dimensionless)},
                    "cp_mol_ig_comp_coeff": {
                       'A': (24.99735, pyunits.J/pyunits.mol/pyunits.K),
                       'B': (55.18696, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**-1),
                       'C': (-33.69137, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**-2),
                       'D': (7.948387, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**-3),
                       'E': (-0.136638, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**2),
                       'F': (-403.6075, pyunits.kJ/pyunits.mol),
                       'G': (228.2431, pyunits.J/pyunits.mol/pyunits.K),
                       'H': (0, pyunits.kJ/pyunits.mol)},
                    "cp_mol_liq_comp_coeff": {
                        '1': (-8.3043E6, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (1.0437E5, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (4.333E2, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (6.0052E-1, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "enth_mol_form_liq_comp_ref": (-285.83, pyunits.kJ/pyunits.mol),
                    "enth_mol_form_vap_comp_ref": (-393.52, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol),
                    "entr_mol_form_vap_comp_ref": (213.6, pyunits.J/pyunits.mol),
                    "pressure_sat_comp_coeff": {
                        'A': (6.81228, None),
                        'B': (1301.679, pyunits.K),
                        'C': (-3.494, pyunits.K)}
                                }
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
        'Na_+': {"type": Cation, "charge": 1,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (22.989769, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.252, pyunits.kmol*pyunits.m**-3),
                        '2': (0.347, pyunits.dimensionless),
                        '3': (1595.8, pyunits.K),
                        '4': (0.6598, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-240.1, pyunits.J/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (167039, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (59, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'Ca_2+': {"type": Cation, "charge": 2,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (40.078, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (13.5, pyunits.kmol*pyunits.m**-3),
                        '2': (1, pyunits.dimensionless),
                        '3': (1, pyunits.K),
                        '4': (1, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-542.83, pyunits.J/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (2.7637E5, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (-2.0901E3, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (8.125, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (-1.4116E-2, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (9.3701E-6, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (-53, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'H2CO3': {"type": Solute, "valid_phase_types": PT.aqueousPhase,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (62.03, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.4495, pyunits.kmol*pyunits.m**-3),
                        '2': (0.427, pyunits.dimensionless),
                        '3': (429.69, pyunits.K),
                        '4': (0.259, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-699.7, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (135749.9, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (187, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'HCO3_-': {"type": Anion, "charge": -1,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (61.0168, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.4495, pyunits.kmol*pyunits.m**-3),
                        '2': (0.427, pyunits.dimensionless),
                        '3': (429.69, pyunits.K),
                        '4': (0.259, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-692, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (135749.9, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (91.2, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'CO3_2-': {"type": Anion, "charge": -2,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (60.01, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.4495, pyunits.kmol*pyunits.m**-3),
                        '2': (0.427, pyunits.dimensionless),
                        '3': (429.69, pyunits.K),
                        '4': (0.259, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-677.1, pyunits.J/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (135749.9, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (-56.9, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'NaHCO3': {"type": Solute, "valid_phase_types": PT.aqueousPhase,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (84.03, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.4495, pyunits.kmol*pyunits.m**-3),
                        '2': (0.427, pyunits.dimensionless),
                        '3': (429.69, pyunits.K),
                        '4': (0.259, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-699.7, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (135749.9, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (187, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'Ca(OH)2': {"type": Solute, "valid_phase_types": PT.aqueousPhase,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (74.03, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.4495, pyunits.kmol*pyunits.m**-3),
                        '2': (0.427, pyunits.dimensionless),
                        '3': (429.69, pyunits.K),
                        '4': (0.259, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-699.7, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (135749.9, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (187, pyunits.J/pyunits.K/pyunits.mol)
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
                         "temperature": (273.15, 300, 500),
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
            "H2CO3_Ka1": {
                    "stoichiometry": {("Liq", "H2CO3"): -1,
                                     ("Liq", "H_+"): 1,
                                     ("Liq", "HCO3_-"): 1},
                   "heat_of_reaction": constant_dh_rxn,
                   "equilibrium_constant": van_t_hoff,
                   "equilibrium_form": log_power_law_equil,
                   "concentration_form": ConcentrationForm.moleFraction,
                   "parameter_data": {
                       "dh_rxn_ref": (7.7, pyunits.kJ/pyunits.mol),
                       "k_eq_ref": (10**-6.35/55.2, pyunits.dimensionless),
                       "T_eq_ref": (298, pyunits.K),

                       # By default, reaction orders follow stoichiometry
                       #    manually set reaction order here to override
                       "reaction_order": {("Liq", "H2CO3"): -1,
                                        ("Liq", "H_+"): 1,
                                        ("Liq", "HCO3_-"): 1}
                        }
                        # End parameter_data
                   },
                   # End R2
            "H2CO3_Ka2": {
                    "stoichiometry": {("Liq", "HCO3_-"): -1,
                                     ("Liq", "H_+"): 1,
                                     ("Liq", "CO3_2-"): 1},
                   "heat_of_reaction": constant_dh_rxn,
                   "equilibrium_constant": van_t_hoff,
                   "equilibrium_form": log_power_law_equil,
                   "concentration_form": ConcentrationForm.moleFraction,
                   "parameter_data": {
                       "dh_rxn_ref": (14.9, pyunits.kJ/pyunits.mol),
                       "k_eq_ref": (10**-10.33/55.2, pyunits.dimensionless),
                       "T_eq_ref": (298, pyunits.K),

                       # By default, reaction orders follow stoichiometry
                       #    manually set reaction order here to override
                       "reaction_order": {("Liq", "HCO3_-"): -1,
                                        ("Liq", "H_+"): 1,
                                        ("Liq", "CO3_2-"): 1}
                        }
                        # End parameter_data
                   },
                   # End R3
            "CO2_to_H2CO3": {
                    "stoichiometry": {("Liq", "H2O"): -1,
                                     ("Liq", "CO2"): -1,
                                     ("Liq", "H2CO3"): 1},
                   "heat_of_reaction": constant_dh_rxn,
                   "equilibrium_constant": van_t_hoff,
                   "equilibrium_form": log_power_law_equil,
                   "concentration_form": ConcentrationForm.moleFraction,
                   "parameter_data": {
                       "dh_rxn_ref": (0, pyunits.kJ / pyunits.mol),
                       "k_eq_ref": (1.7*10**-3, pyunits.dimensionless),
                       "T_eq_ref": (298, pyunits.K),

                       # By default, reaction orders follow stoichiometry
                       #    manually set reaction order here to override
                       "reaction_order": {("Liq", "H2CO3"): 1,
                                        ("Liq", "CO2"): -1,
                                        ("Liq", "H2O"): 0}
                        }
                        # End parameter_data
                   },
             }
             # End equilibrium_reactions
    }
    # End thermo_config definition

reaction_config_cstr = {
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},
    "rate_reactions": {
        "R1": {"stoichiometry": {("Liq", "NaHCO3"): -1,
                                 ("Liq", "HCO3_-"): 1,
                                 ("Liq", "Na_+"): 1},
               "heat_of_reaction": constant_dh_rxn,
               "rate_constant": arrhenius,
               "rate_form": power_law_rate,
               "concentration_form": ConcentrationForm.moleFraction,
               "parameter_data": {
                   "dh_rxn_ref": (0, pyunits.J/pyunits.mol),
                   "arrhenius_const": (1, pyunits.mol/pyunits.m**3/pyunits.s),
                   "energy_activation": (0, pyunits.J/pyunits.mol)}
                   },
        "R2": {"stoichiometry": {("Liq", "Ca(OH)2"): -1,
                                 ("Liq", "OH_-"): 2,
                                 ("Liq", "Ca_2+"): 1},
               "heat_of_reaction": constant_dh_rxn,
               "rate_constant": arrhenius,
               "rate_form": power_law_rate,
               "concentration_form": ConcentrationForm.moleFraction,
               "parameter_data": {
                   "dh_rxn_ref": (0, pyunits.J/pyunits.mol),
                   "arrhenius_const": (1, pyunits.mol/pyunits.m**3/pyunits.s),
                   "energy_activation": (0, pyunits.J/pyunits.mol)}
                   }
         }
         # End rate_reactions
    }
    # End reaction_config definition


# Start test class
class TestRemineralizationCSTR():
    @pytest.fixture(scope="class")
    def remineralization_cstr_kin(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(default={"dynamic": False})
        model.fs.thermo_params = GenericParameterBlock(default=thermo_config_cstr)
        model.fs.rxn_params = GenericReactionParameterBlock(
                default={"property_package": model.fs.thermo_params, **reaction_config_cstr })
        model.fs.unit = CSTR(default={"property_package": model.fs.thermo_params,
                                          "reaction_package": model.fs.rxn_params,
                                          "has_equilibrium_reactions": False,
                                          "has_heat_transfer": False,
                                          "has_heat_of_reaction": False,
                                          "has_pressure_change": False,
                                          "energy_balance_type": EnergyBalanceType.none
                                          })

        model.fs.unit.inlet.pressure.fix(101325.0)
        model.fs.unit.inlet.temperature.fix(298.0)
        model.fs.unit.outlet.temperature.fix(298.0)

        model.fs.unit.volume.fix(100)
        model.fs.unit.inlet.flow_mol.fix(10)

        # Add in our species as if they were Ca(OH)2 and NaHCO3
        #       These are typical addatives for remineralization
        Ca_OH_2 = 6e-4  #mol/L
        NaHCO3 = 0.00206  #mol/L
        total = 55.6 + Ca_OH_2 + NaHCO3      #mol/L
        zero = 1e-20

        model.fs.unit.inlet.mole_frac_comp[0, "NaHCO3"].fix( NaHCO3/total )
        model.fs.unit.inlet.mole_frac_comp[0, "Ca(OH)2"].fix( Ca_OH_2/total )

        model.fs.unit.inlet.mole_frac_comp[0, "H_+"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "CO3_2-"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "H2CO3"].fix( zero )

        model.fs.unit.inlet.mole_frac_comp[0, "Na_+"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "Ca_2+"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "OH_-"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "HCO3_-"].fix( zero )

        model.fs.unit.inlet.mole_frac_comp[0, "CO2"].fix( 0.0005 )

        sum = 0
        for i in model.fs.unit.inlet.mole_frac_comp:
            if i[1] != "H2O":
                sum += value(model.fs.unit.inlet.mole_frac_comp[i[0], i[1]])

        model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix( 1-sum )

        return model

    @pytest.mark.unit
    def test_build_cstr_kin(self, remineralization_cstr_kin):
        model = remineralization_cstr_kin

        assert hasattr(model.fs.thermo_params, 'component_list')
        assert len(model.fs.thermo_params.component_list) == 11
        assert 'H2O' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.H2O, Solvent)
        assert 'H_+' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('H_+'), Cation)
        assert 'OH_-' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('OH_-'), Anion)

        assert 'CO2' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('CO2'), Solute)

        assert 'Na_+' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('Na_+'), Cation)

        assert 'Ca_2+' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('Ca_2+'), Cation)

        assert 'H2CO3' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('H2CO3'), Solute)

        assert 'HCO3_-' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('HCO3_-'), Anion)

        assert 'CO3_2-' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('CO3_2-'), Anion)

        assert 'NaHCO3' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('NaHCO3'), Solute)

        assert 'Ca(OH)2' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('Ca(OH)2'), Solute)

        assert hasattr(model.fs.thermo_params, 'phase_list')
        assert len(model.fs.thermo_params.phase_list) == 1
        assert isinstance(model.fs.thermo_params.Liq, AqueousPhase)

    @pytest.mark.unit
    def test_units_cstr_kin(self, remineralization_cstr_kin):
        model = remineralization_cstr_kin
        assert_units_consistent(model)

    @pytest.mark.unit
    def test_dof_cstr_kin(self, remineralization_cstr_kin):
        model = remineralization_cstr_kin
        assert (degrees_of_freedom(model) == 0)

    @pytest.mark.component
    def test_scaling_cstr_kin(self, remineralization_cstr_kin):
        model = remineralization_cstr_kin

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

        # Scaling for kinetic reactions
        for i in model.fs.rxn_params.rate_reaction_idx:
            scale = value(model.fs.unit.control_volume.reactions[0.0].reaction_rate[i].expr)
            iscale.set_scaling_factor(model.fs.unit.control_volume.rate_reaction_extent[0.0,i], 10/scale)

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

        # CSTR has a volume variable that needs a scaling factor
        #       NOTE: Volume is indexed by time
        iscale.set_scaling_factor(model.fs.unit.control_volume.volume, 10/model.fs.unit.volume[0.0].value)

        iscale.calculate_scaling_factors(model.fs.unit)

        assert hasattr(model.fs.unit.control_volume, 'scaling_factor')
        assert isinstance(model.fs.unit.control_volume.scaling_factor, Suffix)

        assert hasattr(model.fs.unit.control_volume.properties_out[0.0], 'scaling_factor')
        assert isinstance(model.fs.unit.control_volume.properties_out[0.0].scaling_factor, Suffix)

        assert hasattr(model.fs.unit.control_volume.properties_in[0.0], 'scaling_factor')
        assert isinstance(model.fs.unit.control_volume.properties_in[0.0].scaling_factor, Suffix)

    @pytest.mark.component
    def test_initialize_solver_cstr_kin(self, remineralization_cstr_kin):
        model = remineralization_cstr_kin
        solver.options['bound_push'] = 1e-20
        solver.options['mu_init'] = 1e-6
        model.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_solve_cstr_kin(self, remineralization_cstr_kin):
        model = remineralization_cstr_kin
        solver.options['bound_push'] = 1e-20
        solver.options['mu_init'] = 1e-6
        results = solver.solve(model, tee=True, symbolic_solver_labels=True)
        print(results.solver.termination_condition)
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution_cstr_kin(self, remineralization_cstr_kin):
        model = remineralization_cstr_kin

        total_molar_density = value(model.fs.unit.control_volume.properties_out[0.0].dens_mol_phase['Liq'])/1000
        pH = -value(log10(model.fs.unit.outlet.mole_frac_comp[0, "H_+"]*total_molar_density))
        pOH = -value(log10(model.fs.unit.outlet.mole_frac_comp[0, "OH_-"]*total_molar_density))

        assert pytest.approx(8.1593619, rel=1e-4) == pH
        assert pytest.approx(5.8412514, rel=1e-4) == pOH

        # Calculate total hardness
        TH = 2*value(model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp[('Liq', 'Ca_2+')])/1000
        TH = TH*50000
        assert pytest.approx(54.110449, rel=1e-4) == TH

        # Calculating carbonate alkalinity to determine the split of total hardness
        CarbAlk = 2*value(model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp[('Liq', 'CO3_2-')])/1000
        CarbAlk += value(model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp[('Liq', 'HCO3_-')])/1000
        CarbAlk = 50000*CarbAlk
        assert pytest.approx(146.927839, rel=1e-4) == CarbAlk

        # Non-Carbonate Hardness only exists if there is excess hardness above alkalinity
        if TH <= CarbAlk:
            NCH = 0
        else:
            NCH = TH - CarbAlk
        CH = TH - NCH
        assert pytest.approx(TH, rel=1e-5) == CH
