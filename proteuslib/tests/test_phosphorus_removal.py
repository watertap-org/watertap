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
    can solve a pseudo phosphorus removal problem where we assume up front
    that the amount of FeCl3 added is sufficient to reach saturation and
    pseudo-precipitate out phosphorus as FePO4.

    NOTE: We are referring to this as pseudo-precipitation because we are
    not testing the solids phase and solubility product functionality of
    IDAES in this test. Instead, the precipitation is emulated as an
    aqueous equation that is mathematically equivalent to a solubility
    product constraint on the formation of FePO4.

    NOTE: This test does NOT evaluate the precipitation function of IDAES
    We are just trying to establish whether or not IDAES can handle
    these basic problems. Under the condition of saturation, this
    module would be mathematically equivalent to a precipitation
    problem with the same species.

    NOTE: Iron speciation is exceedingly complex and involved. In order to
    try and improve convergence, numerous reactions were left out. These
    reactions are relatively minor under our conditions of interest. The
    removed reactions are noted below.

    Reactions:
        H2O <---> H + OH
        H2CO3 <---> H + HCO3
        HCO3 <---> H + CO3
            #H3PO4 <---> H + H2PO4 (minor reaction, removed)
        H2PO4 <---> H + HPO4
        HPO4 <---> H + PO4
            #FeCl <---> Fe + Cl (minor reaction, removed)
        FeOH <---> Fe + OH (NOTE: This is a minor reaction, but is left in)
        Fe(OH)2 <---> FeOH + OH
        Fe(OH)3 <---> Fe(OH)2 + OH
        Fe(OH)4 <---> Fe(OH)3 + OH
            #FeHPO4 <---> Fe + HPO4 (minor reaction, removed)
            #FeH2PO4 <---> Fe + H2PO4 (minor reaction, removed)
        FePO4 <---> Fe + PO4
    Other species:
        Na
"""
# Importing testing libraries
import pytest

# Importing the object for units from pyomo
from pyomo.environ import units as pyunits

# Imports from idaes core
from idaes.core import AqueousPhase, FlowsheetBlock, EnergyBalanceType
from idaes.core.components import Solvent, Solute, Cation, Anion
from idaes.core.phases import PhaseType as PT

# Imports from idaes generic models
import idaes.generic_models.properties.core.pure.Perrys as Perrys
from idaes.generic_models.properties.core.pure.ConstantProperties import Constant
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

# Import log10 function from pyomo
from pyomo.environ import log10

import idaes.logger as idaeslog

__author__ = "Andres Calderon, Austin Ladshaw"

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
        'Cl_-': {"type": Anion, "charge": -1,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (35.453, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (4.985, pyunits.kmol*pyunits.m**-3),
                        '2': (0.36, pyunits.dimensionless),
                        '3': (1464.06, pyunits.K),
                        '4': (0.739, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-167.2, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (83993.8, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (56.5, pyunits.J/pyunits.K/pyunits.mol)
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
        'H3PO4': {"type": Solute, "valid_phase_types": PT.aqueousPhase,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (97.99518, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.915, pyunits.kmol*pyunits.m**-3),
                        '2': (0.5288, pyunits.dimensionless),
                        '3': (676.2, pyunits.K),
                        '4': (0.31, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-1288.3, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (172459.5, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                     "entr_mol_form_liq_comp_ref": (158, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'H2PO4_-': {"type": Anion, "charge": -1,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (96.98724, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.915, pyunits.kmol*pyunits.m**-3),
                        '2': (0.5288, pyunits.dimensionless),
                        '3': (676.2, pyunits.K),
                        '4': (0.31, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-1296.3, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (172459.5, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                     "entr_mol_form_liq_comp_ref": (90.4, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'HPO4_2-': {"type": Anion, "charge": -2,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (95.9793, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.915, pyunits.kmol*pyunits.m**-3),
                        '2': (0.5288, pyunits.dimensionless),
                        '3': (676.2, pyunits.K),
                        '4': (0.31, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-1292.1, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (172459.5, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                     "entr_mol_form_liq_comp_ref": (-33.4, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'PO4_3-': {"type": Anion, "charge": -3,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (94.97136, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.915, pyunits.kmol*pyunits.m**-3),
                        '2': (0.5288, pyunits.dimensionless),
                        '3': (676.2, pyunits.K),
                        '4': (0.31, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-1277.4, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (172459.5, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (0, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (0, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                     "entr_mol_form_liq_comp_ref": (-222, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'Fe_3+': {"type": Cation, "charge": 3,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (55.845, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (233467, pyunits.J/pyunits.kmol/pyunits.K),
                    "enth_mol_form_liq_comp_ref": (-48.5, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (-315.9, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'FeCl_2+': {"type": Cation, "charge": 2,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (91.3, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (382000, pyunits.J/pyunits.kmol/pyunits.K),
                    # NOTE: these parameters below are unknown
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'FeOH_2+': {"type": Cation, "charge": 2,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (72.8, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (305000, pyunits.J/pyunits.kmol/pyunits.K),
                    # NOTE: these parameters below are not well known
                    "enth_mol_form_liq_comp_ref": (-229.4, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'Fe(OH)2_+': {"type": Cation, "charge": 1,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (89.8, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (375000, pyunits.J/pyunits.kmol/pyunits.K),
                    # NOTE: these parameters below are not well known
                    "enth_mol_form_liq_comp_ref": (-446.7, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'Fe(OH)3': {"type": Solute, "valid_phase_types": PT.aqueousPhase,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (106.8, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (446000, pyunits.J/pyunits.kmol/pyunits.K),
                    # NOTE: these parameters below are not well known
                    "enth_mol_form_liq_comp_ref": (-638.5, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'Fe(OH)4_-': {"type": Anion, "charge": -1,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (123.8, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (518000, pyunits.J/pyunits.kmol/pyunits.K),
                    # NOTE: these parameters below are not well known
                    "enth_mol_form_liq_comp_ref": (-830.0, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'FeHPO4_+': {"type": Cation, "charge": 1,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (151.8, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (635000, pyunits.J/pyunits.kmol/pyunits.K),
                    # NOTE: these parameters below are unknown
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'FeH2PO4_2+': {"type": Cation, "charge": 2,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (152.8, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (635000, pyunits.J/pyunits.kmol/pyunits.K),
                    # NOTE: these parameters below are unknown
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'FePO4(s)': {"type": Solute, "valid_phase_types": PT.aqueousPhase,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (150.8, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (635000, pyunits.J/pyunits.kmol/pyunits.K),
                    # NOTE: these parameters below are unknown
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
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
            "H3PO4_Ka2": {
                        "stoichiometry": {  ("Liq", "H2PO4_-"): -1,
                                            ("Liq", "H_+"): 1,
                                            ("Liq", "HPO4_2-"): 1},
                        "heat_of_reaction": constant_dh_rxn,
                        "equilibrium_constant": van_t_hoff,
                        "equilibrium_form": log_power_law_equil,
                        "concentration_form": ConcentrationForm.moleFraction,
                        "parameter_data": {
                            "dh_rxn_ref": (4.2, pyunits.kJ/pyunits.mol),
                            "k_eq_ref": (10**-5.73/55.2, pyunits.dimensionless),
                            "T_eq_ref": (298, pyunits.K),

                            # By default, reaction orders follow stoichiometry
                            #    manually set reaction order here to override
                            "reaction_order": { ("Liq", "H2PO4_-"): -1,
                                                ("Liq", "H_+"): 1,
                                                ("Liq", "HPO4_2-"): 1}
                            }
                            # End parameter_data
                    },
                    # End R5
            "H3PO4_Ka3": {
                        "stoichiometry": {  ("Liq", "HPO4_2-"): -1,
                                            ("Liq", "H_+"): 1,
                                            ("Liq", "PO4_3-"): 1},
                        "heat_of_reaction": constant_dh_rxn,
                        "equilibrium_constant": van_t_hoff,
                        "equilibrium_form": log_power_law_equil,
                        "concentration_form": ConcentrationForm.moleFraction,
                        "parameter_data": {
                            "dh_rxn_ref": (14.7, pyunits.kJ/pyunits.mol),
                            "k_eq_ref": (10**-7.275/55.2, pyunits.dimensionless),
                            "T_eq_ref": (298, pyunits.K),

                            # By default, reaction orders follow stoichiometry
                            #    manually set reaction order here to override
                            "reaction_order": { ("Liq", "HPO4_2-"): -1,
                                                ("Liq", "H_+"): 1,
                                                ("Liq", "PO4_3-"): 1}
                            }
                            # End parameter_data
                    },
                    # End R6
            "FeOH_K": {
                        "stoichiometry": {  ("Liq", "FeOH_2+"): -1,
                                            ("Liq", "Fe_3+"): 1,
                                            ("Liq", "OH_-"): 1},
                        "heat_of_reaction": constant_dh_rxn,
                        "equilibrium_constant": van_t_hoff,
                        "equilibrium_form": log_power_law_equil,
                        "concentration_form": ConcentrationForm.moleFraction,
                        "parameter_data": {
                            "dh_rxn_ref": (0, pyunits.kJ/pyunits.mol),
                            "k_eq_ref": (1.768e-12/55.2, pyunits.dimensionless),
                            "T_eq_ref": (298, pyunits.K),

                            # By default, reaction orders follow stoichiometry
                            #    manually set reaction order here to override
                            "reaction_order": { ("Liq", "FeOH_2+"): -1,
                                                ("Liq", "Fe_3+"): 1,
                                                ("Liq", "OH_-"): 1}
                            }
                            # End parameter_data
                    },
                    # End R8
            "FeOH2_K": {
                        "stoichiometry": {  ("Liq", "Fe(OH)2_+"): -1,
                                            ("Liq", "FeOH_2+"): 1,
                                            ("Liq", "OH_-"): 1},
                        "heat_of_reaction": constant_dh_rxn,
                        "equilibrium_constant": van_t_hoff,
                        "equilibrium_form": log_power_law_equil,
                        "concentration_form": ConcentrationForm.moleFraction,
                        "parameter_data": {
                            "dh_rxn_ref": (0, pyunits.kJ/pyunits.mol),
                            "k_eq_ref": (3.757e-11/55.2, pyunits.dimensionless),
                            "T_eq_ref": (298, pyunits.K),

                            # By default, reaction orders follow stoichiometry
                            #    manually set reaction order here to override
                            "reaction_order": { ("Liq", "Fe(OH)2_+"): -1,
                                                ("Liq", "FeOH_2+"): 1,
                                                ("Liq", "OH_-"): 1}
                            }
                            # End parameter_data
                    },
                    # End R9
            "FeOH3_K": {
                        "stoichiometry": {  ("Liq", "Fe(OH)3"): -1,
                                            ("Liq", "Fe(OH)2_+"): 1,
                                            ("Liq", "OH_-"): 1},
                        "heat_of_reaction": constant_dh_rxn,
                        "equilibrium_constant": van_t_hoff,
                        "equilibrium_form": log_power_law_equil,
                        "concentration_form": ConcentrationForm.moleFraction,
                        "parameter_data": {
                            "dh_rxn_ref": (0, pyunits.kJ/pyunits.mol),
                            "k_eq_ref": (9.765e-7/55.2, pyunits.dimensionless),
                            "T_eq_ref": (298, pyunits.K),

                            # By default, reaction orders follow stoichiometry
                            #    manually set reaction order here to override
                            "reaction_order": { ("Liq", "Fe(OH)3"): -1,
                                                ("Liq", "Fe(OH)2_+"): 1,
                                                ("Liq", "OH_-"): 1}
                            }
                            # End parameter_data
                    },
                    # End R10
            "FeOH4_K": {
                        "stoichiometry": {  ("Liq", "Fe(OH)4_-"): -1,
                                            ("Liq", "Fe(OH)3"): 1,
                                            ("Liq", "OH_-"): 1},
                        "heat_of_reaction": constant_dh_rxn,
                        "equilibrium_constant": van_t_hoff,
                        "equilibrium_form": log_power_law_equil,
                        "concentration_form": ConcentrationForm.moleFraction,
                        "parameter_data": {
                            "dh_rxn_ref": (0, pyunits.kJ/pyunits.mol),
                            "k_eq_ref": (1.097e-6/55.2, pyunits.dimensionless),
                            "T_eq_ref": (298, pyunits.K),

                            # By default, reaction orders follow stoichiometry
                            #    manually set reaction order here to override
                            "reaction_order": { ("Liq", "Fe(OH)4_-"): -1,
                                                ("Liq", "Fe(OH)3"): 1,
                                                ("Liq", "OH_-"): 1}
                            }
                            # End parameter_data
                    },
                    # End R11
            "FePO4_Ksp": {
                        "stoichiometry": {  ("Liq", "FePO4(s)"): -1,
                                            ("Liq", "Fe_3+"): 1,
                                            ("Liq", "PO4_3-"): 1},
                        "heat_of_reaction": constant_dh_rxn,
                        "equilibrium_constant": van_t_hoff,
                        "equilibrium_form": log_power_law_equil,
                        "concentration_form": ConcentrationForm.moleFraction,
                        "parameter_data": {
                            "dh_rxn_ref": (0, pyunits.kJ/pyunits.mol),
                            "k_eq_ref": (10**-23/55.2/55.2, pyunits.dimensionless),
                            "T_eq_ref": (298, pyunits.K),

                            # By default, reaction orders follow stoichiometry
                            #    manually set reaction order here to override
                            #   NOTE: In a solubility product function, the
                            #       precipitate does not show up in the
                            #       mathematical function.
                            "reaction_order": { ("Liq", "FePO4(s)"): 0,
                                                ("Liq", "Fe_3+"): 1,
                                                ("Liq", "PO4_3-"): 1}
                            }
                            # End parameter_data
                    },
                    # End R14
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

class TestSimplePhosphorusRemoval:
    @pytest.fixture(scope="class")
    def simple_phosphorus_removal(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(default={"dynamic": False})
        model.fs.thermo_params = GenericParameterBlock(default=thermo_config)
        model.fs.rxn_params = GenericReactionParameterBlock(
                default={"property_package": model.fs.thermo_params, **reaction_config})
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

        zero = 1e-20
        model.fs.unit.inlet.mole_frac_comp[0, "H_+"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "CO3_2-"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "PO4_3-"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "H2PO4_-"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "H3PO4"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "FeCl_2+"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "FeOH_2+"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "Fe(OH)2_+"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "Fe(OH)3"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "Fe(OH)4_-"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "FeHPO4_+"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "FeH2PO4_2+"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "FePO4(s)"].fix( zero )

        total_molar_density = 55.2  # mol/L (approximate density of seawater)
        total_nacl_inlet = 0.55 # mol/L (assume seawater salt)
        total_carbonate_inlet = 0.00206 # mol/L (typical value for seawater = 2.06E-3 M)
        frac_CO3_to_NaHCO3 = 0.99
        total_phosphate_inlet = 3.22e-6 # mol/L (typical value for seawater = 3.22E-6 M)
        total_phosphate_inlet += 1e-4 # mol/L (additional phosphorus [what we want to remove])
        total_iron_inlet = 5.38e-8 # mol/L (typical value for seawater = 5.38E-8 M)
        total_iron_inlet += 1e-4 # mol/L (additional iron added for phosphorus removal) [Added as FeCl3]
        NaOH_added = 1e-4 # mol/L (added to raise pH and induce precipitation)

        model.fs.unit.inlet.mole_frac_comp[0, "OH_-"].fix(
                    NaOH_added/total_molar_density )
        model.fs.unit.inlet.mole_frac_comp[0, "Na_+"].fix(
                    (total_nacl_inlet+NaOH_added)/total_molar_density )
        model.fs.unit.inlet.mole_frac_comp[0, "Cl_-"].fix(
                    (total_nacl_inlet+3*total_iron_inlet)/total_molar_density)

        model.fs.unit.inlet.mole_frac_comp[0, "HCO3_-"].fix(
                    (total_carbonate_inlet*frac_CO3_to_NaHCO3)/total_molar_density )
        model.fs.unit.inlet.mole_frac_comp[0, "H2CO3"].fix(
                    (total_carbonate_inlet*(1-frac_CO3_to_NaHCO3))/total_molar_density )

        model.fs.unit.inlet.mole_frac_comp[0, "HPO4_2-"].fix(
                    total_phosphate_inlet/total_molar_density )

        model.fs.unit.inlet.mole_frac_comp[0, "Fe_3+"].fix(
                    total_iron_inlet/total_molar_density )

        # Perform a summation of all non-H2O molefractions to find the H2O molefraction
        sum = 0
        for i in model.fs.unit.inlet.mole_frac_comp:
            # NOTE: i will be a tuple with format (time, component)
            if i[1] != "H2O":
                sum += value(model.fs.unit.inlet.mole_frac_comp[i[0], i[1]])

        model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix( 1-sum )

        model.fs.unit.inlet.pressure.fix(101325.0)
        model.fs.unit.inlet.temperature.fix(298.)
        model.fs.unit.outlet.temperature.fix(298.)
        model.fs.unit.inlet.flow_mol.fix(10)

        return model

    @pytest.mark.unit
    def test_build_simple_phosphorus_removal(self, simple_phosphorus_removal):
        model = simple_phosphorus_removal

        assert hasattr(model.fs.thermo_params, "component_list")
        assert len(model.fs.thermo_params.component_list) == 21
        assert "H2O" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.H2O, Solvent)
        assert "H_+" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("H_+"), Cation)
        assert "OH_-" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("OH_-"), Anion)

        assert "CO3_2-" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("CO3_2-"), Anion)

        assert "HCO3_-" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("HCO3_-"), Anion)

        assert "H2CO3" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("H2CO3"), Solute)

        assert "PO4_3-" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("PO4_3-"), Anion)

        assert "HPO4_2-" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("HPO4_2-"), Anion)

        assert "H2PO4_-" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("H2PO4_-"), Anion)

        assert "H3PO4" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("H3PO4"), Solute)

        assert "FeCl_2+" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("FeCl_2+"), Cation)

        assert "FeOH_2+" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("FeOH_2+"), Cation)

        assert "Fe(OH)2_+" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("Fe(OH)2_+"), Cation)

        assert "Fe(OH)3" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("Fe(OH)3"), Solute)

        assert "Fe(OH)4_-" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("Fe(OH)4_-"), Anion)

        assert "FeHPO4_+" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("FeHPO4_+"), Cation)

        assert "FeH2PO4_2+" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("FeH2PO4_2+"), Cation)

        assert "FePO4(s)" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("FePO4(s)"), Solute)

        assert "Na_+" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("Na_+"), Cation)

        assert "Cl_-" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("Cl_-"), Anion)

        assert "Fe_3+" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("Fe_3+"), Cation)

    @pytest.mark.unit
    def test_units_simple_phosphorus_removal(self, simple_phosphorus_removal):
        model = simple_phosphorus_removal
        assert_units_consistent(model)

    @pytest.mark.unit
    def test_dof_simple_phosphorus_removal(self, simple_phosphorus_removal):
        model = simple_phosphorus_removal
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_scaling_simple_phosphorus_removal(self, simple_phosphorus_removal):
        model = simple_phosphorus_removal

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

        assert isinstance(model.fs.unit.control_volume.scaling_factor, Suffix)
        assert isinstance(model.fs.unit.control_volume.properties_out[0.0].scaling_factor, Suffix)
        assert isinstance(model.fs.unit.control_volume.properties_in[0.0].scaling_factor, Suffix)

    @pytest.mark.component
    def test_initialize_solver(self, simple_phosphorus_removal):
        model = simple_phosphorus_removal
        solver.options["bound_push"] = 1e-20
        solver.options["mu_init"] = 1e-6
        model.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_solve_equilibrium(self, simple_phosphorus_removal):
        model = simple_phosphorus_removal
        solver.options["bound_push"] = 1e-20
        solver.options["mu_init"] = 1e-6
        results = solver.solve(model, tee=True)
        print(results.solver.termination_condition)
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution_simple_phosphorus_removal(self, simple_phosphorus_removal):
        model = simple_phosphorus_removal

        total_molar_density = value(model.fs.unit.control_volume.properties_out[0.0].dens_mol_phase['Liq'])/1000

        carbonate_alk = value(model.fs.unit.outlet.mole_frac_comp[0, "HCO3_-"])*total_molar_density
        carbonate_alk += 2*value(model.fs.unit.outlet.mole_frac_comp[0, "CO3_2-"])*total_molar_density
        carbonate_alk += value(model.fs.unit.outlet.mole_frac_comp[0, "OH_-"])*total_molar_density
        carbonate_alk -= value(model.fs.unit.outlet.mole_frac_comp[0, "H_+"])*total_molar_density
        carbonate_alk = carbonate_alk*50000

        assert pytest.approx(100.7109343, rel=1e-4) == carbonate_alk

        pH = -value(log10(model.fs.unit.outlet.mole_frac_comp[0, "H_+"]*total_molar_density))
        pOH = -value(log10(model.fs.unit.outlet.mole_frac_comp[0, "OH_-"]*total_molar_density))

        assert pytest.approx(8.0910674, rel=1e-4) == pH
        assert pytest.approx(5.9181854, rel=1e-4) == pOH

        total_phosphorus = value(model.fs.unit.outlet.mole_frac_comp[0, "H3PO4"])*total_molar_density
        total_phosphorus += value(model.fs.unit.outlet.mole_frac_comp[0, "H2PO4_-"])*total_molar_density
        total_phosphorus += value(model.fs.unit.outlet.mole_frac_comp[0, "HPO4_2-"])*total_molar_density
        total_phosphorus += value(model.fs.unit.outlet.mole_frac_comp[0, "PO4_3-"])*total_molar_density
        total_phosphorus += value(model.fs.unit.outlet.mole_frac_comp[0, "FeHPO4_+"])*total_molar_density
        total_phosphorus += value(model.fs.unit.outlet.mole_frac_comp[0, "FeH2PO4_2+"])*total_molar_density
        total_phosphorus = total_phosphorus*95000

        phos_precip = value(model.fs.unit.outlet.mole_frac_comp[0, "FePO4(s)"])*total_molar_density
        phos_precip = phos_precip*95000

        assert pytest.approx(0.3235406, rel=1e-4) == total_phosphorus
        assert pytest.approx(9.3784716, rel=1e-4) == phos_precip

        total_iron = value(model.fs.unit.outlet.mole_frac_comp[0, "Fe_3+"])*total_molar_density
        total_iron += value(model.fs.unit.outlet.mole_frac_comp[0, "FeCl_2+"])*total_molar_density
        total_iron += value(model.fs.unit.outlet.mole_frac_comp[0, "FeOH_2+"])*total_molar_density
        total_iron += value(model.fs.unit.outlet.mole_frac_comp[0, "Fe(OH)2_+"])*total_molar_density
        total_iron += value(model.fs.unit.outlet.mole_frac_comp[0, "Fe(OH)3"])*total_molar_density
        total_iron += value(model.fs.unit.outlet.mole_frac_comp[0, "Fe(OH)4_-"])*total_molar_density
        total_iron += value(model.fs.unit.outlet.mole_frac_comp[0, "FeHPO4_+"])*total_molar_density
        total_iron += value(model.fs.unit.outlet.mole_frac_comp[0, "FeH2PO4_2+"])*total_molar_density
        total_iron = total_iron*55800

        iron_precip = value(model.fs.unit.outlet.mole_frac_comp[0, "FePO4(s)"])*total_molar_density
        iron_precip = iron_precip*55800

        assert pytest.approx(0.01523535, rel=1e-4) == total_iron
        assert pytest.approx(5.50861807, rel=1e-4) == iron_precip
