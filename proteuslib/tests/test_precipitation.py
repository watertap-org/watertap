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

    Reactions:
        H2O <---> H + OH
        H2CO3 <---> H + HCO3
        HCO3 <---> H + CO3
        H3PO4 <---> H + H2PO4
        H2PO4 <---> H + HPO4
        HPO4 <---> H + PO4
        FeCl <---> Fe + Cl
        FeOH <---> Fe + OH
        Fe(OH)2 <---> FeOH + OH
        Fe(OH)3 <---> Fe(OH)2 + OH
        Fe(OH)4 <---> Fe(OH)3 + OH
        FeHPO4 <---> Fe + HPO4
        FeH2PO4 <---> Fe + H2PO4
        FePO4 <---> Fe + PO4
    Other species:
        Na
"""
# Importing testing libraries
import pytest

# Importing the object for units from pyomo
from pyomo.environ import units as pyunits
from pyomo.environ import Var
from idaes.core.util.misc import set_param_from_config

# Imports from idaes core
from idaes.core import AqueousPhase, SolidPhase, FlowsheetBlock, EnergyBalanceType
from idaes.core.components import Solvent, Solute, Cation, Anion, Component
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
from idaes.generic_models.properties.core.reactions.equilibrium_forms import log_power_law_equil, power_law_equil

# Import built-in Gibb's Energy function
from idaes.generic_models.properties.core.reactions.equilibrium_constant import gibbs_energy

# Import built-in van't Hoff function
from idaes.generic_models.properties.core.reactions.equilibrium_constant import van_t_hoff

from idaes.generic_models.properties.core.reactions.equilibrium_forms import \
    solubility_product
from idaes.generic_models.properties.core.reactions.equilibrium_constant import \
    ConstantKeq

# Import specific pyomo objects
from pyomo.environ import (ConcreteModel,
                           SolverStatus,
                           TerminationCondition,
                           value,
                           Suffix)

from idaes.core.util import scaling as iscale

# # TODO: Remove this
import idaes.logger as idaeslog

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

__author__ = "Andres Calderon, Austin Ladshaw"

# Heat capacities, enthalpies and entropies
class SolidConstant(object):

    # Ideal liquid properties methods
    class cp_mol_sol_comp(object):

        @staticmethod
        def build_parameters(cobj):
            units = cobj.parent_block().get_metadata().derived_units
            cobj.cp_mol_sol_comp_coeff = Var(
                doc="Parameter for solid phase molar heat capacity",
                units=units["heat_capacity_mole"])
            set_param_from_config(cobj, param="cp_mol_sol_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific heat capacity
            cp = cobj.cp_mol_sol_comp_coeff
            return cp

    class enth_mol_sol_comp(object):

        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_sol_comp_coeff"):
                Constant.cp_mol_sol_comp.build_parameters(cobj)

            if cobj.parent_block().config.include_enthalpy_of_formation:
                units = cobj.parent_block().get_metadata().derived_units

                cobj.enth_mol_form_sol_comp_ref = Var(
                        doc="Solid phase molar heat of formation @ Tref",
                        units=units["energy_mole"])
                set_param_from_config(cobj, param="enth_mol_form_sol_comp_ref")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific enthalpy
            units = b.params.get_metadata().derived_units
            Tr = b.params.temperature_ref

            h_form = (cobj.enth_mol_form_sol_comp_ref if
                    b.params.config.include_enthalpy_of_formation
                    else 0*units["energy_mole"])

            h = cobj.cp_mol_sol_comp_coeff*(T-Tr) + h_form

            return h

    class entr_mol_sol_comp(object):

        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_sol_comp_coeff"):
                Constant.cp_mol_sol_comp.build_parameters(cobj)

            units = cobj.parent_block().get_metadata().derived_units

            cobj.entr_mol_form_sol_comp_ref = Var(
                    doc="Solid phase molar entropy of formation @ Tref",
                    units=units["entropy_mole"])
            set_param_from_config(cobj, param="entr_mol_form_sol_comp_ref")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific entropy
            units = b.params.get_metadata().derived_units
            Tr = b.params.temperature_ref

            s = cobj.cp_mol_sol_comp_coeff*log(T/Tr) + cobj.entr_mol_form_sol_comp_ref

            return s

    class dens_mol_liq_comp(object):

        @staticmethod
        def build_parameters(cobj):
            units = cobj.parent_block().get_metadata().derived_units
            cobj.dens_mol_sol_comp_coeff = Var(
                    doc="Parameter for solid phase molar density",
                    units=units["density_mole"])
            set_param_from_config(cobj, param="dens_mol_sol_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, T):
            # Molar density
            rho = cobj.dens_mol_sol_comp_coeff
            return rho

def dummy_h(b, *args, **kwargs):
    return 0

def dummy_s(b, *args, **kwargs):
    return 0

def dummy_dens(b, *args, **kwargs):
    return 55000 #mol/m**3

def dummy_cp(b, *args, **kwargs):
    return 635 #J/mol/K

# Configuration dictionary
thermo_config = {
    "components": {
        'H2O': {"type": Solvent, "valid_phase_types": PT.aqueousPhase,
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
        # Solid species - Just making this simple change cause major convergence problems
        #'FePO4(s)': {"type": Component, "valid_phase_types": PT.solidPhase,
        #      "dens_mol_sol_comp": dummy_dens,
        #      "enth_mol_sol_comp": dummy_h,
        #      "cp_mol_sol_comp": dummy_cp,
        #      "entr_mol_sol_comp": dummy_s,
        #            },

        'FePO4(s)': {"type": Solute, "valid_phase_types": PT.aqueousPhase,
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              "parameter_data": {
                    "mw": (150.8, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (635000, pyunits.J/pyunits.kmol/pyunits.K),
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    },

              },
              # End Component list
        "phases":  {'Liq': {"type": AqueousPhase,
                            "equation_of_state": Ideal},
                    #'Sol': {"type": SolidPhase,
                    #                    "equation_of_state": Ideal}
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
                   "equilibrium_constant": gibbs_energy,
                   "equilibrium_form": log_power_law_equil,
                   "concentration_form": ConcentrationForm.molarity,
                   "parameter_data": {
                       "dh_rxn_ref": (55.830, pyunits.kJ/pyunits.mol),
                       "ds_rxn_ref": (-80.7, pyunits.J/pyunits.mol/pyunits.K),
                       "T_eq_ref": (300, pyunits.K),

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
                   "equilibrium_constant": gibbs_energy,
                   "equilibrium_form": log_power_law_equil,
                   "concentration_form": ConcentrationForm.molarity,
                   "parameter_data": {
                       "dh_rxn_ref": (7.7, pyunits.kJ/pyunits.mol),
                       "ds_rxn_ref": (-95.8, pyunits.J/pyunits.mol/pyunits.K),
                       "T_eq_ref": (300, pyunits.K),

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
                   "equilibrium_constant": gibbs_energy,
                   "equilibrium_form": log_power_law_equil,
                   "concentration_form": ConcentrationForm.molarity,
                   "parameter_data": {
                       "dh_rxn_ref": (14.9, pyunits.kJ/pyunits.mol),
                       "ds_rxn_ref": (-148.1, pyunits.J/pyunits.mol/pyunits.K),
                       "T_eq_ref": (300, pyunits.K),

                       # By default, reaction orders follow stoichiometry
                       #    manually set reaction order here to override
                       "reaction_order": {("Liq", "HCO3_-"): -1,
                                        ("Liq", "H_+"): 1,
                                        ("Liq", "CO3_2-"): 1}
                        }
                        # End parameter_data
                   },
                   # End R3
            "H3PO4_Ka1": {
                    "stoichiometry": {("Liq", "H3PO4"): -1,
                                     ("Liq", "H_+"): 1,
                                     ("Liq", "H2PO4_-"): 1},
                    "heat_of_reaction": constant_dh_rxn,
                    "equilibrium_constant": gibbs_energy,
                    "equilibrium_form": log_power_law_equil,
                    "concentration_form": ConcentrationForm.molarity,
                    "parameter_data": {
                        "dh_rxn_ref": (-8, pyunits.kJ/pyunits.mol),
                        "ds_rxn_ref": (-67.6, pyunits.J/pyunits.mol/pyunits.K),
                        "T_eq_ref": (300, pyunits.K),

                       # By default, reaction orders follow stoichiometry
                       #    manually set reaction order here to override
                        "reaction_order": { ("Liq", "H3PO4"): -1,
                                            ("Liq", "H_+"): 1,
                                            ("Liq", "H2PO4_-"): 1}
                        }
                        # End parameter_data
                   },
                   # End R4
            "H3PO4_Ka2": {
                        "stoichiometry": {  ("Liq", "H2PO4_-"): -1,
                                            ("Liq", "H_+"): 1,
                                            ("Liq", "HPO4_2-"): 1},
                        "heat_of_reaction": constant_dh_rxn,
                        "equilibrium_constant": gibbs_energy,
                        "equilibrium_form": log_power_law_equil,
                        "concentration_form": ConcentrationForm.molarity,
                        "parameter_data": {
                            "dh_rxn_ref": (4.2, pyunits.kJ/pyunits.mol),
                            "ds_rxn_ref": (-123.8, pyunits.J/pyunits.mol/pyunits.K),
                            "T_eq_ref": (300, pyunits.K),

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
                        "equilibrium_constant": gibbs_energy,
                        "equilibrium_form": log_power_law_equil,
                        "concentration_form": ConcentrationForm.molarity,
                        "parameter_data": {
                            "dh_rxn_ref": (14.7, pyunits.kJ/pyunits.mol),
                            "ds_rxn_ref": (-188.6, pyunits.J/pyunits.mol/pyunits.K),
                            "T_eq_ref": (300, pyunits.K),

                            # By default, reaction orders follow stoichiometry
                            #    manually set reaction order here to override
                            "reaction_order": { ("Liq", "HPO4_2-"): -1,
                                                ("Liq", "H_+"): 1,
                                                ("Liq", "PO4_3-"): 1}
                            }
                            # End parameter_data
                    },
                    # End R6
            "FeCl_K": {
                        "stoichiometry": {  ("Liq", "FeCl_2+"): -1,
                                            ("Liq", "Fe_3+"): 1,
                                            ("Liq", "Cl_-"): 1},
                        "heat_of_reaction": constant_dh_rxn,
                        "equilibrium_constant": van_t_hoff,
                        "equilibrium_form": log_power_law_equil,
                        "concentration_form": ConcentrationForm.molarity,
                        "parameter_data": {
                            "dh_rxn_ref": (0.0, pyunits.J/pyunits.mol),
                            "k_eq_ref": (10**-0.5, pyunits.mol/pyunits.L),
                            "T_eq_ref": (300.0, pyunits.K),

                            # By default, reaction orders follow stoichiometry
                            #    manually set reaction order here to override
                            "reaction_order": { ("Liq", "FeCl_2+"): -1,
                                                ("Liq", "Fe_3+"): 1,
                                                ("Liq", "Cl_-"): 1}
                            }
                            # End parameter_data
                    },
                    # End R7
            "FeOH_K": {
                        "stoichiometry": {  ("Liq", "FeOH_2+"): -1,
                                            ("Liq", "Fe_3+"): 1,
                                            ("Liq", "OH_-"): 1},
                        "heat_of_reaction": constant_dh_rxn,
                        "equilibrium_constant": van_t_hoff,
                        "equilibrium_form": log_power_law_equil,
                        "concentration_form": ConcentrationForm.molarity,
                        "parameter_data": {
                            "dh_rxn_ref": (0.0, pyunits.J/pyunits.mol),
                            "k_eq_ref": (1.768e-12, pyunits.mol/pyunits.L),
                            "T_eq_ref": (300.0, pyunits.K),

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
                        "concentration_form": ConcentrationForm.molarity,
                        "parameter_data": {
                            "dh_rxn_ref": (0.0, pyunits.J/pyunits.mol),
                            "k_eq_ref": (3.757e-11, pyunits.mol/pyunits.L),
                            "T_eq_ref": (300.0, pyunits.K),

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
                        "concentration_form": ConcentrationForm.molarity,
                        "parameter_data": {
                            "dh_rxn_ref": (0.0, pyunits.J/pyunits.mol),
                            "k_eq_ref": (9.765e-7, pyunits.mol/pyunits.L),
                            "T_eq_ref": (300.0, pyunits.K),

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
                        "concentration_form": ConcentrationForm.molarity,
                        "parameter_data": {
                            "dh_rxn_ref": (0.0, pyunits.J/pyunits.mol),
                            "k_eq_ref": (1.097e-6, pyunits.mol/pyunits.L),
                            "T_eq_ref": (300.0, pyunits.K),

                            # By default, reaction orders follow stoichiometry
                            #    manually set reaction order here to override
                            "reaction_order": { ("Liq", "Fe(OH)4_-"): -1,
                                                ("Liq", "Fe(OH)3"): 1,
                                                ("Liq", "OH_-"): 1}
                            }
                            # End parameter_data
                    },
                    # End R11
            "FeHPO4_K": {
                        "stoichiometry": {  ("Liq", "FeHPO4_+"): -1,
                                            ("Liq", "Fe_3+"): 1,
                                            ("Liq", "HPO4_2-"): 1},
                        "heat_of_reaction": constant_dh_rxn,
                        "equilibrium_constant": van_t_hoff,
                        "equilibrium_form": log_power_law_equil,
                        "concentration_form": ConcentrationForm.molarity,
                        "parameter_data": {
                            "dh_rxn_ref": (0.0, pyunits.J/pyunits.mol),
                            "k_eq_ref": (10**-9, pyunits.mol/pyunits.L),
                            "T_eq_ref": (300.0, pyunits.K),

                            # By default, reaction orders follow stoichiometry
                            #    manually set reaction order here to override
                            "reaction_order": { ("Liq", "FeHPO4_+"): -1,
                                                ("Liq", "Fe_3+"): 1,
                                                ("Liq", "HPO4_2-"): 1}
                            }
                            # End parameter_data
                    },
                    # End R12
            "FeH2PO4_K": {
                        "stoichiometry": {  ("Liq", "FeH2PO4_2+"): -1,
                                            ("Liq", "Fe_3+"): 1,
                                            ("Liq", "H2PO4_-"): 1},
                        "heat_of_reaction": constant_dh_rxn,
                        "equilibrium_constant": van_t_hoff,
                        "equilibrium_form": log_power_law_equil,
                        "concentration_form": ConcentrationForm.molarity,
                        "parameter_data": {
                            "dh_rxn_ref": (0.0, pyunits.J/pyunits.mol),
                            "k_eq_ref": (10**-13.4, pyunits.mol/pyunits.L),
                            "T_eq_ref": (300.0, pyunits.K),

                            # By default, reaction orders follow stoichiometry
                            #    manually set reaction order here to override
                            "reaction_order": { ("Liq", "FeH2PO4_2+"): -1,
                                                ("Liq", "Fe_3+"): 1,
                                                ("Liq", "H2PO4_-"): 1}
                            }
                            # End parameter_data
                    },
                    # End R13
            #"FePO4_Ksp": {
            #            "stoichiometry": {  ("Sol", "FePO4(s)"): -1,
            #                                ("Liq", "Fe_3+"): 1,
            #                                ("Liq", "PO4_3-"): 1},
            #            "heat_of_reaction": constant_dh_rxn,
            #            "equilibrium_constant": van_t_hoff,
            #            "equilibrium_form": solubility_product,
            #            "concentration_form": ConcentrationForm.molarity,
            #            "parameter_data": {
            #                "dh_rxn_ref": (0.0, pyunits.J/pyunits.mol),
            #                "k_eq_ref": (10**-23, pyunits.mol**2/pyunits.L**2),
            #                "T_eq_ref": (300.0, pyunits.K),
            #                "reaction_order": { ("Sol", "FePO4(s)"): 0,
            #                                    ("Liq", "Fe_3+"): 1,
            #                                    ("Liq", "PO4_3-"): 1}
            #                }
            #                # End parameter_data
            #        },
            #        # End R14
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
        "FePO4_Ksp": {
                    "stoichiometry": {  ("Sol", "FePO4(s)"): -1,
                                        ("Liq", "Fe_3+"): 1,
                                        ("Liq", "PO4_3-"): 1},
                    "heat_of_reaction": constant_dh_rxn,
                    "equilibrium_constant": van_t_hoff,
                    "equilibrium_form": solubility_product,
                    "concentration_form": ConcentrationForm.molarity,
                    "parameter_data": {
                        "dh_rxn_ref": (0.0, pyunits.J/pyunits.mol),
                        "k_eq_ref": (10**-23, pyunits.mol**2/pyunits.L**2),
                        "T_eq_ref": (300.0, pyunits.K),
                        "reaction_order": { ("Sol", "FePO4(s)"): 0,
                                            ("Liq", "Fe_3+"): 1,
                                            ("Liq", "PO4_3-"): 1}
                        }
                        # End parameter_data
                },
                # End R14
         }
         # End equilibrium_reactions
    }
    # End reaction_config definition

# This config is REQUIRED to use EquilibriumReactor even if we have no equilibrium reactions
reaction_dummy = {
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

# # NOTE: We can now definitively state that there is a fundamental
#           change in how the problem is setup when adding solids

# # TODO: Main
if __name__ == "__main__":
    model = ConcreteModel()
    model.fs = FlowsheetBlock(default={"dynamic": False})
    model.fs.thermo_params = GenericParameterBlock(default=thermo_config)

    model.fs.rxn_params = GenericReactionParameterBlock(
            default={"property_package": model.fs.thermo_params,
                    #**reaction_config
                    **reaction_dummy
                    })

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

    zero = 1e-25
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
    total_nacl_inlet = 0.55 # mol/L (already reduced salt by 4 orders of magnitude)
    total_carbonate_inlet = 0.00206 # mol/L (typical value for seawater = 2.06E-3 M)
    frac_CO3_to_NaHCO3 = 1
    total_phosphate_inlet = 3.22e-6 # mol/L (typical value for seawater = 3.22E-6 M)
    total_phosphate_inlet += 1e-4*0.1 # mol/L (additional phosphorus [e.g., concentrated post-RO])
    #total_phosphate_inlet += 1e-4 # mol/L (additional phosphorus [e.g., concentrated post-RO])
    total_iron_inlet = 5.38e-8 # mol/L (typical value for seawater = 5.38E-8 M)
    #total_iron_inlet += 1e-4*10 # mol/L (additional iron added for phosphorus removal) [Added as FeCl3]
    total_iron_inlet += 1e-4 # mol/L (additional iron added for phosphorus removal) [Added as FeCl3]
    NaOH_added = 1e-4*2 # mol/L (added to raise pH and induce precipitation)

    model.fs.unit.inlet.mole_frac_comp[0, "OH_-"].fix( NaOH_added/total_molar_density )
    model.fs.unit.inlet.mole_frac_comp[0, "Na_+"].fix( (total_nacl_inlet+NaOH_added)/total_molar_density )
    model.fs.unit.inlet.mole_frac_comp[0, "Cl_-"].fix( (total_nacl_inlet+3*total_iron_inlet)/total_molar_density)

    model.fs.unit.inlet.mole_frac_comp[0, "HCO3_-"].fix(
                                (total_carbonate_inlet*frac_CO3_to_NaHCO3)/total_molar_density )
    model.fs.unit.inlet.mole_frac_comp[0, "H2CO3"].fix(
                                (total_carbonate_inlet*(1-frac_CO3_to_NaHCO3))/total_molar_density + zero)

    model.fs.unit.inlet.mole_frac_comp[0, "HPO4_2-"].fix( total_phosphate_inlet/total_molar_density )

    model.fs.unit.inlet.mole_frac_comp[0, "Fe_3+"].fix( total_iron_inlet/total_molar_density )

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


    print("Degrees of freedom = " + str(degrees_of_freedom(model) ) )

    #Custom eps factors for inherent_reactions
    try:
        for rid in model.fs.thermo_params.inherent_reaction_idx:
            scale = value(model.fs.unit.control_volume.properties_out[0.0].k_eq[rid].expr)
            # Want to set eps in some fashion similar to this
            if scale < 1e-16:
                model.fs.thermo_params.component("reaction_"+rid).eps.value = scale*1e-2
            else:
                model.fs.thermo_params.component("reaction_"+rid).eps.value = 1e-16*1e-2

    except:
        pass


    # For equilibrium reactions
    try:
        for rid in model.fs.rxn_params.equilibrium_reaction_idx:
            scale = value(model.fs.unit.control_volume.reactions[0.0].k_eq[rid].expr)
            print(scale)

            # Not all equilibrium reactions will have an eps factor
            try:
                # Want to set eps in some fashion similar to this
                if scale < 1e-16:
                    model.fs.rxn_params.component("reaction_"+rid).eps.value = scale*1e-2
                else:
                    model.fs.rxn_params.component("reaction_"+rid).eps.value = 1e-16*1e-2
            except:
                pass

        #Add scaling factors for reactions
        for i in model.fs.unit.control_volume.equilibrium_reaction_extent_index:
            # i[0] = time, i[1] = reaction
            scale = value(model.fs.unit.control_volume.reactions[0.0].k_eq[i[1]].expr)
            iscale.set_scaling_factor(model.fs.unit.control_volume.equilibrium_reaction_extent[0.0,i[1]], 1)
            iscale.constraint_scaling_transform(model.fs.unit.control_volume.reactions[0.0].
                    equilibrium_constraint[i[1]], 1)
    except:
        pass

    #Add scaling factors for inherent reactions
    for i in model.fs.unit.control_volume.inherent_reaction_extent_index:
        # i[0] = time, i[1] = reaction
        scale = value(model.fs.unit.control_volume.properties_out[0.0].k_eq[i[1]].expr)
        iscale.set_scaling_factor(model.fs.unit.control_volume.inherent_reaction_extent[0.0,i[1]], 10/scale)
        iscale.constraint_scaling_transform(model.fs.unit.control_volume.properties_out[0.0].
                inherent_equilibrium_constraint[i[1]], 0.1)


    # Next, try adding scaling for species
    min = 1e-10
    for i in model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
        # i[0] = phase, i[1] = species
        if model.fs.unit.inlet.mole_frac_comp[0, i[1]].value > min:
            scale = model.fs.unit.inlet.mole_frac_comp[0, i[1]].value
        else:
            scale = min

        #NOTE: Something goes wrong with this scaling for solid species
        if i[0] == 'Liq':
            iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]], 10/scale)
            iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[i], 10/scale)
            iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i], 10/scale)
            iscale.constraint_scaling_transform(
                model.fs.unit.control_volume.properties_out[0.0].component_flow_balances[i[1]], 10/scale)
            iscale.constraint_scaling_transform(model.fs.unit.control_volume.material_balances[0.0,i[1]], 10/scale)



    iscale.calculate_scaling_factors(model.fs.unit)

    #model.fs.unit.control_volume.pprint()
    #exit()


    solver.options['bound_push'] = 1e-20
    solver.options['mu_init'] = 1e-6
    solver.options['max_iter'] = 2000
    model.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)

    solver.options['bound_push'] = 1e-20
    solver.options['mu_init'] = 1e-6
    results = solver.solve(model, tee=True)


    print("comp\toutlet.conc")
    for i in model.fs.unit.inlet.mole_frac_comp:
        print(str(i[1])+"\t"+str(value(model.fs.unit.outlet.mole_frac_comp[i[0], i[1]])*total_molar_density))
    print()

    print("Temperature =\t"+str(value(model.fs.unit.outlet.temperature[0])) )

    print("Flow Mole =\t"+str(value(model.fs.unit.outlet.flow_mol[0])) )

    print("Pressure =\t"+str(value(model.fs.unit.outlet.pressure[0])) )

    carbonate_alk = value(model.fs.unit.outlet.mole_frac_comp[0, "HCO3_-"])*total_molar_density
    carbonate_alk += 2*value(model.fs.unit.outlet.mole_frac_comp[0, "CO3_2-"])*total_molar_density
    carbonate_alk += value(model.fs.unit.outlet.mole_frac_comp[0, "OH_-"])*total_molar_density
    carbonate_alk -= value(model.fs.unit.outlet.mole_frac_comp[0, "H_+"])*total_molar_density
    carbonate_alk = carbonate_alk*50000

    print("Carbonate Alkalinity =\t" + str(carbonate_alk))

    total_carbonate = value(model.fs.unit.outlet.mole_frac_comp[0, "H2CO3"])*total_molar_density
    total_carbonate += value(model.fs.unit.outlet.mole_frac_comp[0, "HCO3_-"])*total_molar_density
    total_carbonate += value(model.fs.unit.outlet.mole_frac_comp[0, "CO3_2-"])*total_molar_density

    print("Total Carbonate =\t" + str(total_carbonate))

    print()

    pH = -value(log10(model.fs.unit.outlet.mole_frac_comp[0, "H_+"]*total_molar_density))
    pOH = -value(log10(model.fs.unit.outlet.mole_frac_comp[0, "OH_-"]*total_molar_density))

    print("pH =\t" + str(pH))
    print("pOH =\t" + str(pOH))

    print()

    total_phosphorus = value(model.fs.unit.outlet.mole_frac_comp[0, "H3PO4"])*total_molar_density
    total_phosphorus += value(model.fs.unit.outlet.mole_frac_comp[0, "H2PO4_-"])*total_molar_density
    total_phosphorus += value(model.fs.unit.outlet.mole_frac_comp[0, "HPO4_2-"])*total_molar_density
    total_phosphorus += value(model.fs.unit.outlet.mole_frac_comp[0, "PO4_3-"])*total_molar_density
    total_phosphorus += value(model.fs.unit.outlet.mole_frac_comp[0, "FeHPO4_+"])*total_molar_density
    total_phosphorus += value(model.fs.unit.outlet.mole_frac_comp[0, "FeH2PO4_2+"])*total_molar_density
    total_phosphorus = total_phosphorus*95000
    total_phosphorus_in = total_phosphate_inlet*95000

    phos_precip = value(model.fs.unit.outlet.mole_frac_comp[0, "FePO4(s)"])*total_molar_density
    phos_precip = phos_precip*95000

    print("Phosphorus in =\t" + str(total_phosphorus_in))
    print("Phosphorus out =\t" + str(total_phosphorus))
    print("Phosphorus precipitated =\t" + str(phos_precip))

    print()

    total_iron = value(model.fs.unit.outlet.mole_frac_comp[0, "Fe_3+"])*total_molar_density
    total_iron += value(model.fs.unit.outlet.mole_frac_comp[0, "FeCl_2+"])*total_molar_density
    total_iron += value(model.fs.unit.outlet.mole_frac_comp[0, "FeOH_2+"])*total_molar_density
    total_iron += value(model.fs.unit.outlet.mole_frac_comp[0, "Fe(OH)2_+"])*total_molar_density
    total_iron += value(model.fs.unit.outlet.mole_frac_comp[0, "Fe(OH)3"])*total_molar_density
    total_iron += value(model.fs.unit.outlet.mole_frac_comp[0, "Fe(OH)4_-"])*total_molar_density
    total_iron += value(model.fs.unit.outlet.mole_frac_comp[0, "FeHPO4_+"])*total_molar_density
    total_iron += value(model.fs.unit.outlet.mole_frac_comp[0, "FeH2PO4_2+"])*total_molar_density
    total_iron = total_iron*55800
    total_iron_in = total_iron_inlet*55800

    iron_precip = value(model.fs.unit.outlet.mole_frac_comp[0, "FePO4(s)"])*total_molar_density
    iron_precip = iron_precip*55800

    print("Iron in =\t" + str(total_iron_in))
    print("Iron out =\t" + str(total_iron))
    print("Iron precipitated =\t" + str(iron_precip))

    print()

    print("Phosphorus Removal Efficiency (%) =\t" + str((total_phosphorus_in-total_phosphorus)/total_phosphorus_in*100))

    try:
        for i in model.fs.unit.control_volume.equilibrium_reaction_extent_index:
            print(value(model.fs.unit.control_volume.reactions[0.0].equilibrium_constraint[i[1]]))
            model.fs.unit.control_volume.reactions[0.0].equilibrium_constraint[i[1]].pprint()
    except:
        pass
