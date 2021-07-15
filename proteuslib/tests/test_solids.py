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
    This test is of the core IDAES components that allow for the declaration
    and usage of solids phases in conjunction with aqueous phases. This test
    is primarily being used to probe for issues that may exist in how IDAES
    handles this new system and properties, since it is brand new to the
    framework.

    Case 1: Declare A and B as aqueous, no reactions (record IDAES debug)

    Case 2a: Declare A and B as solids, no reactions (record IDAES debug) [with bulk H2O]

    Case 2b: Declare A and B as solids, no reactions (record IDAES debug) [without any liq solvent]

    Case 3: [Combine] A and B are aqueous, declare AB a solid, add precipitation reaction (record IDAES debug)

    Solubility Reaction:
        AB <---> A + B


    NOTE: All cases should be tested with and without scaling, and with
            a multitude of inlet concentration ranges.

    NOTE: All cases will involve H2O as the bulk solvent (keeps our problem relevant)
            Exception: Case2b
"""

# Importing testing libraries
import pytest

# Importing the object for units from pyomo
from pyomo.environ import units as pyunits
from pyomo.environ import Var

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

__author__ = "Austin Ladshaw"

# Case 1 Config
case1_thermo_config = {
    "components": {
        'H2O': {"type": Solvent, "valid_phase_types": PT.aqueousPhase,
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (18.0153, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (75.312, pyunits.J/pyunits.mol/pyunits.K),
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'A': {"type": Solute, "valid_phase_types": PT.aqueousPhase,
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              "parameter_data": {
                    "mw": (18, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (75.312, pyunits.J/pyunits.mol/pyunits.K),
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    },
        'B': {"type": Solute, "valid_phase_types": PT.aqueousPhase,
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              "parameter_data": {
                    "mw": (18, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (75.312, pyunits.J/pyunits.mol/pyunits.K),
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    },

              },
              # End Component list
        "phases":  {'Liq': {"type": AqueousPhase,
                            "equation_of_state": Ideal}
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
    # End thermo_config definition

# Case 2a Config
case2a_thermo_config = {
    "components": {
        'H2O': {"type": Solvent, "valid_phase_types": PT.aqueousPhase,
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (18.0153, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (75.312, pyunits.J/pyunits.mol/pyunits.K),
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'A': {"type": Component, "valid_phase_types": PT.solidPhase,
              "dens_mol_sol_comp": Constant,
              "enth_mol_sol_comp": Constant,
              "cp_mol_sol_comp": Constant,
              "entr_mol_sol_comp": Constant,
              "parameter_data": {
                    "mw": (18, pyunits.g/pyunits.mol),
                    "dens_mol_sol_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_sol_comp_coeff": (75.312, pyunits.J/pyunits.mol/pyunits.K),
                    "enth_mol_form_sol_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_sol_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    },
        'B': {"type": Component, "valid_phase_types": PT.solidPhase,
              "dens_mol_sol_comp": Constant,
              "enth_mol_sol_comp": Constant,
              "cp_mol_sol_comp": Constant,
              "entr_mol_sol_comp": Constant,
              "parameter_data": {
                    "mw": (18, pyunits.g/pyunits.mol),
                    "dens_mol_sol_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_sol_comp_coeff": (75.312, pyunits.J/pyunits.mol/pyunits.K),
                    "enth_mol_form_sol_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_sol_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    },

              },
              # End Component list
        "phases":  {'Liq': {"type": AqueousPhase,
                            "equation_of_state": Ideal},
                    'Sol': {"type": SolidPhase,
                                        "equation_of_state": Ideal}
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
    # End thermo_config definition


# Case 2b Config
case2b_thermo_config = {
    "components": {
        'A': {"type": Component, "valid_phase_types": PT.solidPhase,
              "dens_mol_sol_comp": Constant,
              "enth_mol_sol_comp": Constant,
              "cp_mol_sol_comp": Constant,
              "entr_mol_sol_comp": Constant,
              "parameter_data": {
                    "mw": (18, pyunits.g/pyunits.mol),
                    "dens_mol_sol_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_sol_comp_coeff": (75.312, pyunits.J/pyunits.mol/pyunits.K),
                    "enth_mol_form_sol_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_sol_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    },
        'B': {"type": Component, "valid_phase_types": PT.solidPhase,
              "dens_mol_sol_comp": Constant,
              "enth_mol_sol_comp": Constant,
              "cp_mol_sol_comp": Constant,
              "entr_mol_sol_comp": Constant,
              "parameter_data": {
                    "mw": (18, pyunits.g/pyunits.mol),
                    "dens_mol_sol_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_sol_comp_coeff": (75.312, pyunits.J/pyunits.mol/pyunits.K),
                    "enth_mol_form_sol_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_sol_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    },

              },
              # End Component list
        "phases":  {
                    'Sol': {"type": SolidPhase,
                                        "equation_of_state": Ideal}
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
    # End thermo_config definition


# Case 3 Config
case3_thermo_config = {
    "components": {
        'H2O': {"type": Solvent, "valid_phase_types": PT.aqueousPhase,
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (18.0153, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (75.312, pyunits.J/pyunits.mol/pyunits.K),
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'A': {"type": Solute, "valid_phase_types": PT.aqueousPhase,
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              "parameter_data": {
                    "mw": (18, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (75.312, pyunits.J/pyunits.mol/pyunits.K),
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    },
        'B': {"type": Solute, "valid_phase_types": PT.aqueousPhase,
              "dens_mol_liq_comp": Constant,
              "enth_mol_liq_comp": Constant,
              "cp_mol_liq_comp": Constant,
              "entr_mol_liq_comp": Constant,
              "parameter_data": {
                    "mw": (18, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_liq_comp_coeff": (75.312, pyunits.J/pyunits.mol/pyunits.K),
                    "enth_mol_form_liq_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    },
        'AB': {"type": Component, "valid_phase_types": PT.solidPhase,
              "dens_mol_sol_comp": Constant,
              "enth_mol_sol_comp": Constant,
              "cp_mol_sol_comp": Constant,
              "entr_mol_sol_comp": Constant,
              "parameter_data": {
                    "mw": (18, pyunits.g/pyunits.mol),
                    "dens_mol_sol_comp_coeff": (55.2, pyunits.kmol*pyunits.m**-3),
                    "cp_mol_sol_comp_coeff": (75.312, pyunits.J/pyunits.mol/pyunits.K),
                    "enth_mol_form_sol_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    "entr_mol_form_sol_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    },

              },
              # End Component list
        "phases":  {'Liq': {"type": AqueousPhase,
                            "equation_of_state": Ideal},
                    'Sol': {"type": SolidPhase,
                            "equation_of_state": Ideal}
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
    # End thermo_config definition

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

def run_case1(xA, xB, scaling=True):
    print("==========================================================================")
    print("Case 1: A and B are aqueous, no reactions")
    print("xA = "+str(xA))
    print("xB = "+str(xB))
    print("scaling = "+str(scaling))
    print("including water = "+str(True))
    print()
    model = ConcreteModel()
    model.fs = FlowsheetBlock(default={"dynamic": False})
    model.fs.thermo_params = GenericParameterBlock(default=case1_thermo_config)

    model.fs.rxn_params = GenericReactionParameterBlock(
            default={"property_package": model.fs.thermo_params,
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

    model.fs.unit.inlet.mole_frac_comp[0, "A"].fix( xA )
    model.fs.unit.inlet.mole_frac_comp[0, "B"].fix( xB )
    model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix( 1-xA-xB )

    model.fs.unit.inlet.pressure.fix(101325.0)
    model.fs.unit.inlet.temperature.fix(298.)
    model.fs.unit.outlet.temperature.fix(298.)
    model.fs.unit.inlet.flow_mol.fix(10)

    # Scaling
    if scaling == True:

        # For species
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
            #NOTE: trying to scale solids in this way is significantly worse
            else:
                iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]], 10/scale)
                iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[i], 10/scale)
                iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i], 10/scale)
                iscale.constraint_scaling_transform(
                    model.fs.unit.control_volume.properties_out[0.0].component_flow_balances[i[1]], 10/scale)
                iscale.constraint_scaling_transform(model.fs.unit.control_volume.material_balances[0.0,i[1]], 10/scale)

        iscale.calculate_scaling_factors(model.fs.unit)
    #End scaling if statement

    solver.options['bound_push'] = 1e-20
    solver.options['mu_init'] = 1e-6
    solver.options['max_iter'] = 2000
    model.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)

    solver.options['bound_push'] = 1e-20
    solver.options['mu_init'] = 1e-6
    results = solver.solve(model, tee=True)

    print("==========================================================================")

def run_case2a(xA, xB, scaling=True):
    print("==========================================================================")
    print("Case 2a: A and B are solids, no reactions [includes H2O as solvent]")
    print("xA = "+str(xA))
    print("xB = "+str(xB))
    print("scaling = "+str(scaling))
    print("including water = "+str(True))
    print()
    model = ConcreteModel()
    model.fs = FlowsheetBlock(default={"dynamic": False})
    model.fs.thermo_params = GenericParameterBlock(default=case2a_thermo_config)

    model.fs.rxn_params = GenericReactionParameterBlock(
            default={"property_package": model.fs.thermo_params,
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

    model.fs.unit.inlet.mole_frac_comp[0, "A"].fix( xA )
    model.fs.unit.inlet.mole_frac_comp[0, "B"].fix( xB )
    model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix( 1-xA-xB )

    model.fs.unit.inlet.pressure.fix(101325.0)
    model.fs.unit.inlet.temperature.fix(298.)
    model.fs.unit.outlet.temperature.fix(298.)
    model.fs.unit.inlet.flow_mol.fix(10)

    # Scaling
    if scaling == True:

        # For species
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
            #NOTE: trying to scale solids in this way is significantly worse
            else:
                iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]], 10/scale)
                iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[i], 10/scale)
                iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i], 10/scale)
                iscale.constraint_scaling_transform(
                    model.fs.unit.control_volume.properties_out[0.0].component_flow_balances[i[1]], 10/scale)
                iscale.constraint_scaling_transform(model.fs.unit.control_volume.material_balances[0.0,i[1]], 10/scale)


        iscale.calculate_scaling_factors(model.fs.unit)
    #End scaling if statement

    solver.options['bound_push'] = 1e-20
    solver.options['mu_init'] = 1e-6
    solver.options['max_iter'] = 2000
    model.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)

    solver.options['bound_push'] = 1e-20
    solver.options['mu_init'] = 1e-6
    results = solver.solve(model, tee=True)

    print("==========================================================================")

def run_case2b(xA, scaling=True):
    print("==========================================================================")
    print("Case 2a: A and B are solids, no reactions [NO liquids, not even water]")
    print("xA = "+str(xA))
    print("xB = "+str(1-xA))
    print("scaling = "+str(scaling))
    print("including water = "+str(False))
    print()
    model = ConcreteModel()
    model.fs = FlowsheetBlock(default={"dynamic": False})
    model.fs.thermo_params = GenericParameterBlock(default=case2b_thermo_config)

    model.fs.rxn_params = GenericReactionParameterBlock(
            default={"property_package": model.fs.thermo_params,
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

    model.fs.unit.inlet.mole_frac_comp[0, "A"].fix( xA )
    model.fs.unit.inlet.mole_frac_comp[0, "B"].fix( 1-xA )

    model.fs.unit.inlet.pressure.fix(101325.0)
    model.fs.unit.inlet.temperature.fix(298.)
    model.fs.unit.outlet.temperature.fix(298.)
    model.fs.unit.inlet.flow_mol.fix(10)

    # Scaling
    if scaling == True:

        # For species
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
            #NOTE: trying to scale solids in this way is significantly worse (unless there is no liquid phase)
            else:
                iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]], 10/scale)
                iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[i], 10/scale)
                iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i], 10/scale)
                iscale.constraint_scaling_transform(
                    model.fs.unit.control_volume.properties_out[0.0].component_flow_balances[i[1]], 10/scale)
                iscale.constraint_scaling_transform(model.fs.unit.control_volume.material_balances[0.0,i[1]], 10/scale)


        iscale.calculate_scaling_factors(model.fs.unit)
    #End scaling if statement

    solver.options['bound_push'] = 1e-20
    solver.options['mu_init'] = 1e-6
    solver.options['max_iter'] = 2000
    model.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)

    solver.options['bound_push'] = 1e-20
    solver.options['mu_init'] = 1e-6
    results = solver.solve(model, tee=True)

    print("==========================================================================")

def run_case3(xA, xB, xAB=1e-25, scaling=True):
    print("==========================================================================")
    print("Case 2a: A and B are aqueous, AB is solid that forms from reaction")
    print("xA = "+str(xA))
    print("xB = "+str(xB))
    print("xAB = "+str(xAB))
    print("scaling = "+str(scaling))
    print("including water = "+str(True))
    print()
    model = ConcreteModel()
    model.fs = FlowsheetBlock(default={"dynamic": False})
    model.fs.thermo_params = GenericParameterBlock(default=case3_thermo_config)

    model.fs.rxn_params = GenericReactionParameterBlock(
            default={"property_package": model.fs.thermo_params,
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

    model.fs.unit.inlet.mole_frac_comp[0, "A"].fix( xA )
    model.fs.unit.inlet.mole_frac_comp[0, "B"].fix( xB )
    model.fs.unit.inlet.mole_frac_comp[0, "AB"].fix( xAB )
    model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix( 1-xA-xB-xAB )

    model.fs.unit.inlet.pressure.fix(101325.0)
    model.fs.unit.inlet.temperature.fix(298.)
    model.fs.unit.outlet.temperature.fix(298.)
    model.fs.unit.inlet.flow_mol.fix(10)

    # Scaling
    if scaling == True:

        # For species
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
            #NOTE: trying to scale solids in this way is significantly worse
            else:
                iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]], 10/scale)
                iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[i], 10/scale)
                iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i], 10/scale)
                iscale.constraint_scaling_transform(
                    model.fs.unit.control_volume.properties_out[0.0].component_flow_balances[i[1]], 10/scale)
                iscale.constraint_scaling_transform(model.fs.unit.control_volume.material_balances[0.0,i[1]], 10/scale)


        iscale.calculate_scaling_factors(model.fs.unit)
    #End scaling if statement

    solver.options['bound_push'] = 1e-20
    solver.options['mu_init'] = 1e-6
    solver.options['max_iter'] = 2000
    model.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)

    solver.options['bound_push'] = 1e-20
    solver.options['mu_init'] = 1e-6
    results = solver.solve(model, tee=True)

    print("comp\toutlet.tot_molfrac")
    for i in model.fs.unit.inlet.mole_frac_comp:
        print(str(i[1])+"\t"+str(value(model.fs.unit.outlet.mole_frac_comp[i[0], i[1]])))
    print()

    # NOTE: There is a genuine error in IDAES core associated with solid properties
    #       Code throws errors whenever referencing 'conc_mol_phase_comp' and will
    #       crash if I try to report
    for i in model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp:
        print(str(i)+"\t"+str(value(model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp[i])))

    print("==========================================================================")

if __name__ == "__main__":
    #run_case1(xA=1e-9, xB=1e-9, scaling=True)
    #run_case1(xA=1e-2, xB=1e-9, scaling=True)
    #run_case1(xA=1e-9, xB=1e-2, scaling=True)
    #run_case1(xA=1e-2, xB=1e-2, scaling=True)

    #run_case1(xA=1e-9, xB=1e-9, scaling=False)
    #run_case1(xA=1e-2, xB=1e-9, scaling=False)
    #run_case1(xA=1e-9, xB=1e-2, scaling=False)
    #run_case1(xA=1e-2, xB=1e-2, scaling=False)

    #run_case2a(xA=1e-9, xB=1e-9, scaling=True)
    #run_case2a(xA=1e-2, xB=1e-9, scaling=True)
    #run_case2a(xA=1e-9, xB=1e-2, scaling=True)
    #run_case2a(xA=1e-2, xB=1e-2, scaling=True)

    #run_case2a(xA=1e-9, xB=1e-9, scaling=False)
    #run_case2a(xA=1e-2, xB=1e-9, scaling=False)
    #run_case2a(xA=1e-9, xB=1e-2, scaling=False)
    #run_case2a(xA=1e-2, xB=1e-2, scaling=False)

    #run_case2b(xA=1e-9, scaling=True)
    #run_case2b(xA=1e-2, scaling=True)
    #run_case2b(xA=0.5, scaling=True)

    #run_case2b(xA=1e-9, scaling=False)
    #run_case2b(xA=1e-2, scaling=False)
    #run_case2b(xA=0.5, scaling=False)

    run_case3(xA=1e-9, xB=1e-9, xAB=1e-25, scaling=False)
