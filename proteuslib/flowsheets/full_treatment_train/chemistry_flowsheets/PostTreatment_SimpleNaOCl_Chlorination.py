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

""" Simple NaOCl Chlorination process """

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
                           Var,
                           Constraint,
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

# Grab the scaling utilities
from proteuslib.flowsheets.full_treatment_train.electrolyte_scaling_utils import (
    approximate_chemical_state_args,
    calculate_chemical_scaling_factors)

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

    }
    # End simple_naocl_thermo_config definition

# This config is REQUIRED to use EquilibriumReactor even if we have no equilibrium reactions
simple_naocl_reaction_config = {
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

def build_simple_naocl_chlorination_unit(model,
                                mg_per_L_NaOCl_added = 2,
                                inlet_water_density_kg_per_L = 1,
                                inlet_temperature_K = 298,
                                inlet_pressure_Pa = 101325,
                                inlet_flow_mol_per_s = 10):
    model.fs.simple_naocl_thermo_params = GenericParameterBlock(default=simple_naocl_thermo_config)
    model.fs.simple_naocl_rxn_params = GenericReactionParameterBlock(
            default={"property_package": model.fs.simple_naocl_thermo_params, **simple_naocl_reaction_config})
    model.fs.simple_naocl_unit = EquilibriumReactor(default={
            "property_package": model.fs.simple_naocl_thermo_params,
            "reaction_package": model.fs.simple_naocl_rxn_params,
            "has_rate_reactions": False,
            "has_equilibrium_reactions": True,
            "has_heat_transfer": False,
            "has_heat_of_reaction": False,
            "has_pressure_change": False})

    model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "H_+"].fix( 0. )
    model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "OH_-"].fix( 0. )
    model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "HOCl"].fix( 0. )
    model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "Cl_-"].fix( 0. )

    total_molar_density = inlet_water_density_kg_per_L/18*1000 #mol/L

    # Free Chlorine (mg-Cl2/L) = total_chlorine_inlet (mol/L) * 70,900
    #       Assumes chlorine is added as NaOCl
    free_chlorine_added = mg_per_L_NaOCl_added/74.44/1000*70900 #mg/L as NaOCl
    total_chlorine_inlet = free_chlorine_added/70900 # mol/L
    total_molar_density+=total_chlorine_inlet

    model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "OCl_-"].fix( total_chlorine_inlet/total_molar_density )
    model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "Na_+"].fix( total_chlorine_inlet/total_molar_density )

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

    dr = model.fs.simple_naocl_unit.inlet.flow_mol[0].value*model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "OCl_-"].value
    dr = dr*74.44*1000
    model.fs.simple_naocl_unit.dosing_rate = Var(initialize=dr)
    model.fs.simple_naocl_unit.dosing_rate.fix()
    model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "OCl_-"].unfix()

    def _dosing_rate_cons(blk):
        return blk.dosing_rate == blk.inlet.flow_mol[0]*blk.inlet.mole_frac_comp[0, "OCl_-"]*74.44*1000

    model.fs.simple_naocl_unit.dosing_cons = Constraint( rule=_dosing_rate_cons )

    fc = model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "HOCl"].value*total_molar_density
    fc += model.fs.simple_naocl_unit.inlet.mole_frac_comp[0, "OCl_-"].value*total_molar_density
    fc = fc*70900

    model.fs.simple_naocl_unit.free_chlorine = Var(initialize=fc)

    def _free_chlorine_cons(blk):
        return blk.free_chlorine == ((blk.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq","HOCl"]/1000) \
                                    + (blk.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq","OCl_-"]/1000)) \
                                    * 70900
    model.fs.simple_naocl_unit.chlorine_cons = Constraint( rule=_free_chlorine_cons )

def initialize_chlorination_example(unit, state_args, user_scaling=True, debug_out=False):
    solver.options['bound_push'] = 1e-10
    solver.options['mu_init'] = 1e-6

    if user_scaling == True:
        solver.options["nlp_scaling_method"] = "user-scaling"

    unit.inlet.mole_frac_comp[0, "OCl_-"].fix()
    unit.free_chlorine.fix()
    unit.chlorine_cons.deactivate()
    unit.dosing_cons.deactivate()

    if debug_out == True:
        unit.initialize(state_args=state_args, optarg=solver.options, outlvl=idaeslog.DEBUG)
    else:
        unit.initialize(state_args=state_args, optarg=solver.options)

    unit.inlet.mole_frac_comp[0, "OCl_-"].unfix()
    unit.free_chlorine.unfix()
    unit.chlorine_cons.activate()
    unit.dosing_cons.activate()

    iscale.constraint_autoscale_large_jac(unit)

def display_results_of_chlorination(chlorination_unit):
    print()
    print("=========== Chlorination Results ============")
    print("Outlet Temperature:       \t" + str(chlorination_unit.outlet.temperature[0].value))
    print("Outlet Pressure:          \t" + str(chlorination_unit.outlet.pressure[0].value))
    print("Outlet FlowMole:          \t" + str(chlorination_unit.outlet.flow_mol[0].value))
    print()
    total_molar_density = \
        value(chlorination_unit.control_volume.properties_out[0.0].dens_mol_phase['Liq'])/1000
    pH = -value(log10(chlorination_unit.outlet.mole_frac_comp[0, "H_+"]*total_molar_density))
    print("pH at Outlet:             \t" + str(pH))
    total_salt = value(chlorination_unit.outlet.mole_frac_comp[0, "Na_+"])*total_molar_density*23
    total_salt += value(chlorination_unit.outlet.mole_frac_comp[0, "Cl_-"])*total_molar_density*35.4
    psu = total_salt/(total_molar_density*18)*1000
    print("Salinity (PSU):           \t" + str(psu))
    print("NaOCl Dosing Rate (mg/s): \t" + str(chlorination_unit.dosing_rate.value))
    print("Free Chlorine (mg/L):     \t" + str(chlorination_unit.free_chlorine.value))
    print("\tDistribution:")
    hocl = (value(chlorination_unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq","HOCl"])/1000)/(chlorination_unit.free_chlorine.value/70900)
    print("\t % HOCl: \t" + str(hocl*100))
    ocl = (value(chlorination_unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq","OCl_-"])/1000)/(chlorination_unit.free_chlorine.value/70900)
    print("\t % OCl-: \t" + str(ocl*100))
    print("-------------------------------------------")
    print()

def run_chlorination_example():
    model = ConcreteModel()
    model.fs = FlowsheetBlock(default={"dynamic": False})

    build_simple_naocl_chlorination_unit(model, mg_per_L_NaOCl_added = 2)
    state_args, stoich_extents = approximate_chemical_state_args(model.fs.simple_naocl_unit,
                                model.fs.simple_naocl_rxn_params, simple_naocl_reaction_config)

    calculate_chemical_scaling_factors(model.fs.simple_naocl_unit,
                                model.fs.simple_naocl_thermo_params,
                                model.fs.simple_naocl_rxn_params, state_args)

    initialize_chlorination_example(model.fs.simple_naocl_unit, state_args)

    solve_with_user_scaling(model, tee=True, bound_push=1e-10, mu_init=1e-6)

    display_results_of_chlorination(model.fs.simple_naocl_unit)

    return model

# This example will try to find the dosing rate needed to yield 2mg/L free chlorine
def run_chlorination_constrained_outlet_example():
    model = ConcreteModel()
    model.fs = FlowsheetBlock(default={"dynamic": False})

    # Give bad initial guess for the dosing rate
    build_simple_naocl_chlorination_unit(model, mg_per_L_NaOCl_added = 1)
    state_args, stoich_extents = approximate_chemical_state_args(model.fs.simple_naocl_unit,
                                model.fs.simple_naocl_rxn_params, simple_naocl_reaction_config)

    calculate_chemical_scaling_factors(model.fs.simple_naocl_unit,
                                model.fs.simple_naocl_thermo_params,
                                model.fs.simple_naocl_rxn_params, state_args)

    initialize_chlorination_example(model.fs.simple_naocl_unit, state_args)

    #Redefine the constraints and fixed vars here
    model.fs.simple_naocl_unit.dosing_rate.unfix()
    model.fs.simple_naocl_unit.free_chlorine.fix(2)

    solve_with_user_scaling(model, tee=True, bound_push=1e-10, mu_init=1e-6)

    display_results_of_chlorination(model.fs.simple_naocl_unit)

    return model

if __name__ == "__main__":
    model = run_chlorination_example()
    model = run_chlorination_constrained_outlet_example()
