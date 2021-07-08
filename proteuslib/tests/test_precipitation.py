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
    This test is to establish...

    Reactions:
        H2O <---> H + OH
    Precipitation:
        Ca + 2 OH --> Ca(OH)2
"""
# Importing testing libraries
import pytest

# Importing the object for units from pyomo
from pyomo.environ import units as pyunits

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

__author__ = "Austin Ladshaw"

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
        'OH_-': {"type": Anion, "charge": -1,
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
        'Ca(OH)2': {"type": Component, "valid_phase_types": PT.solidPhase,
              # Define the methods used to calculate the following properties
              #     Nothing currently exists for solids in IDAES, so we
              #     are passing a dummy function that does nothing
              "dens_mol_sol_comp": dummy_dens,
              "enth_mol_sol_comp": dummy_h,
              "cp_mol_sol_comp": dummy_cp,
              "entr_mol_sol_comp": dummy_s,
                    },
              },
              # End Component list
        "phases":  {'Liq': {"type": AqueousPhase,
                            "equation_of_state": Ideal},
                    'Sol': {"type": SolidPhase,
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

if __name__ == "__main__":
    model = ConcreteModel()
    model.fs = FlowsheetBlock(default={"dynamic": False})
    model.fs.thermo_params = GenericParameterBlock(default=thermo_config)
    model.fs.rxn_params = GenericReactionParameterBlock(
            default={"property_package": model.fs.thermo_params, **reaction_config})

    # Cannot build unit model without methods for solids
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
    model.fs.unit.inlet.mole_frac_comp[0, "Ca(OH)2"].fix( zero )

    total_molar_density = 55.2  # mol/L (approximate density of seawater)
    total_base = 1e-10
    total_ca = 1e-4

    model.fs.unit.inlet.mole_frac_comp[0, "OH_-"].fix( total_base/total_molar_density )
    model.fs.unit.inlet.mole_frac_comp[0, "Ca_2+"].fix( total_ca/total_molar_density )

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

    # This is in a try block because a unit model may not have inherent reactions
    '''try:
        for rid in model.fs.thermo_params.inherent_reaction_idx:
            scale = value(model.fs.unit.control_volume.properties_out[0.0].k_eq[rid].expr)
            # Want to set eps in some fashion similar to this
            if scale < 1e-16:
                model.fs.thermo_params.component("reaction_"+rid).eps.value = scale*1e-2
            else:
                model.fs.thermo_params.component("reaction_"+rid).eps.value = 1e-16*1e-2

        #Add scaling factors for reactions
        for i in model.fs.unit.control_volume.inherent_reaction_extent_index:
            # i[0] = time, i[1] = reaction
            scale = value(model.fs.unit.control_volume.properties_out[0.0].k_eq[i[1]].expr)
            iscale.set_scaling_factor(model.fs.unit.control_volume.inherent_reaction_extent[0.0,i[1]], 10/scale)
            iscale.constraint_scaling_transform(model.fs.unit.control_volume.properties_out[0.0].
                    inherent_equilibrium_constraint[i[1]], 0.1)
    except:
        pass

    #exit()

    # Next, try adding scaling for species
    min = 1e-10
    for i in model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
        # i[0] = phase, i[1] = species
        if model.fs.unit.inlet.mole_frac_comp[0, i[1]].value > min:
            scale = model.fs.unit.inlet.mole_frac_comp[0, i[1]].value
        else:
            scale = min
        iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]],
            10/scale)
        iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[i],
            10/scale)
        iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i],
            10/scale)
        iscale.constraint_scaling_transform(
            model.fs.unit.control_volume.properties_out[0.0].component_flow_balances[i[1]], 10/scale)
        iscale.constraint_scaling_transform(model.fs.unit.control_volume.material_balances[0.0,i[1]], 10/scale)'''

    #iscale.calculate_scaling_factors(model.fs.unit)

    # NOTE: Right now this works better without scaling???

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

    print()

    pH = -value(log10(model.fs.unit.outlet.mole_frac_comp[0, "H_+"]*total_molar_density))
    pOH = -value(log10(model.fs.unit.outlet.mole_frac_comp[0, "OH_-"]*total_molar_density))

    print("pH =\t" + str(pH))
    print("pOH =\t" + str(pOH))
