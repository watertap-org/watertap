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
    This test is to establish that the core chemistry packages in IDAES solve
    a simple water dissociation problem and return the correct pH value, as well
    as excerising the implementation of the ENRTL model within that same context.

    This test also includes checking the pH value of a typical acid and checking
    the calculated activity coefficients for a relatively dilute system of ions.
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
from idaes.generic_models.properties.core.pure.electrolyte import relative_permittivity_constant
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.enrtl import ENRTL
from idaes.generic_models.properties.core.eos.enrtl_reference_states import Unsymmetric

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

# Import the scaling methods
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
water_thermo_config = {
    "components": {
        'H2O': {"type": Solvent,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              "relative_permittivity_liq_comp": relative_permittivity_constant,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (18.0153, pyunits.g/pyunits.mol),
                    "relative_permittivity_liq_comp": 78.54,
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
                    }
              },
              # End Component list
        "phases":  {'Liq': {"type": AqueousPhase,
                            "equation_of_state": ENRTL,
                            "equation_of_state_options": {
                                    "reference_state": Unsymmetric}
                            },
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
                   "concentration_form": ConcentrationForm.activity,
                   "parameter_data": {
                       #NOTE: The k value on the activity basis is UNITLESS
                       #        based on a standard molar concentration of 1 mol/L
                       #        HOWEVER, the typical Kw dissociation constant of
                       #        1e-14 is defined on a molar basis. Thus, we must
                       #        divide by the total (~55.2 M) concentration raised to the
                       #        net reaction order (i.e., 2 in this case).
                       "dh_rxn_ref": (0, pyunits.kJ/pyunits.mol),
                       "k_eq_ref": (10**-14/55.2**2,pyunits.dimensionless),
                       "T_eq_ref": (298, pyunits.K),

                       # By default, reaction orders follow stoichiometry
                       #    manually set reaction order here to override
                       "reaction_order": {("Liq", "H2O"): 0,
                                        ("Liq", "H_+"): 1,
                                        ("Liq", "OH_-"): 1}
                        }
                        # End parameter_data
                   }
                   # End R1
             }
             # End equilibrium_reactions
    }
    # End thermo_config definition

# Define the reaction_config for water dissociation
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
class TestENRTLwater():
    @pytest.fixture(scope="class")
    def water_model(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(default={"dynamic": False})
        model.fs.thermo_params = GenericParameterBlock(default=water_thermo_config)
        model.fs.rxn_params = GenericReactionParameterBlock(
                default={"property_package": model.fs.thermo_params, **reaction_config})
        model.fs.unit = EquilibriumReactor(default={
                "property_package": model.fs.thermo_params,
                "reaction_package": model.fs.rxn_params,
                "has_rate_reactions": False,
                "has_equilibrium_reactions": False,
                "has_heat_transfer": False,
                "has_heat_of_reaction": False,
                "has_pressure_change": False})

        #NOTE: ENRTL model cannot initialize if the inlet values are 0
        zero = 1e-20
        model.fs.unit.inlet.mole_frac_comp[0, "H_+"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "OH_-"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix( 1.-2*zero )
        model.fs.unit.inlet.pressure.fix(101325.0)
        model.fs.unit.inlet.temperature.fix(298.)
        model.fs.unit.inlet.flow_mol.fix(10)

        return model

    @pytest.mark.unit
    def test_build_water(self, water_model):
        model = water_model

        assert hasattr(model.fs.thermo_params, 'component_list')
        assert len(model.fs.thermo_params.component_list) == 3
        assert 'H2O' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.H2O, Solvent)
        assert 'H_+' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('H_+'), Cation)
        assert 'OH_-' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('OH_-'), Anion)

    @pytest.mark.unit
    def test_units_water(self, water_model):
        model = water_model
        assert_units_consistent(model)

    @pytest.mark.unit
    def test_dof_water(self, water_model):
        model = water_model
        assert (degrees_of_freedom(model) == 0)

    @pytest.mark.component
    def test_scaling_water(self, water_model):
        model = water_model

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
    def test_initialize_solver_water(self, water_model):
        model = water_model
        solver.options['bound_push'] = 1e-20
        solver.options['mu_init'] = 1e-6
        model.fs.unit.initialize(optarg=solver.options)
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_solve_water(self, water_model):
        model = water_model
        solver.options['max_iter'] = 2
        results = solver.solve(model)
        print(results.solver.termination_condition)
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution_water(self, water_model):
        model = water_model

        pH = -log10(value(model.fs.unit.control_volume.properties_out[0.0].act_phase_comp["Liq","H_+"])*55.2)
        pOH = -log10(value(model.fs.unit.control_volume.properties_out[0.0].act_phase_comp["Liq","OH_-"])*55.2)

        assert pytest.approx(7.00000, rel=1e-5) == pH
        assert pytest.approx(7.00000, rel=1e-5) == pOH

# Configuration dictionary
carbonic_thermo_config = {
    "components": {
        'H2O': {"type": Solvent,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              "relative_permittivity_liq_comp": relative_permittivity_constant,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (18.0153, pyunits.g/pyunits.mol),
                    "relative_permittivity_liq_comp": 78.54,
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
                    }
              },
              # End Component list
        "phases":  {'Liq': {"type": AqueousPhase,
                            "equation_of_state": ENRTL,
                            "equation_of_state_options": {
                                    "reference_state": Unsymmetric}
                            },
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
                   "concentration_form": ConcentrationForm.activity,
                   "parameter_data": {
                       #NOTE: The k value on the activity basis is UNITLESS
                       #        based on a standard molar concentration of 1 mol/L
                       #        HOWEVER, the typical Kw dissociation constant of
                       #        1e-14 is defined on a molar basis. Thus, we must
                       #        divide by the total (~55.2 M) concentration raised to the
                       #        net reaction order (i.e., 2 in this case).
                       "dh_rxn_ref": (0, pyunits.kJ/pyunits.mol),
                       "k_eq_ref": (10**-14/55.2**2, pyunits.dimensionless),
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
                   "concentration_form": ConcentrationForm.activity,
                   "parameter_data": {
                       #NOTE: The k value on the activity basis is UNITLESS
                       #        based on a standard molar concentration of 1 mol/L
                       #        HOWEVER, the typical Ka1 dissociation constant of
                       #        10**-6.33 is defined on a molar basis. Thus, we must
                       #        divide by the total (~55.2 M) concentration raised to the
                       #        net reaction order (i.e., 1 in this case). 
                       "dh_rxn_ref": (0, pyunits.kJ/pyunits.mol),
                       "k_eq_ref": (10**-6.33/55.2,pyunits.dimensionless),
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
                   "equilibrium_constant": van_t_hoff,
                   "equilibrium_form": log_power_law_equil,
                   "concentration_form": ConcentrationForm.activity,
                   "parameter_data": {
                       #NOTE: The k value on the activity basis is UNITLESS
                       #        based on a standard molar concentration of 1 mol/L
                       #        HOWEVER, the typical Ka1 dissociation constant of
                       #        10**-10.35 is defined on a molar basis. Thus, we must
                       #        divide by the total (~55.2 M) concentration raised to the
                       #        net reaction order (i.e., 1 in this case).                       
                       "dh_rxn_ref": (0, pyunits.kJ/pyunits.mol),
                       "k_eq_ref": (10**-10.35/55.2,pyunits.dimensionless),
                       "T_eq_ref": (300, pyunits.K),

                       # By default, reaction orders follow stoichiometry
                       #    manually set reaction order here to override
                       "reaction_order": {("Liq", "HCO3_-"): -1,
                                        ("Liq", "H_+"): 1,
                                        ("Liq", "CO3_2-"): 1}
                        }
                        # End parameter_data
                   }
                   # End R3
             }
             # End equilibrium_reactions
    }
    # End thermo_config definition

# Start test class
class TestENRTLcarbonicAcid():
    @pytest.fixture(scope="class")
    def carbonic_acid_model(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(default={"dynamic": False})
        model.fs.thermo_params = GenericParameterBlock(default=carbonic_thermo_config)
        model.fs.rxn_params = GenericReactionParameterBlock(
                default={"property_package": model.fs.thermo_params, **reaction_config})
        model.fs.unit = EquilibriumReactor(default={
                "property_package": model.fs.thermo_params,
                "reaction_package": model.fs.rxn_params,
                "has_rate_reactions": False,
                "has_equilibrium_reactions": False,
                "has_heat_transfer": False,
                "has_heat_of_reaction": False,
                "has_pressure_change": False})

        #NOTE: ENRTL model cannot initialize if the inlet values are 0
        zero = 1e-20
        acid = 0.00206/(55.2+0.00206)
        model.fs.unit.inlet.mole_frac_comp[0, "H_+"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "OH_-"].fix( zero )

        # Added as conjugate base form
        model.fs.unit.inlet.mole_frac_comp[0, "CO3_2-"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "HCO3_-"].fix( acid )
        model.fs.unit.inlet.mole_frac_comp[0, "H2CO3"].fix( zero )
        model.fs.unit.inlet.mole_frac_comp[0, "Na_+"].fix( acid )

        model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix( 1.-4*zero-acid- \
                        value(model.fs.unit.inlet.mole_frac_comp[0, "Na_+"]) )
        model.fs.unit.inlet.pressure.fix(101325.0)
        model.fs.unit.inlet.temperature.fix(298.)
        model.fs.unit.inlet.flow_mol.fix(10)

        return model

    @pytest.mark.unit
    def test_build_carbonic_acid(self, carbonic_acid_model):
        model = carbonic_acid_model

        assert hasattr(model.fs.thermo_params, 'component_list')
        assert len(model.fs.thermo_params.component_list) == 7
        assert 'H2O' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.H2O, Solvent)
        assert 'H_+' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('H_+'), Cation)
        assert 'OH_-' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('OH_-'), Anion)

        assert 'Na_+' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('Na_+'), Cation)
        assert 'HCO3_-' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('HCO3_-'), Anion)
        assert 'CO3_2-' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('CO3_2-'), Anion)
        assert 'H2CO3' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.H2CO3, Solute)

    #NOTE: This test takes ~8s and I have no idea why...
    @pytest.mark.unit
    def test_units_carbonic_acid(self, carbonic_acid_model):
        model = carbonic_acid_model
        assert_units_consistent(model)

    @pytest.mark.unit
    def test_dof_carbonic_acid(self, carbonic_acid_model):
        model = carbonic_acid_model
        assert (degrees_of_freedom(model) == 0)

    @pytest.mark.component
    def test_scaling_carbonic_acid(self, carbonic_acid_model):
        model = carbonic_acid_model

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
    def test_initialize_solver_carbonic_acid(self, carbonic_acid_model):
        model = carbonic_acid_model
        solver.options['bound_push'] = 1e-20
        solver.options['mu_init'] = 1e-6
        solver.options['max_iter'] = 1000
        model.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_solve_carbonic_acid(self, carbonic_acid_model):
        model = carbonic_acid_model
        solver.options['max_iter'] = 2
        results = solver.solve(model)
        print(results.solver.termination_condition)
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution_carbonic_acid(self, carbonic_acid_model):
        model = carbonic_acid_model

        pH = -log10(value(model.fs.unit.control_volume.properties_out[0.0].act_phase_comp["Liq","H_+"])*55.2)
        pOH = -log10(value(model.fs.unit.control_volume.properties_out[0.0].act_phase_comp["Liq","OH_-"])*55.2)

        assert pytest.approx(8.28, rel=1e-2) == pH
        assert pytest.approx(5.72, rel=1e-2) == pOH
        assert pytest.approx(14.00, rel=1e-2) == pH + pOH

        gamma = {}
        for index in model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp:
            gamma[index] = value(model.fs.unit.control_volume.properties_out[0.0].act_phase_comp[index])/value(model.fs.unit.outlet.mole_frac_comp[0,index[1]])

        assert pytest.approx(0.95075, rel=1e-5) == gamma[('Liq', 'OH_-')]
        assert pytest.approx(0.95075, rel=1e-5) == gamma[('Liq', 'H_+')]
        assert pytest.approx(0.95075, rel=1e-5) == gamma[('Liq', 'Na_+')]
        assert pytest.approx(0.95075, rel=1e-5) == gamma[('Liq', 'HCO3_-')]
        assert pytest.approx(0.81710, rel=1e-5) == gamma[('Liq', 'CO3_2-')]
        assert pytest.approx(1.00000, rel=1e-5) == gamma[('Liq', 'H2O')]
        assert pytest.approx(1.00000, rel=1e-5) == gamma[('Liq', 'H2CO3')]

        total_acid = value(model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq","H2CO3"])/1000
        total_acid += value(model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq","HCO3_-"])/1000
        total_acid += value(model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp["Liq","CO3_2-"])/1000

        assert pytest.approx(0.002061178349769601, rel=1e-5) == total_acid
