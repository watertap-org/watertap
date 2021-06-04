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

'''
    This test is to establish that the core chemistry packages in IDAES solve
    a simple water dissociation problem and return the correct pH value.
'''
# Importing testing libraries
import unittest
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

# Import built-in Gibb's Energy function
from idaes.generic_models.properties.core.reactions.equilibrium_constant import gibbs_energy

# Import built-in van't Hoff function
from idaes.generic_models.properties.core.reactions.equilibrium_constant import van_t_hoff

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

# Import the core idaes objects for Flowsheets and types of balances
from idaes.core import FlowsheetBlock

# Import log10 function from pyomo
from pyomo.environ import log10

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
                   }
                   # End R1
             }
             # End equilibrium_reactions
    }
    # End thermo_config definition

# Define the reaction_config for water dissociation
water_reaction_config = {
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
               }
               # End R1
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
class TestPureWater():
    @pytest.fixture(scope="class")
    def inherent_reactions_config(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(default={"dynamic": False})
        model.fs.thermo_params = GenericParameterBlock(default=thermo_config)
        model.fs.rxn_params = GenericReactionParameterBlock(
                default={"property_package": model.fs.thermo_params, **water_reaction_config})
        model.fs.unit = EquilibriumReactor(default={
                "property_package": model.fs.thermo_params,
                "reaction_package": model.fs.rxn_params,
                "has_rate_reactions": False,
                "has_equilibrium_reactions": False,
                "has_heat_transfer": False,
                "has_heat_of_reaction": False,
                "has_pressure_change": False})

        model.fs.unit.inlet.mole_frac_comp[0, "H_+"].fix( 0. )
        model.fs.unit.inlet.mole_frac_comp[0, "OH_-"].fix( 0. )
        model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix( 1. )
        model.fs.unit.inlet.pressure.fix(101325.0)
        model.fs.unit.inlet.temperature.fix(298.)
        model.fs.unit.inlet.flow_mol.fix(10)

        return model

    @pytest.fixture(scope="class")
    def equilibrium_reactions_config(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(default={"dynamic": False})
        model.fs.thermo_params = GenericParameterBlock(default=thermo_only_config)
        model.fs.rxn_params = GenericReactionParameterBlock(
                default={"property_package": model.fs.thermo_params, **water_reaction_config})
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
        model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix( 1. )
        model.fs.unit.inlet.pressure.fix(101325.0)
        model.fs.unit.inlet.temperature.fix(298.)
        model.fs.unit.inlet.flow_mol.fix(10)

        return model

    @pytest.mark.unit
    def test_build_model_inherent(self, inherent_reactions_config):
        model = inherent_reactions_config

        assert hasattr(model.fs.thermo_params, 'component_list')
        assert len(model.fs.thermo_params.component_list) == 3
        assert 'H2O' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.H2O, Solvent)
        assert 'H_+' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('H_+'), Cation)
        assert 'OH_-' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('OH_-'), Anion)

        assert hasattr(model.fs.thermo_params, 'phase_list')
        assert len(model.fs.thermo_params.phase_list) == 1
        assert isinstance(model.fs.thermo_params.Liq, AqueousPhase)

    @pytest.mark.unit
    def test_build_model_equilibrium(self, equilibrium_reactions_config):
        model = equilibrium_reactions_config

        assert hasattr(model.fs.thermo_params, 'component_list')
        assert len(model.fs.thermo_params.component_list) == 3
        assert 'H2O' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.H2O, Solvent)
        assert 'H_+' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('H_+'), Cation)
        assert 'OH_-' in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component('OH_-'), Anion)

        assert hasattr(model.fs.thermo_params, 'phase_list')
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
        assert (degrees_of_freedom(model) == 0)

    @pytest.mark.unit
    def test_dof_equilibrium(self, equilibrium_reactions_config):
        model = equilibrium_reactions_config
        assert (degrees_of_freedom(model) == 0)

    @pytest.mark.unit
    def test_stats_inherent(self, inherent_reactions_config):
        model = inherent_reactions_config
        assert (number_variables(model) == 77)
        assert (number_total_constraints(model) == 24)
        assert (number_unused_variables(model) == 13)

    @pytest.mark.unit
    def test_stats_equilibrium(self, equilibrium_reactions_config):
        model = equilibrium_reactions_config
        assert (number_variables(model) == 71)
        assert (number_total_constraints(model) == 24)
        assert (number_unused_variables(model) == 7)

    @pytest.mark.component
    def test_scaling_inherent(self, inherent_reactions_config):
        model = inherent_reactions_config
        iscale.calculate_scaling_factors(model.fs.unit)

        assert isinstance(model.fs.unit.control_volume.scaling_factor, Suffix)

        assert isinstance(model.fs.unit.control_volume.properties_out[0.0].scaling_factor, Suffix)

        assert isinstance(model.fs.unit.control_volume.properties_in[0.0].scaling_factor, Suffix)

    @pytest.mark.component
    def test_scaling_equilibrium(self, equilibrium_reactions_config):
        model = equilibrium_reactions_config
        iscale.calculate_scaling_factors(model.fs.unit)

        assert isinstance(model.fs.unit.control_volume.scaling_factor, Suffix)

        assert isinstance(model.fs.unit.control_volume.properties_out[0.0].scaling_factor, Suffix)

        assert isinstance(model.fs.unit.control_volume.properties_in[0.0].scaling_factor, Suffix)

        # When using equilibrium reactions, there are another set of scaling factors calculated
        assert isinstance(model.fs.unit.control_volume.reactions[0.0].scaling_factor, Suffix)

    @pytest.mark.component
    def test_initialize_inherent(self, inherent_reactions_config):
        model = inherent_reactions_config
        state_args = {'mole_frac_comp':
                        {   'H2O': 1,
                            'H_+': 10**-7/55.6,
                            'OH_-': 10**-7/55.6
                        },
                        'pressure': 101325,
                        'temperature': 298,
                        'flow_mol': 10
                    }
        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        solver.options['bound_push'] = 1e-10
        solver.options['mu_init'] = 1e-6
        model.fs.unit.initialize(state_args=state_args, optarg=solver.options)

        fin_fixed_vars = fixed_variables_set(model)
        fin_act_consts = activated_constraints_set(model)

        assert degrees_of_freedom(model) == 0

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

    @pytest.mark.component
    def test_initialize_equilibrium(self, equilibrium_reactions_config):
        model = equilibrium_reactions_config
        state_args = {'mole_frac_comp':
                        {   'H2O': 1,
                            'H_+': 10**-7/55.6,
                            'OH_-': 10**-7/55.6
                        },
                        'pressure': 101325,
                        'temperature': 298,
                        'flow_mol': 10
                    }
        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        solver.options['bound_push'] = 1e-10
        solver.options['mu_init'] = 1e-6
        model.fs.unit.initialize(state_args=state_args, optarg=solver.options)

        fin_fixed_vars = fixed_variables_set(model)
        fin_act_consts = activated_constraints_set(model)

        assert degrees_of_freedom(model) == 0

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

    @pytest.mark.component
    def test_solve_inherent(self, inherent_reactions_config):
        model = inherent_reactions_config
        solver.options['max_iter'] = 2
        results = solver.solve(model)
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solve_equilibrium(self, equilibrium_reactions_config):
        model = equilibrium_reactions_config
        solver.options['max_iter'] = 2
        results = solver.solve(model)
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution_inherent(self, inherent_reactions_config):
        model = inherent_reactions_config

        assert pytest.approx(298, rel=1e-5) == value(model.fs.unit.outlet.temperature[0])
        assert pytest.approx(10, rel=1e-5) == value(model.fs.unit.outlet.flow_mol[0])
        assert pytest.approx(101325, rel=1e-5) == value(model.fs.unit.outlet.pressure[0])

        total_molar_density = \
            value(model.fs.unit.control_volume.properties_out[0.0].dens_mol_phase['Liq'])/1000
        assert pytest.approx(55.2336, rel=1e-5) == total_molar_density
        pH = -value(log10(model.fs.unit.outlet.mole_frac_comp[0, "H_+"]*total_molar_density))
        pOH = -value(log10(model.fs.unit.outlet.mole_frac_comp[0, "OH_-"]*total_molar_density))
        assert pytest.approx(7.00059, rel=1e-5) == pH
        assert pytest.approx(7.00059, rel=1e-5) == pOH
        assert pytest.approx(0.99999999, rel=1e-5) == value(model.fs.unit.outlet.mole_frac_comp[0.0, 'H2O'])

    @pytest.mark.component
    def test_solution_equilibrium(self, equilibrium_reactions_config):
        model = equilibrium_reactions_config

        assert pytest.approx(298, rel=1e-5) == value(model.fs.unit.outlet.temperature[0])
        assert pytest.approx(10, rel=1e-5) == value(model.fs.unit.outlet.flow_mol[0])
        assert pytest.approx(101325, rel=1e-5) == value(model.fs.unit.outlet.pressure[0])

        total_molar_density = \
            value(model.fs.unit.control_volume.properties_out[0.0].dens_mol_phase['Liq'])/1000
        assert pytest.approx(55.2336, rel=1e-5) == total_molar_density
        pH = -value(log10(model.fs.unit.outlet.mole_frac_comp[0, "H_+"]*total_molar_density))
        pOH = -value(log10(model.fs.unit.outlet.mole_frac_comp[0, "OH_-"]*total_molar_density))
        assert pytest.approx(7.00059, rel=1e-5) == pH
        assert pytest.approx(7.00059, rel=1e-5) == pOH
        assert pytest.approx(0.99999999, rel=1e-5) == value(model.fs.unit.outlet.mole_frac_comp[0.0, 'H2O'])
