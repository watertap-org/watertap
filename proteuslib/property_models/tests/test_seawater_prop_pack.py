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
import pytest
import proteuslib.property_models.seawater_prop_pack as props
from idaes.generic_models.properties.tests.test_harness import \
    PropertyTestHarness as PropertyTestHarness_idaes
from proteuslib.property_models.tests.property_test_harness import (PropertyTestHarness,
                                                                    # PropertyComponentTestHarness,
                                                                    # PropertyRegressionTest
                                                                    )
from pyomo.environ import (ConcreteModel,
                           Param,
                           Expression,
                           Var,
                           units as pyunits)
from idaes.core import (FlowsheetBlock,
                        MaterialFlowBasis,
                        MaterialBalanceType,
                        EnergyBalanceType)
from idaes.core.components import Solute, Solvent
from idaes.core.phases import LiquidPhase

from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables,
                                              unused_variables_set)
from idaes.core.util import get_solver
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     set_scaling_factor,
                                     unscaled_variables_generator,
                                     unscaled_constraints_generator,
                                     badly_scaled_var_generator)
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()
solver.options["nlp_scaling_method"] = "user-scaling"

# -----------------------------------------------------------------------------
class TestSeawaterProperty_idaes(PropertyTestHarness_idaes):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}
        self.prop_args = {}
        self.has_density_terms = False


class TestSeawaterProperty(PropertyTestHarness):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}
        self.standard_scaling = {('flow_mass_phase_comp', ('Liq', 'H2O')): 1,
                                 ('flow_mass_phase_comp', ('Liq', 'TDS')): 1e2}
        self.stateblock_statistics = {'number_variables': 24,
                                      'number_total_constraints': 20,
                                      'number_unused_variables': 1,  # pressure is unused
                                      'default_degrees_of_freedom': 3}  # 4 state vars, but pressure is not active
        self.default_solution = {('mass_frac_phase_comp', ('Liq', 'H2O')): 0.965,
                                 ('mass_frac_phase_comp', ('Liq', 'TDS')): 0.035,
                                 ('dens_mass_phase', 'Liq'): 1023.5,
                                 ('dens_mass_solvent', None): 996.9,
                                 ('flow_vol_phase', 'Liq'): 9.770e-4,
                                 ('conc_mass_phase_comp', ('Liq', 'H2O')): 987.7,
                                 ('conc_mass_phase_comp', ('Liq', 'TDS')): 35.82,
                                 ('flow_mol_phase_comp', ('Liq', 'H2O')): 53.57,
                                 ('flow_mol_phase_comp', ('Liq', 'TDS')): 1.1145,
                                 ('mole_frac_phase_comp', ('Liq', 'H2O')): 0.9796,
                                 ('mole_frac_phase_comp', ('Liq', 'TDS')): 2.038e-2,
                                 ('molality_comp', 'TDS'): 1.155,
                                 ('visc_d_phase', 'Liq'): 9.588e-4,
                                 ('osm_coeff', None): 0.9068,
                                 ('pressure_osm', None): 2.588e6,
                                 ('enth_mass_phase', 'Liq'): 9.9765e4,
                                 ('pressure_sat', None): 3111,
                                 ('cp_phase', 'Liq'): 4001,
                                 ('therm_cond_phase', 'Liq'): 0.6086,
                                 ('dh_vap', None): 2.356e6}
#
#
# @pytest.mark.component
# class TestSeawaterPropertyComponent(PropertyComponentTestHarness):
#     def configure(self):
#         self.prop_pack = props.SeawaterParameterBlock
#         self.param_args = {}
#         self.set_default_scaling_dict = {('flow_mass_phase_comp', ('Liq', 'H2O')): 1,
#                                          ('flow_mass_phase_comp', ('Liq', 'TDS')): 1e2}
#
#         self.solver = 'ipopt'
#         self.optarg = {'nlp_scaling_method': 'user-scaling'}
#
#         self.default_state_values_dict = {('flow_mass_phase_comp', ('Liq', 'H2O')): 0.965,
#                                           ('flow_mass_phase_comp', ('Liq', 'TDS')): 0.035,
#                                           ('temperature', None): 298.15,
#                                           ('pressure', None): 101325}
#         self.default_initialize_results_dict = {('mass_frac_phase_comp', ('Liq', 'H2O')): 0.965,
#                                                 ('mass_frac_phase_comp', ('Liq', 'TDS')): 0.035,
#                                                 ('dens_mass_phase', 'Liq'): 1023.5,
#                                                 ('dens_mass_solvent', None): 996.9,
#                                                 ('flow_vol_phase', 'Liq'): 9.770e-4,
#                                                 ('conc_mass_phase_comp', ('Liq', 'H2O')): 987.7,
#                                                 ('conc_mass_phase_comp', ('Liq', 'TDS')): 35.82,
#                                                 ('flow_mol_phase_comp', ('Liq', 'H2O')): 53.57,
#                                                 ('flow_mol_phase_comp', ('Liq', 'TDS')): 1.1145,
#                                                 ('mole_frac_phase_comp', ('Liq', 'H2O')): 0.9796,
#                                                 ('mole_frac_phase_comp', ('Liq', 'TDS')): 2.038e-2,
#                                                 ('molality_comp', 'TDS'): 1.155,
#                                                 ('visc_d_phase', 'Liq'): 9.588e-4,
#                                                 ('osm_coeff', None): 0.9068,
#                                                 ('pressure_osm', None): 2.588e6,
#                                                 ('enth_mass_phase', 'Liq'): 9.9765e4,
#                                                 ('pressure_sat', None): 3111,
#                                                 ('cp_phase', 'Liq'): 4001,
#                                                 ('therm_cond_phase', 'Liq'): 0.6086,
#                                                 ('dh_vap', None): 2.356e6}
#
#
# @pytest.mark.component
# class TestSeawaterPropertySolution(PropertyRegressionTest):
#     def configure(self):
#         self.prop_pack = props.SeawaterParameterBlock
#         self.param_args = {}
#
#         self.solver = 'ipopt'
#         self.optarg = {'nlp_scaling_method': 'user-scaling'}
#
#         self.set_default_scaling_dict = {('flow_mass_phase_comp', ('Liq', 'H2O')): 1,
#                                          ('flow_mass_phase_comp', ('Liq', 'TDS')): 1e2}
#         self.set_state_dict = {('flow_mass_phase_comp', ('Liq', 'H2O')): 0.965,
#                                ('flow_mass_phase_comp', ('Liq', 'TDS')): 0.035,
#                                ('temperature', None): 298.15,
#                                ('pressure', None): 101325}
#         self.results_dict = {('mass_frac_phase_comp', ('Liq', 'H2O')): 0.965,
#                              ('mass_frac_phase_comp', ('Liq', 'TDS')): 0.035,
#                              ('dens_mass_phase', 'Liq'): 1023.5,
#                              ('dens_mass_solvent', None): 996.9,
#                              ('flow_vol_phase', 'Liq'): 9.770e-4,
#                              ('conc_mass_phase_comp', ('Liq', 'H2O')): 987.7,
#                              ('conc_mass_phase_comp', ('Liq', 'TDS')): 35.82,
#                              ('flow_mol_phase_comp', ('Liq', 'H2O')): 53.57,
#                              ('flow_mol_phase_comp', ('Liq', 'TDS')): 1.115,
#                              ('mole_frac_phase_comp', ('Liq', 'H2O')): 0.9796,
#                              ('mole_frac_phase_comp', ('Liq', 'TDS')): 2.038e-2,
#                              ('molality_comp', 'TDS'): 1.155,
#                              ('visc_d_phase', 'Liq'): 9.588e-4,
#                              ('osm_coeff', None): 0.9068,
#                              ('pressure_osm', None): 2.588e6,
#                              ('enth_mass_phase', 'Liq'): 9.9765e4,
#                              ('pressure_sat', None): 3111,
#                              ('cp_phase', 'Liq'): 4001,
#                              ('therm_cond_phase', 'Liq'): 0.6086,
#                              ('dh_vap', None): 2.356e6}
#
#
# @pytest.mark.component
# class TestSeawaterPropertySolution_1(PropertyRegressionTest):
#     def configure(self):
#         self.prop_pack = props.SeawaterParameterBlock
#         self.param_args = {}
#
#         self.solver = 'ipopt'
#         self.optarg = {'nlp_scaling_method': 'user-scaling'}
#
#         self.set_default_scaling_dict = {('flow_mass_phase_comp', ('Liq', 'H2O')): 1,
#                                          ('flow_mass_phase_comp', ('Liq', 'TDS')): 1e2}
#         self.set_state_dict = {('flow_mass_phase_comp', ('Liq', 'H2O')): 0.95,
#                                ('flow_mass_phase_comp', ('Liq', 'TDS')): 0.05,
#                                ('temperature', None): 273.15 + 50,
#                                ('pressure', None): 5e5}
#         self.results_dict = {('mass_frac_phase_comp', ('Liq', 'H2O')): 0.95,
#                              ('mass_frac_phase_comp', ('Liq', 'TDS')): 0.05,
#                              ('dens_mass_phase', 'Liq'): 1025.0,
#                              ('dens_mass_solvent', None): 988.0,
#                              ('flow_vol_phase', 'Liq'): 9.755e-4,
#                              ('conc_mass_phase_comp', ('Liq', 'H2O')): 973.8,
#                              ('conc_mass_phase_comp', ('Liq', 'TDS')): 51.25,
#                              ('flow_mol_phase_comp', ('Liq', 'H2O')): 52.73,
#                              ('flow_mol_phase_comp', ('Liq', 'TDS')): 1.592,
#                              ('mole_frac_phase_comp', ('Liq', 'H2O')): 0.9707,
#                              ('mole_frac_phase_comp', ('Liq', 'TDS')): 2.931e-2,
#                              ('molality_comp', 'TDS'): 1.676,
#                              ('visc_d_phase', 'Liq'): 6.169e-4,
#                              ('osm_coeff', None): 0.9172,
#                              ('pressure_osm', None): 4.081e6,
#                              ('enth_mass_phase', 'Liq'): 1.959e5,
#                              ('pressure_sat', None): 1.201e4,
#                              ('cp_phase', 'Liq'): 3941,
#                              ('therm_cond_phase', 'Liq'): 0.6381,
#                              ('dh_vap', None): 2.263e6}
