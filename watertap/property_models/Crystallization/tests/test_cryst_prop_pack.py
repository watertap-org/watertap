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
import pytest
import cryst_prop_pack as props
from idaes.generic_models.properties.tests.test_harness import \
    PropertyTestHarness as PropertyTestHarness_idaes
from watertap.property_models.tests.property_test_harness import \
    (PropertyTestHarness, PropertyRegressionTest, PropertyCalculateStateTest)

# -----------------------------------------------------------------------------
class TestNaClProperty_idaes(PropertyTestHarness_idaes):
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}
        self.prop_args = {}
        self.has_density_terms = False


eps = 1e-8 # lower bound on flows in property package

# class TestSeawaterProperty(PropertyTestHarness):
#     def configure(self):
#         self.prop_pack = props.SeawaterParameterBlock
#         self.param_args = {}
#         self.scaling_args = {('flow_mass_phase_comp', ('Liq', 'H2O')): 1,
#                              ('flow_mass_phase_comp', ('Liq', 'TDS')): 1e2}
#         self.stateblock_statistics = {'number_variables': 25,
#                                       'number_total_constraints': 21,
#                                       'number_unused_variables': 1,  # pressure is unused
#                                       'default_degrees_of_freedom': 3}  # 4 state vars, but pressure is not active
#         self.default_solution = {('mass_frac_phase_comp', ('Liq', 'H2O')): 0.965,
#                                  ('mass_frac_phase_comp', ('Liq', 'TDS')): 0.035,
#                                  ('dens_mass_phase', 'Liq'): 1023.5,
#                                  ('dens_mass_solvent', None): 996.9,
#                                  ('flow_vol_phase', 'Liq'): 9.770e-4,
#                                  ('conc_mass_phase_comp', ('Liq', 'H2O')): 987.7,
#                                  ('conc_mass_phase_comp', ('Liq', 'TDS')): 35.82,
#                                  ('flow_mol_phase_comp', ('Liq', 'H2O')): 53.57,
#                                  ('flow_mol_phase_comp', ('Liq', 'TDS')): 1.1145,
#                                  ('mole_frac_phase_comp', ('Liq', 'H2O')): 0.9796,
#                                  ('mole_frac_phase_comp', ('Liq', 'TDS')): 2.038e-2,
#                                  ('molality_comp', 'TDS'): 1.155,
#                                  ('visc_d_phase', 'Liq'): 9.588e-4,
#                                  ('osm_coeff', None): 0.9068,
#                                  ('pressure_osm', None): 2.588e6,
#                                  ('enth_mass_phase', 'Liq'): 9.9765e4,
#                                  ('pressure_sat', None): 3111,
#                                  ('cp_phase', 'Liq'): 4001,
#                                  ('therm_cond_phase', 'Liq'): 0.6086,
#                                  ('dh_vap', None): 2.356e6,
#                                  ('diffus_phase', 'Liq'): 1.471e-9}


# @pytest.mark.component
# class TestSeawaterPropertySolution_1(PropertyRegressionTest):
#     def configure(self):
#         self.prop_pack = props.SeawaterParameterBlock
#         self.param_args = {}

#         self.solver = 'ipopt'
#         self.optarg = {'nlp_scaling_method': 'user-scaling'}

#         self.scaling_args = {('flow_mass_phase_comp', ('Liq', 'H2O')): 1,
#                              ('flow_mass_phase_comp', ('Liq', 'TDS')): 1e2}
#         self.state_args = {('flow_mass_phase_comp', ('Liq', 'H2O')): 0.99,
#                            ('flow_mass_phase_comp', ('Liq', 'TDS')): 0.01,
#                            ('temperature', None): 273.15 + 50,
#                            ('pressure', None): 2e5}
#         self.regression_solution = {('mass_frac_phase_comp', ('Liq', 'H2O')): 0.99,
#                                     ('mass_frac_phase_comp', ('Liq', 'TDS')): 0.01,
#                                     ('dens_mass_phase', 'Liq'): 995.4,
#                                     ('dens_mass_solvent', None): 988.0,
#                                     ('flow_vol_phase', 'Liq'): 1.005e-3,
#                                     ('conc_mass_phase_comp', ('Liq', 'H2O')): 985.5,
#                                     ('conc_mass_phase_comp', ('Liq', 'TDS')): 9.954,
#                                     ('flow_mol_phase_comp', ('Liq', 'H2O')): 54.95,
#                                     ('flow_mol_phase_comp', ('Liq', 'TDS')): 0.3184,
#                                     ('mole_frac_phase_comp', ('Liq', 'H2O')): 0.9942,
#                                     ('mole_frac_phase_comp', ('Liq', 'TDS')): 5.761e-3,
#                                     ('molality_comp', 'TDS'): 0.3216,
#                                     ('visc_d_phase', 'Liq'): 5.596e-4,
#                                     ('osm_coeff', None): 0.9029,
#                                     ('pressure_osm', None): 7.710e5,
#                                     ('enth_mass_phase', 'Liq'): 2.066e5,
#                                     ('pressure_sat', None): 1.229e4,
#                                     ('cp_phase', 'Liq'): 4.130e3,
#                                     ('therm_cond_phase', 'Liq'): 0.6400,
#                                     ('dh_vap', None): 2.358e6,
#                                     ('diffus_phase', 'Liq'): 1.493e-9}


# @pytest.mark.component
# class TestSeawaterPropertySolution_2(PropertyRegressionTest):
#     def configure(self):
#         self.prop_pack = props.SeawaterParameterBlock
#         self.param_args = {}

#         self.solver = 'ipopt'
#         self.optarg = {'nlp_scaling_method': 'user-scaling'}

#         self.scaling_args = {('flow_mass_phase_comp', ('Liq', 'H2O')): 1,
#                              ('flow_mass_phase_comp', ('Liq', 'TDS')): 1e2}
#         self.state_args = {('flow_mass_phase_comp', ('Liq', 'H2O')): 0.95,
#                            ('flow_mass_phase_comp', ('Liq', 'TDS')): 0.05,
#                            ('temperature', None): 273.15 + 10,
#                            ('pressure', None): 100e5}
#         self.regression_solution = {('mass_frac_phase_comp', ('Liq', 'H2O')): 0.95,
#                                     ('mass_frac_phase_comp', ('Liq', 'TDS')): 0.05,
#                                     ('dens_mass_phase', 'Liq'): 1.039e3,
#                                     ('dens_mass_solvent', None): 999.5,
#                                     ('flow_vol_phase', 'Liq'): 9.628e-4,
#                                     ('conc_mass_phase_comp', ('Liq', 'H2O')): 986.8,
#                                     ('conc_mass_phase_comp', ('Liq', 'TDS')): 51.93,
#                                     ('flow_mol_phase_comp', ('Liq', 'H2O')): 52.73,
#                                     ('flow_mol_phase_comp', ('Liq', 'TDS')): 1.592,
#                                     ('mole_frac_phase_comp', ('Liq', 'H2O')): 0.9707,
#                                     ('mole_frac_phase_comp', ('Liq', 'TDS')): 2.931e-2,
#                                     ('molality_comp', 'TDS'): 1.676,
#                                     ('visc_d_phase', 'Liq'): 1.443e-3,
#                                     ('osm_coeff', None): 0.9106,
#                                     ('pressure_osm', None): 3.591e6,
#                                     ('enth_mass_phase', 'Liq'): 3.897e4,
#                                     ('pressure_sat', None): 1.194e3,
#                                     ('cp_phase', 'Liq'): 3.916e3,
#                                     ('therm_cond_phase', 'Liq'): 0.5854,
#                                     ('dh_vap', None): 2.353e6,
#                                     ('diffus_phase', 'Liq'): 1.471e-9}


@pytest.mark.component
class TestNaClCalculateState_1(PropertyCalculateStateTest):
    # Test pure liquid solution with mass fractions
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.solver = 'ipopt'
        self.optarg = {'nlp_scaling_method': 'user-scaling'}

        self.scaling_args = {('flow_mass_phase_comp', ('Liq', 'H2O')): 1e-1,
                             ('flow_mass_phase_comp', ('Liq', 'NaCl')): 1e1,
                             }

        self.var_args = {('flow_vol_phase', 'Liq'): 2e-2,
                         ('flow_vol_phase', 'Sol'): 0,
                         ('flow_vol_phase', 'Vap'): 0,
                         ('mass_frac_phase_comp', ('Liq', 'NaCl')): 0.05,
                         ('mass_frac_phase_comp', ('Sol', 'NaCl')): 0,
                         ('mass_frac_phase_comp', ('Vap', 'H2O')): 0,
                         ('temperature', None): 273.15 + 25,
                         ('pressure', None): 5e5
                         }

        self.state_solution = {('flow_mass_phase_comp', ('Liq', 'H2O')): 19.6,
                               ('flow_mass_phase_comp', ('Liq', 'NaCl')): 1.032,
                               ('flow_mass_phase_comp', ('Sol', 'NaCl')): 0,
                               ('flow_mass_phase_comp', ('Sol', 'H2O')): 0,
                               ('flow_mass_phase_comp', ('Vap', 'NaCl')): 0,
                               ('flow_mass_phase_comp', ('Vap', 'H2O')): 0
                               }


@pytest.mark.component
class TestNaClCalculateState_2(PropertyCalculateStateTest):
    # Test pure liquid with mole fractions
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.solver = 'ipopt'
        self.optarg = {'nlp_scaling_method': 'user-scaling'}

        self.scaling_args = {('flow_mass_phase_comp', ('Liq', 'H2O')): 1e1,
                             ('flow_mass_phase_comp', ('Liq', 'NaCl')): 1e2,
                             # The rest are expected to be zero
                             ('flow_mass_phase_comp', ('Sol', 'NaCl')): 1/eps,
                             ('flow_mass_phase_comp', ('Sol', 'H2O')): 1/eps,
                             ('flow_mass_phase_comp', ('Vap', 'NaCl')): 1/eps,
                             ('flow_mass_phase_comp', ('Vap', 'H2O')): 1/eps,
                             }
        self.var_args = {('flow_vol_phase', 'Liq'): 2e-4,
                         ('flow_vol_phase', 'Sol'): 0,
                         ('flow_vol_phase', 'Vap'): 0,
                         ('mole_frac_phase_comp', ('Liq', 'NaCl')): 0.05,
                         ('mole_frac_phase_comp', ('Sol', 'H2O')): 0,
                         ('mole_frac_phase_comp', ('Vap', 'NaCl')): 0,
                         ('temperature', None): 273.15 + 25,
                         ('pressure', None): 5e5}
        self.state_solution = {('flow_mass_phase_comp', ('Liq', 'H2O')): 18.84e-2,
                               ('flow_mass_phase_comp', ('Liq', 'NaCl')): 3.215e-2,
                               ('flow_mass_phase_comp', ('Sol', 'NaCl')): 0,
                               ('flow_mass_phase_comp', ('Sol', 'H2O')): 0,
                               ('flow_mass_phase_comp', ('Vap', 'NaCl')): 0,
                               ('flow_mass_phase_comp', ('Vap', 'H2O')): 0
                               }


@pytest.mark.component
class TestNaClCalculateState_3(PropertyCalculateStateTest):
    # Test pure liquid solution with pressure_sat defined instead of temperature
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.solver = 'ipopt'
        self.optarg = {'nlp_scaling_method': 'user-scaling'}

        self.scaling_args = {('flow_mass_phase_comp', ('Liq', 'H2O')): 1e-1,
                             ('flow_mass_phase_comp', ('Liq', 'NaCl')): 1e1,
                             }

        self.var_args = {('flow_vol_phase', 'Liq'): 2e-2,
                         ('mass_frac_phase_comp', ('Liq', 'NaCl')): 0.05,
                         ('flow_vol_phase', 'Sol'): 0,
                         ('mass_frac_phase_comp', ('Sol', 'NaCl')): 0,
                         ('flow_vol_phase', 'Vap'): 0,
                         ('mass_frac_phase_comp', ('Vap', 'H2O')): 0,
                         ('pressure_sat', None):  2905,
                         ('pressure', None): 5e5
                         }

        self.state_solution = {('flow_mass_phase_comp', ('Liq', 'H2O')): 19.6,
                               ('flow_mass_phase_comp', ('Liq', 'NaCl')): 1.032,
                               ('flow_mass_phase_comp', ('Sol', 'NaCl')): 0,
                               ('flow_mass_phase_comp', ('Sol', 'H2O')): 0,
                               ('flow_mass_phase_comp', ('Vap', 'NaCl')): 0,
                               ('flow_mass_phase_comp', ('Vap', 'H2O')): 0,
                               ('temperature', None): 273.15 + 25
                               }


@pytest.mark.component
class TestNaClCalculateState_4(PropertyCalculateStateTest):
    # Test pure solid solution with mass fractions
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.solver = 'ipopt'
        self.optarg = {'nlp_scaling_method': 'user-scaling'}

        self.scaling_args = {('flow_mass_phase_comp', ('Liq', 'H2O')): 1/eps,
                             ('flow_mass_phase_comp', ('Sol', 'NaCl')): 1e-2,
                             ('flow_mass_phase_comp', ('Sol', 'H2O')): 1/eps,
                             ('flow_mass_phase_comp', ('Vap', 'H2O')): 1/eps
                             }

        self.var_args = {('flow_vol_phase', 'Liq'): 0,
                         ('mass_frac_phase_comp', ('Liq', 'NaCl')): 0,
                         ('flow_vol_phase', 'Sol'): 2e-2,
                         ('mass_frac_phase_comp', ('Sol', 'H2O')): 0,
                         ('flow_vol_phase', 'Vap'): 0,
                         ('mass_frac_phase_comp', ('Vap', 'NaCl')): 0,
                         ('temperature', None): 273.15 + 25,
                         ('pressure', None): 5e5
                         }

        self.state_solution = {('flow_mass_phase_comp', ('Sol', 'NaCl')): 2115 * 2e-2, # solid density is constant
                               ('flow_mass_phase_comp', ('Sol', 'H2O')):  0, 
                               ('flow_mass_phase_comp', ('Liq', 'H2O')): 9.968e-09,
                               ('flow_mass_phase_comp', ('Liq', 'NaCl')): 0,
                               ('flow_mass_phase_comp', ('Vap', 'NaCl')): 0,
                               ('flow_mass_phase_comp', ('Vap', 'H2O')):  2.708e-09
                               }


@pytest.mark.component
class TestNaClCalculateState_5(PropertyCalculateStateTest):
    # Test solid-liquid-vapor mixture solution with mass fractions
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.solver = 'ipopt'
        self.optarg = {'nlp_scaling_method': 'user-scaling'}

        self.scaling_args = {('flow_mass_phase_comp', ('Liq', 'H2O')): 1e-1,
                             ('flow_mass_phase_comp', ('Liq', 'NaCl')): 1e1,
                             ('flow_mass_phase_comp', ('Sol', 'NaCl')): 1e-2,
                             ('flow_mass_phase_comp', ('Sol', 'H2O')): 1/eps}

        self.var_args = {('flow_vol_phase', 'Liq'): 2e-2,
                         ('mass_frac_phase_comp', ('Liq', 'NaCl')): 0.05,
                         ('flow_vol_phase', 'Sol'): 2e-2,
                         ('mass_frac_phase_comp', ('Sol', 'H2O')): 0,
                         ('flow_vol_phase', 'Vap'): 2e-2,
                         ('mass_frac_phase_comp', ('Vap', 'NaCl')): 0,
                         ('temperature', None): 273.15 + 25,
                         ('pressure', None): 5e5
                         }

        self.state_solution = {('flow_mass_phase_comp', ('Liq', 'H2O')):  19.6,
                               ('flow_mass_phase_comp', ('Liq', 'NaCl')): 1.032,
                               ('flow_mass_phase_comp', ('Vap', 'NaCl')): 0,
                               ('flow_mass_phase_comp', ('Vap', 'H2O')):  0.07265, 
                               ('flow_mass_phase_comp', ('Sol', 'NaCl')): 42.3, 
                               ('flow_mass_phase_comp', ('Sol', 'H2O')):  0, 
                               }


@pytest.mark.component
class TestNaClCalculateState_6(PropertyCalculateStateTest):
    # Test liquid-solid-vapor mixture with mole fractions
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}

        self.solver = 'ipopt'
        self.optarg = {'tol': 1e-8, 'nlp_scaling_method': 'user-scaling'}

        self.scaling_args = {('flow_mass_phase_comp', ('Liq', 'H2O')): 1e1,
                             ('flow_mass_phase_comp', ('Liq', 'NaCl')): 1e2,
                             ('flow_mass_phase_comp', ('Sol', 'NaCl')): 1e2,
                             ('flow_mass_phase_comp', ('Vap', 'H2O')): 1e3,
                             ('flow_mass_phase_comp', ('Sol', 'H2O')): 1/eps, # expected to be zero
                             ('flow_mass_phase_comp', ('Vap', 'NaCl')): 1/eps, # expected to be zero
                             }
        self.var_args = {('flow_vol_phase', 'Liq'): 2e-4,
                         ('flow_vol_phase', 'Sol'): 2e-4,
                         ('flow_vol_phase', 'Vap'): 2e-4,
                         ('mole_frac_phase_comp', ('Liq', 'NaCl')): 0.05,
                         ('mole_frac_phase_comp', ('Sol', 'H2O')): 0,
                         ('mole_frac_phase_comp', ('Vap', 'NaCl')): 0,
                         ('temperature', None): 273.15 + 25,
                         ('pressure', None): 5e5}
        self.state_solution = {('flow_mass_phase_comp', ('Liq', 'H2O')): 18.84e-2,
                               ('flow_mass_phase_comp', ('Liq', 'NaCl')): 3.215e-2,
                               ('flow_mass_phase_comp', ('Sol', 'NaCl')): 0.423,
                               ('flow_mass_phase_comp', ('Sol', 'H2O')): 0,
                               ('flow_mass_phase_comp', ('Vap', 'NaCl')): 0,
                               ('flow_mass_phase_comp', ('Vap', 'H2O')): 3.632 * 2e-4 # Density from ideal gas law * vol. flow
                               }
