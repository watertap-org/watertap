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
from proteuslib.property_models.tests.property_test_harness import (PropertyUnitTestHarness,
                                                                    PropertyComponentTestHarness,
                                                                    PropertyRegressionTest)
from idaes.generic_models.properties.tests.test_harness import \
    PropertyTestHarness as PropertyTestHarness_idaes
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
from idaes.core.util.scaling import calculate_scaling_factors
from idaes.core.util import get_solver

# -----------------------------------------------------------------------------
@pytest.mark.unit
class TestSeawaterProperty(PropertyTestHarness_idaes):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}
        self.prop_args = {}
        self.has_density_terms = False


@pytest.mark.unit
class TestSeawaterPropertyUnit(PropertyUnitTestHarness):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}
        self.set_default_scaling_dict = {('flow_mass_phase_comp', ('Liq', 'H2O')): 1,
                                         ('flow_mass_phase_comp', ('Liq', 'TDS')): 1e2}

        self.component_dict = {'H2O': Solvent,
                               'TDS': Solute}
        self.phase_dict = {'Liq': LiquidPhase}

        dens_units = pyunits.kg/pyunits.m**3
        t_inv_units = pyunits.K**-1
        visc_d_units = pyunits.Pa * pyunits.s
        s_inv_units = pyunits.kg/pyunits.g
        enth_mass_units = pyunits.J/pyunits.kg
        cp_units = pyunits.J/(pyunits.kg * pyunits.K)
        self.param_obj_dict = {'mw_comp': (Param, pyunits.kg/pyunits.mol),
                               'dens_mass_param_A1': (Var, dens_units),
                               'dens_mass_param_A2': (Var, dens_units * t_inv_units),
                               'dens_mass_param_A3': (Var, dens_units * t_inv_units**2),
                               'dens_mass_param_A4': (Var, dens_units * t_inv_units**3),
                               'dens_mass_param_A5': (Var, dens_units * t_inv_units**4),
                               'dens_mass_param_B1': (Var, dens_units),
                               'dens_mass_param_B2': (Var, dens_units * t_inv_units),
                               'dens_mass_param_B3': (Var, dens_units * t_inv_units**2),
                               'dens_mass_param_B4': (Var, dens_units * t_inv_units**3),
                               'dens_mass_param_B5': (Var, dens_units * t_inv_units**2),
                               'visc_d_param_muw_A': (Var, visc_d_units),
                               'visc_d_param_muw_B': (Var, visc_d_units**-1 * t_inv_units**2),
                               'visc_d_param_muw_C': (Var, pyunits.K),
                               'visc_d_param_muw_D': (Var, visc_d_units**-1),
                               'visc_d_param_A_1': (Var, pyunits.dimensionless),
                               'visc_d_param_A_2': (Var, t_inv_units),
                               'visc_d_param_A_3': (Var, t_inv_units**2),
                               'visc_d_param_B_1': (Var, pyunits.dimensionless),
                               'visc_d_param_B_2': (Var, t_inv_units),
                               'visc_d_param_B_3': (Var, t_inv_units**2),
                               'osm_coeff_param_1': (Var, pyunits.dimensionless),
                               'osm_coeff_param_2': (Var, t_inv_units),
                               'osm_coeff_param_3': (Var, t_inv_units**2),
                               'osm_coeff_param_4': (Var, t_inv_units**4),
                               'osm_coeff_param_5': (Var, pyunits.dimensionless),
                               'osm_coeff_param_6': (Var, t_inv_units),
                               'osm_coeff_param_7': (Var, t_inv_units**3),
                               'osm_coeff_param_8': (Var, pyunits.dimensionless),
                               'osm_coeff_param_9': (Var, t_inv_units),
                               'osm_coeff_param_10': (Var, t_inv_units**2),
                               'enth_mass_param_A1': (Var, enth_mass_units),
                               'enth_mass_param_A2': (Var, enth_mass_units * t_inv_units),
                               'enth_mass_param_A3': (Var, enth_mass_units * t_inv_units**2),
                               'enth_mass_param_A4': (Var, enth_mass_units * t_inv_units**3),
                               'enth_mass_param_B1': (Var, enth_mass_units),
                               'enth_mass_param_B2': (Var, enth_mass_units),
                               'enth_mass_param_B3': (Var, enth_mass_units),
                               'enth_mass_param_B4': (Var, enth_mass_units),
                               'enth_mass_param_B5': (Var, enth_mass_units * t_inv_units),
                               'enth_mass_param_B6': (Var, enth_mass_units * t_inv_units**2),
                               'enth_mass_param_B7': (Var, enth_mass_units * t_inv_units**3),
                               'enth_mass_param_B8': (Var, enth_mass_units * t_inv_units),
                               'enth_mass_param_B9': (Var, enth_mass_units * t_inv_units),
                               'enth_mass_param_B10': (Var, enth_mass_units * t_inv_units**2),
                               'pressure_sat_param_psatw_A1': (Var, pyunits.K),
                               'pressure_sat_param_psatw_A2': (Var, pyunits.dimensionless),
                               'pressure_sat_param_psatw_A3': (Var, t_inv_units),
                               'pressure_sat_param_psatw_A4': (Var, t_inv_units**2),
                               'pressure_sat_param_psatw_A5': (Var, t_inv_units**3),
                               'pressure_sat_param_psatw_A6': (Var, pyunits.dimensionless),
                               'pressure_sat_param_B1': (Var, s_inv_units),
                               'pressure_sat_param_B2': (Var, s_inv_units**2),
                               'cp_phase_param_A1': (Var, cp_units),
                               'cp_phase_param_A2': (Var, cp_units * s_inv_units),
                               'cp_phase_param_A3': (Var, cp_units * s_inv_units**2),
                               'cp_phase_param_B1': (Var, cp_units * t_inv_units),
                               'cp_phase_param_B2': (Var, cp_units * s_inv_units * t_inv_units),
                               'cp_phase_param_B3': (Var, cp_units * s_inv_units**2 * t_inv_units),
                               'cp_phase_param_C1': (Var, cp_units * t_inv_units**2),
                               'cp_phase_param_C2': (Var, cp_units * s_inv_units * t_inv_units**2),
                               'cp_phase_param_C3': (Var, cp_units * s_inv_units**2 * t_inv_units**2),
                               'cp_phase_param_D1': (Var, cp_units * t_inv_units**3),
                               'cp_phase_param_D2': (Var, cp_units * s_inv_units * t_inv_units**3),
                               'cp_phase_param_D3': (Var, cp_units * s_inv_units**2 * t_inv_units**3),
                               'therm_cond_phase_param_1': (Var, pyunits.dimensionless),
                               'therm_cond_phase_param_2': (Var, s_inv_units),
                               'therm_cond_phase_param_3': (Var, pyunits.dimensionless),
                               'therm_cond_phase_param_4': (Var, pyunits.dimensionless),
                               'therm_cond_phase_param_5': (Var, pyunits.K),
                               'therm_cond_phase_param_6': (Var, s_inv_units * pyunits.K),
                               'therm_cond_phase_param_7': (Var, pyunits.K),
                               'therm_cond_phase_param_8': (Var, s_inv_units * pyunits.K),
                               'dh_vap_w_param_0': (Var, enth_mass_units),
                               'dh_vap_w_param_1': (Var, cp_units),
                               'dh_vap_w_param_2': (Var, enth_mass_units * t_inv_units**2),
                               'dh_vap_w_param_3': (Var, enth_mass_units * t_inv_units**3),
                               'dh_vap_w_param_4': (Var, enth_mass_units * t_inv_units**4)}
        self.default_scaling_factor_dict = {('temperature', None): 1e-2,
                                            ('pressure', None): 1e-6,
                                            ('dens_mass_phase', 'Liq'): 1e-3,
                                            ('dens_mass_solvent', None): 1e-3,
                                            ('visc_d_phase', 'Liq'): 1e3,
                                            ('osm_coeff', None): 1e0,
                                            ('enth_mass_phase', 'Liq'): 1e-5,
                                            ('pressure_sat', None): 1e-5,
                                            ('cp_phase', 'Liq'): 1e-3,
                                            ('therm_cond_phase', 'Liq'): 1e0,
                                            ('dh_vap', None): 1e-6,
                                            ('flow_mass_phase_comp', ('Liq', 'H2O')): 1,
                                            ('flow_mass_phase_comp', ('Liq', 'TDS')): 1e2}

        self.state_vars_list = ['flow_mass_phase_comp', 'temperature', 'pressure']
        self.stateblock_obj_dict = {'flow_mass_phase_comp': (Var, pyunits.kg/pyunits.s),
                                    'temperature': (Var, pyunits.K),
                                    'pressure': (Var, pyunits.Pa),
                                    'mass_frac_phase_comp': (Var, pyunits.dimensionless),
                                    'dens_mass_phase': (Var, pyunits.kg * pyunits.m**-3),
                                    'dens_mass_solvent': (Var, pyunits.kg * pyunits.m**-3),
                                    'flow_vol_phase': (Var, pyunits.m**3 / pyunits.s),
                                    'flow_vol': (Expression, pyunits.m**3 / pyunits.s),
                                    'conc_mass_phase_comp': (Var, pyunits.kg * pyunits.m**-3),
                                    'flow_mol_phase_comp': (Var, pyunits.mol / pyunits.s),
                                    'mole_frac_phase_comp': (Var, pyunits.dimensionless),
                                    'molality_comp': (Var, pyunits.mole / pyunits.kg),
                                    'visc_d_phase': (Var, pyunits.Pa * pyunits.s),
                                    'osm_coeff': (Var, pyunits.dimensionless),
                                    'pressure_osm': (Var, pyunits.Pa),
                                    'enth_mass_phase': (Var, enth_mass_units),
                                    'enth_flow': (Expression, pyunits.J / pyunits.s),
                                    'pressure_sat': (Var, pyunits.Pa),
                                    'cp_phase': (Var, cp_units),
                                    'therm_cond_phase': (Var, pyunits.W/pyunits.m/pyunits.K),
                                    'dh_vap': (Var, enth_mass_units)
                                    }

        self.statistics_dict = {'number_variables': 101,
                                'number_total_constraints': 20,
                                'number_unused_variables': 1,  # pressure is unused
                                'default_degrees_of_freedom': 3}  # 4 state vars, but pressure is not active
        self.general_methods_dict = {
            ('get_material_flow_terms', ('Liq', 'H2O')): ('flow_mass_phase_comp', ('Liq', 'H2O')),
            ('get_material_flow_terms', ('Liq', 'TDS')): ('flow_mass_phase_comp', ('Liq', 'TDS')),
            ('get_enthalpy_flow_terms', None): ('enth_flow', None)}
        self.default_methods_dict = {'default_material_balance_type': MaterialBalanceType.componentTotal,
                                     'default_energy_balance_type': EnergyBalanceType.enthalpyTotal,
                                     'get_material_flow_basis': MaterialFlowBasis.mass}


@pytest.mark.component
class TestSeawaterPropertyComponent(PropertyComponentTestHarness):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}
        self.set_default_scaling_dict = {('flow_mass_phase_comp', ('Liq', 'H2O')): 1,
                                         ('flow_mass_phase_comp', ('Liq', 'TDS')): 1e2}

        self.solver = 'ipopt'
        self.optarg = {'nlp_scaling_method': 'user-scaling'}

        self.default_state_values_dict = {('flow_mass_phase_comp', ('Liq', 'H2O')): 0.965,
                                          ('flow_mass_phase_comp', ('Liq', 'TDS')): 0.035,
                                          ('temperature', None): 298.15,
                                          ('pressure', None): 101325}
        self.default_initialize_results_dict = {('mass_frac_phase_comp', ('Liq', 'H2O')): 0.965,
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


@pytest.mark.component
class TestSeawaterPropertySolution(PropertyRegressionTest):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}

        self.solver = 'ipopt'
        self.optarg = {'nlp_scaling_method': 'user-scaling'}

        self.set_default_scaling_dict = {('flow_mass_phase_comp', ('Liq', 'H2O')): 1,
                                         ('flow_mass_phase_comp', ('Liq', 'TDS')): 1e2}
        self.set_state_dict = {('flow_mass_phase_comp', ('Liq', 'H2O')): 0.965,
                               ('flow_mass_phase_comp', ('Liq', 'TDS')): 0.035,
                               ('temperature', None): 298.15,
                               ('pressure', None): 101325}
        self.results_dict = {('mass_frac_phase_comp', ('Liq', 'H2O')): 0.965,
                             ('mass_frac_phase_comp', ('Liq', 'TDS')): 0.035,
                             ('dens_mass_phase', 'Liq'): 1023.5,
                             ('dens_mass_solvent', None): 996.9,
                             ('flow_vol_phase', 'Liq'): 9.770e-4,
                             ('conc_mass_phase_comp', ('Liq', 'H2O')): 987.7,
                             ('conc_mass_phase_comp', ('Liq', 'TDS')): 35.82,
                             ('flow_mol_phase_comp', ('Liq', 'H2O')): 53.57,
                             ('flow_mol_phase_comp', ('Liq', 'TDS')): 1.115,
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


@pytest.mark.component
class TestSeawaterPropertySolution_1(PropertyRegressionTest):
    def configure(self):
        self.prop_pack = props.SeawaterParameterBlock
        self.param_args = {}

        self.solver = 'ipopt'
        self.optarg = {'nlp_scaling_method': 'user-scaling'}

        self.set_default_scaling_dict = {('flow_mass_phase_comp', ('Liq', 'H2O')): 1,
                                         ('flow_mass_phase_comp', ('Liq', 'TDS')): 1e2}
        self.set_state_dict = {('flow_mass_phase_comp', ('Liq', 'H2O')): 0.95,
                               ('flow_mass_phase_comp', ('Liq', 'TDS')): 0.05,
                               ('temperature', None): 273.15 + 50,
                               ('pressure', None): 5e5}
        self.results_dict = {('mass_frac_phase_comp', ('Liq', 'H2O')): 0.95,
                             ('mass_frac_phase_comp', ('Liq', 'TDS')): 0.05,
                             ('dens_mass_phase', 'Liq'): 1025.0,
                             ('dens_mass_solvent', None): 988.0,
                             ('flow_vol_phase', 'Liq'): 9.755e-4,
                             ('conc_mass_phase_comp', ('Liq', 'H2O')): 973.8,
                             ('conc_mass_phase_comp', ('Liq', 'TDS')): 51.25,
                             ('flow_mol_phase_comp', ('Liq', 'H2O')): 52.73,
                             ('flow_mol_phase_comp', ('Liq', 'TDS')): 1.592,
                             ('mole_frac_phase_comp', ('Liq', 'H2O')): 0.9707,
                             ('mole_frac_phase_comp', ('Liq', 'TDS')): 2.931e-2,
                             ('molality_comp', 'TDS'): 1.676,
                             ('visc_d_phase', 'Liq'): 6.169e-4,
                             ('osm_coeff', None): 0.9172,
                             ('pressure_osm', None): 4.081e6,
                             ('enth_mass_phase', 'Liq'): 1.959e5,
                             ('pressure_sat', None): 1.201e4,
                             ('cp_phase', 'Liq'): 3941,
                             ('therm_cond_phase', 'Liq'): 0.6381,
                             ('dh_vap', None): 2.263e6}

#
# @pytest.mark.component
# def test_solution_seawater_alternative_condition_1():
#     m = ConcreteModel()
#     m.fs = FlowsheetBlock(default={"dynamic": False})
#     m.fs.properties = props.SeawaterParameterBlock()
#     m.fs.stream = m.fs.properties.build_state_block([0], default={})
#
#     # specify stream
#     mass_flow = 1
#     x_TDS = 0.05
#     m.fs.stream[0].flow_mass_phase_comp['Liq', 'TDS'].fix(x_TDS * mass_flow)
#     m.fs.stream[0].flow_mass_phase_comp['Liq', 'H2O'].fix((1 - x_TDS) * mass_flow)
#     m.fs.stream[0].temperature.fix(273.15 + 50)
#     m.fs.stream[0].pressure.fix(5e5)
#
#     # touch all properties
#     metadata = m.fs.properties.get_metadata().properties
#     for v_str in metadata.keys():
#         getattr(m.fs.stream[0], v_str)
#
#     # scale model
#     m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
#     m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))
#     calculate_scaling_factors(m.fs)
#
#     # specify results
    results_dict = {('mass_frac_phase_comp', ('Liq', 'H2O')): 0.95,
                    ('mass_frac_phase_comp', ('Liq', 'TDS')): 0.05,
                    ('dens_mass_phase', 'Liq'): 1025.0,
                    ('dens_mass_solvent', None): 988.0,
                    ('flow_vol_phase', 'Liq'): 9.755e-4,
                    ('conc_mass_phase_comp', ('Liq', 'H2O')): 973.8,
                    ('conc_mass_phase_comp', ('Liq', 'TDS')): 51.25,
                    ('flow_mol_phase_comp', ('Liq', 'H2O')): 52.73,
                    ('flow_mol_phase_comp', ('Liq', 'TDS')): 1.592,
                    ('mole_frac_phase_comp', ('Liq', 'H2O')): 0.9707,
                    ('mole_frac_phase_comp', ('Liq', 'TDS')): 2.931e-2,
                    ('molality_comp', 'TDS'): 1.676,
                    ('visc_d_phase', 'Liq'): 6.169e-4,
                    ('osm_coeff', None): 0.9172,
                    ('pressure_osm', None): 4.081e6,
                    ('enth_mass_phase', 'Liq'): 1.959e5,
                    ('pressure_sat', None): 1.201e4,
                    ('cp_phase', 'Liq'): 3941,
                    ('therm_cond_phase', 'Liq'): 0.6381,
                    ('dh_vap', None): 2.263e6}
#
#     # check solve with no initialization
#     property_regression_test(m.fs.stream[0], results_dict)
