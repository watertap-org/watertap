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
import proteuslib.property_models.NaCl_prop_pack as props
from idaes.generic_models.properties.tests.test_harness import \
    PropertyTestHarness
from pyomo.environ import (ConcreteModel,
                           Constraint,
                           TerminationCondition,
                           SolverStatus,
                           value,
                           Suffix,
                           Param,
                           Expression,
                           Var,
                           Set,
                           units)
from idaes.core import (FlowsheetBlock,
                        MaterialFlowBasis,
                        MaterialBalanceType,
                        EnergyBalanceType)
from idaes.core.components import Solute, Solvent
from idaes.core.phases import LiquidPhase
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables)
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

# -----------------------------------------------------------------------------
@pytest.mark.unit
class TestBasicMix(PropertyTestHarness):
    def configure(self):
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}
        self.prop_args = {}
        self.has_density_terms = False

class TestNaClPropPack():
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.properties = props.NaClParameterBlock()
        m.fs.stream = m.fs.properties.build_state_block([0], default={})

        # specify conditions
        mass_flow = 1
        x_NaCl = 0.035
        m.fs.stream[0].flow_mass_phase_comp['Liq', 'NaCl'].fix(x_NaCl * mass_flow)
        m.fs.stream[0].flow_mass_phase_comp['Liq', 'H2O'].fix((1 - x_NaCl) * mass_flow)
        m.fs.stream[0].temperature.fix(273.15 + 25)
        m.fs.stream[0].pressure.fix(101325)
        return m

    @pytest.mark.unit
    def test_components_phases(self, frame):
        m = frame

        # test components, package only supports H2O and NaCl
        assert isinstance(m.fs.properties.component_list, Set)
        assert len(m.fs.properties.component_list) == 2
        assert 'H2O' in m.fs.properties.component_list
        assert isinstance(m.fs.properties.H2O, Solvent)
        assert 'NaCl' in m.fs.properties.component_list
        assert isinstance(m.fs.properties.NaCl, Solute)

        # test phase, package only supports Liquid
        assert isinstance(m.fs.properties.phase_list, Set)
        assert len(m.fs.properties.phase_list) == 1
        assert isinstance(m.fs.properties.Liq, LiquidPhase)

    @staticmethod
    def check_block_obj(blk, obj_type_dict):
        # check that all objects are their expected type
        for (obj_str, obj_type) in obj_type_dict.items():
            obj = getattr(blk, obj_str)
            assert isinstance(obj, obj_type)

    @staticmethod
    def check_block_obj_coverage(blk, obj_type_dict,
                         type_tpl= (Param, Var, Expression, Constraint)):
        # check that all added objects are tested
        for obj in blk.component_objects(type_tpl, descend_into=False):
            obj_str = obj.local_name
            assert obj_str in obj_type_dict

    @pytest.mark.unit
    def test_parameters(self, frame):
        m = frame

        # create a dictionary with all pyomo objects and their type on the parameter block
        param_obj_type_dict = {'mw_comp': Param,
                                'dens_mass_param': Var,
                                'visc_d_param': Var,
                                'diffus_param': Var,
                                'osm_coeff_param': Var}

        enth_mass_param_end_list = ['A1', 'A2', 'A3', 'A4', 'B1', 'B2']
        for e in enth_mass_param_end_list:
            param_obj_type_dict['enth_mass_param_' + e] = Var

        self.check_block_obj(m.fs.properties, param_obj_type_dict)
        self.check_block_obj_coverage(m.fs.properties, param_obj_type_dict)

        # test that the parameter variables are fixed
        for var in m.fs.properties.component_objects(Var):
            for ind, v in var.items():
                assert v.is_fixed()

    @pytest.mark.unit
    def test_default_scaling(self, frame):
        m = frame

        assert isinstance(m.fs.properties.default_scaling_factor, dict)
        default_scaling_factor_test = {('temperature', None): 1e-2,
                                       ('pressure', None): 1e-6,
                                       ('dens_mass_phase', 'Liq'): 1e-3,
                                       ('visc_d_phase', 'Liq'): 1e3,
                                       ('diffus_phase', 'Liq'): 1e9,
                                       ('osm_coeff', None): 1e0,
                                       ('enth_mass_phase', 'Liq'): 1e-5}

        assert len(default_scaling_factor_test) == len(m.fs.properties.default_scaling_factor)
        for t, sf in m.fs.properties.default_scaling_factor.items():
            assert t in default_scaling_factor_test.keys()
            assert default_scaling_factor_test[t] == sf

    @pytest.mark.unit
    def test_metadata_exists(self, frame):
        m = frame
        assert hasattr(m.fs.properties, 'define_metadata')

    @pytest.mark.unit
    def test_build(self, frame):
        m = frame

        # create a dictionary with all pyomo objects and their type on the state block
        sb_obj_type_dict = {}

        def add_obj_type_from_list(obj_str_list=None, obj_type=None):
            for obj_str in obj_str_list:
                sb_obj_type_dict[obj_str] = obj_type

        # scaling factor
        sb_obj_type_dict['scaling_factor'] = Suffix

        # state variables
        state_vars_list = ['flow_mass_phase_comp', 'temperature', 'pressure']
        state_vars_dict = m.fs.stream[0].define_state_vars()
        assert len(state_vars_dict) == len(state_vars_list)
        add_obj_type_from_list(obj_str_list=state_vars_list, obj_type=Var)

        # on demand properties
        var_str_list = ['mass_frac_phase_comp', 'dens_mass_phase', 'flow_vol_phase',
                        'conc_mass_phase_comp', 'flow_mol_phase_comp', 'mole_frac_phase_comp',
                        'molality_comp', 'visc_d_phase', 'diffus_phase', 'osm_coeff',
                        'pressure_osm', 'enth_mass_phase']
        add_obj_type_from_list(obj_str_list=var_str_list, obj_type=Var)

        exp_str_list = ['flow_vol', 'enth_flow']
        add_obj_type_from_list(obj_str_list=exp_str_list, obj_type=Expression)

        # test that properties are not built if not demanded
        for obj_str in var_str_list + exp_str_list:
            assert not m.fs.stream[0].is_property_constructed(obj_str)

        # test that properties are built on demand
        for obj_str in var_str_list + exp_str_list:
            assert hasattr(m.fs.stream[0], obj_str)

        # constraints for on demand variables
        con_str_list = []
        for obj_str in var_str_list:
            con_str_list.append('eq_' + obj_str)
        add_obj_type_from_list(obj_str_list=con_str_list, obj_type=Constraint)

        self.check_block_obj(m.fs.stream[0], sb_obj_type_dict)
        self.check_block_obj_coverage(m.fs.stream[0], sb_obj_type_dict)

    @pytest.mark.unit
    def test_statistics(self, frame):
        m = frame

        assert number_variables(m) == 38
        assert number_total_constraints(m) == 16
        assert number_unused_variables(m) == 1  # pressure is unused

    @pytest.mark.unit
    def test_general_methods(self, frame):
        m = frame

        assert hasattr(m.fs.stream[0], 'get_material_flow_terms')
        assert (m.fs.stream[0].get_material_flow_terms('Liq', 'H2O')
                is m.fs.stream[0].flow_mass_phase_comp['Liq', 'H2O'])
        assert (m.fs.stream[0].get_material_flow_terms('Liq', 'NaCl')
                is m.fs.stream[0].flow_mass_phase_comp['Liq', 'NaCl'])

        assert hasattr(m.fs.stream[0], 'get_enthalpy_flow_terms')
        assert (m.fs.stream[0].get_enthalpy_flow_terms(None)
                is m.fs.stream[0].enth_flow)

        assert hasattr(m.fs.stream[0], 'default_material_balance_type')
        assert (m.fs.stream[0].default_material_balance_type()
                is MaterialBalanceType.componentTotal)

        assert hasattr(m.fs.stream[0], 'default_energy_balance_type')
        assert (m.fs.stream[0].default_energy_balance_type()
                is EnergyBalanceType.enthalpyTotal)

        assert hasattr(m.fs.stream[0], 'get_material_flow_basis')
        assert (m.fs.stream[0].get_material_flow_basis()
                is MaterialFlowBasis.mass)

    @pytest.mark.unit
    def test_units(self, frame):
        m = frame

        # # check all variables have their assigned units - parameter block
        var_unit_dict = {
            'dens_mass_param': units.kg / units.m**3,
            'visc_d_param': units.Pa * units.s,
            'diffus_param': units.m**2 / units.s,
            'osm_coeff_param': units.dimensionless,
            'enth_mass_param_A1': units.J / units.kg,
            'enth_mass_param_A2': (units.J / units.kg) * units.K**-1,
            'enth_mass_param_A3': (units.J / units.kg) * units.K**-2,
            'enth_mass_param_A4': (units.J / units.kg) * units.K**-3,
            'enth_mass_param_B1': units.dimensionless,
            'enth_mass_param_B2': units.dimensionless}

        for (v, u) in var_unit_dict.items():
            var = getattr(m.fs.properties, v)
            assert_units_equivalent(var, u)

        # check all variables have their assigned units - state block
        var_unit_dict = {
            'flow_mass_phase_comp': units.kg / units.s,
            'temperature': units.K,
            'pressure': units.Pa,
            'mass_frac_phase_comp': units.dimensionless,
            'dens_mass_phase': units.kg * units.m**-3,
            'flow_vol_phase': units.m ** 3 / units.s,
            'flow_vol': units.m ** 3 / units.s,
            'conc_mass_phase_comp': units.kg * units.m**-3,
            'flow_mol_phase_comp': units.mol / units.s,
            'mole_frac_phase_comp': units.dimensionless,
            'molality_comp': units.mole / units.kg,
            'visc_d_phase': units.Pa * units.s,
            'diffus_phase': units.m**2 * units.s**-1,
            'osm_coeff': units.dimensionless,
            'pressure_osm': units.Pa,
            'enth_mass_phase': units.J * units.kg**-1,
            'enth_flow': units.J / units.s}

        for (v, u) in var_unit_dict.items():
            var = getattr(m.fs.stream[0], v)
            assert_units_equivalent(var, u)

        assert_units_consistent(m)

    @pytest.mark.unit
    def test_scaling(self, frame):
        m = frame

        set_scaling_factor(m.fs.stream[0].flow_mass_phase_comp['Liq', 'H2O'], 1)
        set_scaling_factor(m.fs.stream[0].flow_mass_phase_comp['Liq', 'NaCl'], 1e2)
        calculate_scaling_factors(m.fs)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        for v in unscaled_var_list:
            print(v)
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

        # check if any variables are badly scaled
        badly_scaled_var_list = list(badly_scaled_var_generator(m))
        assert len(badly_scaled_var_list) == 0

    @pytest.mark.unit
    def test_dof(self, frame):
        m = frame
        assert (degrees_of_freedom(m) == 0)

    @pytest.mark.component
    def test_solve(self, frame):
        m = frame

        m.fs.stream.initialize()

        results = solver.solve(m)

        # Check for optimal solution
        assert results.solver.termination_condition == \
               TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution(self, frame):
        m = frame

        # specified conditions
        assert pytest.approx(0.035, rel=1e-3) == value(m.fs.stream[0].flow_mass_phase_comp['Liq', 'NaCl'])
        assert pytest.approx(0.965, rel=1e-3) == value(m.fs.stream[0].flow_mass_phase_comp['Liq', 'H2O'])
        assert pytest.approx(273.15 + 25, rel=1e-3) == value(m.fs.stream[0].temperature)
        assert pytest.approx(101325, rel=1e-3) == value(m.fs.stream[0].pressure)

        # calculated properties
        assert pytest.approx(0.035, rel=1e-3) == value(m.fs.stream[0].mass_frac_phase_comp['Liq', 'NaCl'])
        assert pytest.approx(0.965, rel=1e-3) == value(m.fs.stream[0].mass_frac_phase_comp['Liq', 'H2O'])
        assert pytest.approx(1021.5, rel=1e-3) == value(m.fs.stream[0].dens_mass_phase['Liq'])
        assert pytest.approx(9.790e-4, rel=1e-3) == value(m.fs.stream[0].flow_vol_phase['Liq'])
        assert pytest.approx(9.790e-4, rel=1e-3) == value(m.fs.stream[0].flow_vol)
        assert pytest.approx(985.7, rel=1e-3) == value(m.fs.stream[0].conc_mass_phase_comp['Liq', 'H2O'])
        assert pytest.approx(35.75, rel=1e-3) == value(m.fs.stream[0].conc_mass_phase_comp['Liq', 'NaCl'])
        assert pytest.approx(53.57, rel=1e-3) == value(m.fs.stream[0].flow_mol_phase_comp['Liq', 'H2O'])
        assert pytest.approx(0.5989, rel=1e-3) == value(m.fs.stream[0].flow_mol_phase_comp['Liq', 'NaCl'])
        assert pytest.approx(0.9889, rel=1e-3) == value(m.fs.stream[0].mole_frac_phase_comp['Liq', 'H2O'])
        assert pytest.approx(1.106e-2, rel=1e-3) == value(m.fs.stream[0].mole_frac_phase_comp['Liq', 'NaCl'])
        assert pytest.approx(0.6206, rel=1e-3) == value(m.fs.stream[0].molality_comp['NaCl'])
        assert pytest.approx(1.055e-3, rel=1e-3) == value(m.fs.stream[0].visc_d_phase['Liq'])
        assert pytest.approx(1.472e-9, rel=1e-3) == value(m.fs.stream[0].diffus_phase['Liq'])
        assert pytest.approx(0.9271, rel=1e-3) == value(m.fs.stream[0].osm_coeff)
        assert pytest.approx(2.853e6, rel=1e-3) == value(m.fs.stream[0].pressure_osm)
        assert pytest.approx(9.974e4, rel=1e-3) == value(m.fs.stream[0].enth_mass_phase['Liq'])
        assert pytest.approx(9.974e4, rel=1e-3) == value(m.fs.stream[0].enth_flow)
