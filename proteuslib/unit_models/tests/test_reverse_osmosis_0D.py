import pytest
from pyomo.environ import (ConcreteModel,
                           Constraint,
                           TerminationCondition,
                           SolverStatus,
                           value,
                           Var)
from pyomo.network import Port
from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType)
from proteuslib.unit_models.reverse_osmosis_0D import RO_0D
import proteuslib.property_models.NaCl_prop_pack as props

from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              fixed_variables_set,
                                              activated_constraints_set,
                                              number_unused_variables)
from idaes.core.util.testing import get_default_solver
from idaes.core.util.exceptions import BalanceTypeNotSupportedError
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     unscaled_variables_generator,
                                     unscaled_constraints_generator,
                                     badly_scaled_var_generator)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()


# -----------------------------------------------------------------------------

class TestReverseOsmosis():
    @pytest.fixture(scope="class")
    def RO_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = RO_0D(default={
            "property_package": m.fs.properties,
            "has_pressure_change": True, })

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_pressure = 50e5
        feed_temperature = 273.15 + 25
        membrane_pressure_drop = 3e5
        membrane_area = 50
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 101325

        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
        m.fs.unit.inlet.flow_mass_comp[0, 'NaCl'].fix(
            feed_flow_mass * feed_mass_frac_NaCl)
        m.fs.unit.inlet.flow_mass_comp[0, 'H2O'].fix(
            feed_flow_mass * feed_mass_frac_H2O)
        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.deltaP.fix(-membrane_pressure_drop)
        m.fs.unit.area.fix(membrane_area)
        m.fs.unit.A.fix(A)
        m.fs.unit.B.fix(B)
        m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
        return m

    @pytest.mark.unit
    def test_config(self, RO_frame):
        m = RO_frame
        # check unit config arguments
        assert len(m.fs.unit.config) == 8

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.material_balance_type == \
               MaterialBalanceType.useDefault
        assert m.fs.unit.config.energy_balance_type == \
               EnergyBalanceType.useDefault
        assert m.fs.unit.config.momentum_balance_type == \
               MomentumBalanceType.pressureTotal
        assert m.fs.unit.config.has_pressure_change
        assert m.fs.unit.config.property_package is \
               m.fs.properties

    @pytest.mark.unit
    def test_build(self, RO_frame):
        m = RO_frame

        # test ports and variables
        port_lst = ['inlet', 'retentate', 'permeate']
        port_vars_lst = ['flow_mass_comp', 'pressure', 'temperature']
        for port_str in port_lst:
            assert hasattr(m.fs.unit, port_str)
            port = getattr(m.fs.unit, port_str)
            assert len(port.vars) == 3
            assert isinstance(port, Port)
            for var_str in port_vars_lst:
                assert hasattr(port, var_str)
                var = getattr(port, var_str)
                assert isinstance(var, Var)

        # test unit objects (including parameters, variables, and constraints)
        unit_objs_lst = ['A', 'B', 'dens_H2O',
                         'flux_mass_comp_in', 'flux_mass_comp_out', 'area',
                         'deltaP', 'mass_transfer_comp', 'flux_mass_comp_avg',
                         'eq_mass_transfer_term', 'eq_permeate_production',
                         'eq_flux_in', 'eq_flux_out',
                         'eq_connect_mass_transfer', 'eq_connect_enthalpy_transfer',
                         'eq_permeate_isothermal']
        for obj_str in unit_objs_lst:
            assert hasattr(m.fs.unit, obj_str)


        # test state block objects
        cv_name = 'feed_side'
        cv_stateblock_lst = ['properties_in', 'properties_out']
        stateblock_objs_lst = \
            ['flow_mass_comp', 'pressure', 'temperature', 'pressure_osm',
             'osm_coeff', 'mass_frac_comp', 'conc_mass_comp', 'dens_mass',
             'enth_mass',
             'eq_pressure_osm', 'eq_osm_coeff', 'eq_mass_frac_comp',
             'eq_conc_mass_comp', 'eq_dens_mass', 'eq_enth_mass'
             ]
        # control volume
        assert hasattr(m.fs.unit, cv_name)
        cv_blk = getattr(m.fs.unit, cv_name)
        for blk_str in cv_stateblock_lst:
            assert hasattr(cv_blk, blk_str)
            blk = getattr(cv_blk, blk_str)
            for obj_str in stateblock_objs_lst:
                assert hasattr(blk[0], obj_str)
        # permeate
        assert hasattr(m.fs.unit, 'properties_permeate')
        blk = getattr(m.fs.unit, 'properties_permeate')
        for var_str in stateblock_objs_lst:
            assert hasattr(blk[0], var_str)

        # test statistics
        assert number_variables(m) == 104
        assert number_total_constraints(m) == 40
        assert number_unused_variables(m) == 13  # vars from property package parameters

    @pytest.mark.unit
    def test_dof(self, RO_frame):
        m = RO_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, RO_frame):
        m = RO_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_initialize(self, RO_frame):
        m = RO_frame

        orig_fixed_vars = fixed_variables_set(m)
        orig_act_consts = activated_constraints_set(m)

        m.fs.unit.initialize(
            optarg={'nlp_scaling_method': 'user-scaling'})

        assert degrees_of_freedom(m) == 0

        fin_fixed_vars = fixed_variables_set(m)
        fin_act_consts = activated_constraints_set(m)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    @pytest.mark.component
    def test_var_scaling(self, RO_frame):
        m = RO_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, RO_frame):
        m = RO_frame
        solver.options = {'nlp_scaling_method': 'user-scaling'}
        results = solver.solve(m)

        # Check for optimal solution
        assert results.solver.termination_condition == \
               TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_conservation(self, RO_frame):
        m = RO_frame
        b = m.fs.unit
        comp_lst = ['NaCl', 'H2O']

        flow_mass_inlet = sum(
            b.feed_side.properties_in[0].flow_mass_comp[j] for j in comp_lst)
        flow_mass_retentate = sum(
            b.feed_side.properties_out[0].flow_mass_comp[j] for j in comp_lst)
        flow_mass_permeate = sum(
            b.properties_permeate[0].flow_mass_comp[j] for j in comp_lst)

        assert (abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate
                          )) <= 1e-6)

        assert (abs(value(
            flow_mass_inlet * b.feed_side.properties_in[0].enth_mass
            - flow_mass_retentate * b.feed_side.properties_out[0].enth_mass
            - flow_mass_permeate * b.properties_permeate[0].enth_mass
        )) <= 1e-6)

    @pytest.mark.component
    def test_solution(self, RO_frame):
        m = RO_frame
        assert (pytest.approx(5.90e-3, abs=1e-5) ==
                value(m.fs.unit.flux_mass_comp_avg[0, 'H2O']))
        assert (pytest.approx(1.51e-6, abs=1e-8) ==
                value(m.fs.unit.flux_mass_comp_avg[0, 'NaCl']))
        assert (pytest.approx(0.295, abs=1e-3) ==
                value(m.fs.unit.properties_permeate[0].
                      flow_mass_comp['H2O']))
        assert (pytest.approx(7.57e-5, abs=1e-7) ==
                value(m.fs.unit.properties_permeate[0].
                      flow_mass_comp['NaCl']))