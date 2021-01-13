import pytest
from pyomo.environ import (ConcreteModel,
                           Constraint,
                           TerminationCondition,
                           SolverStatus,
                           value,
                           Var)
from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType)
from idaes.water_treatment.unit_models.reverse_osmosis import RO
import idaes.water_treatment.properties.NaCl_prop_pack as props

from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              fixed_variables_set,
                                              activated_constraints_set,
                                              number_unused_variables)
from idaes.core.util.testing import (get_default_solver)
from idaes.core.util.exceptions import BalanceTypeNotSupportedError
from idaes.core.util.scaling import badly_scaled_var_generator

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

        m.fs.unit = RO(default={
            "property_package": m.fs.properties,
            "has_pressure_change": True, })
        return m

    def test_config(self, RO_frame):
        # check unit config arguments
        assert len(RO_frame.fs.unit.config) == 16

        assert not RO_frame.fs.unit.config.dynamic
        assert not RO_frame.fs.unit.config.has_holdup
        assert RO_frame.fs.unit.config.material_balance_type == \
               MaterialBalanceType.useDefault
        assert RO_frame.fs.unit.config.energy_balance_type == \
               EnergyBalanceType.useDefault
        assert RO_frame.fs.unit.config.momentum_balance_type == \
               MomentumBalanceType.pressureTotal
        assert not RO_frame.fs.unit.config.has_rate_reactions
        assert not RO_frame.fs.unit.config.has_equilibrium_reactions
        assert not RO_frame.fs.unit.config.has_phase_equilibrium
        assert RO_frame.fs.unit.config.has_mass_transfer
        assert RO_frame.fs.unit.config.has_heat_transfer
        assert not RO_frame.fs.unit.config.has_heat_of_reaction
        assert RO_frame.fs.unit.config.property_package is \
               RO_frame.fs.properties
        assert RO_frame.fs.unit.config.reaction_package is None

    @pytest.mark.build
    def test_build(self, RO_frame):
        # test ports and variables
        port_lst = ['inlet', 'retentate', 'permeate']
        port_vars_lst = ['flow_mass_comp', 'pressure', 'temperature']
        for port_str in port_lst:
            assert hasattr(RO_frame.fs.unit, port_str)
            port = getattr(RO_frame.fs.unit, port_str)
            assert len(port.vars) == 3
            for var_str in port_vars_lst:
                assert hasattr(port, var_str)

        # test unit objects (including parameters, variables, and constraints)
        unit_objs_lst = ['A', 'B', 'dens_H2O',
                         'flux_mass_comp', 'Area', 'mass_transfer', 'heat_duty',
                         'deltaP',
                         'eq_mass_transfer', 'eq_flux',
                         'eq_permeate_production',
                         'eq_permeate_isothermal', 'eq_enthalpy_and_heat']
        for obj_str in unit_objs_lst:
            assert hasattr(RO_frame.fs.unit, obj_str)

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
        assert hasattr(RO_frame.fs.unit, cv_name)
        cv_blk = getattr(RO_frame.fs.unit, cv_name)
        for blk_str in cv_stateblock_lst:
            assert hasattr(cv_blk, blk_str)
            blk = getattr(cv_blk, blk_str)
            for obj_str in stateblock_objs_lst:
                assert hasattr(blk[0], obj_str)
        # permeate
        assert hasattr(RO_frame.fs.unit, 'properties_permeate')
        blk = getattr(RO_frame.fs.unit, 'properties_permeate')
        for var_str in stateblock_objs_lst:
            assert hasattr(blk[0], var_str)

        # test statistics
        assert number_variables(RO_frame) == 47
        assert number_total_constraints(RO_frame) == 40
        assert number_unused_variables(RO_frame) == 0

    def test_dof(self, RO_frame):
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_pressure = 50e5
        feed_temperature = 273.15 + 25
        membrane_pressure_drop = 3e5
        membrane_area = 50

        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
        RO_frame.fs.unit.inlet.flow_mass_comp[0, 'NaCl'].fix(
            feed_flow_mass * feed_mass_frac_NaCl)
        RO_frame.fs.unit.inlet.flow_mass_comp[0, 'H2O'].fix(
            feed_flow_mass * feed_mass_frac_H2O)
        RO_frame.fs.unit.inlet.pressure[0].fix(feed_pressure)
        RO_frame.fs.unit.inlet.temperature[0].fix(feed_temperature)
        RO_frame.fs.unit.deltaP.fix(-membrane_pressure_drop)
        RO_frame.fs.unit.Area.fix(membrane_area)

        assert degrees_of_freedom(RO_frame) == 0

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_initialize(self, RO_frame):
        orig_fixed_vars = fixed_variables_set(RO_frame)
        orig_act_consts = activated_constraints_set(RO_frame)

        RO_frame.fs.unit.initialize(
            optarg={'tol': 1e-6, 'max_iter': 5000,
                    'nlp_scaling_method': 'user-scaling'})

        assert degrees_of_freedom(RO_frame) == 0

        fin_fixed_vars = fixed_variables_set(RO_frame)
        fin_act_consts = activated_constraints_set(RO_frame)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_scaling(self, RO_frame):
        # assert RO_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(RO_frame))
        assert badly_scaled_var_lst == []

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solve(self, RO_frame):
        solver.options = {'tol': 1e-6, 'max_iter': 5000,
                          'nlp_scaling_method': 'user-scaling'}
        results = solver.solve(RO_frame)

        # Check for optimal solution
        assert results.solver.termination_condition == \
               TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_conservation(self, RO_frame):
        b = RO_frame.fs.unit
        comp_lst = ['NaCl', 'H2O']

        flow_mass_inlet = sum(
            b.feed_side.properties_in[0].flow_mass_comp[j] for j in comp_lst)
        flow_mass_retentate = sum(
            b.feed_side.properties_out[0].flow_mass_comp[j] for j in comp_lst)
        flow_mass_permeate = sum(
            b.properties_permeate[0].flow_mass_comp[j] for j in comp_lst)

        assert (abs(value(
            flow_mass_inlet - flow_mass_retentate - flow_mass_permeate
        )) <= 1e-6)

        assert (abs(value(
            flow_mass_inlet * b.feed_side.properties_in[0].enth_mass
            - flow_mass_retentate * b.feed_side.properties_out[0].enth_mass
            - flow_mass_permeate * b.properties_permeate[0].enth_mass
        )) <= 1e-6)

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solution(self, RO_frame):
        assert (pytest.approx(5.90e-3, abs=1e-5) ==
                value(RO_frame.fs.unit.flux_mass_comp[0, 'H2O', 'Average']))
        assert (pytest.approx(1.52e-6, abs=1e-8) ==
                value(RO_frame.fs.unit.flux_mass_comp[0, 'NaCl', 'Average']))
        assert (pytest.approx(0.295, abs=1e-3) ==
                value(RO_frame.fs.unit.properties_permeate[0].
                      flow_mass_comp['H2O']))
        assert (pytest.approx(7.61e-5, abs=1e-7) ==
                value(RO_frame.fs.unit.properties_permeate[0].
                      flow_mass_comp['NaCl']))