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
from pyomo.environ import (ConcreteModel,
                           Var,
                           Expression,
                           Constraint,
                           TerminationCondition,
                           SolverStatus,
                           value)
from pyomo.network import Port
from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType)
from proteuslib.unit_models.pressure_exchanger import PressureExchanger
import proteuslib.property_models.seawater_prop_pack as props

from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              fixed_variables_set,
                                              activated_constraints_set,
                                              number_unused_variables)
from idaes.core.util import get_solver
from idaes.core.util.testing import initialization_tester
from idaes.core.util.initialization import solve_indexed_blocks
from idaes.core.util.exceptions import BalanceTypeNotSupportedError
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     unscaled_variables_generator,
                                     unscaled_constraints_generator,
                                     badly_scaled_var_generator)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# # -----------------------------------------------------------------------------

@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={'dynamic': False})
    m.fs.properties = props.SeawaterParameterBlock()
    m.fs.unit = PressureExchanger(default={'property_package': m.fs.properties})

    # check unit config arguments
    assert len(m.fs.unit.config) == 7

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == \
           MaterialBalanceType.useDefault
    assert m.fs.unit.config.energy_balance_type == \
           EnergyBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == \
           MomentumBalanceType.pressureTotal
    assert m.fs.unit.config.property_package is m.fs.properties

@pytest.mark.unit
def test_build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={'dynamic': False})
    m.fs.properties = props.SeawaterParameterBlock()
    m.fs.unit = PressureExchanger(default={'property_package': m.fs.properties})

    # test ports and state variables
    port_lst = ['low_pressure_inlet', 'low_pressure_outlet',
                'high_pressure_inlet', 'high_pressure_outlet']
    port_vars_lst = ['flow_mass_phase_comp', 'pressure', 'temperature']
    for port_str in port_lst:
        assert hasattr(m.fs.unit, port_str)
        port = getattr(m.fs.unit, port_str)
        assert len(port.vars) == 3
        assert isinstance(port, Port)
        for var_str in port_vars_lst:
            assert hasattr(port, var_str)
            var = getattr(port, var_str)
            assert isinstance(var, Var)

    # test unit variables
    assert hasattr(m.fs.unit, 'efficiency_pressure_exchanger')
    assert isinstance(m.fs.unit.efficiency_pressure_exchanger, Var)

    # test unit constraints
    unit_cons_lst = ['eq_pressure_transfer', 'eq_equal_flow_vol', 'eq_equal_low_pressure']
    for c in unit_cons_lst:
        assert hasattr(m.fs.unit, c)
        con = getattr(m.fs.unit, c)
        assert isinstance(con, Constraint)

    # test control volumes, only terms directly used by pressure exchanger
    cv_list = ['low_pressure_side', 'high_pressure_side']
    cv_var_lst = ['deltaP']
    cv_con_lst = ['eq_isothermal_temperature']
    cv_exp_lst = ['work']
    for cv_str in cv_list:
        assert hasattr(m.fs.unit, cv_str)
        cv = getattr(m.fs.unit, cv_str)
        for cv_var_str in cv_var_lst:
            cv_var = getattr(cv, cv_var_str)
            assert isinstance(cv_var, Var)
        for cv_con_str in cv_con_lst:
            cv_con = getattr(cv, cv_con_str)
            assert isinstance(cv_con, Constraint)
        for cv_exp_str in cv_exp_lst:
            cv_exp = getattr(cv, cv_exp_str)
            assert isinstance(cv_exp, Expression)

    # test state blocks, only terms directly used by pressure exchanger
    cv_stateblock_lst = ['properties_in', 'properties_out']
    stateblock_var_lst = ['pressure', 'temperature']
    stateblock_exp_lst = ['flow_vol']
    for cv_str in cv_list:
        cv = getattr(m.fs.unit, cv_str)
        for cv_sb_str in cv_stateblock_lst:
            assert hasattr(cv, cv_sb_str)
            cv_sb = getattr(cv, cv_sb_str)
            for sb_var_str in stateblock_var_lst:
                assert hasattr(cv_sb[0], sb_var_str)
                sb_var = getattr(cv_sb[0], sb_var_str)
                assert isinstance(sb_var, Var)
            for sb_exp_str in stateblock_exp_lst:
                assert hasattr(cv_sb[0], sb_exp_str)
                sb_exp = getattr(cv_sb[0], sb_exp_str)
                assert isinstance(sb_exp, Expression)

    # test statistics
    assert number_variables(m.fs.unit) == 39
    assert number_total_constraints(m.fs.unit) == 31
    assert number_unused_variables(m.fs.unit) == 0

class TestPressureExchanger():
    @pytest.fixture(scope="class")
    def unit_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={'dynamic': False})
        m.fs.properties = props.SeawaterParameterBlock()
        m.fs.unit = PressureExchanger(default={'property_package': m.fs.properties})

        # Specify inlet conditions
        temperature = 25 + 273.15
        flow_vol = 1e-3
        lowP_mass_frac_TDS = 0.035
        lowP_pressure = 101325
        highP_mass_frac_TDS = 0.07
        highP_pressure = 50e5

        m.fs.unit.low_pressure_side.properties_in[0].flow_vol_phase['Liq'].fix(flow_vol)
        m.fs.unit.low_pressure_side.properties_in[0].mass_frac_phase_comp['Liq', 'TDS'].fix(lowP_mass_frac_TDS)
        m.fs.unit.low_pressure_side.properties_in[0].pressure.fix(lowP_pressure)
        m.fs.unit.low_pressure_side.properties_in[0].temperature.fix(temperature)

        m.fs.unit.high_pressure_side.properties_in[0].flow_vol_phase['Liq'].fix(flow_vol)
        m.fs.unit.high_pressure_side.properties_in[0].mass_frac_phase_comp['Liq', 'TDS'].fix(highP_mass_frac_TDS)
        m.fs.unit.high_pressure_side.properties_in[0].pressure.fix(highP_pressure)
        m.fs.unit.high_pressure_side.properties_in[0].temperature.fix(temperature)

        # solve inlet conditions and only fix state variables (i.e. unfix flow_vol and mass_frac_phase)
        solver.options['nlp_scaling_method'] = 'user-scaling'
        results = solver.solve(m.fs.unit.low_pressure_side.properties_in[0])
        assert results.solver.termination_condition == TerminationCondition.optimal
        m.fs.unit.low_pressure_side.properties_in[0].flow_mass_phase_comp['Liq', 'TDS'].fix()
        m.fs.unit.low_pressure_side.properties_in[0].flow_vol_phase['Liq'].unfix()
        m.fs.unit.low_pressure_side.properties_in[0].mass_frac_phase_comp['Liq', 'TDS'].unfix()

        results = solver.solve(m.fs.unit.high_pressure_side.properties_in[0])
        assert results.solver.termination_condition == TerminationCondition.optimal
        m.fs.unit.high_pressure_side.properties_in[0].flow_mass_phase_comp['Liq', 'H2O'].fix()
        m.fs.unit.high_pressure_side.properties_in[0].flow_mass_phase_comp['Liq', 'TDS'].fix()
        m.fs.unit.high_pressure_side.properties_in[0].flow_vol_phase['Liq'].unfix()
        m.fs.unit.high_pressure_side.properties_in[0].mass_frac_phase_comp['Liq', 'TDS'].unfix()

        # Specify unit
        efficiency = 0.95
        m.fs.unit.efficiency_pressure_exchanger.fix(efficiency)
        return m

    @pytest.mark.unit
    def test_dof(self, unit_frame):
        assert degrees_of_freedom(unit_frame) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, unit_frame):
        m = unit_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m.fs.unit, include_fixed=True))
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_initialize(self, unit_frame):
        m = unit_frame
        kwargs = {'solver': 'ipopt',
                  'optarg': {'nlp_scaling_method': 'user-scaling'}}
        initialization_tester(unit_frame, **kwargs)


    @pytest.mark.component
    def test_var_scaling(self, unit_frame):
        m = unit_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, unit_frame):
        m = unit_frame
        solver.options = {'nlp_scaling_method': 'user-scaling'}
        results = solver.solve(m)

        # Check for optimal solution
        assert results.solver.termination_condition == \
               TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution(self, unit_frame):
        m = unit_frame
        assert (pytest.approx(0.9877, rel=1e-3) ==
                value(m.fs.unit.low_pressure_inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']))
        assert (pytest.approx(3.582e-2, rel=1e-3) ==
                value(m.fs.unit.low_pressure_inlet.flow_mass_phase_comp[0, 'Liq', 'TDS']))
        assert (pytest.approx(4.755e6, rel=1e-3) ==
                value(m.fs.unit.low_pressure_outlet.pressure[0]))
        assert (pytest.approx(0.9767, rel=1e-3) ==
                value(m.fs.unit.high_pressure_inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']))
        assert (pytest.approx(7.352e-2, rel=1e-3) ==
                value(m.fs.unit.high_pressure_inlet.flow_mass_phase_comp[0, 'Liq', 'TDS']))
        assert (pytest.approx(4.654e6, rel=1e-3) ==
                value(m.fs.unit.low_pressure_side.deltaP[0]))
        assert (pytest.approx(4.654e3, rel=1e-3) ==
                value(m.fs.unit.low_pressure_side.work[0]))
        assert (pytest.approx(-4.899e6, rel=1e-3) ==
                value(m.fs.unit.high_pressure_side.deltaP[0]))
        assert (pytest.approx(-4.899e3, rel=1e-3) ==
                value(m.fs.unit.high_pressure_side.work[0]))

    @pytest.mark.unit
    def test_report(self, unit_frame):
        unit_frame.fs.unit.report()
