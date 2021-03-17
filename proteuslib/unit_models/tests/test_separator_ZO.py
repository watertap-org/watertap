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
from proteuslib.unit_models.separator_ZO import SeparatorZO
import proteuslib.property_models.seawater_prop_pack as props

from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              fixed_variables_set,
                                              activated_constraints_set,
                                              number_unused_variables)
from idaes.core.util.testing import (get_default_solver,
                                     initialization_tester)
from idaes.core.util.exceptions import BalanceTypeNotSupportedError
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     unscaled_variables_generator,
                                     unscaled_constraints_generator,
                                     badly_scaled_var_generator)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()


# -----------------------------------------------------------------------------

@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={'dynamic': False})
    m.fs.properties = props.SeawaterParameterBlock()
    m.fs.unit = SeparatorZO(default={'property_package': m.fs.properties})

    # check unit config arguments
    assert len(m.fs.unit.config) == 5

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert not m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.property_package is m.fs.properties

@pytest.mark.unit
def test_default_build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={'dynamic': False})
    m.fs.properties = props.SeawaterParameterBlock()
    m.fs.unit = SeparatorZO(default={'property_package': m.fs.properties})

    # test ports and state variables
    port_lst = ['inlet', 'outlet', 'waste']
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
    assert hasattr(m.fs.unit, 'recovery_frac_phase_comp')
    assert isinstance(m.fs.unit.recovery_frac_phase_comp, Var)
    assert hasattr(m.fs.unit, 'removal_frac_phase_comp')
    assert isinstance(m.fs.unit.removal_frac_phase_comp, Var)

    # test constraints
    unit_cons_lst = ['eq_component_mass_balance', 'eq_component_removal', 'eq_removal_to_recovery',
                     'eq_outlet_temperature', 'eq_waste_temperature',
                     'eq_outlet_pressure', 'eq_waste_pressure']
    for c in unit_cons_lst:
        assert hasattr(m.fs.unit, c)
        con = getattr(m.fs.unit, c)
        assert isinstance(con, Constraint)

    # test state block objects
    stateblock_lst = ['inlet_state', 'outlet_state', 'waste_state']
    stateblock_objs_lst = ['flow_mass_phase_comp', 'pressure', 'temperature']
    for blk_str in stateblock_lst:
        assert hasattr(m.fs.unit, blk_str)
        blk = getattr(m.fs.unit, blk_str)
        for obj_str in stateblock_objs_lst:
            assert hasattr(blk[0], obj_str)

    # test statistics
    assert number_variables(m) == 52
    assert number_total_constraints(m) == 10
    assert number_unused_variables(m) == 36  # vars from property package parameters

class TestSeparatorZO():
    @pytest.mark.unit
    def test_has_pressure_change_build(self, unit_frame):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={'dynamic': False})
        m.fs.properties = props.SeawaterParameterBlock()
        m.fs.unit = SeparatorZO(default={'property_package': m.fs.properties,
                                         'has_pressure_change': True})

        assert hasattr(m.fs.unit, 'deltaP_outlet')
        assert isinstance(m.fs.unit.deltaP_outlet, Var)
        assert hasattr(m.fs.unit, 'deltaP_waste')
        assert isinstance(m.fs.unit.deltaP_waste, Var)

        # test statistics
        assert number_variables(m) == 54
        assert number_total_constraints(m) == 10
        assert number_unused_variables(m) == 36  # vars from property package parameters

    @pytest.fixture(scope="class")
    def unit_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={'dynamic': False})
        m.fs.properties = props.SeawaterParameterBlock()
        m.fs.unit = SeparatorZO(default={'property_package': m.fs.properties,
                                         'has_pressure_change': True})

        # Specify feed conditions
        feed_flow_mass = 1
        feed_mass_frac_TDS = 0.035
        feed_pressure = 50e5
        feed_temperature = 273.15 + 25

        feed_mass_frac_H2O = 1 - feed_mass_frac_TDS
        m.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'TDS'].fix(
            feed_flow_mass * feed_mass_frac_TDS)
        m.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(
            feed_flow_mass * feed_mass_frac_H2O)
        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)

        # Specify unit
        water_recovery = 0.5
        TDS_removal = 0.95
        pressure_drop = 3e5
        pressure_atmospheric = 101325

        m.fs.unit.recovery_frac_phase_comp[:, :, 'H2O'].fix(water_recovery)
        m.fs.unit.removal_frac_phase_comp[:, :, 'TDS'].fix(TDS_removal)
        m.fs.unit.deltaP_waste.fix(-pressure_drop)
        m.fs.unit.outlet.pressure[0].fix(pressure_atmospheric)
        return m

    @pytest.mark.unit
    def test_dof(self, unit_frame):
        assert degrees_of_freedom(unit_frame) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, unit_frame):
        m = unit_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_initialize(self, unit_frame):
        initialization_tester(unit_frame)

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
    def test_conservation(self, unit_frame):
        m = unit_frame
        b = m.fs.unit
        comp_lst = ['TDS', 'H2O']


        # mass balance
        flow_mass_inlet = sum(
            b.inlet_state[0].flow_mass_phase_comp['Liq', j] for j in comp_lst)
        flow_mass_outlet = sum(
            b.outlet_state[0].flow_mass_phase_comp['Liq', j] for j in comp_lst)
        flow_mass_waste = sum(
            b.waste_state[0].flow_mass_phase_comp['Liq', j] for j in comp_lst)

        assert (abs(value(flow_mass_inlet - flow_mass_outlet - flow_mass_waste)) <= 1e-6)

        # energy balance, isothermal
        assert (abs(value(b.inlet_state[0].temperature - b.outlet_state[0].temperature) <= 1e-6))
        assert (abs(value(b.inlet_state[0].temperature - b.waste_state[0].temperature) <= 1e-6))

    @pytest.mark.component
    def test_solution(self, unit_frame):
        m = unit_frame
        assert (pytest.approx(0.4825, rel=1e-3) ==
                value(m.fs.unit.outlet_state[0].flow_mass_phase_comp['Liq', 'H2O']))
        assert (pytest.approx(1.75e-3, rel=1e-3) ==
                value(m.fs.unit.outlet_state[0].flow_mass_phase_comp['Liq', 'TDS']))
        assert (pytest.approx(-4898675, rel=1e-3) == value(m.fs.unit.deltaP_outlet[0]))
        assert (pytest.approx(0.4825, rel=1e-3) ==
                value(m.fs.unit.waste_state[0].flow_mass_phase_comp['Liq', 'H2O']))
        assert (pytest.approx(3.325e-2, rel=1e-3) ==
                value(m.fs.unit.waste_state[0].flow_mass_phase_comp['Liq', 'TDS']))
        assert (pytest.approx(4.7e6, rel=1e-3) == value(m.fs.unit.waste_state[0].pressure))