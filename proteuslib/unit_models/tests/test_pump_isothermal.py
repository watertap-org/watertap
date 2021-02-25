import pytest
from pyomo.environ import (ConcreteModel,
                           TerminationCondition,
                           SolverStatus,
                           value,
                           Constraint,
                           Var)
from pyomo.network import Port
from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType)
from proteuslib.unit_models.pump_isothermal import Pump
import proteuslib.property_models.seawater_prop_pack as props

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import get_default_solver
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     unscaled_variables_generator,
                                     unscaled_constraints_generator,
                                     get_scaling_factor,
                                     set_scaling_factor)


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()


# -----------------------------------------------------------------------------

class TestReverseOsmosis():
    @pytest.fixture(scope="class")
    def Pump_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = props.SeawaterParameterBlock()

        m.fs.unit = Pump(default={"property_package": m.fs.properties})

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_TDS = 0.035
        feed_pressure_in = 1e5
        feed_pressure_out = 5e5
        feed_temperature = 273.15 + 25
        efi_pump = 0.75

        feed_mass_frac_H2O = 1 - feed_mass_frac_TDS
        m.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'TDS'].fix(
            feed_flow_mass * feed_mass_frac_TDS)
        m.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(
            feed_flow_mass * feed_mass_frac_H2O)
        m.fs.unit.inlet.pressure[0].fix(feed_pressure_in)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.outlet.pressure[0].fix(feed_pressure_out)
        m.fs.unit.efficiency_pump.fix(efi_pump)
        return m

    @pytest.mark.unit
    def test_build(self, Pump_frame):
        m = Pump_frame

        # check that IDAES pump model has not changed variable and constraint names
        var_list = ['work_mechanical', 'deltaP', 'ratioP', 'work_fluid', 'efficiency_pump']
        for v in var_list:
            assert hasattr(m.fs.unit, v)

        con_list = ['ratioP_calculation', 'fluid_work_calculation', 'actual_work']
        for c in con_list:
            assert hasattr(m.fs.unit, c)

        # check that IDAES control volume has not changed variable and constraint names
        assert hasattr(m.fs.unit.control_volume, 'work')
        assert hasattr(m.fs.unit.control_volume, 'deltaP')

        assert hasattr(m.fs.unit.control_volume, 'material_balances')
        assert hasattr(m.fs.unit.control_volume, 'pressure_balance')

        # check that energy balance was removed
        assert not hasattr(m.fs.unit.control_volume, 'enthalpy_balances')

        # check that isothermal constraint was added
        assert hasattr(m.fs.unit.control_volume, 'isothermal_balance')
        assert isinstance(m.fs.unit.control_volume.isothermal_balance, Constraint)

    @pytest.mark.unit
    def test_dof(self, Pump_frame):
        m = Pump_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, Pump_frame):
        m = Pump_frame

        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'TDS'))
        calculate_scaling_factors(m)

        if get_scaling_factor(m.fs.unit.ratioP) is None:  # if IDAES hasn't specified a scaling factor
            set_scaling_factor(m.fs.unit.ratioP, 1)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_initialize(self, Pump_frame):
        m = Pump_frame

        # orig_fixed_vars = fixed_variables_set(m)
        # orig_act_consts = activated_constraints_set(m)
        m.fs.unit.initialize()

        assert degrees_of_freedom(m) == 0
        #
        # fin_fixed_vars = fixed_variables_set(m)
        # fin_act_consts = activated_constraints_set(m)
        #
        # assert len(fin_act_consts) == len(orig_act_consts)
        # assert len(fin_fixed_vars) == len(orig_fixed_vars)
        #
        # for c in fin_act_consts:
        #     assert c in orig_act_consts
        # for v in fin_fixed_vars:
        #     assert v in orig_fixed_vars
