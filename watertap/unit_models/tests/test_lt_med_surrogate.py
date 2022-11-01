import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    Constraint,
    assert_optimal_termination,
)
from pyomo.network import Port
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from watertap.unit_models.lt_med_surrogate import LT_MED_surrogate
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from watertap.property_models.seawater_ion_generic import configuration
import watertap.property_models.seawater_prop_pack as sw_props
import watertap.property_models.water_prop_pack as w_props
from watertap.core.util.initialization import assert_no_degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent

from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    constraint_scaling_transform,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)


class TestLTMED:
    @pytest.fixture(scope="class")
    def LT_MED_frame(self):
        # create model, flowsheet
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties1 = sw_props.SeawaterParameterBlock()
        m.fs.properties2 = w_props.WaterParameterBlock()
        m.fs.unit = LT_MED_surrogate(
            property_package=m.fs.properties1, property_package2=m.fs.properties2
        )

        # System specification

        feed_salinity = 35  # g/L
        feed_temperature = 25  # degC
        steam_temperature = 80  # deg C
        sys_capacity = 2000  # m3/day
        recovery_rate = 0.5  # dimentionless

        m.fs.lt_med.feed_props[0].conc_mass_phase_comp["Liq", "TDS"].fix(feed_salinity)
        m.fs.lt_med.feed_props[0].temperature.fix(feed_temperature + 273.15)
        m.fs.lt_med.steam_props[0].temperature.fix(steam_temperature + 273.15)
        m.fs.lt_med.Capacity.fix(sys_capacity)
        m.fs.lt_med.RR.fix(recovery_rate)

        return m

    @pytest.mark.unit
    def test_config(self, LT_MED_frame):
        m = LT_MED_frame
        # check unit config arguments
        assert len(m.fs.unit.config) == 5

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.property_package is m.fs.properties1
        assert m.fs.unit.config.property_package2 is m.fs.properties2

    @pytest.mark.unit
    def test_build(self, LT_MED_frame):
        m = LT_MED_frame

        # test ports
        port_lst = ["feed", "distillate", "brine"]
        for port_str in port_lst:
            port = getattr(m.fs.unit, port_str)
            assert isinstance(port, Port)
            # number of state variables for seawater property package
            assert len(port.vars) == 0

        # test statistics
        assert number_variables(m) == 52
        assert number_total_constraints(m) == 52
        assert number_unused_variables(m) == 7  # vars from property package parameters

    @pytest.mark.unit
    def test_dof(self, LT_MED_frame):
        m = LT_MED_frame
        assert degrees_of_freedom(m) == 0
