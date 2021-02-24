import pytest
from pyomo.environ import (ConcreteModel,
                           TerminationCondition,
                           SolverStatus,
                           value,
                           Var)
from pyomo.network import Port
from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType)
from proteuslib.unit_models.pump_isothermal import Pump
import proteuslib.property_models.seawater_prop_pack as props

from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              fixed_variables_set,
                                              activated_constraints_set,
                                              number_unused_variables)
from idaes.core.util.testing import get_default_solver

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

        feed_mass_frac_H2O = 1 - feed_mass_frac_TDS
        m.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'TDS'].fix(
            feed_flow_mass * feed_mass_frac_TDS)
        m.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(
            feed_flow_mass * feed_mass_frac_H2O)
        m.fs.unit.inlet.pressure[0].fix(feed_pressure_in)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.outlet.pressure[0].fix(feed_pressure_out)
        return m