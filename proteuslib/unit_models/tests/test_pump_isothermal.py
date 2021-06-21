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
                           TerminationCondition,
                           SolverStatus,
                           value,
                           Constraint,
                           Var)
from idaes.core import FlowsheetBlock
from proteuslib.unit_models.pump_isothermal import Pump
import proteuslib.property_models.seawater_prop_pack as props

from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     unscaled_variables_generator,
                                     unscaled_constraints_generator,
                                     get_scaling_factor,
                                     set_scaling_factor)


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------

class TestPumpIsothermal():
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
            var = getattr(m.fs.unit, v)
            assert isinstance(var, Var)

        con_list = ['ratioP_calculation', 'fluid_work_calculation', 'actual_work']
        for c in con_list:
            assert hasattr(m.fs.unit, c)
            con = getattr(m.fs.unit, c)
            assert isinstance(con, Constraint)

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
        initialization_tester(Pump_frame)

    @pytest.mark.component
    def test_solve(self, Pump_frame):
        m = Pump_frame
        solver.options = {'nlp_scaling_method': 'user-scaling'}
        results = solver.solve(m)

        # Check for optimal solution
        assert results.solver.termination_condition == \
               TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_conservation(self, Pump_frame):
        m = Pump_frame
        b = m.fs.unit.control_volume
        comp_lst = ['TDS', 'H2O']

        for t in m.fs.config.time:
            # mass balance
            for j in comp_lst:
                assert (abs(
                    value(b.properties_in[t].flow_mass_phase_comp['Liq', j]
                          - b.properties_out[t].flow_mass_phase_comp['Liq', j]))
                        <= 1e-6)
            # energy balance
            assert (abs(
                value(b.properties_in[t].enth_flow
                      - b.properties_out[t].enth_flow))
                    <= 1e-6)

    @pytest.mark.component
    def test_solution(self, Pump_frame):
        m = Pump_frame

        # pump variables
        assert (pytest.approx(521.06, rel=1e-3) == value(m.fs.unit.work_mechanical[0]))
        assert (pytest.approx(4e5, rel=1e-3) == value(m.fs.unit.deltaP[0]))
        assert (pytest.approx(5, rel=1e-3) == value(m.fs.unit.ratioP[0]))
        assert (pytest.approx(390.79, rel=1e-3) == value(m.fs.unit.work_fluid[0]))

        # outlet state variables
        assert (pytest.approx(0.965, rel=1e-3) ==
                value(m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp['Liq', 'H2O']))
        assert (pytest.approx(0.035, rel=1e-3) ==
                value(m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp['Liq', 'TDS']))
        assert (pytest.approx(298.15, rel=1e-5) ==
                value(m.fs.unit.control_volume.properties_out[0].temperature))
