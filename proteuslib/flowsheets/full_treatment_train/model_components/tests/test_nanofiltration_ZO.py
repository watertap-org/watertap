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
                           value,
                           Param,
                           Var,
                           Expression,
                           Constraint,
                           assert_optimal_termination)
from pyomo.network import Port
from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        ControlVolume0DBlock)
from proteuslib.flowsheets.full_treatment_train.model_components.nanofiltration_ZO import NanofiltrationZO
import proteuslib.flowsheets.full_treatment_train.model_components.seawater_ion_prop_pack as props

from idaes.core.util import get_solver
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     constraint_scaling_transform,
                                     unscaled_variables_generator,
                                     unscaled_constraints_generator,
                                     badly_scaled_var_generator)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.PropParameterBlock()
    m.fs.unit = NanofiltrationZO(default={"property_package": m.fs.properties})

    assert len(m.fs.unit.config) == 8

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == \
           MaterialBalanceType.useDefault
    assert m.fs.unit.config.energy_balance_type == \
           EnergyBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == \
           MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.property_package is \
           m.fs.properties

@pytest.mark.unit
def test_option_has_pressure_change():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.PropParameterBlock()
    m.fs.unit = NanofiltrationZO(
        default={"property_package": m.fs.properties,
                 "has_pressure_change": True})

    assert isinstance(m.fs.unit.feed_side.deltaP, Var)
    assert isinstance(m.fs.unit.deltaP, Var)

class TestNanofiltration():
    @pytest.fixture(scope="class")
    def unit_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = props.PropParameterBlock()

        m.fs.unit = NanofiltrationZO(default={"property_package": m.fs.properties})

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac = {'Na': 11122e-6,
                          'Ca': 382e-6,
                          'Mg': 1394e-6,
                          'SO4': 2136e-6,
                          'Cl': 20316.88e-6}
        m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp['Liq', 'H2O'].fix(
            feed_flow_mass * (1 - sum(x for x in feed_mass_frac.values())))
        for j in feed_mass_frac:
            m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp['Liq', j].fix(feed_flow_mass * feed_mass_frac[j])
        m.fs.unit.feed_side.properties_in[0].pressure.fix(101325)
        m.fs.unit.feed_side.properties_in[0].temperature.fix(298.15)

        m.fs.unit.flux_vol_solvent.fix(1.67e-6)
        m.fs.unit.area.fix(500)
        m.fs.unit.properties_permeate[0].pressure.fix(101325)

        m.fs.unit.rejection_phase_comp[0, 'Liq', 'Na'].fix(0.01)
        m.fs.unit.rejection_phase_comp[0, 'Liq', 'Ca'].fix(0.79)
        m.fs.unit.rejection_phase_comp[0, 'Liq', 'Mg'].fix(0.94)
        m.fs.unit.rejection_phase_comp[0, 'Liq', 'SO4'].fix(0.87)
        m.fs.unit.rejection_phase_comp[0, 'Liq', 'Cl'] = 0.15  # guess, but electroneutrality enforced below
        charge_comp = {'Na': 1, 'Ca': 2, 'Mg': 2, 'SO4': -2, 'Cl': -1}
        m.fs.unit.eq_electroneutrality = Constraint(
            expr=0 ==
                 sum(charge_comp[j]
                     * m.fs.unit.feed_side.properties_out[0].conc_mol_phase_comp['Liq', j]
                     for j in charge_comp))
        constraint_scaling_transform(m.fs.unit.eq_electroneutrality, 1)

        return m

    @pytest.mark.unit
    def test_build(self, unit_frame):
        m = unit_frame

        # test ports
        port_lst = ['inlet', 'retentate', 'permeate']
        for port_str in port_lst:
            port = getattr(m.fs.unit, port_str)
            assert len(port.vars) == 3  # number of state variables for NaCl property package
            assert isinstance(port, Port)

        # test statistics
        assert number_variables(m) == 89
        assert number_total_constraints(m) == 74
        assert number_unused_variables(m) == 1

    @pytest.mark.unit
    def test_dof(self, unit_frame):
        m = unit_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, unit_frame):
        m = unit_frame

        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'Na'))
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e4, index=('Liq', 'Ca'))
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e3, index=('Liq', 'Mg'))
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e3, index=('Liq', 'SO4'))
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'Cl'))
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
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_conservation(self, unit_frame):
        m = unit_frame
        b = m.fs.unit
        comp_lst = m.fs.properties.component_list

        flow_mass_inlet = sum(
            b.feed_side.properties_in[0].flow_mass_phase_comp['Liq', j] for j in comp_lst)
        flow_mass_retentate = sum(
            b.feed_side.properties_out[0].flow_mass_phase_comp['Liq', j] for j in comp_lst)
        flow_mass_permeate = sum(
            b.properties_permeate[0].flow_mass_phase_comp['Liq', j] for j in comp_lst)

        assert (abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate
                          )) <= 1e-5)

    @pytest.mark.component
    def test_solution(self, unit_frame):
        m = unit_frame
        assert (pytest.approx(0.8350, rel=1e-3) ==
                value(m.fs.unit.properties_permeate[0].flow_mass_phase_comp['Liq', 'H2O']))
        assert (pytest.approx(2.002, rel=1e-3) ==
                value(m.fs.unit.properties_permeate[0].conc_mol_phase_comp['Liq', 'Ca']))
        assert (pytest.approx(484.0, rel=1e-3) ==
                value(m.fs.unit.properties_permeate[0].conc_mol_phase_comp['Liq', 'Cl']))
        assert (pytest.approx(3.441, rel=1e-3) ==
                value(m.fs.unit.properties_permeate[0].conc_mol_phase_comp['Liq', 'Mg']))
        assert (pytest.approx(478.9, rel=1e-3) ==
                value(m.fs.unit.properties_permeate[0].conc_mol_phase_comp['Liq', 'Na']))
        assert (pytest.approx(2.891, rel=1e-3) ==
                value(m.fs.unit.properties_permeate[0].conc_mol_phase_comp['Liq', 'SO4']))

    @pytest.mark.unit
    def test_report(self, unit_frame):
        unit_frame.fs.unit.report()
