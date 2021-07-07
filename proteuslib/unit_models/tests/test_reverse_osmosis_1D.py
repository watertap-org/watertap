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
                           Param,
                           Var,
                           Expression,
                           Constraint)
from pyomo.util.check_units import (assert_units_consistent,
                                    assert_units_equivalent)
from pyomo.network import Port
from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        ControlVolume1DBlock,
                        StateBlock)
from proteuslib.unit_models.reverse_osmosis_1D import ReverseOsmosis1D
# ConcentrationPolarizationType,
# MassTransferCoefficient,
# PressureChangeType)
import proteuslib.property_models.NaCl_prop_pack \
    as props

from idaes.core.util import get_solver
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_statistics import *
#                                               (degrees_of_freedom,
#                                               number_variables,
#                                               number_total_constraints,
#                                               number_unused_variables,
#                                               unused_variables_set,
#                                               unfixed_variables_in_activated_equalities_set,
#                                               variables_set,
#                                               fixed_variables_in_activated_equalities_set)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     unscaled_variables_generator,
                                     unscaled_constraints_generator,
                                     badly_scaled_var_generator,
                                     get_scaling_factor,
                                     set_scaling_factor)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()
# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = ReverseOsmosis1D(default={"property_package": m.fs.properties})

    assert len(m.fs.unit.config) == 13

    #
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
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = ReverseOsmosis1D(default={
        "property_package": m.fs.properties,
        "has_pressure_change": True})

    assert isinstance(m.fs.unit.feed_side.deltaP, Var)
    assert isinstance(m.fs.unit.deltaP, Var)

class TestReverseOsmosis():
    @pytest.fixture(scope="class")
    def RO_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis1D(default={
            "property_package": m.fs.properties,
            "has_pressure_change": True})

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_pressure = 50e5
        feed_temperature = 273.15 + 25
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 101325
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].fix(
            feed_flow_mass * feed_mass_frac_NaCl)

        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(
            feed_flow_mass * feed_mass_frac_H2O)

        m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.deltaP.fix(0)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate_outlet.pressure[0].fix(pressure_atmospheric)
        m.fs.unit.permeate_out[0].flow_mass_phase_comp['Liq', 'H2O'].setub(0.5)
        m.fs.unit.permeate_out[0].flow_mass_phase_comp['Liq', 'H2O'].setlb(0.25)
        m.fs.unit.length.fix(8)
        m.fs.unit.recovery_vol_phase[0, 'Liq'].fix(0.4)

        return m

    @pytest.mark.unit
    def test_build(self, RO_frame):
        m = RO_frame

        # test ports
        port_lst = ['feed_inlet', 'feed_outlet', 'permeate_outlet']
        for port_str in port_lst:
            port = getattr(m.fs.unit, port_str)
            assert len(port.vars) == 3  # number of state variables for NaCl property package
            assert isinstance(port, Port)

        # test pyomo objects on unit
        #TODO: uncomment relevant vars/constraints once added into model and delete remaining
        unit_objs_type_dict = {'dens_solvent': Param,
                               'A_comp': Var,
                               'B_comp': Var,
                               'flux_mass_phase_comp': Var,
                               'area': Var,
                               'feed_area_cross': Var,
                               'width': Var,
                               'length': Var,
                               'recovery_vol_phase': Var,
                               # 'recovery_mass_phase_comp': Var,
                               # 'rejection_phase_comp': Var,
                               # 'over_pressure_ratio': Var,
                               # 'deltaP': Var,
                               # 'cp_modulus': Var,
                               'mass_transfer_phase_comp': Var,
                               # 'flux_mass_phase_comp_sum': Expression,
                               'eq_mass_transfer_term': Constraint,
                               'eq_permeate_production': Constraint,
                               'eq_flux_mass': Constraint,
                               'eq_connect_mass_transfer': Constraint,
                               # 'eq_connect_enthalpy_transfer': Constraint,
                               'eq_permeate_isothermal': Constraint,
                               'eq_recovery_vol_phase': Constraint,
                               # 'eq_recovery_mass_phase_comp': Constraint,
                               # 'eq_rejection_phase_comp': Constraint,
                               # 'eq_over_pressure_ratio': Constraint,
                               'eq_area': Constraint,
                               'eq_mass_flux_equal_mass_transfer': Constraint,
                               'eq_mass_frac_permeate': Constraint,
                               # 'eq_feed_area_cross': Constraint,
                               # 'eq_permeate_area_cross': Constraint,
                               'eq_permeate_outlet_isothermal': Constraint,
                               'eq_permeate_outlet_isobaric': Constraint
                               }
        for (obj_str, obj_type) in unit_objs_type_dict.items():
            obj = getattr(m.fs.unit, obj_str)
            assert isinstance(obj, obj_type)
        # check that all added unit objects are tested
        for obj in m.fs.unit.component_objects(
                [Param, Var, Expression, Constraint], descend_into=False):
            obj_str = obj.local_name
            if obj_str[0] == '_':
                continue  # do not test hidden references
            assert obj_str in unit_objs_type_dict

        # test feed-side control volume and associated stateblocks
        assert isinstance(m.fs.unit.feed_side, ControlVolume1DBlock)
        cv_stateblock_lst = ['properties']
        for sb_str in cv_stateblock_lst:
            sb = getattr(m.fs.unit.feed_side, sb_str)
            assert isinstance(sb, props.NaClStateBlock)

        stateblock_lst = ['permeate_side', 'permeate_out']
        for sb_str in stateblock_lst:
            sb = getattr(m.fs.unit, sb_str)
            assert isinstance(sb, StateBlock)
            assert isinstance(sb, props.NaClStateBlock)

        # test statistics
        assert number_variables(m) == 777
        assert number_total_constraints(m) == 717
        unused_list = unused_variables_set(m)
        [print(i) for i in unused_list]
        assert number_unused_variables(m) == 26  # TODO: vars from property package parameters (old message from 0dRO)
        # unused areas,  (return later)

    @pytest.mark.integration
    def test_units(self, RO_frame):
        m = RO_frame
        #TODO: add assert_units_equivalent for several variables
        assert_units_consistent(m.fs.unit)

    @pytest.mark.unit
    def test_dof(self, RO_frame):
        m = RO_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, RO_frame):
        m = RO_frame

        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))

        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_initialize(self, RO_frame):
        initialization_tester(RO_frame)

    @pytest.mark.component
    def test_var_scaling(self, RO_frame):
        m = RO_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []
#
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
            b.feed_side.properties[0, 0].flow_mass_phase_comp['Liq', j] for j in comp_lst)
        flow_mass_retentate = sum(
            b.feed_side.properties[0, 1].flow_mass_phase_comp['Liq', j] for j in comp_lst)
        flow_mass_permeate = sum(
            b.permeate_out[0].flow_mass_phase_comp['Liq', j] for j in comp_lst)

        assert value(flow_mass_inlet) == pytest.approx(1.0, rel=1e-3)
        assert value(flow_mass_retentate) == pytest.approx(0.6103, rel=1e-3)
        assert value(flow_mass_permeate) == pytest.approx(0.3898, rel=1e-3)

        assert (abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate
                          )) <= 1e-2)

    @pytest.mark.component
    def test_solution(self, RO_frame):
        m = RO_frame
        assert (pytest.approx(71.54, rel=1e-3) ==
                value(m.fs.unit.area))
        assert (pytest.approx(2.262e-3, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, 1, 'Liq', 'H2O']))
        assert (pytest.approx(2.045e-6, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, 1, 'Liq', 'NaCl']))
        assert (pytest.approx(0.3896, rel=1e-3) ==
                value(m.fs.unit.permeate_out[0].flow_mass_phase_comp['Liq', 'H2O']))
        assert (pytest.approx(1.2159e-4, rel=1e-3) ==
                value(m.fs.unit.permeate_out[0].flow_mass_phase_comp['Liq', 'NaCl']))
        assert (pytest.approx(0, rel=1e-3) == value(m.fs.unit.deltaP[0, .5]))
