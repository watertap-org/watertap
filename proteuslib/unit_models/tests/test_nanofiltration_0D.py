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
                           Constraint,
                           TerminationCondition,
                           SolverStatus,
                           value,
                           Var,
                           Expression,
                           Param)
from pyomo.network import Port
from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        ControlVolume0DBlock)
from proteuslib.unit_models.nanofiltration_0D import NanoFiltration0D, ConcentrationPolarizationType
import proteuslib.property_models.NaCl_prop_pack as props

from idaes.core.util import get_solver
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     unscaled_variables_generator,
                                     unscaled_constraints_generator,
                                     badly_scaled_var_generator)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

class TestNanoFiltration():
    @pytest.fixture(scope="class")
    def NF_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = NanoFiltration0D(default={
            "property_package": m.fs.properties,
            "has_pressure_change": True,
            "concentration_polarization_type": ConcentrationPolarizationType.fixed})

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_pressure = 6e5
        feed_temperature = 273.15 + 25
        membrane_pressure_drop = 1e5
        membrane_area = 50 * feed_flow_mass
        A = 3.77e-11
        B = 4.724e-5
        sigma = 0.28
        pressure_atmospheric = 101325

        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
        m.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].fix(feed_flow_mass * feed_mass_frac_NaCl)
        m.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(feed_flow_mass * feed_mass_frac_H2O)
        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.deltaP.fix(-membrane_pressure_drop)
        m.fs.unit.area.fix(membrane_area)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.sigma.fix(sigma)
        m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
        m.fs.unit.cp_modulus.fix(1.1)
        return m

    @pytest.mark.unit
    def test_config(self, NF_frame):
        m = NF_frame
        # check unit config arguments
        assert len(m.fs.unit.config) == 11

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
        assert m.fs.unit.config.concentration_polarization_type == \
               ConcentrationPolarizationType.fixed
        assert isinstance(m.fs.unit.cp_modulus, Var)

    @pytest.mark.unit
    def test_build(self, NF_frame):
        m = NF_frame

        # test ports and variables
        port_lst = ['inlet', 'retentate', 'permeate']
        port_vars_lst = ['flow_mass_phase_comp', 'pressure', 'temperature']
        for port_str in port_lst:
            port = getattr(m.fs.unit, port_str)
            assert len(port.vars) == 3
            assert isinstance(port, Port)

        # test pyomo objects on unit
        unit_objs_type_dict = {'dens_solvent': Param,
                               'A_comp': Var,
                               'B_comp': Var,
                               'sigma': Var,
                               'flux_mass_io_phase_comp': Var,
                               'avg_conc_mass_io_phase_comp': Var,
                               'area': Var,
                               'cp_modulus': Var,
                               'recovery_vol_phase': Var,
                               'recovery_mass_phase_comp': Var,
                               'rejection_phase_comp': Var,
                               'over_pressure_ratio': Var,
                               'deltaP': Var,
                               'mass_transfer_phase_comp': Var,
                               'flux_mass_phase_comp_avg': Expression,
                               'eq_mass_transfer_term': Constraint,
                               'eq_permeate_production': Constraint,
                               'eq_flux_io': Constraint,
                               'eq_avg_conc_io': Constraint,
                               'eq_connect_mass_transfer': Constraint,
                               'eq_connect_enthalpy_transfer': Constraint,
                               'eq_permeate_isothermal': Constraint,
                               'eq_recovery_vol_phase': Constraint,
                               'eq_recovery_mass_phase_comp': Constraint,
                               'eq_rejection_phase_comp': Constraint,
                               'eq_over_pressure_ratio': Constraint}
        for obj_str, obj_type in unit_objs_type_dict.items():
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
        assert isinstance(m.fs.unit.feed_side, ControlVolume0DBlock)
        cv_stateblock_lst = ['properties_in', 'properties_out',
                             'properties_interface_in', 'properties_interface_out']
        for sb_str in cv_stateblock_lst:
            sb = getattr(m.fs.unit.feed_side, sb_str)
            assert isinstance(sb, props.NaClStateBlock)
        # test objects added to control volume
        cv_objs_type_dict = {'eq_concentration_polarization_io': Constraint,
                             'eq_equal_temp_interface_io': Constraint,
                             'eq_equal_pressure_interface_io': Constraint,
                             'eq_equal_flow_vol_interface_io': Constraint}
        for (obj_str, obj_type) in cv_objs_type_dict.items():
            obj = getattr(m.fs.unit.feed_side, obj_str)
            assert isinstance(obj, obj_type)

        # test permeate-side control volume and associated stateblocks
        assert isinstance(m.fs.unit.permeate_side, ControlVolume0DBlock)
        cv_stateblock_lst = ['properties_in', 'properties_out', 'properties_mixed']
        for sb_str in cv_stateblock_lst:
            sb = getattr(m.fs.unit.permeate_side, sb_str)
            assert isinstance(sb, props.NaClStateBlock)
        # test objects added to control volume
        cv_objs_type_dict = {'eq_mass_frac_permeate_io': Constraint,
                             'eq_temperature_permeate_io': Constraint,
                             'eq_pressure_permeate_io': Constraint,
                             'eq_flow_vol_permeate_io': Constraint}
        for (obj_str, obj_type) in cv_objs_type_dict.items():
            obj = getattr(m.fs.unit.permeate_side, obj_str)
            assert isinstance(obj, obj_type)

        # test statistics
        assert number_variables(m) == 128
        assert number_total_constraints(m) == 99
        assert number_unused_variables(m) == 7  # vars from property package parameters

    @pytest.mark.unit
    def test_dof(self, NF_frame):
        m = NF_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, NF_frame):
        m = NF_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_initialize(self, NF_frame):
        initialization_tester(NF_frame)

    @pytest.mark.component
    def test_var_scaling(self, NF_frame):
        m = NF_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, NF_frame):
        m = NF_frame
        solver.options = {'nlp_scaling_method': 'user-scaling'}
        results = solver.solve(m)

        # Check for optimal solution
        assert results.solver.termination_condition == \
               TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_conservation(self, NF_frame):
        m = NF_frame
        b = m.fs.unit
        comp_lst = ['NaCl', 'H2O']

        flow_mass_inlet = sum(
            b.feed_side.properties_in[0].flow_mass_phase_comp['Liq', j] for j in comp_lst)
        flow_mass_retentate = sum(
            b.feed_side.properties_out[0].flow_mass_phase_comp['Liq', j] for j in comp_lst)
        flow_mass_permeate = sum(
            b.permeate_side.properties_mixed[0].flow_mass_phase_comp['Liq', j] for j in comp_lst)

        assert (abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate
                          )) <= 1e-6)

        assert (abs(value(
            flow_mass_inlet * b.feed_side.properties_in[0].enth_mass_phase['Liq']
            - flow_mass_retentate * b.feed_side.properties_out[0].enth_mass_phase['Liq']
            - flow_mass_permeate * b.permeate_side.properties_mixed[0].enth_mass_phase['Liq']
        )) <= 1e-6)

    @pytest.mark.component
    def test_solution(self, NF_frame):
        m = NF_frame
        assert (pytest.approx(1.055e-2, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp_avg[0, 'Liq', 'H2O']))
        assert (pytest.approx(3.563e-4, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp_avg[0, 'Liq', 'NaCl']))
        assert (pytest.approx(0.5276, rel=1e-3) ==
                value(m.fs.unit.permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']))
        assert (pytest.approx(1.782e-2, rel=1e-3) ==
                value(m.fs.unit.permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl']))
