##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################

import pytest
from pyomo.environ import (ConcreteModel,
                           TerminationCondition,
                           SolverStatus,
                           value,
                           Var,
                           Constraint)
from pyomo.network import Port
from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType)
from proteuslib.unit_models.reverse_osmosis_0D import ReverseOsmosis0D
import proteuslib.property_models.NaCl_prop_pack as props

from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables)
from idaes.core.util.testing import (get_default_solver,
                                     initialization_tester)
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     unscaled_variables_generator,
                                     unscaled_constraints_generator,
                                     badly_scaled_var_generator)

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

        m.fs.unit = ReverseOsmosis0D(default={
            "property_package": m.fs.properties,
            "has_pressure_change": True})

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_pressure = 50e5
        feed_temperature = 273.15 + 25
        membrane_pressure_drop = 3e5
        membrane_area = 50
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 101325

        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
        m.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].fix(
            feed_flow_mass * feed_mass_frac_NaCl)
        m.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(
            feed_flow_mass * feed_mass_frac_H2O)
        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.deltaP.fix(-membrane_pressure_drop)
        m.fs.unit.area.fix(membrane_area)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
        return m

    @pytest.mark.unit
    def test_config(self, RO_frame):
        m = RO_frame
        # check unit config arguments
        assert len(m.fs.unit.config) == 9

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
        assert not m.fs.unit.config.has_concentration_polarization

    @pytest.mark.unit
    def test_build(self, RO_frame):
        m = RO_frame

        # test ports and state variables
        port_lst = ['inlet', 'retentate', 'permeate']
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

        # test unit objects (including parameters, variables, and constraints)
        unit_objs_lst = ['A_comp', 'B_comp', 'dens_solvent',
                         'flux_mass_phase_comp_in', 'flux_mass_phase_comp_out', 'area',
                         'deltaP', 'mass_transfer_phase_comp', 'flux_mass_phase_comp_avg'
                         ]

        unit_constr_lst = ['eq_mass_transfer_term', 'eq_permeate_production',
                           'eq_flux_in', 'eq_flux_out',
                           'eq_connect_mass_transfer', 'eq_connect_enthalpy_transfer',
                           'eq_permeate_isothermal']
        for obj_str in unit_objs_lst:
            assert hasattr(m.fs.unit, obj_str)
            
        # TODO: test variables, constraints, expression separately (like done in prop pack); make var lists for
        # interface and bulk; implement is_property_constructed for variables RO model needs?

        cv_name = 'feed_side'

        # test feed-side/membrane-interface relational constraints
        cv_constraint_lst = ['eq_concentration_polarization_in', 'eq_concentration_polarization_out',
                             'eq_equal_temp_inter_in', 'eq_equal_temp_inter_out',
                             'eq_equal_pressure_inter_in', 'eq_equal_pressure_inter_out',
                             'eq_equal_flow_vol_inter_in', 'eq_equal_flow_vol_inter_out']

        # test state block objects for feed-side inlet/outlet, bulk and interface
        cv_stateblock_lst = ['properties_in', 'properties_out', 'properties_inter_in', 'properties_inter_out']
        stateblock_objs_lst = \
            ['flow_mass_phase_comp', 'pressure', 'temperature', 'pressure_osm',
             'osm_coeff', 'mass_frac_phase_comp', 'conc_mass_phase_comp',
             'dens_mass_phase', 'enth_mass_phase',
             'eq_pressure_osm', 'eq_osm_coeff', 'eq_mass_frac_phase_comp',
             'eq_conc_mass_phase_comp', 'eq_dens_mass_phase', 'eq_enth_mass_phase']

        # control volume
        assert hasattr(m.fs.unit, cv_name)
        cv_blk = getattr(m.fs.unit, cv_name)

        # feed-side constraints
        for constr_str in cv_constraint_lst:
            # check that constraints exist on control volume
            assert hasattr(cv_blk, constr_str)
            blk_constr = getattr(cv_blk, constr_str)
            # check that listed constraints on control volume are constraints
            assert isinstance(blk_constr, Constraint)

        # stateblocks on feed-side
        for blk_str in cv_stateblock_lst:
            # check that stateblocks exist for inlet and outlet for both the bulk and membrane interface
            assert hasattr(cv_blk, blk_str)
            blk = getattr(cv_blk, blk_str)
            # check that listed attributes exist on each stateblock
            for obj_str in stateblock_objs_lst:
                assert hasattr(blk[0], obj_str)

        # permeate
        assert hasattr(m.fs.unit, 'properties_permeate')
        blk = getattr(m.fs.unit, 'properties_permeate')
        for var_str in stateblock_objs_lst:
            assert hasattr(blk[0], var_str)

        # test statistics
        assert number_variables(m) == 100
        assert number_total_constraints(m) == 73
        assert number_unused_variables(m) == 7  # vars from property package parameters

    @pytest.mark.unit
    def test_dof(self, RO_frame):
        m = RO_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, RO_frame):
        m = RO_frame
        #TODO: seems that I can put any string argument as the first arg in set_default_scaling
        # and nothing goes wrong. Perhaps there should be a test to check that the arg is an existing
        # variable in the property package
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index='H2O')
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index='NaCl')
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
            b.feed_side.properties_in[0].flow_mass_phase_comp['Liq', j] for j in comp_lst)
        flow_mass_retentate = sum(
            b.feed_side.properties_out[0].flow_mass_phase_comp['Liq', j] for j in comp_lst)
        flow_mass_permeate = sum(
            b.properties_permeate[0].flow_mass_phase_comp['Liq', j] for j in comp_lst)

        assert (abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate
                          )) <= 1e-6)

        assert (abs(value(
            flow_mass_inlet * b.feed_side.properties_in[0].enth_mass_phase['Liq']
            - flow_mass_retentate * b.feed_side.properties_out[0].enth_mass_phase['Liq']
            - flow_mass_permeate * b.properties_permeate[0].enth_mass_phase['Liq']
        )) <= 1e-6)

    @pytest.mark.component
    def test_solution(self, RO_frame):
        m = RO_frame
        assert (pytest.approx(5.574e-3, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp_avg[0, 'Liq', 'H2O']))
        assert (pytest.approx(1.491e-6, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp_avg[0, 'Liq', 'NaCl']))
        assert (pytest.approx(0.2787, rel=1e-3) ==
                value(m.fs.unit.properties_permeate[0].flow_mass_phase_comp['Liq', 'H2O']))
        assert (pytest.approx(7.453e-5, rel=1e-3) ==
                value(m.fs.unit.properties_permeate[0].flow_mass_phase_comp['Liq', 'NaCl']))
