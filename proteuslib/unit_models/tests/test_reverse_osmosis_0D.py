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
                           Param,
                           Var,
                           Expression,
                           Constraint)
from pyomo.network import Port
from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        ControlVolume0DBlock)
from proteuslib.unit_models.reverse_osmosis_0D import (ReverseOsmosis0D,
                                                       ConcentrationPolarizationType)
import proteuslib.property_models.NaCl_prop_pack\
    as props

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
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = ReverseOsmosis0D(default={"property_package": m.fs.properties})

    assert len(m.fs.unit.config) == 9

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
    assert m.fs.unit.config.concentration_polarization_type == \
           ConcentrationPolarizationType.none

@pytest.mark.unit
def test_option_has_pressure_change():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = ReverseOsmosis0D(default={
        "property_package": m.fs.properties,
        "has_pressure_change": True})

    assert isinstance(m.fs.unit.feed_side.deltaP, Var)
    assert isinstance(m.fs.unit.deltaP, Var)

@pytest.mark.unit
def test_option_concentration_polarization_type_fixed():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = ReverseOsmosis0D(default={
        "property_package": m.fs.properties,
        "has_pressure_change": True,
        "concentration_polarization_type": ConcentrationPolarizationType.fixed})

    assert m.fs.unit.config.concentration_polarization_type == \
           ConcentrationPolarizationType.fixed
    assert isinstance(m.fs.unit.cp_modulus, Var)

@pytest.mark.unit
def test_option_concentration_polarization_type_calculated():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = ReverseOsmosis0D(default={
        "property_package": m.fs.properties,
        "has_pressure_change": True,
        "concentration_polarization_type": ConcentrationPolarizationType.calculated})

    assert m.fs.unit.config.concentration_polarization_type == \
           ConcentrationPolarizationType.calculated
    assert isinstance(m.fs.unit.Kf, Var)
    assert isinstance(m.fs.unit.Kf, Var)

class TestReverseOsmosis():
    @pytest.fixture(scope="class")
    def RO_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis0D(default={
            "property_package": m.fs.properties,
            "has_pressure_change": True,
            "concentration_polarization_type": ConcentrationPolarizationType.fixed})

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
        concentration_polarization_modulus = 1.1

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
        m.fs.unit.cp_modulus.fix(concentration_polarization_modulus)
        return m

    @pytest.mark.unit
    def test_build(self, RO_frame):
        m = RO_frame

        # test ports
        port_lst = ['inlet', 'retentate', 'permeate']
        for port_str in port_lst:
            port = getattr(m.fs.unit, port_str)
            assert len(port.vars) == 3  # number of state variables for NaCl property package
            assert isinstance(port, Port)

        # test pyomo objects on unit
        unit_objs_type_dict = {'dens_solvent': Param,
                               'A_comp': Var,
                               'B_comp': Var,
                               'flux_mass_phase_comp': Var,
                               'area': Var,
                               'deltaP': Var,
                               'cp_modulus': Var,
                               'mass_transfer_phase_comp': Var,
                               'flux_mass_phase_comp_avg': Expression,
                               'eq_mass_transfer_term': Constraint,
                               'eq_permeate_production': Constraint,
                               'eq_flux': Constraint,
                               'eq_connect_mass_transfer': Constraint,
                               'eq_connect_enthalpy_transfer': Constraint,
                               'eq_permeate_isothermal': Constraint}
        for (obj_str, obj_type) in unit_objs_type_dict.items():
            obj = getattr(m.fs.unit, obj_str)
            assert isinstance(obj, obj_type)
        # check that all added unit objects are tested
        for obj in m.fs.unit.component_objects(
                [Param, Var, Expression, Constraint], descend_into=False):
            obj_str = obj.local_name
            assert obj_str in unit_objs_type_dict

        # test control volume and associated stateblocks
        assert isinstance(m.fs.unit.feed_side, ControlVolume0DBlock)
        cv_stateblock_lst = ['properties_in', 'properties_out',
                             'properties_interface_in', 'properties_interface_out']
        for sb_str in cv_stateblock_lst:
            sb = getattr(m.fs.unit.feed_side, sb_str)
            assert isinstance(sb, props.NaClStateBlock)
        # test objects added to control volume
        cv_objs_type_dict = {'eq_concentration_polarization': Constraint,
                             'eq_equal_temp_interface': Constraint,
                             'eq_equal_pressure_interface': Constraint,
                             'eq_equal_flow_vol_interface': Constraint}
        for (obj_str, obj_type) in cv_objs_type_dict.items():
            obj = getattr(m.fs.unit.feed_side, obj_str)
            assert isinstance(obj, obj_type)

        # test permeate stateblock
        assert isinstance(m.fs.unit.properties_permeate, props.NaClStateBlock)

        # test statistics
        assert number_variables(m) == 93
        assert number_total_constraints(m) == 65
        assert number_unused_variables(m) == 7  # vars from property package parameters

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
        assert (pytest.approx(4.682e-3, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp_avg[0, 'Liq', 'H2O']))
        assert (pytest.approx(1.580e-6, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp_avg[0, 'Liq', 'NaCl']))
        assert (pytest.approx(0.2341, rel=1e-3) ==
                value(m.fs.unit.properties_permeate[0].flow_mass_phase_comp['Liq', 'H2O']))
        assert (pytest.approx(7.901e-5, rel=1e-3) ==
                value(m.fs.unit.properties_permeate[0].flow_mass_phase_comp['Liq', 'NaCl']))
        assert (pytest.approx(value(m.fs.unit.cp_modulus[0, 'NaCl']), rel=1e-3) ==
                value(m.fs.unit.feed_side.properties_interface_in[0].conc_mass_phase_comp['Liq', 'NaCl'])
                / value(m.fs.unit.feed_side.properties_in[0].conc_mass_phase_comp['Liq', 'NaCl']))
        assert (pytest.approx(value(m.fs.unit.cp_modulus[0, 'NaCl']), rel=1e-3) ==
                value(m.fs.unit.feed_side.properties_interface_out[0].conc_mass_phase_comp['Liq', 'NaCl'])
                / value(m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp['Liq', 'NaCl']))

    @pytest.mark.unit
    def test_CP_calculation(self):
        """ Testing 0D-RO with ConcentrationPolarizationType.calculated option enabled.
        This option makes use of an alternative constraint for the feed-side, membrane-interface concentration.
        Additionally, two more variables are created when this option is enabled: Kf_in and Kf_out- feed-channel
        mass transfer coefficients at the channel inlet and outlet, respectively.
        """
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis0D(default={
            "property_package": m.fs.properties,
            "has_pressure_change": True,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated})

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
        kf = 2e-5

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
        m.fs.unit.Kf.fix(kf)

        # test statistics
        assert number_variables(m) == 94
        assert number_total_constraints(m) == 65
        assert number_unused_variables(m) == 7  # vars from property package parameters

        # test degrees of freedom
        assert degrees_of_freedom(m) == 0

        # test scaling
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))
        calculate_scaling_factors(m)

        # check that all variables have scaling factors.
        # TODO: Setting the "include_fixed" arg as True reveals
        #  unscaled vars that weren't being accounted for previously. However, calling the whole block (i.e.,
        #  m) shows that several NaCl property parameters are unscaled. For now, we are just interested in ensuring
        #  unit variables are scaled (hence, calling m.fs.unit) but might need to revisit scaling and associated
        #  testing for property models.

        unscaled_var_list = list(unscaled_variables_generator(m.fs.unit, include_fixed=True))

        # Uncomment the line below to see which vars are unscaled:
        # [print(unscaled_var_list[i]) for i in range(len(unscaled_var_list))]
        assert len(unscaled_var_list) == 0

        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

        # test initialization
        initialization_tester(m)

        # test variable scaling
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

        # test solve
        solver.options = {'nlp_scaling_method': 'user-scaling'}
        results = solver.solve(m)

        # Check for optimal solution
        assert results.solver.termination_condition == \
               TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

        # test solution
        assert (pytest.approx(3.807e-3, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp_avg[0, 'Liq', 'H2O']))
        assert (pytest.approx(1.668e-6, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp_avg[0, 'Liq', 'NaCl']))
        assert (pytest.approx(0.1904, rel=1e-3) ==
                value(m.fs.unit.properties_permeate[0].flow_mass_phase_comp['Liq', 'H2O']))
        assert (pytest.approx(8.342e-5, rel=1e-3) ==
                value(m.fs.unit.properties_permeate[0].flow_mass_phase_comp['Liq', 'NaCl']))
        assert (pytest.approx(35.751, rel=1e-3) ==
                value(m.fs.unit.feed_side.properties_in[0].conc_mass_phase_comp['Liq', 'NaCl']))
        assert (pytest.approx(46.123, rel=1e-3) ==
                value(m.fs.unit.feed_side.properties_interface_in[0].conc_mass_phase_comp['Liq', 'NaCl']))
        assert (pytest.approx(44.321, rel=1e-3) ==
                value(m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp['Liq', 'NaCl']))
        assert (pytest.approx(50.081, rel=1e-3) ==
                value(m.fs.unit.feed_side.properties_interface_out[0].conc_mass_phase_comp['Liq', 'NaCl']))
