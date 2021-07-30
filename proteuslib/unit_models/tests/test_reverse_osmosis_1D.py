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
                           assert_optimal_termination)
from pyomo.util.check_units import (assert_units_consistent,
                                    assert_units_equivalent)
from pyomo.network import Port
from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        ControlVolume1DBlock,
                        StateBlock)
from proteuslib.unit_models.reverse_osmosis_1D import (ReverseOsmosis1D,
                                                       ConcentrationPolarizationType,
                                                       MassTransferCoefficient,
                                                       PressureChangeType)
import proteuslib.property_models.NaCl_prop_pack \
    as props

from idaes.core.util import get_solver
from idaes.core.util.model_statistics import *

from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     unscaled_variables_generator,
                                     unscaled_constraints_generator,
                                     badly_scaled_var_generator,
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

    assert len(m.fs.unit.config) == 17


    #
    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == \
           MaterialBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == \
           MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.property_package is \
           m.fs.properties
    assert m.fs.unit.config.pressure_change_type is PressureChangeType.fixed_per_stage
    assert m.fs.unit.config.concentration_polarization_type is ConcentrationPolarizationType.calculated
    assert m.fs.unit.config.mass_transfer_coefficient is MassTransferCoefficient.calculated
    assert not m.fs.unit.config.has_full_reporting

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
    assert isinstance(m.fs.unit.deltaP_stage, Var)

@pytest.mark.unit
def test_option_concentration_polarization_type_fixed():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = ReverseOsmosis1D(default={
        "property_package": m.fs.properties,
        "has_pressure_change": False,
        "concentration_polarization_type": ConcentrationPolarizationType.fixed,
        "mass_transfer_coefficient": MassTransferCoefficient.none})

    assert m.fs.unit.config.concentration_polarization_type == \
           ConcentrationPolarizationType.fixed
    assert isinstance(m.fs.unit.cp_modulus, Var)

@pytest.mark.unit
def test_option_concentration_polarization_type_calculated_kf_fixed():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = ReverseOsmosis1D(default={
        "property_package": m.fs.properties,
        "has_pressure_change": False,
        "concentration_polarization_type": ConcentrationPolarizationType.calculated,
        "mass_transfer_coefficient": MassTransferCoefficient.fixed})

    assert m.fs.unit.config.concentration_polarization_type == \
           ConcentrationPolarizationType.calculated
    assert m.fs.unit.config.mass_transfer_coefficient == \
           MassTransferCoefficient.fixed
    assert isinstance(m.fs.unit.Kf, Var)

@pytest.mark.unit
def test_option_concentration_polarization_type_calculated_kf_calculated():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = ReverseOsmosis1D(default={
        "property_package": m.fs.properties,
        "has_pressure_change": False,
        "concentration_polarization_type": ConcentrationPolarizationType.calculated,
        "mass_transfer_coefficient": MassTransferCoefficient.calculated})

    assert m.fs.unit.config.concentration_polarization_type == \
           ConcentrationPolarizationType.calculated
    assert m.fs.unit.config.mass_transfer_coefficient == \
           MassTransferCoefficient.calculated
    assert isinstance(m.fs.unit.Kf, Var)
    assert isinstance(m.fs.unit.channel_height, Var)
    assert isinstance(m.fs.unit.dh, Var)
    assert isinstance(m.fs.unit.spacer_porosity, Var)
    assert isinstance(m.fs.unit.N_Sc, Var)
    assert isinstance(m.fs.unit.N_Sh, Var)
    assert isinstance(m.fs.unit.N_Re, Var)


@pytest.mark.unit
def test_option_pressure_change_calculated():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = ReverseOsmosis1D(default={
        "property_package": m.fs.properties,
        "has_pressure_change": True,
        "concentration_polarization_type": ConcentrationPolarizationType.none,
        "mass_transfer_coefficient": MassTransferCoefficient.none,
        "pressure_change_type": PressureChangeType.calculated})

    assert m.fs.unit.config.concentration_polarization_type == \
           ConcentrationPolarizationType.none
    assert m.fs.unit.config.mass_transfer_coefficient == \
           MassTransferCoefficient.none
    assert m.fs.unit.config.pressure_change_type == \
           PressureChangeType.calculated
    assert isinstance(m.fs.unit.feed_side.deltaP, Var)
    assert isinstance(m.fs.unit.deltaP, Var)
    assert isinstance(m.fs.unit.deltaP_stage, Var)
    assert isinstance(m.fs.unit.channel_height, Var)
    assert isinstance(m.fs.unit.dh, Var)
    assert isinstance(m.fs.unit.spacer_porosity, Var)
    assert isinstance(m.fs.unit.N_Re, Var)

class Test_initialization_args():
    @pytest.fixture(scope="class")
    def m(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis1D(default={
            "property_package": m.fs.properties,
            "has_pressure_change": True,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.calculated,
            "pressure_change_type": PressureChangeType.calculated,
            "transformation_scheme": "BACKWARD",
            "transformation_method": "dae.finite_difference",
            "finite_elements": 5,
            "has_full_reporting": True
        })

        # fully specify system
        feed_flow_mass = 1000 / 3600
        feed_mass_frac_NaCl = 0.034283
        feed_pressure = 70e5

        feed_temperature = 273.15 + 25
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 1e5
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].fix(
            feed_flow_mass * feed_mass_frac_NaCl)

        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(
            feed_flow_mass * feed_mass_frac_H2O)

        m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate_outlet.pressure[0].fix(pressure_atmospheric)
        m.fs.unit.N_Re[0, 0].fix(400)
        m.fs.unit.recovery_mass_phase_comp[0, 'Liq', 'H2O'].fix(0.5)
        m.fs.unit.spacer_porosity.fix(0.97)
        m.fs.unit.channel_height.fix(0.001)

        return m

    @pytest.mark.component
    def test_initialize_ignore_dof(self, m):
         # Unfix a variable and check that ignore_dof arg works
        temp = m.fs.unit.N_Re[0, 0].value
        m.fs.unit.N_Re[0, 0].unfix()
        # Next line should pass since we are ignoring non-zero DOF
        m.fs.unit.initialize(ignore_dof=True, fail_on_warning=False)
        # Check that ValueError is thrown when ignore_dof and fail_on_warning are set to do so
        with pytest.raises(ValueError, match="Non-zero degrees of freedom: Degrees of freedom on fs.unit = 1. "
                                             "Fix 1 more variable or set keyword arg to ignore_dof=True"):
            m.fs.unit.initialize(ignore_dof=False, fail_on_warning=True)
        # Refix that variable
        m.fs.unit.N_Re[0, 0].fix(temp)

    @pytest.mark.component
    def test_initialize_fail_on_warning(self, m):
        # Cause a failed initialization and check that ValueError is thrown
        temp = m.fs.unit.feed_inlet.pressure[0].value
        m.fs.unit.feed_inlet.pressure[0].fix(0.1 * temp)
        # Next line should pass since only a warning is given
        m.fs.unit.initialize(ignore_dof=True, fail_on_warning=False)
        # Set fail_on_warning to True
        with pytest.raises(ValueError, match="Initialization Step 2: solve indexed blocks failed. The solver failed to "
                                             "converge to an optimal solution. This suggests that the user provided "
                                             "infeasible inputs or that the model is poorly scaled."):
            m.fs.unit.initialize(ignore_dof=True, fail_on_warning=True)

class TestReverseOsmosis():
    @pytest.fixture(scope="class")
    def RO_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis1D(default={
            "property_package": m.fs.properties,
            "has_pressure_change": True,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.calculated,
            "pressure_change_type": PressureChangeType.calculated,
            "transformation_scheme": "BACKWARD",
            "transformation_method": "dae.finite_difference",
            "finite_elements": 100,
            "has_full_reporting": True
        })

        # fully specify system
        feed_flow_mass = 1000 / 3600
        feed_mass_frac_NaCl = 0.034283
        feed_pressure = 70e5

        feed_temperature = 273.15 + 25
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 1e5
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].fix(
            feed_flow_mass * feed_mass_frac_NaCl)

        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(
            feed_flow_mass * feed_mass_frac_H2O)

        m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate_outlet.pressure[0].fix(pressure_atmospheric)
        m.fs.unit.N_Re[0, 0].fix(400)
        m.fs.unit.recovery_mass_phase_comp[0, 'Liq', 'H2O'].fix(0.5)
        m.fs.unit.spacer_porosity.fix(0.97)
        m.fs.unit.channel_height.fix(0.001)

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
        unit_objs_type_dict = {'dens_solvent': Param,
                               'A_comp': Var,
                               'B_comp': Var,
                               'flux_mass_phase_comp': Var,
                               'area': Var,
                               'area_cross': Var,
                               'width': Var,
                               'length': Var,
                               'recovery_vol_phase': Var,
                               'recovery_mass_phase_comp': Var,
                               'Kf': Var,
                               'channel_height': Var,
                               'spacer_porosity': Var,
                               'dh': Var,
                               'N_Re': Var,
                               'N_Sc': Var,
                               'N_Sh': Var,
                               'deltaP': Var,
                               'deltaP_stage': Var,
                               'velocity': Var,
                               'friction_factor_darcy': Var,
                               'mass_transfer_phase_comp': Var,
                               'nfe': Param,
                               'eq_mass_transfer_term': Constraint,
                               'eq_permeate_production': Constraint,
                               'eq_flux_mass': Constraint,
                               'eq_connect_mass_transfer': Constraint,
                               'eq_feed_isothermal': Constraint,
                               'eq_permeate_isothermal': Constraint,
                               'eq_recovery_vol_phase': Constraint,
                               'eq_recovery_mass_phase_comp': Constraint,
                               'eq_area': Constraint,
                               'eq_mass_flux_equal_mass_transfer': Constraint,
                               'eq_permeate_outlet_isothermal': Constraint,
                               'eq_permeate_outlet_isobaric': Constraint,
                               'eq_Kf': Constraint,
                               'eq_N_Re': Constraint,
                               'eq_N_Sc': Constraint,
                               'eq_N_Sh': Constraint,
                               'eq_area_cross': Constraint,
                               'eq_dh': Constraint,
                               'eq_pressure_drop': Constraint,
                               'eq_velocity': Constraint,
                               'eq_friction_factor_darcy': Constraint,
                               'eq_dP_dx': Constraint,
                               'N_Re_avg': Expression,
                               'Kf_avg': Expression,
                               'flux_mass_phase_comp_avg': Expression
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
        assert number_variables(m) == 5583
        assert number_total_constraints(m) == 5540
        assert number_unused_variables(m) == 20

    @pytest.mark.integration
    def test_units(self, RO_frame):
        m = RO_frame
        assert_units_consistent(m.fs.unit)

    @pytest.mark.unit
    def test_dof(self, RO_frame):
        m = RO_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, RO_frame):
        m = RO_frame

        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e1, index=('Liq', 'H2O'))
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e4, index=('Liq', 'NaCl'))

        set_scaling_factor(m.fs.unit.permeate_side[0, 0].flow_mass_phase_comp['Liq', 'H2O'], 0)
        set_scaling_factor(m.fs.unit.permeate_side[0, 0].flow_mass_phase_comp['Liq', 'NaCl'], 0)

        for x in m.fs.unit.feed_side.length_domain:
            set_scaling_factor(m.fs.unit.mass_transfer_phase_comp[0, x, 'Liq', 'NaCl'], 1e3)

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
        assert_optimal_termination(results)

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

        assert value(flow_mass_inlet) == pytest.approx(1/3.6, rel=1e-3)
        assert value(flow_mass_retentate) == pytest.approx(0.1437, rel=1e-3)
        assert value(flow_mass_permeate) == pytest.approx(0.1341, rel=1e-3)

        assert (abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate
                          )) <= 1e-2)

    @pytest.mark.component
    def test_solution(self, RO_frame):
        m = RO_frame
        x_interface_in = value(1 / m.fs.unit.nfe)
        assert (pytest.approx(-1.5596e5, rel=1e-3) == value(m.fs.unit.deltaP_stage[0]))
        assert (pytest.approx(1.219e-2, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, x_interface_in, 'Liq', 'H2O']))
        assert (pytest.approx(1.733e-6, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, x_interface_in, 'Liq', 'NaCl']))
        assert (pytest.approx(2.500e-3, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, 1, 'Liq', 'H2O']))
        assert (pytest.approx(2.599e-6, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, 1, 'Liq', 'NaCl']))
        assert (pytest.approx(24.001, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp_avg[0, 'Liq', 'H2O'] * 3.6e3))
        assert (pytest.approx(8.0338, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp_avg[0, 'Liq', 'NaCl'] * 3.6e6))
        assert (pytest.approx(0.1340, rel=1e-3) ==
                value(m.fs.unit.permeate_out[0].flow_mass_phase_comp['Liq', 'H2O']))
        assert (pytest.approx(4.4897e-5, rel=1e-3) ==
                value(m.fs.unit.permeate_out[0].flow_mass_phase_comp['Liq', 'NaCl']))
        assert (pytest.approx(274.24, rel=1e-3) == value(m.fs.unit.N_Re_avg[0]))
        assert (pytest.approx(112.00, rel=1e-3) == value(m.fs.unit.Kf_avg[0, 'NaCl'] * 3.6e6))
        assert (pytest.approx(20.1100, rel=1e-3) == value(m.fs.unit.area))
    @pytest.mark.component
    def testReverseOsmosis_basic(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis1D(default={
            "property_package": m.fs.properties,
            "has_pressure_change": False,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "transformation_scheme": "BACKWARD",
            "transformation_method": "dae.finite_difference",
            "finite_elements": 20
            })

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
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate_outlet.pressure[0].fix(pressure_atmospheric)
        m.fs.unit.length.fix(8)
        m.fs.unit.recovery_vol_phase[0, 'Liq'].fix(0.4)

        # test pyomo objects on unit
        unit_objs_type_dict = {'dens_solvent': Param,
                               'A_comp': Var,
                               'B_comp': Var,
                               'flux_mass_phase_comp': Var,
                               'area': Var,
                               'area_cross': Var,
                               'width': Var,
                               'length': Var,
                               'recovery_vol_phase': Var,
                               'recovery_mass_phase_comp': Var,
                               'mass_transfer_phase_comp': Var,
                               'nfe': Param,
                               'eq_mass_transfer_term': Constraint,
                               'eq_permeate_production': Constraint,
                               'eq_flux_mass': Constraint,
                               'eq_connect_mass_transfer': Constraint,
                               'eq_feed_isothermal': Constraint,
                               'eq_permeate_isothermal': Constraint,
                               'eq_recovery_vol_phase': Constraint,
                               'eq_recovery_mass_phase_comp': Constraint,
                               'eq_area': Constraint,
                               'eq_mass_flux_equal_mass_transfer': Constraint,
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

        # test statistics
        assert number_variables(m) == 982
        assert number_total_constraints(m) == 941
        unused_list = unused_variables_set(m)
        assert number_unused_variables(m) == 27

        # Test units
        assert_units_consistent(m.fs.unit)

        # Check degrees of freedom = 0
        assert degrees_of_freedom(m) == 0

        # Scaling
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e4, index=('Liq', 'NaCl'))

        set_scaling_factor(m.fs.unit.mass_transfer_phase_comp[0, 0, 'Liq', 'NaCl'], 1e2)

        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))

        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

        # Test initialization
        initialization_tester(m)

        # Test variable scaling
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

        # Solve
        solver.options = {'nlp_scaling_method': 'user-scaling'}
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        # Check mass conservation
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

        assert (abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
                <= 1e-2)

        # Test solution
        x_interface_in = value(1 / m.fs.unit.nfe)
        assert (pytest.approx(8.089e-3, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, x_interface_in, 'Liq', 'H2O']))
        assert (pytest.approx(1.302e-6, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, x_interface_in, 'Liq', 'NaCl']))
        assert (pytest.approx(104.70, rel=1e-3) ==
                value(m.fs.unit.area))
        assert (pytest.approx(9.891e-4, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, 1, 'Liq', 'H2O']))
        assert (pytest.approx(2.003e-6, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, 1, 'Liq', 'NaCl']))
        assert (pytest.approx(0.3896, rel=1e-3) ==
                value(m.fs.unit.permeate_out[0].flow_mass_phase_comp['Liq', 'H2O']))
        assert (pytest.approx(1.8186e-4, rel=1e-3) ==
                value(m.fs.unit.permeate_out[0].flow_mass_phase_comp['Liq', 'NaCl']))

    @pytest.mark.component
    def testReverseOsmosis_cp_mod_fixed(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis1D(default={
            "property_package": m.fs.properties,
            "has_pressure_change": False,
            "concentration_polarization_type": ConcentrationPolarizationType.fixed,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "transformation_scheme": "BACKWARD",
            "transformation_method": "dae.finite_difference",
            "finite_elements": 20
            })

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
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate_outlet.pressure[0].fix(pressure_atmospheric)
        m.fs.unit.length.fix(8)
        m.fs.unit.recovery_vol_phase[0, 'Liq'].fix(0.4)
        m.fs.unit.cp_modulus.fix(1.1)

        # test pyomo objects on unit
        unit_objs_type_dict = {'dens_solvent': Param,
                               'A_comp': Var,
                               'B_comp': Var,
                               'flux_mass_phase_comp': Var,
                               'area': Var,
                               'area_cross': Var,
                               'width': Var,
                               'length': Var,
                               'recovery_vol_phase': Var,
                               'recovery_mass_phase_comp': Var,
                               'cp_modulus': Var,
                               'mass_transfer_phase_comp': Var,
                               'nfe': Param,
                               'eq_mass_transfer_term': Constraint,
                               'eq_permeate_production': Constraint,
                               'eq_flux_mass': Constraint,
                               'eq_connect_mass_transfer': Constraint,
                               'eq_feed_isothermal': Constraint,
                               'eq_permeate_isothermal': Constraint,
                               'eq_recovery_vol_phase': Constraint,
                               'eq_recovery_mass_phase_comp': Constraint,
                               'eq_area': Constraint,
                               'eq_mass_flux_equal_mass_transfer': Constraint,
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

       # test statistics
        assert number_variables(m) == 1003
        assert number_total_constraints(m) == 941
        assert number_unused_variables(m) == 28

        assert_units_consistent(m.fs.unit)

        assert degrees_of_freedom(m) == 0

        # Scaling
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))

        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

        # Test initialization
        initialization_tester(m)
        #Check for poorly scaled variables
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

        # Solve
        solver.options = {'nlp_scaling_method': 'user-scaling'}
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        # Check mass conservation
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

        assert (abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
                <= 1e-2)

        # Test solution
        x_interface_in = value(1 / m.fs.unit.nfe)
        assert (pytest.approx(6.374e-3, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, x_interface_in, 'Liq', 'H2O']))
        assert (pytest.approx(1.475e-6, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, x_interface_in, 'Liq', 'NaCl']))
        assert (pytest.approx(221.92, rel=1e-3) ==
                value(m.fs.unit.area))
        assert (pytest.approx(3.366e-4, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, 1, 'Liq', 'H2O']))
        assert (pytest.approx(2.052e-6, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, 1, 'Liq', 'NaCl']))
        assert (pytest.approx(0.3896, rel=1e-3) ==
                value(m.fs.unit.permeate_out[0].flow_mass_phase_comp['Liq', 'H2O']))
        assert (pytest.approx(4.271e-4, rel=1e-3) ==
                value(m.fs.unit.permeate_out[0].flow_mass_phase_comp['Liq', 'NaCl']))

    @pytest.mark.component
    def testReverseOsmosis_cp_calculated_kf_fixed(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis1D(default={
            "property_package": m.fs.properties,
            "has_pressure_change": False,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.fixed,
            "transformation_scheme": "BACKWARD",
            "transformation_method": "dae.finite_difference",
            "finite_elements": 20
            })

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
        # m.fs.unit.deltaP.fix(0)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate_outlet.pressure[0].fix(pressure_atmospheric)
        m.fs.unit.length.fix(8)
        m.fs.unit.recovery_vol_phase[0, 'Liq'].fix(0.4)
        m.fs.unit.Kf.fix(2e-5)

        # test pyomo objects on unit
        unit_objs_type_dict = {'dens_solvent': Param,
                               'A_comp': Var,
                               'B_comp': Var,
                               'flux_mass_phase_comp': Var,
                               'area': Var,
                               'area_cross': Var,
                               'width': Var,
                               'length': Var,
                               'recovery_vol_phase': Var,
                               'recovery_mass_phase_comp': Var,
                               'Kf': Var,
                               'mass_transfer_phase_comp': Var,
                               'nfe': Param,
                               'eq_mass_transfer_term': Constraint,
                               'eq_permeate_production': Constraint,
                               'eq_flux_mass': Constraint,
                               'eq_connect_mass_transfer': Constraint,
                               'eq_feed_isothermal': Constraint,
                               'eq_permeate_isothermal': Constraint,
                               'eq_recovery_vol_phase': Constraint,
                               'eq_recovery_mass_phase_comp': Constraint,
                               'eq_area': Constraint,
                               'eq_mass_flux_equal_mass_transfer': Constraint,
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

        # test statistics
        assert number_variables(m) == 1003
        assert number_total_constraints(m) == 941
        assert number_unused_variables(m) == 28

        assert_units_consistent(m.fs.unit)

        assert degrees_of_freedom(m) == 0

        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))
        calculate_scaling_factors(m)

        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

        initialization_tester(m)

        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

        solver.options = {'nlp_scaling_method': 'user-scaling'}
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        # Check mass conservation
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

        assert (abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
                <= 1e-2)

        # Test solution
        x_interface_in = value(1 / m.fs.unit.nfe)
        assert (pytest.approx(4.721e-3, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, x_interface_in, 'Liq', 'H2O']))
        assert (pytest.approx(1.640e-6, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, x_interface_in, 'Liq', 'NaCl']))
        assert (pytest.approx(179.596, rel=1e-3) ==
                value(m.fs.unit.area))
        assert (pytest.approx(6.855e-4, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, 1, 'Liq', 'H2O']))
        assert (pytest.approx(2.029e-6, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, 1, 'Liq', 'NaCl']))
        assert (pytest.approx(0.3896, rel=1e-3) ==
                value(m.fs.unit.permeate_out[0].flow_mass_phase_comp['Liq', 'H2O']))
        assert (pytest.approx(3.393e-4, rel=1e-3) ==
                value(m.fs.unit.permeate_out[0].flow_mass_phase_comp['Liq', 'NaCl']))

    @pytest.mark.component
    def testReverseOsmosis_cp_calculated_kf_calculated(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis1D(default={
            "property_package": m.fs.properties,
            "has_pressure_change": False,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.calculated,
            "transformation_scheme": "BACKWARD",
            "transformation_method": "dae.finite_difference",
            "finite_elements": 20
        })

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
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate_outlet.pressure[0].fix(pressure_atmospheric)
        m.fs.unit.length.fix(8)
        m.fs.unit.recovery_vol_phase[0, 'Liq'].fix(0.4)
        m.fs.unit.spacer_porosity.fix(0.75)
        m.fs.unit.channel_height.fix(0.002)

         # test pyomo objects on unit
        unit_objs_type_dict = {'dens_solvent': Param,
                               'A_comp': Var,
                               'B_comp': Var,
                               'flux_mass_phase_comp': Var,
                               'area': Var,
                               'area_cross': Var,
                               'width': Var,
                               'length': Var,
                               'recovery_vol_phase': Var,
                               'recovery_mass_phase_comp': Var,
                               'Kf': Var,
                               'channel_height': Var,
                               'spacer_porosity': Var,
                               'dh': Var,
                               'N_Re': Var,
                               'N_Sc': Var,
                               'N_Sh': Var,
                               'mass_transfer_phase_comp': Var,
                               'nfe': Param,
                               'eq_mass_transfer_term': Constraint,
                               'eq_permeate_production': Constraint,
                               'eq_flux_mass': Constraint,
                               'eq_connect_mass_transfer': Constraint,
                               'eq_feed_isothermal': Constraint,
                               'eq_permeate_isothermal': Constraint,
                               'eq_recovery_vol_phase': Constraint,
                               'eq_recovery_mass_phase_comp': Constraint,
                               'eq_area': Constraint,
                               'eq_mass_flux_equal_mass_transfer': Constraint,
                               'eq_permeate_outlet_isothermal': Constraint,
                               'eq_permeate_outlet_isobaric': Constraint,
                               'eq_Kf': Constraint,
                               'eq_N_Re': Constraint,
                               'eq_N_Sc': Constraint,
                               'eq_N_Sh': Constraint,
                               'eq_area_cross': Constraint,
                               'eq_dh': Constraint
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

        # test statistics
        assert number_variables(m) == 1111
        assert number_total_constraints(m) == 1068
        assert number_unused_variables(m) == 20

        assert_units_consistent(m.fs.unit)

        assert degrees_of_freedom(m) == 0

        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))

        set_scaling_factor(m.fs.unit.feed_side.properties[0, 1].flow_mass_phase_comp['Liq', 'NaCl'], 1e4)
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

        initialization_tester(m)

        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

        solver.options = {'nlp_scaling_method': 'user-scaling'}
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        # Check mass conservation
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

        assert (abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
                <= 1e-2)

        # Test solution
        x_interface_in = value(1 / m.fs.unit.nfe)
        assert (pytest.approx(4.419e-3, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, x_interface_in, 'Liq', 'H2O']))
        assert (pytest.approx(1.670e-6, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, x_interface_in, 'Liq', 'NaCl']))
        assert (pytest.approx(198.413, rel=1e-3) ==
                value(m.fs.unit.area))
        assert (pytest.approx(6.289e-4, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, 1, 'Liq', 'H2O']))
        assert (pytest.approx(2.034e-6, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, 1, 'Liq', 'NaCl']))
        assert (pytest.approx(0.3896, rel=1e-3) ==
                value(m.fs.unit.permeate_out[0].flow_mass_phase_comp['Liq', 'H2O']))
        assert (pytest.approx(3.788e-4, rel=1e-3) ==
                value(m.fs.unit.permeate_out[0].flow_mass_phase_comp['Liq', 'NaCl']))

    @pytest.mark.component
    def testRO_cp_calculated_kf_calculated_pdrop_fixed_by_dx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis1D(default={
            "property_package": m.fs.properties,
            "has_pressure_change": True,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.calculated,
            "pressure_change_type": PressureChangeType.fixed_per_unit_length,
            "transformation_scheme": "BACKWARD",
            "transformation_method": "dae.finite_difference",
            "finite_elements": 20
        })

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
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate_outlet.pressure[0].fix(pressure_atmospheric)
        m.fs.unit.length.fix(8)
        m.fs.unit.recovery_vol_phase[0, 'Liq'].fix(0.4)
        m.fs.unit.spacer_porosity.fix(0.75)
        m.fs.unit.channel_height.fix(0.002)
        m.fs.unit.deltaP.fix(-0.1e5)

        # test pyomo objects on unit
        unit_objs_type_dict = {'dens_solvent': Param,
                               'A_comp': Var,
                               'B_comp': Var,
                               'flux_mass_phase_comp': Var,
                               'area': Var,
                               'area_cross': Var,
                               'width': Var,
                               'length': Var,
                               'recovery_vol_phase': Var,
                               'recovery_mass_phase_comp': Var,
                               'Kf': Var,
                               'channel_height': Var,
                               'spacer_porosity': Var,
                               'dh': Var,
                               'N_Re': Var,
                               'N_Sc': Var,
                               'N_Sh': Var,
                               'deltaP': Var,
                               'deltaP_stage': Var,
                               'mass_transfer_phase_comp': Var,
                               'nfe': Param,
                               'eq_mass_transfer_term': Constraint,
                               'eq_permeate_production': Constraint,
                               'eq_flux_mass': Constraint,
                               'eq_connect_mass_transfer': Constraint,
                               'eq_feed_isothermal': Constraint,
                               'eq_permeate_isothermal': Constraint,
                               'eq_recovery_vol_phase': Constraint,
                               'eq_recovery_mass_phase_comp': Constraint,
                               'eq_area': Constraint,
                               'eq_mass_flux_equal_mass_transfer': Constraint,
                               'eq_permeate_outlet_isothermal': Constraint,
                               'eq_permeate_outlet_isobaric': Constraint,
                               'eq_Kf': Constraint,
                               'eq_N_Re': Constraint,
                               'eq_N_Sc': Constraint,
                               'eq_N_Sh': Constraint,
                               'eq_area_cross': Constraint,
                               'eq_dh': Constraint,
                               'eq_pressure_drop': Constraint
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

        # test statistics
        assert number_variables(m) == 1133
        assert number_total_constraints(m) == 1069
        unused_list = unused_variables_set(m)
        assert number_unused_variables(m) == 21

        assert_units_consistent(m.fs.unit)

        assert degrees_of_freedom(m) == 0

        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))

        set_scaling_factor(m.fs.unit.feed_side.properties[0, 1].flow_mass_phase_comp['Liq', 'NaCl'], 1e4)
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

        initialization_tester(m)

        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

        solver.options = {'nlp_scaling_method': 'user-scaling'}
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        # Check mass conservation
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

        assert (abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate)) <= 1e-2)

        # Test solution
        x_interface_in = value(1 / m.fs.unit.nfe)
        assert (pytest.approx(-8.000e4, rel=1e-3) == value(m.fs.unit.deltaP_stage[0]))
        assert (pytest.approx(4.343e-3, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, x_interface_in, 'Liq', 'H2O']))
        assert (pytest.approx(1.676e-6, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, x_interface_in, 'Liq', 'NaCl']))
        assert (pytest.approx(211.30, rel=1e-3) ==
                value(m.fs.unit.area))
        assert (pytest.approx(5.488e-4, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, 1, 'Liq', 'H2O']))
        assert (pytest.approx(2.008e-6, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, 1, 'Liq', 'NaCl']))
        assert (pytest.approx(0.3896, rel=1e-3) ==
                value(m.fs.unit.permeate_out[0].flow_mass_phase_comp['Liq', 'H2O']))
        assert (pytest.approx(4.021e-4, rel=1e-3) ==
                value(m.fs.unit.permeate_out[0].flow_mass_phase_comp['Liq', 'NaCl']))

    @pytest.mark.component
    def testRO_cp_calculated_kf_calculated_pdrop_fixed_by_stage(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis1D(default={
            "property_package": m.fs.properties,
            "has_pressure_change": True,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.calculated,
            "pressure_change_type": PressureChangeType.fixed_per_stage,
            "transformation_scheme": "BACKWARD",
            "transformation_method": "dae.finite_difference",
            "finite_elements": 20
        })

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
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate_outlet.pressure[0].fix(pressure_atmospheric)
        m.fs.unit.length.fix(8)
        m.fs.unit.recovery_vol_phase[0, 'Liq'].fix(0.4)
        m.fs.unit.spacer_porosity.fix(0.75)
        m.fs.unit.channel_height.fix(0.002)
        m.fs.unit.deltaP_stage.fix(-62435.6)

        # test pyomo objects on unit
        unit_objs_type_dict = {'dens_solvent': Param,
                               'A_comp': Var,
                               'B_comp': Var,
                               'flux_mass_phase_comp': Var,
                               'area': Var,
                               'area_cross': Var,
                               'width': Var,
                               'length': Var,
                               'recovery_vol_phase': Var,
                               'recovery_mass_phase_comp': Var,
                               'Kf': Var,
                               'channel_height': Var,
                               'spacer_porosity': Var,
                               'dh': Var,
                               'N_Re': Var,
                               'N_Sc': Var,
                               'N_Sh': Var,
                               'deltaP': Var,
                               'deltaP_stage': Var,
                               'mass_transfer_phase_comp': Var,
                               'nfe': Param,
                               'eq_mass_transfer_term': Constraint,
                               'eq_permeate_production': Constraint,
                               'eq_flux_mass': Constraint,
                               'eq_connect_mass_transfer': Constraint,
                               'eq_feed_isothermal': Constraint,
                               'eq_permeate_isothermal': Constraint,
                               'eq_recovery_vol_phase': Constraint,
                               'eq_recovery_mass_phase_comp': Constraint,
                               'eq_area': Constraint,
                               'eq_mass_flux_equal_mass_transfer': Constraint,
                               'eq_permeate_outlet_isothermal': Constraint,
                               'eq_permeate_outlet_isobaric': Constraint,
                               'eq_Kf': Constraint,
                               'eq_N_Re': Constraint,
                               'eq_N_Sc': Constraint,
                               'eq_N_Sh': Constraint,
                               'eq_area_cross': Constraint,
                               'eq_dh': Constraint,
                               'eq_pressure_drop': Constraint
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

        # test statistics
        assert number_variables(m) == 1133
        assert number_total_constraints(m) == 1089
        unused_list = unused_variables_set(m)
        assert number_unused_variables(m) == 20

        assert_units_consistent(m.fs.unit)

        assert degrees_of_freedom(m) == 0

        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))

        set_scaling_factor(m.fs.unit.feed_side.properties[0, 1].flow_mass_phase_comp['Liq', 'NaCl'], 1e4)
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

        initialization_tester(m)

        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

        solver.options = {'nlp_scaling_method': 'user-scaling'}
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        # Check mass conservation
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

        assert (abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
                <= 1e-2)

        # Test solution
        x_interface_in = value(1 / m.fs.unit.nfe)
        assert (pytest.approx(-6.2436e4, rel=1e-3) == value(m.fs.unit.deltaP_stage[0]))
        assert (pytest.approx(4.360e-3, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, x_interface_in, 'Liq', 'H2O']))
        assert (pytest.approx(1.675e-6, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, x_interface_in, 'Liq', 'NaCl']))
        assert (pytest.approx(208.378, rel=1e-3) ==
                value(m.fs.unit.area))
        assert (pytest.approx(5.656e-4, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, 1, 'Liq', 'H2O']))
        assert (pytest.approx(2.013e-6, rel=1e-3) ==
                value(m.fs.unit.flux_mass_phase_comp[0, 1, 'Liq', 'NaCl']))
        assert (pytest.approx(0.3896, rel=1e-3) ==
                value(m.fs.unit.permeate_out[0].flow_mass_phase_comp['Liq', 'H2O']))
        assert (pytest.approx(3.969e-4, rel=1e-3) ==
                value(m.fs.unit.permeate_out[0].flow_mass_phase_comp['Liq', 'NaCl']))

