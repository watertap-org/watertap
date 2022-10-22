###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    Param,
    Var,
    Expression,
    Constraint,
    assert_optimal_termination,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.network import Port
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    ControlVolume0DBlock,
    FlowDirection,
)
from watertap.unit_models.osmotically_assisted_reverse_osmosis_0D import (
    OsmoticallyAssistedReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
import watertap.property_models.NaCl_prop_pack as props

from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)

from watertap.core import MembraneChannel0DBlock

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = OsmoticallyAssistedReverseOsmosis0D(property_package=m.fs.properties)

    assert len(m.fs.unit.config) == 12

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.property_package is m.fs.properties
    assert (
        m.fs.unit.config.concentration_polarization_type
        == ConcentrationPolarizationType.calculated
    )
    assert (
        m.fs.unit.config.mass_transfer_coefficient == MassTransferCoefficient.calculated
    )
    assert m.fs.unit.config.pressure_change_type == PressureChangeType.fixed_per_stage
    assert m.fs.unit.feed_side._flow_direction == FlowDirection.forward
    assert m.fs.unit.permeate_side._flow_direction == FlowDirection.backward


@pytest.mark.unit
def test_option_has_pressure_change():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = OsmoticallyAssistedReverseOsmosis0D(
        property_package=m.fs.properties, has_pressure_change=True
    )

    assert isinstance(m.fs.unit.feed_side.deltaP, Var)


@pytest.mark.unit
def test_option_concentration_polarization_type_fixed():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = OsmoticallyAssistedReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.fixed,
        mass_transfer_coefficient=MassTransferCoefficient.none,
    )

    assert (
        m.fs.unit.config.concentration_polarization_type
        == ConcentrationPolarizationType.fixed
    )
    assert isinstance(m.fs.unit.feed_side.cp_modulus, Var)


@pytest.mark.unit
def test_option_concentration_polarization_type_calculated_kf_fixed():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = OsmoticallyAssistedReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.fixed,
    )

    assert (
        m.fs.unit.config.concentration_polarization_type
        == ConcentrationPolarizationType.calculated
    )
    assert m.fs.unit.config.mass_transfer_coefficient == MassTransferCoefficient.fixed
    assert isinstance(m.fs.unit.feed_side.K, Var)


@pytest.mark.unit
def test_option_concentration_polarization_type_calculated_kf_calculated():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = OsmoticallyAssistedReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
    )

    assert (
        m.fs.unit.config.concentration_polarization_type
        == ConcentrationPolarizationType.calculated
    )
    assert (
        m.fs.unit.config.mass_transfer_coefficient == MassTransferCoefficient.calculated
    )
    assert isinstance(m.fs.unit.feed_side.K, Var)
    assert isinstance(m.fs.unit.permeate_side.K, Var)
    assert isinstance(m.fs.unit.feed_side.channel_height, Var)
    assert isinstance(m.fs.unit.permeate_side.channel_height, Var)
    assert isinstance(m.fs.unit.width, Var)
    assert isinstance(m.fs.unit.length, Var)
    assert isinstance(m.fs.unit.feed_side.dh, Var)
    assert isinstance(m.fs.unit.feed_side.spacer_porosity, Var)
    assert isinstance(m.fs.unit.permeate_side.dh, Var)
    assert isinstance(m.fs.unit.permeate_side.spacer_porosity, Var)
    assert isinstance(m.fs.unit.feed_side.N_Sc, Var)
    assert isinstance(m.fs.unit.feed_side.N_Sh, Var)
    assert isinstance(m.fs.unit.feed_side.N_Re, Var)
    assert isinstance(m.fs.unit.permeate_side.N_Sc, Var)
    assert isinstance(m.fs.unit.permeate_side.N_Sh, Var)
    assert isinstance(m.fs.unit.permeate_side.N_Re, Var)


@pytest.mark.unit
def test_option_pressure_change_calculated():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = OsmoticallyAssistedReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.none,
        mass_transfer_coefficient=MassTransferCoefficient.none,
        pressure_change_type=PressureChangeType.calculated,
    )

    assert (
        m.fs.unit.config.concentration_polarization_type
        == ConcentrationPolarizationType.none
    )
    assert m.fs.unit.config.mass_transfer_coefficient == MassTransferCoefficient.none
    assert m.fs.unit.config.pressure_change_type == PressureChangeType.calculated
    assert isinstance(m.fs.unit.feed_side.deltaP, Var)
    assert isinstance(m.fs.unit.feed_side.channel_height, Var)
    assert isinstance(m.fs.unit.width, Var)
    assert isinstance(m.fs.unit.length, Var)
    assert isinstance(m.fs.unit.feed_side.dh, Var)
    assert isinstance(m.fs.unit.feed_side.spacer_porosity, Var)
    assert isinstance(m.fs.unit.feed_side.N_Re, Var)

    assert isinstance(m.fs.unit.permeate_side.deltaP, Var)
    assert isinstance(m.fs.unit.permeate_side.channel_height, Var)
    assert isinstance(m.fs.unit.permeate_side.dh, Var)
    assert isinstance(m.fs.unit.permeate_side.spacer_porosity, Var)
    assert isinstance(m.fs.unit.permeate_side.N_Re, Var)
    assert isinstance(m.fs.unit.eq_area, Constraint)


class TestOsmoticallyAssistedReverseOsmosis:
    @pytest.fixture(scope="class")
    def RO_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = OsmoticallyAssistedReverseOsmosis0D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.fixed,
            mass_transfer_coefficient=MassTransferCoefficient.none,
        )

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_pressure = 50e5
        feed_temperature = 273.15 + 25
        membrane_pressure_drop = -0.5e5
        membrane_area = 50
        A = 4.2e-12
        B = 1.3e-8
        pressure_atmospheric = 101325
        feed_cp_mod = 1.1
        permeate_cp_mod = 0.9
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )
        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )

        m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.area.fix(membrane_area)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.feed_side.cp_modulus.fix(feed_cp_mod)

        perm_flow_mass = 0.5
        perm_mass_frac_NaCl = 0.005
        perm_mass_frac_H2O = 1 - perm_mass_frac_NaCl
        m.fs.unit.permeate_outlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            perm_flow_mass * perm_mass_frac_NaCl
        )
        m.fs.unit.permeate_outlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            perm_flow_mass * perm_mass_frac_H2O
        )
        m.fs.unit.permeate_outlet.pressure[0].fix(pressure_atmospheric)

        m.fs.unit.permeate_side.cp_modulus.fix(permeate_cp_mod)
        m.fs.unit.feed_side.deltaP.fix(0)
        m.fs.unit.permeate_side.deltaP.fix(0)

        return m

    @pytest.mark.unit
    def test_build(self, RO_frame):
        m = RO_frame

        # test ports
        port_lst = ["feed_inlet", "feed_outlet", "permeate_inlet", "permeate_outlet"]
        for port_str in port_lst:
            port = getattr(m.fs.unit, port_str)
            assert isinstance(port, Port)
            # number of state variables for NaCl property package
            assert len(port.vars) == 3

        # test feed-side control volume and associated stateblocks
        assert isinstance(m.fs.unit.feed_side, MembraneChannel0DBlock)
        assert isinstance(m.fs.unit.permeate_side, MembraneChannel0DBlock)

        # test statistics
        assert number_variables(m) == 139
        assert number_total_constraints(m) == 105
        from idaes.core.util.model_statistics import unused_variables_set

        [print(i) for i in unused_variables_set(m)]
        assert number_unused_variables(m) == 7  # vars from property package parameters

    @pytest.mark.unit
    def test_dof(self, RO_frame):
        m = RO_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, RO_frame):
        m = RO_frame

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
        )
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        for i in badly_scaled_var_generator(m):
            print(i[0].name, i[1])

    @pytest.mark.component
    def test_initialize(self, RO_frame):
        initialization_tester(RO_frame)

    @pytest.mark.component
    def test_var_scaling(self, RO_frame):
        m = RO_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        [print(i[0], i[1]) for i in badly_scaled_var_lst]
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, RO_frame):
        m = RO_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_conservation(self, RO_frame):
        m = RO_frame
        b = m.fs.unit
        comp_lst = ["NaCl", "H2O"]

        feed_flow_mass_inlet = sum(
            b.feed_side.properties_in[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        feed_flow_mass_outlet = sum(
            b.feed_side.properties_out[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        perm_flow_mass_inlet = sum(
            b.permeate_side.properties_in[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        perm_flow_mass_outlet = sum(
            b.permeate_side.properties_out[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )

        assert (
            abs(
                value(
                    feed_flow_mass_inlet
                    + perm_flow_mass_inlet
                    - feed_flow_mass_outlet
                    - perm_flow_mass_outlet
                )
            )
            <= 1e-5
        )

        assert (
            abs(
                value(
                    feed_flow_mass_inlet
                    * b.feed_side.properties_in[0].enth_mass_phase["Liq"]
                    - feed_flow_mass_outlet
                    * b.feed_side.properties_out[0].enth_mass_phase["Liq"]
                    + perm_flow_mass_inlet
                    * b.permeate_side.properties_in[0].enth_mass_phase["Liq"]
                    - perm_flow_mass_outlet
                    * b.permeate_side.properties_out[0].enth_mass_phase["Liq"]
                )
            )
            <= 1e-5
        )

    @pytest.mark.component
    def test_solution(self, RO_frame):
        m = RO_frame
        assert pytest.approx(6.712e-03, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        )
        assert pytest.approx(5.2729e-7, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        )
        assert pytest.approx(0.62939, rel=1e-3) == value(
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        )
        assert pytest.approx(0.03497, rel=1e-3) == value(
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
        )
        assert pytest.approx(
            value(m.fs.unit.feed_side.cp_modulus[0, 0.0, "NaCl"]), rel=1e-3
        ) == value(
            m.fs.unit.feed_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ) / value(
            m.fs.unit.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(
            value(m.fs.unit.feed_side.cp_modulus[0, 1.0, "NaCl"]), rel=1e-3
        ) == value(
            m.fs.unit.feed_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ) / value(
            m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"]
        )

        assert pytest.approx(
            value(m.fs.unit.permeate_side.cp_modulus[0, 0.0, "NaCl"]), rel=1e-3
        ) == value(
            m.fs.unit.permeate_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ) / value(
            m.fs.unit.permeate_side.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(
            value(m.fs.unit.permeate_side.cp_modulus[0, 1.0, "NaCl"]), rel=1e-3
        ) == value(
            m.fs.unit.permeate_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ) / value(
            m.fs.unit.permeate_side.properties_out[0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(0, abs=1e-3) == value(m.fs.unit.feed_side.deltaP[0])
        assert pytest.approx(0, abs=1e-3) == value(m.fs.unit.permeate_side.deltaP[0])

    @pytest.mark.component
    def test_CP_calculation_with_kf_fixed(self):
        """Testing 0D-RO with ConcentrationPolarizationType.calculated option enabled.
        This option makes use of an alternative constraint for the feed-side, membrane-interface concentration.
        Additionally, two more variables are created when this option is enabled: Kf - feed-channel
        mass transfer coefficients at the channel inlet and outlet.
        """
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = OsmoticallyAssistedReverseOsmosis0D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.fixed,
        )

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_pressure = 50e5
        feed_temperature = 273.15 + 25
        membrane_pressure_drop = -0.5e5
        membrane_area = 50
        A = 4.2e-12
        B = 1.3e-8
        pressure_atmospheric = 101325
        kf = 3.15e-5
        kp = 8.15e-5

        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )
        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)

        m.fs.unit.area.fix(membrane_area)

        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.structural_parameter.fix(300e-6)

        perm_flow_mass = 0.5
        perm_mass_frac_NaCl = 0.005
        perm_mass_frac_H2O = 1 - perm_mass_frac_NaCl

        m.fs.unit.permeate_outlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            perm_flow_mass * perm_mass_frac_NaCl
        )
        m.fs.unit.permeate_outlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            perm_flow_mass * perm_mass_frac_H2O
        )
        m.fs.unit.permeate_outlet.pressure[0].fix(pressure_atmospheric)

        m.fs.unit.feed_side.deltaP.fix(membrane_pressure_drop)
        m.fs.unit.permeate_side.deltaP.fix(membrane_pressure_drop)

        m.fs.unit.feed_side.K[0, 0.0, "NaCl"].fix(kf)
        m.fs.unit.feed_side.K[0, 1.0, "NaCl"].fix(kf)
        m.fs.unit.permeate_side.K[0, 0.0, "NaCl"].fix(kp)
        m.fs.unit.permeate_side.K[0, 1.0, "NaCl"].fix(kp)

        # test statistics
        assert number_variables(m) == 146
        assert number_total_constraints(m) == 111
        assert number_unused_variables(m) == 2

        # Test units
        assert_units_consistent(m.fs.unit)

        # test degrees of freedom
        assert degrees_of_freedom(m) == 0

        # test scaling
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
        )
        calculate_scaling_factors(m)

        # check that all variables have scaling factors.
        unscaled_var_list = list(
            unscaled_variables_generator(m.fs.unit, include_fixed=True)
        )
        [print(i) for i in unscaled_var_list]
        assert len(unscaled_var_list) == 0

        # # test initialization
        initialization_tester(m)

        # test variable scaling
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

        # test solve
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        # test solution
        assert pytest.approx(4.8519e-3, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        )
        assert pytest.approx(5.9458e-7, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        )
        assert pytest.approx(0.2549, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(0.00247, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(0.4975, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties_out[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(0.0025, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties_out[0].flow_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(35.7511, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(47.553, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(43.810, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(52.7622, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(9.6197, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(4.9939, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties_out[0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        # TODO: It seems the properties reference is not doing what was intended
        assert value(
            m.fs.unit.permeate_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        ) == value(
            m.fs.unit.permeate_side.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert value(
            m.fs.unit.permeate_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        ) == value(
            m.fs.unit.permeate_side.properties_out[0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(2.51225, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(2.58608, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(1.225, rel=1e-3) == value(
            m.fs.unit.feed_side.cp_modulus[0, 0, "NaCl"]
        )
        assert pytest.approx(1.1095, rel=1e-3) == value(
            m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"]
        )
        assert pytest.approx(0.2611, rel=1e-3) == value(
            m.fs.unit.permeate_side.cp_modulus[0, 0, "NaCl"]
        )
        assert pytest.approx(0.517, rel=1e-3) == value(
            m.fs.unit.permeate_side.cp_modulus[0, 1, "NaCl"]
        )

    @pytest.mark.component
    def test_CP_calculation_with_kf_calculation(self):
        """Testing 0D-RO with ConcentrationPolarizationType.calculated option and MassTransferCoefficient.calculated
        option enabled.
        """
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = OsmoticallyAssistedReverseOsmosis0D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
        )

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
        feed_pressure = 50e5
        feed_temperature = 273.15 + 25
        membrane_pressure_drop = 3e5
        length = 20
        membrane_area = 50
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 101325

        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )
        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.area.fix(membrane_area)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate_inlet.pressure[0].fix(pressure_atmospheric)

        m.fs.unit.feed_side.channel_height.fix(0.002)
        m.fs.unit.feed_side.spacer_porosity.fix(0.75)
        m.fs.unit.length.fix(length)

        # test statistics
        assert number_variables(m) == 143
        assert number_total_constraints(m) == 113
        assert number_unused_variables(m) == 0  # vars from property package parameters

        # Test units
        assert_units_consistent(m.fs.unit)

        # test degrees of freedom
        assert degrees_of_freedom(m) == 0

        # test scaling
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors.
        # TODO: see aforementioned TODO on revisiting scaling and associated testing for property models.
        unscaled_var_list = list(
            unscaled_variables_generator(m.fs.unit, include_fixed=True)
        )
        assert len(unscaled_var_list) == 0

        # test initialization
        initialization_tester(m)

        # test variable scaling
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

        # test solve
        results = solver.solve(m, tee=True)

        # Check for optimal solution
        assert_optimal_termination(results)

        # test solution
        assert pytest.approx(4.562e-3, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        )
        assert pytest.approx(1.593e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        )
        assert pytest.approx(0.2281, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(7.963e-5, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(35.75, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(41.96, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(46.57, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(49.94, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )

    @pytest.mark.component
    def test_Pdrop_calculation(self):
        """Testing 0D-RO with PressureChangeType.calculated option."""
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = OsmoticallyAssistedReverseOsmosis0D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            pressure_change_type=PressureChangeType.calculated,
        )

        # fully specify system
        feed_flow_mass = 1 / 3.6
        feed_mass_frac_NaCl = 0.035
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
        feed_pressure = 70e5
        feed_temperature = 273.15 + 25
        membrane_area = 19
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 101325

        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )
        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.area.fix(membrane_area)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate_inlet.pressure[0].fix(pressure_atmospheric)
        m.fs.unit.feed_side.channel_height.fix(0.001)
        m.fs.unit.feed_side.spacer_porosity.fix(0.97)
        m.fs.unit.length.fix(16)

        # test statistics
        assert number_variables(m) == 149
        assert number_total_constraints(m) == 120
        assert number_unused_variables(m) == 0  # vars from property package parameters

        # Test units
        assert_units_consistent(m.fs.unit)

        # test degrees of freedom
        assert degrees_of_freedom(m) == 0

        # test scaling
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors.
        # TODO: see aforementioned TODO on revisiting scaling and associated testing for property models.
        unscaled_var_list = list(
            unscaled_variables_generator(m.fs.unit, include_fixed=True)
        )
        assert len(unscaled_var_list) == 0

        # test initialization
        initialization_tester(m)

        # test variable scaling
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

        # test solve
        results = solver.solve(m, tee=True)

        # Check for optimal solution
        assert_optimal_termination(results)

        # test solution
        assert pytest.approx(-1.661e5, rel=1e-3) == value(m.fs.unit.deltaP[0])
        assert pytest.approx(-1.038e4, rel=1e-3) == value(
            m.fs.unit.deltaP[0] / m.fs.unit.length
        )
        assert pytest.approx(395.8, rel=1e-3) == value(m.fs.unit.feed_side.N_Re[0, 0.0])
        assert pytest.approx(0.2361, rel=1e-3) == value(
            m.fs.unit.feed_side.velocity[0, 0.0]
        )
        assert pytest.approx(191.1, rel=1e-3) == value(m.fs.unit.feed_side.N_Re[0, 1.0])
        assert pytest.approx(0.1187, rel=1e-3) == value(
            m.fs.unit.feed_side.velocity[0, 1.0]
        )
        assert pytest.approx(7.089e-3, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        )
        assert pytest.approx(2.188e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        )
        assert pytest.approx(0.1347, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(4.157e-5, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(50.08, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(70.80, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(76.32, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )

    @pytest.mark.component
    def test_Pdrop_fixed_per_unit_length(self):
        """Testing 0D-RO with PressureChangeType.fixed_per_unit_length option."""
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = OsmoticallyAssistedReverseOsmosis0D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            pressure_change_type=PressureChangeType.fixed_per_unit_length,
        )

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
        feed_pressure = 50e5
        feed_temperature = 273.15 + 25
        membrane_area = 50
        length = 20
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 101325
        membrane_pressure_drop = 3e5

        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )
        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.area.fix(membrane_area)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate_inlet.pressure[0].fix(pressure_atmospheric)

        m.fs.unit.feed_side.channel_height.fix(0.002)
        m.fs.unit.feed_side.spacer_porosity.fix(0.75)
        m.fs.unit.length.fix(length)
        m.fs.unit.feed_side.dP_dx.fix(-membrane_pressure_drop / length)

        # test statistics
        assert number_variables(m) == 144
        assert number_total_constraints(m) == 114
        assert number_unused_variables(m) == 0

        # Test units
        assert_units_consistent(m.fs.unit)

        # test degrees of freedom
        assert degrees_of_freedom(m) == 0

        # test scaling
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors.
        # TODO: see aforementioned TODO on revisiting scaling and associated testing for property models.
        unscaled_var_list = list(
            unscaled_variables_generator(m.fs.unit, include_fixed=True)
        )
        assert len(unscaled_var_list) == 0

        # test initialization
        initialization_tester(m)

        # test variable scaling
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

        # test solve
        results = solver.solve(m, tee=True)

        # Check for optimal solution
        assert_optimal_termination(results)

        # test solution
        assert pytest.approx(-3.000e5, rel=1e-3) == value(m.fs.unit.deltaP[0])
        assert pytest.approx(4.562e-3, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        )
        assert pytest.approx(1.593e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        )
        assert pytest.approx(0.2281, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(7.963e-5, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(41.96, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(46.57, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(49.94, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )

    @pytest.mark.unit
    def test_report(self, RO_frame):
        RO_frame.fs.unit.report()

    @pytest.mark.unit
    def test_repeated_scaling(self, RO_frame):
        # check repeated scaling does not create badly scaled vars
        calculate_scaling_factors(RO_frame)
        for _ in badly_scaled_var_generator(RO_frame):
            assert False
