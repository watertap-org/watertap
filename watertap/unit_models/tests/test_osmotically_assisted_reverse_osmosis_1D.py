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
    Var,
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
    FlowDirection,
)
from watertap.unit_models.osmotically_assisted_reverse_osmosis_1D import (
    OsmoticallyAssistedReverseOsmosis1D,
)
from watertap.unit_models.reverse_osmosis_base import TransportModel
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
    badly_scaled_var_generator,
)

from watertap.core import (
    MembraneChannel1DBlock,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
    FrictionFactor,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(property_package=m.fs.properties)

    assert len(m.fs.unit.config) == 19

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
    assert m.fs.unit.config.friction_factor == FrictionFactor.flat_sheet
    assert not m.fs.unit.config.has_full_reporting


@pytest.mark.unit
def test_option_has_pressure_change():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
        property_package=m.fs.properties, has_pressure_change=True
    )

    assert isinstance(m.fs.unit.feed_side.deltaP, Var)
    assert isinstance(m.fs.unit.permeate_side.deltaP, Var)
    assert isinstance(m.fs.unit.feed_side.deltaP_stage, Var)
    assert isinstance(m.fs.unit.permeate_side.deltaP_stage, Var)


@pytest.mark.unit
def test_option_concentration_polarization_type_fixed():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
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
    assert isinstance(m.fs.unit.permeate_side.cp_modulus, Var)


@pytest.mark.unit
def test_option_concentration_polarization_type_calculated_kf_fixed():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
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
    assert isinstance(m.fs.unit.permeate_side.K, Var)


@pytest.mark.unit
def test_option_concentration_polarization_type_calculated_kf_calculated():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
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
    assert isinstance(m.fs.unit.feed_side.N_Sc_comp, Var)
    assert isinstance(m.fs.unit.feed_side.N_Sh_comp, Var)
    assert isinstance(m.fs.unit.feed_side.N_Re, Var)
    assert isinstance(m.fs.unit.permeate_side.N_Sc_comp, Var)
    assert isinstance(m.fs.unit.permeate_side.N_Sh_comp, Var)
    assert isinstance(m.fs.unit.permeate_side.N_Re, Var)


@pytest.mark.unit
def test_option_pressure_change_calculated():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
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
    assert isinstance(m.fs.unit.feed_side.deltaP_stage, Var)

    assert isinstance(m.fs.unit.feed_side.channel_height, Var)
    assert isinstance(m.fs.unit.width, Var)
    assert isinstance(m.fs.unit.length, Var)
    assert isinstance(m.fs.unit.feed_side.dh, Var)
    assert isinstance(m.fs.unit.feed_side.spacer_porosity, Var)
    assert isinstance(m.fs.unit.feed_side.N_Re, Var)

    assert isinstance(m.fs.unit.permeate_side.deltaP, Var)
    assert isinstance(m.fs.unit.permeate_side.deltaP_stage, Var)

    assert isinstance(m.fs.unit.permeate_side.channel_height, Var)
    assert isinstance(m.fs.unit.permeate_side.dh, Var)
    assert isinstance(m.fs.unit.permeate_side.spacer_porosity, Var)
    assert isinstance(m.fs.unit.permeate_side.N_Re, Var)
    assert isinstance(m.fs.unit.eq_area, Constraint)


@pytest.mark.unit
def test_option_has_mass_transfer_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
        property_package=m.fs.properties, transport_model=TransportModel.SKK
    )

    assert isinstance(m.fs.unit.reflect_coeff, Var)
    assert isinstance(m.fs.unit.alpha, Var)

    assert value(m.fs.unit.reflect_coeff) == 0.9

    assert pytest.approx(0.9, rel=1e-3) == value(m.fs.unit.reflect_coeff)
    assert pytest.approx(1e8, rel=1e-3) == value(m.fs.unit.alpha)


class TestOsmoticallyAssistedReverseOsmosis:
    @pytest.fixture(scope="class")
    def RO_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.fixed,
            mass_transfer_coefficient=MassTransferCoefficient.none,
            has_full_reporting=True,
        )

        # fully specify system
        feed_flow_mass = 5 / 18
        feed_mass_frac_NaCl = 0.075
        feed_pressure = 65e5
        feed_temperature = 273.15 + 25
        membrane_area = 155
        width = 1.1
        A = 1e-12
        B = 7.7e-8

        feed_cp_mod = 1.05
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
        m.fs.unit.feed_side.cp_modulus.fix(feed_cp_mod)
        m.fs.unit.feed_side.deltaP_stage.fix(0)

        permeate_flow_mass = 0.33 * feed_flow_mass
        permeate_mass_frac_NaCl = 0.1
        permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            permeate_flow_mass * permeate_mass_frac_H2O
        )
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            permeate_flow_mass * permeate_mass_frac_NaCl
        )
        m.fs.unit.permeate_inlet.pressure[0].fix(5e5)
        m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.permeate_side.cp_modulus.fix(permeate_cp_mod)
        m.fs.unit.permeate_side.deltaP_stage.fix(0)

        m.fs.unit.area.fix(membrane_area)
        m.fs.unit.width.fix(width)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)

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
        assert isinstance(m.fs.unit.feed_side, MembraneChannel1DBlock)
        assert isinstance(m.fs.unit.permeate_side, MembraneChannel1DBlock)

        # test statistics
        assert number_variables(m) == 752
        assert number_total_constraints(m) == 683
        assert number_unused_variables(m) == 30

    @pytest.mark.unit
    def test_dof(self, RO_frame):
        m = RO_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, RO_frame):
        m = RO_frame

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
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
            b.feed_side.properties[0, 0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        feed_flow_mass_outlet = sum(
            b.feed_side.properties[0, 1].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        perm_flow_mass_inlet = sum(
            b.permeate_side.properties[0, 1].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        perm_flow_mass_outlet = sum(
            b.permeate_side.properties[0, 0].flow_mass_phase_comp["Liq", j]
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
                    * b.feed_side.properties[0, 0].enth_mass_phase["Liq"]
                    - feed_flow_mass_outlet
                    * b.feed_side.properties[0, 1].enth_mass_phase["Liq"]
                    + perm_flow_mass_inlet
                    * b.permeate_side.properties[0, 1].enth_mass_phase["Liq"]
                    - perm_flow_mass_outlet
                    * b.permeate_side.properties[0, 0].enth_mass_phase["Liq"]
                )
            )
            <= 1e-5
        )

    @pytest.mark.component
    def test_solution(self, RO_frame):
        m = RO_frame
        assert pytest.approx(8.0587e-4, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        )
        assert pytest.approx(4.0855e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        )
        assert pytest.approx(0.1320, rel=1e-3) == value(
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        )
        assert pytest.approx(0.0202, rel=1e-3) == value(
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
        )
        assert pytest.approx(
            value(
                m.fs.unit.feed_side.cp_modulus[
                    0, m.fs.unit.difference_elements.first(), "NaCl"
                ]
            ),
            rel=1e-3,
        ) == value(
            m.fs.unit.feed_side.properties_interface[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ) / value(
            m.fs.unit.feed_side.properties[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(
            value(m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"]), rel=1e-3
        ) == value(
            m.fs.unit.feed_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ) / value(
            m.fs.unit.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        )

        assert pytest.approx(
            value(
                m.fs.unit.permeate_side.cp_modulus[
                    0, m.fs.unit.difference_elements.first(), "NaCl"
                ]
            ),
            rel=1e-3,
        ) == value(
            m.fs.unit.permeate_side.properties_interface[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ) / value(
            m.fs.unit.permeate_side.properties[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(
            value(m.fs.unit.permeate_side.cp_modulus[0, 1.0, "NaCl"]), rel=1e-3
        ) == value(
            m.fs.unit.permeate_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ) / value(
            m.fs.unit.permeate_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(0, abs=1e-3) == value(m.fs.unit.feed_side.deltaP_stage[0])
        assert pytest.approx(0, abs=1e-3) == value(
            m.fs.unit.permeate_side.deltaP_stage[0]
        )

    @pytest.fixture(scope="class")
    def RO_SKK_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.fixed,
            mass_transfer_coefficient=MassTransferCoefficient.none,
            transport_model=TransportModel.SKK,
            has_full_reporting=True,
        )

        # fully specify system
        feed_flow_mass = 5 / 18
        feed_mass_frac_NaCl = 0.075
        feed_pressure = 65e5
        feed_temperature = 273.15 + 25
        membrane_area = 155
        width = 1.1
        A = 1e-12
        B = 7.7e-8

        feed_cp_mod = 1.05
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
        m.fs.unit.feed_side.cp_modulus.fix(feed_cp_mod)
        m.fs.unit.feed_side.deltaP_stage.fix(0)

        permeate_flow_mass = 0.33 * feed_flow_mass
        permeate_mass_frac_NaCl = 0.1
        permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            permeate_flow_mass * permeate_mass_frac_H2O
        )
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            permeate_flow_mass * permeate_mass_frac_NaCl
        )
        m.fs.unit.permeate_inlet.pressure[0].fix(5e5)
        m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.permeate_side.cp_modulus.fix(permeate_cp_mod)
        m.fs.unit.permeate_side.deltaP_stage.fix(0)

        m.fs.unit.area.fix(membrane_area)
        m.fs.unit.width.fix(width)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.reflect_coeff.fix(0.95)

        return m

    @pytest.fixture(scope="class")
    def RO_SKK_calculated_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            pressure_change_type=PressureChangeType.calculated,
            transport_model=TransportModel.SKK,
            has_full_reporting=True,
        )

        # fully specify system
        feed_flow_mass = 5 / 18
        feed_mass_frac_NaCl = 0.075
        feed_pressure = 65e5
        feed_temperature = 273.15 + 25
        membrane_area = 155
        width = 1.1
        A = 1e-12
        B = 7.7e-8

        feed_cp_mod = 1.05
        permeate_cp_mod = 0.9

        membrane_pressure_drop = -0.5e-5

        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )
        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)
        # m.fs.unit.feed_side.deltaP.fix(membrane_pressure_drop)

        permeate_flow_mass = 0.33 * feed_flow_mass
        permeate_mass_frac_NaCl = 0.1
        permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            permeate_flow_mass * permeate_mass_frac_H2O
        )
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            permeate_flow_mass * permeate_mass_frac_NaCl
        )
        m.fs.unit.permeate_inlet.pressure[0].fix(5e5)
        m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)

        # m.fs.unit.permeate_side.deltaP.fix(membrane_pressure_drop)

        m.fs.unit.area.fix(membrane_area)
        m.fs.unit.width.fix(width)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.reflect_coeff.fix(0.95)
        m.fs.unit.structural_parameter.fix(1200e-6)

        m.fs.unit.permeate_side.channel_height.fix(0.002)
        m.fs.unit.permeate_side.spacer_porosity.fix(0.97)
        m.fs.unit.feed_side.channel_height.fix(0.002)
        m.fs.unit.feed_side.spacer_porosity.fix(0.97)
        return m

    @pytest.mark.unit
    def test_skk_build(self, RO_SKK_frame):
        m = RO_SKK_frame

        # test ports
        port_lst = ["feed_inlet", "feed_outlet", "permeate_inlet", "permeate_outlet"]
        for port_str in port_lst:
            port = getattr(m.fs.unit, port_str)
            assert isinstance(port, Port)
            # number of state variables for NaCl property package
            assert len(port.vars) == 3

        # test feed-side control volume and associated stateblocks
        assert isinstance(m.fs.unit.feed_side, MembraneChannel1DBlock)
        assert isinstance(m.fs.unit.permeate_side, MembraneChannel1DBlock)

        # test statistics
        assert number_variables(m) == 754
        assert number_total_constraints(m) == 684
        assert number_unused_variables(m) == 30

    @pytest.mark.unit
    def test_skk_dof(self, RO_SKK_frame):
        m = RO_SKK_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_skk_general_dof(self, RO_SKK_calculated_frame):
        m = RO_SKK_calculated_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_skk_calculate_scaling(self, RO_SKK_frame):
        m = RO_SKK_frame

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
        )
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        for i in badly_scaled_var_generator(m):
            print(i[0].name, i[1])

    @pytest.mark.component
    def test_skk_initialize(self, RO_SKK_frame):
        initialization_tester(RO_SKK_frame)

    @pytest.mark.component
    def test_skk_var_scaling(self, RO_SKK_frame):
        m = RO_SKK_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        [print(i[0], i[1]) for i in badly_scaled_var_lst]
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_skk_solve(self, RO_SKK_frame):
        m = RO_SKK_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_skk_conservation(self, RO_SKK_frame):
        m = RO_SKK_frame
        b = m.fs.unit
        comp_lst = ["NaCl", "H2O"]

        feed_flow_mass_inlet = sum(
            b.feed_side.properties[0, 0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        feed_flow_mass_outlet = sum(
            b.feed_side.properties[0, 1].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        perm_flow_mass_inlet = sum(
            b.permeate_side.properties[0, 1].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        perm_flow_mass_outlet = sum(
            b.permeate_side.properties[0, 0].flow_mass_phase_comp["Liq", j]
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
                    * b.feed_side.properties[0, 0].enth_mass_phase["Liq"]
                    - feed_flow_mass_outlet
                    * b.feed_side.properties[0, 1].enth_mass_phase["Liq"]
                    + perm_flow_mass_inlet
                    * b.permeate_side.properties[0, 1].enth_mass_phase["Liq"]
                    - perm_flow_mass_outlet
                    * b.permeate_side.properties[0, 0].enth_mass_phase["Liq"]
                )
            )
            <= 1e-5
        )

    @pytest.mark.component
    def test_skk_solution(self, RO_SKK_frame):
        m = RO_SKK_frame
        assert pytest.approx(8.5717e-4, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        )
        assert pytest.approx(9.242e-06, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        )
        assert pytest.approx(1.241e-01, rel=1e-3) == value(
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        )
        assert pytest.approx(1.940e-02, rel=1e-3) == value(
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
        )
        assert pytest.approx(
            value(
                m.fs.unit.feed_side.cp_modulus[
                    0, m.fs.unit.difference_elements.first(), "NaCl"
                ]
            ),
            rel=1e-3,
        ) == value(
            m.fs.unit.feed_side.properties_interface[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ) / value(
            m.fs.unit.feed_side.properties[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(
            value(m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"]), rel=1e-3
        ) == value(
            m.fs.unit.feed_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ) / value(
            m.fs.unit.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        )

        assert pytest.approx(
            value(
                m.fs.unit.permeate_side.cp_modulus[
                    0, m.fs.unit.difference_elements.first(), "NaCl"
                ]
            ),
            rel=1e-3,
        ) == value(
            m.fs.unit.permeate_side.properties_interface[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ) / value(
            m.fs.unit.permeate_side.properties[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(
            value(m.fs.unit.permeate_side.cp_modulus[0, 1.0, "NaCl"]), rel=1e-3
        ) == value(
            m.fs.unit.permeate_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ) / value(
            m.fs.unit.permeate_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(0, abs=1e-3) == value(m.fs.unit.feed_side.deltaP_stage[0])
        assert pytest.approx(0, abs=1e-3) == value(
            m.fs.unit.permeate_side.deltaP_stage[0]
        )

    @pytest.mark.component
    def test_skk_CP_calculation_with_kf_calculation(self, RO_SKK_calculated_frame):
        m = RO_SKK_calculated_frame

        # Test units
        assert_units_consistent(m.fs.unit)

        # test degrees of freedom
        assert degrees_of_freedom(m) == 0

        # test scaling
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
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

        assert (
            m.fs.unit.config.concentration_polarization_type
            == ConcentrationPolarizationType.calculated
        )
        assert (
            m.fs.unit.config.mass_transfer_coefficient
            == MassTransferCoefficient.calculated
        )

        assert number_variables(m) == 937
        assert number_total_constraints(m) == 884
        assert number_unused_variables(m) == 23

        assert pytest.approx(7.215e-4, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        )
        assert pytest.approx(8.463e-06, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        )
        assert pytest.approx(0.1451, rel=1e-3) == value(
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        )
        assert pytest.approx(0.01952, rel=1e-3) == value(
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
        )

        assert pytest.approx(78.8775, rel=1e-3) == value(
            m.fs.unit.feed_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        x_interface_in = m.fs.unit.length_domain.at(2)
        assert pytest.approx(88.360, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(128.619, rel=1e-3) == value(
            m.fs.unit.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(132.351, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(107.06, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(74.889, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(52.882, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(27.392, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(
            value(
                m.fs.unit.feed_side.cp_modulus[
                    0, m.fs.unit.difference_elements.first(), "NaCl"
                ]
            ),
            rel=1e-3,
        ) == value(
            m.fs.unit.feed_side.properties_interface[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ) / value(
            m.fs.unit.feed_side.properties[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(
            value(m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"]), rel=1e-3
        ) == value(
            m.fs.unit.feed_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ) / value(
            m.fs.unit.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        )

        assert pytest.approx(
            value(
                m.fs.unit.permeate_side.cp_modulus[
                    0, m.fs.unit.difference_elements.first(), "NaCl"
                ]
            ),
            rel=1e-3,
        ) == value(
            m.fs.unit.permeate_side.properties_interface[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ) / value(
            m.fs.unit.permeate_side.properties[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(
            value(m.fs.unit.permeate_side.cp_modulus[0, 1.0, "NaCl"]), rel=1e-3
        ) == value(
            m.fs.unit.permeate_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ) / value(
            m.fs.unit.permeate_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(-197270.178, abs=1e-1) == value(
            m.fs.unit.feed_side.deltaP_stage[0]
        )
        assert pytest.approx(-109664.455, abs=1e-1) == value(
            m.fs.unit.permeate_side.deltaP_stage[0]
        )

    @pytest.mark.component
    def test_CP_calculation_with_kf_fixed(self):
        """Testing 1D-OARO with ConcentrationPolarizationType.calculated option enabled.
        This option makes use of an alternative constraint for the feed-side, membrane-interface concentration.
        Additionally, feed and permeate-side mass trasnfer coefficients are fixed.
        """
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.fixed,
        )

        # fully specify system
        feed_flow_mass = 5 / 18
        feed_mass_frac_NaCl = 0.075
        feed_pressure = 65e5
        feed_temperature = 273.15 + 25
        membrane_area = 150
        membrane_pressure_drop = -0.5e5
        A = 1e-12
        B = 7.7e-8
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

        permeate_flow_mass = 0.33 * feed_flow_mass
        permeate_mass_frac_NaCl = 0.1
        permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            permeate_flow_mass * permeate_mass_frac_H2O
        )
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            permeate_flow_mass * permeate_mass_frac_NaCl
        )
        m.fs.unit.permeate_inlet.pressure[0].fix(5e5)
        m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)

        m.fs.unit.feed_side.deltaP.fix(membrane_pressure_drop)
        m.fs.unit.permeate_side.deltaP.fix(membrane_pressure_drop)

        m.fs.unit.feed_side.K[0, 0.0, "NaCl"].fix(kf)
        m.fs.unit.feed_side.K[0, 1.0, "NaCl"].fix(kf)
        m.fs.unit.permeate_side.K[0, 0.0, "NaCl"].fix(kp)
        m.fs.unit.permeate_side.K[0, 1.0, "NaCl"].fix(kp)

        # test statistics
        assert number_variables(m) == 785
        assert number_total_constraints(m) == 713
        assert number_unused_variables(m) == 27

        # Test units
        assert_units_consistent(m.fs.unit)

        # test degrees of freedom
        assert degrees_of_freedom(m) == 0

        # test scaling
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
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
        results = solver.solve(m, tee=True)

        # Check for optimal solution
        assert_optimal_termination(results)

        # test solution
        assert pytest.approx(9.1458e-4, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        )
        assert pytest.approx(3.9757e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        )
        assert pytest.approx(0.0825, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties[0, 1].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(0.009167, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties[0, 1].flow_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(0.2197, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties[0, 0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(0.009763, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties[0, 0].flow_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(78.8775, rel=1e-3) == value(
            m.fs.unit.feed_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(159.6307, rel=1e-3) == value(
            m.fs.unit.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        )
        x_interface_in = m.fs.unit.length_domain.at(2)
        assert pytest.approx(89.700, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(160.4275, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(107.06, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(43.7056, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(45.111, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(106.0555, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(1.1, rel=1e-3) == value(
            m.fs.unit.feed_side.cp_modulus[0, 0, "NaCl"]
        )
        assert pytest.approx(1.005, rel=1e-3) == value(
            m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"]
        )
        assert pytest.approx(1.1, rel=1e-3) == value(
            m.fs.unit.permeate_side.cp_modulus[0, 0, "NaCl"]
        )
        assert pytest.approx(0.9906, rel=1e-3) == value(
            m.fs.unit.permeate_side.cp_modulus[0, 1, "NaCl"]
        )

    @pytest.mark.component
    def test_CP_calculation_with_kf_calculation(self):
        """Testing 1D-OARO with ConcentrationPolarizationType.calculated option and MassTransferCoefficient.calculated
        option enabled.
        """
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=False,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
        )

        # fully specify system
        feed_flow_mass = 5 / 18
        feed_mass_frac_NaCl = 0.075
        feed_pressure = 65e5
        feed_temperature = 273.15 + 25
        membrane_area = 150
        width = 1.1
        A = 1e-12
        B = 7.7e-8

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
        m.fs.unit.width.fix(width)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)

        permeate_flow_mass = 0.33 * feed_flow_mass
        permeate_mass_frac_NaCl = 0.1
        permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            permeate_flow_mass * permeate_mass_frac_H2O
        )
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            permeate_flow_mass * permeate_mass_frac_NaCl
        )
        m.fs.unit.permeate_inlet.pressure[0].fix(5e5)
        m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)

        m.fs.unit.structural_parameter.fix(1200e-6)

        m.fs.unit.permeate_side.channel_height.fix(0.002)
        m.fs.unit.permeate_side.spacer_porosity.fix(0.97)
        m.fs.unit.feed_side.channel_height.fix(0.002)
        m.fs.unit.feed_side.spacer_porosity.fix(0.97)

        # # test statistics
        assert number_variables(m) == 867
        assert number_total_constraints(m) == 815
        assert number_unused_variables(m) == 23

        # Test units
        assert_units_consistent(m.fs.unit)

        # test degrees of freedom
        assert degrees_of_freedom(m) == 0

        # test scaling
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors.
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
        assert pytest.approx(6.756e-4, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        )
        assert pytest.approx(4.538e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        )
        assert pytest.approx(78.88, rel=1e-3) == value(
            m.fs.unit.feed_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        x_interface_in = m.fs.unit.length_domain.at(2)
        assert pytest.approx(87.056, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(124.028, rel=1e-3) == value(
            m.fs.unit.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(127.961, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(107.06, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(70.902, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(52.542, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(27.905, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        )

    @pytest.mark.component
    def test_Pdrop_fixed_per_unit_length(self):
        """Testing 1D-OARO with PressureChangeType.fixed_per_unit_length option."""
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            pressure_change_type=PressureChangeType.fixed_per_unit_length,
        )

        # fully specify system
        feed_flow_mass = 5 / 18
        feed_mass_frac_NaCl = 0.075
        feed_pressure = 65e5
        feed_temperature = 273.15 + 25
        membrane_area = 150
        A = 1e-12
        B = 7.7e-8
        feed_pressure_drop = -1.9e5
        perm_pressure_drop = -1.7e5
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

        permeate_flow_mass = 0.33 * feed_flow_mass
        permeate_mass_frac_NaCl = 0.1
        permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            permeate_flow_mass * permeate_mass_frac_H2O
        )
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            permeate_flow_mass * permeate_mass_frac_NaCl
        )
        m.fs.unit.permeate_inlet.pressure[0].fix(5e5)
        m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.structural_parameter.fix(1200e-6)

        m.fs.unit.permeate_side.channel_height.fix(0.002)
        m.fs.unit.permeate_side.spacer_porosity.fix(0.97)
        m.fs.unit.feed_side.channel_height.fix(0.002)
        m.fs.unit.feed_side.spacer_porosity.fix(0.97)

        length = 141
        m.fs.unit.length.fix(length)
        m.fs.unit.feed_side.dP_dx.fix(feed_pressure_drop / length)
        m.fs.unit.permeate_side.dP_dx.fix(perm_pressure_drop / length)

        # test statistics
        assert number_variables(m) == 891
        assert number_total_constraints(m) == 817
        assert number_unused_variables(m) == 25

        # Test units
        assert_units_consistent(m.fs.unit)

        # test degrees of freedom
        assert degrees_of_freedom(m) == 0

        # test scaling
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors.
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
        assert pytest.approx(-1.9e5, rel=1e-3) == value(
            m.fs.unit.feed_side.deltaP_stage[0]
        )
        assert pytest.approx(-1.7e5, rel=1e-3) == value(
            m.fs.unit.permeate_side.deltaP_stage[0]
        )
        assert pytest.approx(408.561, rel=1e-3) == value(
            m.fs.unit.feed_side.N_Re[0, 0.0]
        )
        assert pytest.approx(240.933, rel=1e-3) == value(
            m.fs.unit.feed_side.N_Re[0, 1.0]
        )
        assert pytest.approx(298.053, rel=1e-3) == value(
            m.fs.unit.permeate_side.N_Re[0, 0.0]
        )
        assert pytest.approx(128.761, rel=1e-3) == value(
            m.fs.unit.permeate_side.N_Re[0, 1.0]
        )
        assert pytest.approx(6.740e-4, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        )
        assert pytest.approx(4.513e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        )
        assert pytest.approx(78.878, rel=1e-3) == value(
            m.fs.unit.feed_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        x_interface_in = m.fs.unit.length_domain.at(2)
        assert pytest.approx(87.409, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(123.864, rel=1e-3) == value(
            m.fs.unit.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(127.557, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(107.06, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(72.301, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(52.590, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(27.153, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        )

    @pytest.mark.component
    def test_Pdrop_calculation(self):
        """Testing 1D-OARO with PressureChangeType.calculated option."""
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            pressure_change_type=PressureChangeType.calculated,
        )

        # fully specify system
        feed_flow_mass = 5 / 18
        feed_mass_frac_NaCl = 0.075
        feed_pressure = 65e5
        feed_temperature = 273.15 + 25
        membrane_area = 150
        A = 1e-12
        B = 7.7e-8

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

        permeate_flow_mass = 0.33 * feed_flow_mass
        permeate_mass_frac_NaCl = 0.1
        permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            permeate_flow_mass * permeate_mass_frac_H2O
        )
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            permeate_flow_mass * permeate_mass_frac_NaCl
        )
        m.fs.unit.permeate_inlet.pressure[0].fix(5e5)
        m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.structural_parameter.fix(1200e-6)

        m.fs.unit.permeate_side.channel_height.fix(0.002)
        m.fs.unit.permeate_side.spacer_porosity.fix(0.97)
        m.fs.unit.feed_side.channel_height.fix(0.002)
        m.fs.unit.feed_side.spacer_porosity.fix(0.97)
        m.fs.unit.feed_side.N_Re[0, 0].fix(400)

        # test statistics
        assert number_variables(m) == 935
        assert number_total_constraints(m) == 883
        assert number_unused_variables(m) == 23

        # Test units
        assert_units_consistent(m.fs.unit)

        # test degrees of freedom
        assert degrees_of_freedom(m) == 0

        # test scaling
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors.
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
        assert pytest.approx(400, rel=1e-3) == value(m.fs.unit.feed_side.N_Re[0, 0.0])
        assert pytest.approx(0.1253, rel=1e-3) == value(
            m.fs.unit.feed_side.velocity[0, 0.0]
        )
        assert pytest.approx(237.296, rel=1e-3) == value(
            m.fs.unit.feed_side.N_Re[0, 1.0]
        )
        assert pytest.approx(0.0776, rel=1e-3) == value(
            m.fs.unit.feed_side.velocity[0, 1.0]
        )
        assert pytest.approx(290.356, rel=1e-3) == value(
            m.fs.unit.permeate_side.N_Re[0, 0.0]
        )
        assert pytest.approx(0.0884, rel=1e-3) == value(
            m.fs.unit.permeate_side.velocity[0, 0.0]
        )
        assert pytest.approx(126.063, rel=1e-3) == value(
            m.fs.unit.permeate_side.N_Re[0, 1.0]
        )
        assert pytest.approx(0.04062, rel=1e-3) == value(
            m.fs.unit.permeate_side.velocity[0, 1.0]
        )
        assert pytest.approx(6.681e-4, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        )
        assert pytest.approx(4.471e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        )
        assert pytest.approx(78.878, rel=1e-3) == value(
            m.fs.unit.feed_side.properties[0, 0.0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        x_interface_in = m.fs.unit.length_domain.at(2)
        assert pytest.approx(87.252, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(123.235, rel=1e-3) == value(
            m.fs.unit.feed_side.properties[0, 1.0].conc_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(126.981, rel=1e-3) == value(
            m.fs.unit.feed_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(107.06, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(71.894, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(52.807, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        )
        assert pytest.approx(27.618, rel=1e-3) == value(
            m.fs.unit.permeate_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        )

    water_recovery_list = [0.15, 0.2, 0.3, 0.4, 0.5, 0.55]

    @pytest.mark.parametrize("water_recovery", water_recovery_list)
    @pytest.mark.component
    def test_water_recovery_sweep(self, water_recovery):
        """Testing 1D-OARO with PressureChangeType.calculated option."""
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            pressure_change_type=PressureChangeType.calculated,
        )

        # fully specify system
        feed_flow_mass = 5 / 18
        feed_mass_frac_NaCl = 0.075
        feed_pressure = 65e5
        feed_temperature = 273.15 + 25
        membrane_area = 150
        A = 1e-12
        B = 7.7e-8

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

        permeate_flow_mass = 0.33 * feed_flow_mass
        permeate_mass_frac_NaCl = 0.1
        permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            permeate_flow_mass * permeate_mass_frac_H2O
        )
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            permeate_flow_mass * permeate_mass_frac_NaCl
        )
        m.fs.unit.permeate_inlet.pressure[0].fix(5e5)
        m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.structural_parameter.fix(1200e-6)

        m.fs.unit.permeate_side.channel_height.fix(0.002)
        m.fs.unit.permeate_side.spacer_porosity.fix(0.97)
        m.fs.unit.feed_side.channel_height.fix(0.002)
        m.fs.unit.feed_side.spacer_porosity.fix(0.97)
        m.fs.unit.feed_side.N_Re[0, 0].fix(400)

        # test degrees of freedom
        assert degrees_of_freedom(m) == 0

        # test scaling
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
        )

        calculate_scaling_factors(m)

        # test initialization
        initialization_tester(m)

        m.fs.unit.permeate_inlet.pressure[0].unfix()
        m.fs.unit.area.unfix()
        m.fs.unit.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(water_recovery)
        m.fs.unit.permeate_outlet.pressure[0].fix(1e5)

        # test solve
        results = solver.solve(m, tee=True)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.unit
    def test_report(self, RO_frame):
        RO_frame.fs.unit.report()

    @pytest.mark.unit
    def test_repeated_scaling(self, RO_frame):
        # check repeated scaling does not create badly scaled vars
        calculate_scaling_factors(RO_frame)
        for _ in badly_scaled_var_generator(RO_frame):
            assert False
