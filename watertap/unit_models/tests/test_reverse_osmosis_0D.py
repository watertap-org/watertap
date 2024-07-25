#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
from pyomo.environ import ConcreteModel, units as pyunits
from pyomo.network import Port
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    StateBlock,
)
from idaes.core.util.model_statistics import (
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
import idaes.core.util.scaling as iscale
from watertap.core import (
    MembraneChannel0DBlock,
    FrictionFactor,
    ModuleType,
)
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.core.solvers import get_solver
from watertap.unit_models.reverse_osmosis_base import TransportModel

import watertap.property_models.NaCl_prop_pack as props

from watertap.unit_models.tests.unit_test_harness import UnitTestHarness
import pytest

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------


@pytest.mark.unit
def test_default_config_and_build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = ReverseOsmosis0D(property_package=m.fs.properties)

    assert len(m.fs.unit.config) == 15
    config_keys = [
        "dynamic",
        "has_holdup",
        "property_package",
        "property_package_args",
        "material_balance_type",
        "energy_balance_type",
        "momentum_balance_type",
        "concentration_polarization_type",
        "mass_transfer_coefficient",
        "transport_model",
        "module_type",
        "has_pressure_change",
        "pressure_change_type",
        "friction_factor",
        "has_full_reporting",
    ]

    for key in m.fs.unit.config:
        assert key in config_keys

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
    assert m.fs.unit.config.transport_model == TransportModel.SD
    assert m.fs.unit.config.friction_factor == FrictionFactor.default_by_module_type
    assert m.fs.unit.config.module_type == ModuleType.flat_sheet
    assert not m.fs.unit.config.has_full_reporting
    # test ports
    port_lst = ["inlet", "retentate", "permeate"]
    for port_str in port_lst:
        port = getattr(m.fs.unit, port_str)
        assert isinstance(port, Port)
        # number of state variables for NaCl property package
        assert len(port.vars) == 3

    # test feed-side control volume and associated stateblocks
    assert isinstance(m.fs.unit.feed_side, MembraneChannel0DBlock)
    assert isinstance(m.fs.unit.permeate_side, StateBlock)
    assert isinstance(m.fs.unit.mixed_permeate, StateBlock)
    assert (
        str(m.fs.unit.permeate_side.index_set())
        == "fs._time*fs.unit.feed_side.length_domain"
    )
    # test statistics
    assert number_variables(m) == 135
    assert number_total_constraints(m) == 110
    assert number_unused_variables(m) == 2


def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = ReverseOsmosis0D(
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
    membrane_pressure_drop = 3e5
    membrane_area = 50
    A = 4.2e-12
    B = 3.5e-8
    pressure_atmospheric = 101325
    concentration_polarization_modulus = 1.1

    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.deltaP.fix(-membrane_pressure_drop)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
    m.fs.unit.feed_side.cp_modulus.fix(concentration_polarization_modulus)

    # Set scaling factors for badly scaled variables
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestReverseOsmosis0D(UnitTestHarness):
    def configure(self):
        m = build()

        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            0.004721771
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            1.5757670e-6
        )
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.23608853
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 7.8788350e-5
        self.unit_solutions[m.fs.unit.feed_side.cp_modulus[0, 0, "NaCl"]] = 1.1
        self.unit_solutions[m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"]] = 1.1
        self.unit_solutions[m.fs.unit.deltaP[0]] = -3e5

        comp_lst = ["NaCl", "H2O"]

        flow_mass_inlet = sum(
            m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            m.fs.unit.feed_side.properties_out[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.feed_side.properties_interface[
                    0, 0.0
                ].conc_mass_phase_comp["Liq", "NaCl"]
                / m.fs.unit.feed_side.properties_in[0].conc_mass_phase_comp[
                    "Liq", "NaCl"
                ],
                "out": m.fs.unit.feed_side.cp_modulus[0, 0.0, "NaCl"],
            },
            "Check 2": {
                "in": m.fs.unit.feed_side.properties_interface[
                    0, 1
                ].conc_mass_phase_comp["Liq", "NaCl"]
                / m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp[
                    "Liq", "NaCl"
                ],
                "out": m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"],
            },
            "Check 3": {
                "in": flow_mass_inlet,
                "out": flow_mass_retentate + flow_mass_permeate,
            },
        }

        return m


def build_SKK():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = ReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.fixed,
        mass_transfer_coefficient=MassTransferCoefficient.none,
        transport_model=TransportModel.SKK,
    )

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
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.deltaP.fix(-membrane_pressure_drop)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.reflect_coeff.fix(0.9)
    m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
    m.fs.unit.feed_side.cp_modulus.fix(concentration_polarization_modulus)

    # Set scaling factors for badly scaled variables
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestReverseOsmosis0D_SKK(UnitTestHarness):
    def configure(self):
        m = build_SKK()

        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            0.006710409
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            3.1094137e-5
        )
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.33552046
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 0.001554707
        self.unit_solutions[m.fs.unit.feed_side.cp_modulus[0, 0, "NaCl"]] = 1.1
        self.unit_solutions[m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"]] = 1.1
        self.unit_solutions[m.fs.unit.deltaP[0]] = -3e5

        comp_lst = ["NaCl", "H2O"]

        flow_mass_inlet = sum(
            m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            m.fs.unit.feed_side.properties_out[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.feed_side.properties_interface[
                    0, 0.0
                ].conc_mass_phase_comp["Liq", "NaCl"]
                / m.fs.unit.feed_side.properties_in[0].conc_mass_phase_comp[
                    "Liq", "NaCl"
                ],
                "out": m.fs.unit.feed_side.cp_modulus[0, 0.0, "NaCl"],
            },
            "Check 2": {
                "in": m.fs.unit.feed_side.properties_interface[
                    0, 1
                ].conc_mass_phase_comp["Liq", "NaCl"]
                / m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp[
                    "Liq", "NaCl"
                ],
                "out": m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"],
            },
            "Check 3": {
                "in": flow_mass_inlet,
                "out": flow_mass_retentate + flow_mass_permeate,
            },
        }
        return m


def build_kf_fixed():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = ReverseOsmosis0D(
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
    membrane_pressure_drop = 3e5
    membrane_area = 50
    A = 4.2e-12
    B = 3.5e-8
    pressure_atmospheric = 101325
    kf = 2e-5

    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.deltaP.fix(-membrane_pressure_drop)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
    m.fs.unit.feed_side.K[0, 0.0, "NaCl"].fix(kf)
    m.fs.unit.feed_side.K[0, 1.0, "NaCl"].fix(kf)

    # Set scaling factors for badly scaled variables
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestReverseOsmosis0D_kf_fixed(UnitTestHarness):
    def configure(self):
        m = build_kf_fixed()

        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            0.0038152554
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            1.6673684e-6
        )
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.19076277
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 8.3368422e-5
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 46.07121447
        self.unit_solutions[
            m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 44.344015627
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 50.20324317

        comp_lst = ["NaCl", "H2O"]

        flow_mass_inlet = sum(
            m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            m.fs.unit.feed_side.properties_out[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.feed_side.properties_interface[
                    0, 0.0
                ].conc_mass_phase_comp["Liq", "NaCl"]
                / m.fs.unit.feed_side.properties_in[0].conc_mass_phase_comp[
                    "Liq", "NaCl"
                ],
                "out": m.fs.unit.feed_side.cp_modulus[0, 0.0, "NaCl"],
            },
            "Check 2": {
                "in": m.fs.unit.feed_side.properties_interface[
                    0, 1
                ].conc_mass_phase_comp["Liq", "NaCl"]
                / m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp[
                    "Liq", "NaCl"
                ],
                "out": m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"],
            },
            "Check 3": {
                "in": flow_mass_inlet,
                "out": flow_mass_retentate + flow_mass_permeate,
            },
        }

        return m


def build_kf_calculated():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = ReverseOsmosis0D(
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

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.deltaP.fix(-membrane_pressure_drop)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)

    m.fs.unit.feed_side.channel_height.fix(0.002)
    m.fs.unit.feed_side.spacer_porosity.fix(0.75)
    m.fs.unit.length.fix(length)

    # Set scaling factors for badly scaled variables
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestReverseOsmosis0D_kf_calculated(UnitTestHarness):
    def configure(self):
        m = build_kf_calculated()

        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            0.00456244
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            1.5926761e-6
        )
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.22812202
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 7.96338071e-5
        self.unit_solutions[
            m.fs.unit.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 35.7511
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 41.95562266
        self.unit_solutions[
            m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 46.56687242
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 49.9360425

        comp_lst = ["NaCl", "H2O"]

        flow_mass_inlet = sum(
            m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            m.fs.unit.feed_side.properties_out[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.feed_side.properties_interface[
                    0, 0.0
                ].conc_mass_phase_comp["Liq", "NaCl"]
                / m.fs.unit.feed_side.properties_in[0].conc_mass_phase_comp[
                    "Liq", "NaCl"
                ],
                "out": m.fs.unit.feed_side.cp_modulus[0, 0.0, "NaCl"],
            },
            "Check 2": {
                "in": m.fs.unit.feed_side.properties_interface[
                    0, 1
                ].conc_mass_phase_comp["Liq", "NaCl"]
                / m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp[
                    "Liq", "NaCl"
                ],
                "out": m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"],
            },
            "Check 3": {
                "in": flow_mass_inlet,
                "out": flow_mass_retentate + flow_mass_permeate,
            },
        }
        return m


def build_p_drop_calculation():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = ReverseOsmosis0D(
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

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
    m.fs.unit.feed_side.channel_height.fix(0.001)
    m.fs.unit.feed_side.spacer_porosity.fix(0.97)
    m.fs.unit.length.fix(16)

    # Set scaling factors for badly scaled variables
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestReverseOsmosis0D_p_drop_calculation(UnitTestHarness):
    def configure(self):
        m = build_p_drop_calculation()

        self.unit_solutions[m.fs.unit.deltaP[0]] = -1.661080199e5
        self.unit_solutions[m.fs.unit.feed_side.N_Re[0, 0]] = 395.840743
        self.unit_solutions[m.fs.unit.feed_side.velocity[0, 0]] = 0.2360863
        self.unit_solutions[m.fs.unit.feed_side.N_Re[0, 1]] = 191.1195577
        self.unit_solutions[m.fs.unit.feed_side.velocity[0, 1]] = 0.1187048845
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            0.007089
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            2.1880044e-6
        )
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.134691
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 4.157208e-5
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 50.075075
        self.unit_solutions[
            m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 70.79956388
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 76.3246904

        comp_lst = ["NaCl", "H2O"]

        flow_mass_inlet = sum(
            m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            m.fs.unit.feed_side.properties_out[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.feed_side.properties_interface[
                    0, 0.0
                ].conc_mass_phase_comp["Liq", "NaCl"]
                / m.fs.unit.feed_side.properties_in[0].conc_mass_phase_comp[
                    "Liq", "NaCl"
                ],
                "out": m.fs.unit.feed_side.cp_modulus[0, 0.0, "NaCl"],
            },
            "Check 2": {
                "in": m.fs.unit.feed_side.properties_interface[
                    0, 1
                ].conc_mass_phase_comp["Liq", "NaCl"]
                / m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp[
                    "Liq", "NaCl"
                ],
                "out": m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"],
            },
            "Check 3": {
                "in": flow_mass_inlet,
                "out": flow_mass_retentate + flow_mass_permeate,
            },
        }

        return m


def build_p_drop_calculation_fixed_per_unit_length():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = ReverseOsmosis0D(
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

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)

    m.fs.unit.feed_side.channel_height.fix(0.002)
    m.fs.unit.feed_side.spacer_porosity.fix(0.75)
    m.fs.unit.length.fix(length)
    m.fs.unit.feed_side.dP_dx.fix(-membrane_pressure_drop / length)

    # Set scaling factors for badly scaled variables
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestReverseOsmosis0D_p_drop_fixed_per_unit_length(UnitTestHarness):
    def configure(self):
        m = build_p_drop_calculation_fixed_per_unit_length()

        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            0.0045624403
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            1.5926761e-6
        )
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.22812202
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 7.9633807e-5
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 41.955623
        self.unit_solutions[
            m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 46.566872
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 49.936042
        self.unit_solutions[m.fs.unit.deltaP[0]] = -3e5

        comp_lst = ["NaCl", "H2O"]

        flow_mass_inlet = sum(
            m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            m.fs.unit.feed_side.properties_out[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.feed_side.properties_interface[
                    0, 0.0
                ].conc_mass_phase_comp["Liq", "NaCl"]
                / m.fs.unit.feed_side.properties_in[0].conc_mass_phase_comp[
                    "Liq", "NaCl"
                ],
                "out": m.fs.unit.feed_side.cp_modulus[0, 0.0, "NaCl"],
            },
            "Check 2": {
                "in": m.fs.unit.feed_side.properties_interface[
                    0, 1
                ].conc_mass_phase_comp["Liq", "NaCl"]
                / m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp[
                    "Liq", "NaCl"
                ],
                "out": m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"],
            },
            "Check 3": {
                "in": flow_mass_inlet,
                "out": flow_mass_retentate + flow_mass_permeate,
            },
        }

        return m


def build_friction_factor_spiral_wound():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = ReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        pressure_change_type=PressureChangeType.calculated,
        module_type=ModuleType.spiral_wound,
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

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
    m.fs.unit.feed_side.channel_height.fix(0.001)
    m.fs.unit.feed_side.spacer_porosity.fix(0.97)
    m.fs.unit.length.fix(16)

    # Set scaling factors for badly scaled variables
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestReverseOsmosis0D_friction_factor_spiral_wound(UnitTestHarness):
    def configure(self):
        m = build_friction_factor_spiral_wound()

        self.unit_solutions[m.fs.unit.deltaP[0]] = -4.68073933e5
        self.unit_solutions[m.fs.unit.feed_side.N_Re[0, 0]] = 791.6814859
        self.unit_solutions[m.fs.unit.feed_side.velocity[0, 0]] = 0.472172596
        self.unit_solutions[m.fs.unit.feed_side.N_Re[0, 1]] = 374.5383911
        self.unit_solutions[m.fs.unit.feed_side.velocity[0, 1]] = 0.2329687
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            0.0072232991
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            2.11255771e-6
        )
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.13724268
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 4.01385965e-5
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 47.435640
        self.unit_solutions[
            m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 72.1598903
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 75.1415897247

        comp_lst = ["NaCl", "H2O"]

        flow_mass_inlet = sum(
            m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            m.fs.unit.feed_side.properties_out[0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.feed_side.properties_interface[
                    0, 0.0
                ].conc_mass_phase_comp["Liq", "NaCl"]
                / m.fs.unit.feed_side.properties_in[0].conc_mass_phase_comp[
                    "Liq", "NaCl"
                ],
                "out": m.fs.unit.feed_side.cp_modulus[0, 0.0, "NaCl"],
            },
            "Check 2": {
                "in": m.fs.unit.feed_side.properties_interface[
                    0, 1
                ].conc_mass_phase_comp["Liq", "NaCl"]
                / m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp[
                    "Liq", "NaCl"
                ],
                "out": m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"],
            },
            "Check 3": {
                "in": flow_mass_inlet,
                "out": flow_mass_retentate + flow_mass_permeate,
            },
        }

        return m


@pytest.mark.unit
def test_RO_dynamic_instantiation():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_set=[0, 1], time_units=pyunits.minute)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = ReverseOsmosis0D(
        dynamic=True,
        has_holdup=True,
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        pressure_change_type=PressureChangeType.calculated,
        module_type=ModuleType.spiral_wound,
    )

    m.fs.unit2 = ReverseOsmosis0D(
        dynamic=True,
        has_holdup=True,
        property_package=m.fs.properties,
        has_pressure_change=False,
        concentration_polarization_type=ConcentrationPolarizationType.none,
        mass_transfer_coefficient=MassTransferCoefficient.calculated.none,
    )
