#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
from pyomo.environ import ConcreteModel

from idaes.core.solvers import get_solver

from idaes.core import FlowsheetBlock

import idaes.core.util.scaling as iscale

from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)

from watertap.unit_models.reverse_osmosis_base import TransportModel

import watertap.property_models.NaCl_prop_pack as props

from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

from watertap.core import FrictionFactor

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------


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

        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        ] = 0.004721771
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        ] = 1.5757670e-6
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.23608853
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 7.8788350e-5
        self.unit_solutions[m.fs.unit.feed_side.cp_modulus[0, 0, "NaCl"]] = 1.1
        self.unit_solutions[m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"]] = 1.1
        self.unit_solutions[m.fs.unit.deltaP[0]] = -3e5

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

        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        ] = 0.006710409
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        ] = 3.1094137e-5
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.33552046
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 0.001554707
        self.unit_solutions[m.fs.unit.feed_side.cp_modulus[0, 0, "NaCl"]] = 1.1
        self.unit_solutions[m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"]] = 1.1
        self.unit_solutions[m.fs.unit.deltaP[0]] = -3e5

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

        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        ] = 0.0038152554
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        ] = 1.6673684e-6
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

        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        ] = 0.00456244
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        ] = 1.5926761e-6
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
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        ] = 0.007089
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        ] = 2.1880044e-6
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

        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        ] = 0.0045624403
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        ] = 1.5926761e-6
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
        friction_factor=FrictionFactor.spiral_wound,
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

        self.unit_solutions[m.fs.unit.deltaP[0]] = -1.801003299e5
        self.unit_solutions[m.fs.unit.feed_side.N_Re[0, 0]] = 395.840743
        self.unit_solutions[m.fs.unit.feed_side.velocity[0, 0]] = 0.2360863
        self.unit_solutions[m.fs.unit.feed_side.N_Re[0, 1]] = 191.3192807
        self.unit_solutions[m.fs.unit.feed_side.velocity[0, 1]] = 0.1188201
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        ] = 0.00708203
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        ] = 2.18590835e-6
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.13455865
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 4.15322586e-5
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 50.075075
        self.unit_solutions[
            m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 70.73118
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 76.211017

        return m
