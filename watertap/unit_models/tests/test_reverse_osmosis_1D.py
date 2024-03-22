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
from pyomo.environ import ConcreteModel

from idaes.core.solvers import get_solver

from idaes.core import FlowsheetBlock

import idaes.core.util.scaling as iscale

from watertap.unit_models.reverse_osmosis_base import TransportModel

import watertap.property_models.NaCl_prop_pack as props

from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

from watertap.core import FrictionFactor


from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()
# -----------------------------------------------------------------------------


def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        pressure_change_type=PressureChangeType.calculated,
        finite_elements=3,
        has_full_reporting=True,
    )

    # fully specify system
    feed_flow_mass = 1000 / 3600
    feed_mass_frac_NaCl = 0.034283
    feed_pressure = 70e5

    feed_temperature = 273.15 + 25
    A = 4.2e-12
    B = 3.5e-8
    pressure_atmospheric = 1e5
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )

    m.fs.unit.inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
    m.fs.unit.feed_side.N_Re[0, 0].fix(400)
    m.fs.unit.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(0.5)
    m.fs.unit.feed_side.spacer_porosity.fix(0.97)
    m.fs.unit.feed_side.channel_height.fix(0.001)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )

    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestReverseOsmosis1D(UnitTestHarness):
    def configure(self):
        m = build()
        x_interface_in = m.fs.unit.feed_side.length_domain.at(2)

        self.unit_solutions[m.fs.unit.deltaP[0]] = -1.75535e5
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "H2O"]
        ] = 0.00799178
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "NaCl"]
        ] = 2.113811e-6
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        ] = 0.00248645
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "NaCl"]
        ] = 2.593074e-6
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        ] = 0.005036466
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        ] = 2.372961e-6
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.1341274
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 6.319489e-5
        self.unit_solutions[m.fs.unit.feed_side.N_Re_avg[0]] = 371.01255
        self.unit_solutions[m.fs.unit.feed_side.K_avg[0, "NaCl"]] = 2.985544e-5
        self.unit_solutions[m.fs.unit.area] = 26.63124

        return m


def build_basic():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=False,
        concentration_polarization_type=ConcentrationPolarizationType.none,
        mass_transfer_coefficient=MassTransferCoefficient.none,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        finite_elements=3,
    )

    # fully specify system
    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.035
    feed_pressure = 50e5
    feed_temperature = 273.15 + 25
    A = 4.2e-12
    B = 3.5e-8
    pressure_atmospheric = 101325
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )

    m.fs.unit.inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
    m.fs.unit.length.fix(8)
    m.fs.unit.recovery_vol_phase[0, "Liq"].fix(0.4)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )

    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestReverseOsmosis1D_basic(UnitTestHarness):
    def configure(self):
        m = build_basic()
        x_interface_in = m.fs.unit.feed_side.length_domain.at(2)

        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "H2O"]
        ] = 0.00484064
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "NaCl"]
        ] = 1.628559e-6
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        ] = 0.001019483
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "NaCl"]
        ] = 1.99996e-6
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.3895746
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 0.000265183
        self.unit_solutions[m.fs.unit.area] = 144.307206

        return m


def build_SKK():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        pressure_change_type=PressureChangeType.calculated,
        transport_model=TransportModel.SKK,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        finite_elements=3,
        has_full_reporting=True,
    )

    # fully specify system
    feed_flow_mass = 1000 / 3600
    feed_mass_frac_NaCl = 0.034283
    feed_pressure = 70e5

    feed_temperature = 273.15 + 25
    A = 4.2e-12
    B = 3.5e-8
    pressure_atmospheric = 1e5
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )

    m.fs.unit.inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.reflect_coeff.fix(0.9)
    m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
    m.fs.unit.feed_side.N_Re[0, 0].fix(400)
    m.fs.unit.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(0.5)
    m.fs.unit.feed_side.spacer_porosity.fix(0.97)
    m.fs.unit.feed_side.channel_height.fix(0.001)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )

    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestReverseOsmosis1D_SKK(UnitTestHarness):
    def configure(self):
        m = build_SKK()
        x_interface_in = m.fs.unit.feed_side.length_domain.at(2)

        self.unit_solutions[m.fs.unit.deltaP[0]] = -98760.812
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "H2O"]
        ] = 0.01225735
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "NaCl"]
        ] = 7.4932174e-5
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        ] = 0.006690554
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "NaCl"]
        ] = 5.422075e-5
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
        ] = 0.0094392
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]
        ] = 6.5286415e-5
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.1341274
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 0.0009276945
        self.unit_solutions[m.fs.unit.feed_side.N_Re_avg[0]] = 383.229748
        self.unit_solutions[m.fs.unit.feed_side.K_avg[0, "NaCl"]] = 3.0276974e-5
        self.unit_solutions[m.fs.unit.area] = 14.20961

        return m


def build_cp_mod_fixed():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=False,
        concentration_polarization_type=ConcentrationPolarizationType.fixed,
        mass_transfer_coefficient=MassTransferCoefficient.none,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        finite_elements=3,
    )

    # fully specify system
    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.035
    feed_pressure = 50e5
    feed_temperature = 273.15 + 25
    A = 4.2e-12
    B = 3.5e-8
    pressure_atmospheric = 101325
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )

    m.fs.unit.inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
    m.fs.unit.length.fix(8)
    m.fs.unit.recovery_vol_phase[0, "Liq"].fix(0.4)
    m.fs.unit.feed_side.cp_modulus.fix(1.1)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )

    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestReverseOsmosis1D_cp_mod_fixed(UnitTestHarness):
    def configure(self):
        m = build_cp_mod_fixed()
        x_interface_in = m.fs.unit.feed_side.length_domain.at(2)

        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "H2O"]
        ] = 0.0024494142
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "NaCl"]
        ] = 1.863938e-6
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        ] = 0.000359974
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "NaCl"]
        ] = 2.051563e-6
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.3894809
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 0.000652886
        self.unit_solutions[m.fs.unit.area] = 329.708005

        return m


def build_cp_calculated_kf_fixed():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=False,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.fixed,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        finite_elements=3,
    )

    # fully specify system
    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.035
    feed_pressure = 50e5
    feed_temperature = 273.15 + 25
    A = 4.2e-12
    B = 3.5e-8
    pressure_atmospheric = 101325
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )

    m.fs.unit.inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    # m.fs.unit.dP_dx.fix(0)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
    m.fs.unit.length.fix(8)
    m.fs.unit.recovery_vol_phase[0, "Liq"].fix(0.4)
    m.fs.unit.feed_side.K.fix(2e-5)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )

    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestReverseOsmosis1D_cp_calculated_kf_fixed(UnitTestHarness):
    def configure(self):
        m = build_cp_calculated_kf_fixed()
        x_interface_in = m.fs.unit.feed_side.length_domain.at(2)

        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "H2O"]
        ] = 0.002791588
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "NaCl"]
        ] = 1.83061e-6
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        ] = 0.0007081512
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "NaCl"]
        ] = 2.027294e-6
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.3895265
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 0.000464535

        return m


def build_cp_calculated_kf_calculated():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=False,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        finite_elements=3,
    )

    # fully specify system
    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.035
    feed_pressure = 50e5
    feed_temperature = 273.15 + 25
    A = 4.2e-12
    B = 3.5e-8
    pressure_atmospheric = 101325
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )

    m.fs.unit.inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
    m.fs.unit.length.fix(8)
    m.fs.unit.recovery_vol_phase[0, "Liq"].fix(0.4)
    m.fs.unit.feed_side.spacer_porosity.fix(0.75)
    m.fs.unit.feed_side.channel_height.fix(0.002)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )

    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestReverseOsmosis1D_cp_calculated_kf_calculated(UnitTestHarness):
    def configure(self):
        m = build_cp_calculated_kf_calculated()
        x_interface_in = m.fs.unit.feed_side.length_domain.at(2)

        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "H2O"]
        ] = 0.002382908
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "NaCl"]
        ] = 1.8703962e-6
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        ] = 0.000626977
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "NaCl"]
        ] = 2.0339201e-6
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.3895066
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 0.000546681

        return m


def build_friction_factor_spiral_wound():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        pressure_change_type=PressureChangeType.calculated,
        friction_factor=FrictionFactor.spiral_wound,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        finite_elements=3,
    )

    # fully specify system
    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.035
    feed_pressure = 50e5
    feed_temperature = 273.15 + 25
    membrane_area = 19
    A = 4.2e-12
    B = 3.5e-8
    pressure_atmospheric = 101325
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

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
    m.fs.unit.length.fix(8)
    m.fs.unit.feed_side.spacer_porosity.fix(0.75)
    m.fs.unit.feed_side.channel_height.fix(0.002)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )
    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestReverseOsmosis1D_friction_factor_spiral_wound(UnitTestHarness):
    def configure(self):
        m = build_friction_factor_spiral_wound()
        x_interface_in = m.fs.unit.feed_side.length_domain.at(2)

        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "H2O"]
        ] = 0.00590823
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "NaCl"]
        ] = 1.494063e-6
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        ] = 0.00474506
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "NaCl"]
        ] = 1.559319e-6
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.1010762
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 2.9015866e-5

        return m


def build_cp_calculated_kf_calculated_pdrop_fixed_by_dx():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        pressure_change_type=PressureChangeType.fixed_per_unit_length,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        finite_elements=3,
    )

    # fully specify system
    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.035
    feed_pressure = 50e5
    feed_temperature = 273.15 + 25
    A = 4.2e-12
    B = 3.5e-8
    pressure_atmospheric = 101325
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )

    m.fs.unit.inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
    m.fs.unit.length.fix(8)
    m.fs.unit.recovery_vol_phase[0, "Liq"].fix(0.4)
    m.fs.unit.feed_side.spacer_porosity.fix(0.75)
    m.fs.unit.feed_side.channel_height.fix(0.002)
    m.fs.unit.feed_side.dP_dx.fix(-0.1e5)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )

    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestReverseOsmosis1D_cp_calculated_kf_calculated_pdrop_fixed_by_dx(
    UnitTestHarness
):
    def configure(self):
        m = build_cp_calculated_kf_calculated_pdrop_fixed_by_dx()
        x_interface_in = m.fs.unit.feed_side.length_domain.at(2)

        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "H2O"]
        ] = 0.00224943
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "NaCl"]
        ] = 1.872372e-6
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        ] = 0.000551092
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "NaCl"]
        ] = 2.0075516e-6
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.3894964
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 0.000588787

        return m


def build_cp_calculated_kf_calculated_pdrop_fixed_by_stage():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        pressure_change_type=PressureChangeType.fixed_per_stage,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        finite_elements=3,
    )

    # fully specify system
    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.035
    feed_pressure = 50e5
    feed_temperature = 273.15 + 25
    A = 4.2e-12
    B = 3.5e-8
    pressure_atmospheric = 101325
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )

    m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )

    m.fs.unit.inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
    m.fs.unit.length.fix(8)
    m.fs.unit.recovery_vol_phase[0, "Liq"].fix(0.4)
    m.fs.unit.feed_side.spacer_porosity.fix(0.75)
    m.fs.unit.feed_side.channel_height.fix(0.002)
    m.fs.unit.deltaP.fix(-62435.6)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )

    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestReverseOsmosis1D_cp_calculated_kf_calculated_pdrop_fixed_by_stage(
    UnitTestHarness
):
    def configure(self):
        m = build_cp_calculated_kf_calculated_pdrop_fixed_by_stage()

        x_interface_in = m.fs.unit.feed_side.length_domain.at(2)

        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "H2O"]
        ] = 0.002278454
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "NaCl"]
        ] = 1.871967e-6
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        ] = 0.000566998
        self.unit_solutions[
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "NaCl"]
        ] = 2.0134278e-6
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.3894987
        self.unit_solutions[
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 0.0005792222

        return m
