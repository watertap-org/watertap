#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
from pyomo.environ import ConcreteModel

from watertap.core.solvers import get_solver


from idaes.core import (
    FlowsheetBlock,
    FlowDirection,
)
import idaes.core.util.scaling as iscale
from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
import watertap.property_models.seawater_prop_pack as props_sw
import watertap.property_models.water_prop_pack as props_w
from watertap.unit_models.MD.membrane_distillation_1D import MembraneDistillation1D
from watertap.unit_models.MD.MD_channel_base import (
    ConcentrationPolarizationType,
    TemperaturePolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.MD.membrane_distillation_base import MDconfigurationType

from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

solver = get_solver()


def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_hot_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_cold_ch = props_w.WaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation1D(
        hot_ch={
            "property_package": m.fs.properties_hot_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.fixed,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "pressure_change_type": PressureChangeType.fixed_per_stage,
            "flow_direction": FlowDirection.forward,
        },
        cold_ch={
            "property_package": m.fs.properties_cold_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.fixed,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "pressure_change_type": PressureChangeType.fixed_per_stage,
            "flow_direction": FlowDirection.backward,
        },
    )

    # fully specify system
    hot_ch_flow_mass = 1
    hot_ch_mass_frac_TDS = 0.035
    hot_ch_pressure = 7e5
    membrane_pressure_drop = -0.5e5
    membrane_area = 12
    hot_ch_mass_frac_H2O = 1 - hot_ch_mass_frac_TDS
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_TDS
    )
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_H2O
    )

    m.fs.unit.hot_ch_inlet.pressure[0].fix(hot_ch_pressure)
    m.fs.unit.hot_ch_inlet.temperature[0].fix(273.15 + 90)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.length.fix(3)
    m.fs.unit.permeability_coef.fix(1e-10)
    m.fs.unit.membrane_thickness.fix(1e-4)
    m.fs.unit.membrane_thermal_conductivity.fix(0.2)

    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(hot_ch_flow_mass)
    m.fs.unit.cold_ch_inlet.pressure[0].fix(101325)
    m.fs.unit.cold_ch_inlet.temperature[0].fix(273.15 + 25)

    m.fs.unit.hot_ch.deltaP_channel.fix(membrane_pressure_drop)
    m.fs.unit.cold_ch.deltaP_channel.fix(membrane_pressure_drop)

    m.fs.unit.hot_ch.h_conv.fix(2400)
    m.fs.unit.cold_ch.h_conv.fix(2400)

    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestMembraneDisillation1D(UnitTestHarness):
    def configure(self):
        m = build()

        self.unit_solutions[
            m.fs.unit.cold_ch_outlet.flow_mass_phase_comp[0, "Vap", "H2O"]
        ] = 0
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9282725
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 314.5297587
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 343.790759
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.0030606
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.4751264
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.702165528

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"],
                "out": m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"],
            },
        }

        return m


def build_temperature_polarization_none():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_hot_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_cold_ch = props_w.WaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation1D(
        hot_ch={
            "property_package": m.fs.properties_hot_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.none,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "pressure_change_type": PressureChangeType.fixed_per_stage,
            "flow_direction": FlowDirection.forward,
        },
        cold_ch={
            "property_package": m.fs.properties_cold_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.none,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "pressure_change_type": PressureChangeType.fixed_per_stage,
            "flow_direction": FlowDirection.backward,
        },
    )

    # fully specify system
    hot_ch_flow_mass = 1
    hot_ch_mass_frac_TDS = 0.035
    hot_ch_pressure = 7e5
    membrane_pressure_drop = -0.5e5
    membrane_area = 5
    hot_ch_mass_frac_H2O = 1 - hot_ch_mass_frac_TDS
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_TDS
    )
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_H2O
    )

    m.fs.unit.hot_ch_inlet.pressure[0].fix(hot_ch_pressure)
    m.fs.unit.hot_ch_inlet.temperature[0].fix(273.15 + 90)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.length.fix(3)
    m.fs.unit.permeability_coef.fix(1e-10)
    m.fs.unit.membrane_thickness.fix(1e-4)
    m.fs.unit.membrane_thermal_conductivity.fix(0.2)

    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(hot_ch_flow_mass)
    m.fs.unit.cold_ch_inlet.pressure[0].fix(7e5)
    m.fs.unit.cold_ch_inlet.temperature[0].fix(273.15 + 25)

    m.fs.unit.hot_ch.deltaP_channel.fix(0)
    m.fs.unit.cold_ch.deltaP_channel.fix(0)

    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestMembraneDisillation1D_temperature_polarization_none(UnitTestHarness):
    def configure(self):
        m = build_temperature_polarization_none()

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.922776
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 308.1848897
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 349.046787
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.008444798
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.4870051
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.7830275

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"],
                "out": m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"],
            },
        }

        return m


def build_temperature_polarization_fixed():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_hot_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_cold_ch = props_w.WaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation1D(
        hot_ch={
            "property_package": m.fs.properties_hot_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.fixed,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "pressure_change_type": PressureChangeType.fixed_per_stage,
            "flow_direction": FlowDirection.forward,
        },
        cold_ch={
            "property_package": m.fs.properties_cold_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.fixed,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "pressure_change_type": PressureChangeType.fixed_per_stage,
            "flow_direction": FlowDirection.backward,
        },
    )

    # fully specify system
    hot_ch_flow_mass = 1
    hot_ch_mass_frac_TDS = 0.035
    hot_ch_pressure = 7e5
    membrane_pressure_drop = -0.5e5
    membrane_area = 12
    hot_ch_mass_frac_H2O = 1 - hot_ch_mass_frac_TDS
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_TDS
    )
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_H2O
    )

    m.fs.unit.hot_ch_inlet.pressure[0].fix(hot_ch_pressure)
    m.fs.unit.hot_ch_inlet.temperature[0].fix(273.15 + 90)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.length.fix(3)
    m.fs.unit.permeability_coef.fix(1e-10)
    m.fs.unit.membrane_thickness.fix(1e-4)
    m.fs.unit.membrane_thermal_conductivity.fix(0.2)

    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(hot_ch_flow_mass)
    m.fs.unit.cold_ch_inlet.pressure[0].fix(7e5)
    m.fs.unit.cold_ch_inlet.temperature[0].fix(273.15 + 25)

    m.fs.unit.hot_ch.deltaP_channel.fix(0)
    m.fs.unit.cold_ch.deltaP_channel.fix(0)

    m.fs.unit.hot_ch.h_conv.fix(2400)
    m.fs.unit.cold_ch.h_conv.fix(2400)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation1D_temperature_polarization_fixed(UnitTestHarness):
    def configure(self):
        m = build_temperature_polarization_fixed()

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.928775
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 314.523531
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 343.78807
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.0030602
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.4751043
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.7021242

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"],
                "out": m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"],
            },
        }

        return m


def build_temperature_polarization_calculated():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_hot_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_cold_ch = props_w.WaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation1D(
        hot_ch={
            "property_package": m.fs.properties_hot_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.calculated,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "pressure_change_type": PressureChangeType.fixed_per_stage,
            "flow_direction": FlowDirection.forward,
        },
        cold_ch={
            "property_package": m.fs.properties_cold_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "pressure_change_type": PressureChangeType.fixed_per_stage,
            "flow_direction": FlowDirection.backward,
        },
    )

    # fully specify system
    hot_ch_flow_mass = 1
    hot_ch_mass_frac_TDS = 0.035
    hot_ch_pressure = 7e5
    membrane_pressure_drop = -0.5e5
    membrane_area = 12
    hot_ch_mass_frac_H2O = 1 - hot_ch_mass_frac_TDS
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_TDS
    )
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_H2O
    )

    m.fs.unit.hot_ch_inlet.pressure[0].fix(hot_ch_pressure)
    m.fs.unit.hot_ch_inlet.temperature[0].fix(273.15 + 90)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.permeability_coef.fix(1e-10)
    m.fs.unit.membrane_thickness.fix(1e-4)
    m.fs.unit.membrane_thermal_conductivity.fix(0.2)

    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(hot_ch_flow_mass)
    m.fs.unit.cold_ch_inlet.pressure[0].fix(7e5)
    m.fs.unit.cold_ch_inlet.temperature[0].fix(273.15 + 25)

    m.fs.unit.hot_ch.deltaP_channel.fix(0)
    m.fs.unit.cold_ch.deltaP_channel.fix(0)

    m.fs.unit.length.fix(8)
    m.fs.unit.hot_ch.channel_height.fix(0.0019)
    m.fs.unit.hot_ch.spacer_porosity.fix(0.77)
    m.fs.unit.cold_ch.channel_height.fix(0.0019)
    m.fs.unit.cold_ch.spacer_porosity.fix(0.77)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation1D_temperature_polarization_calculated(UnitTestHarness):
    def configure(self):
        m = build_temperature_polarization_calculated()

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.92163799
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 305.028018
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 351.87209
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.0036135
        self.unit_solutions[m.fs.unit.recovery_mass[0]] = 0.0449347
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.47367225
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.8264937

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"],
                "out": m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"],
            },
        }

        return m


def build_temperature_polarization_calculated_concentration_polarization_fixed():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_hot_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_cold_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation1D(
        hot_ch={
            "property_package": m.fs.properties_hot_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.calculated,
            "concentration_polarization_type": ConcentrationPolarizationType.fixed,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "pressure_change_type": PressureChangeType.fixed_per_stage,
            "flow_direction": FlowDirection.forward,
        },
        cold_ch={
            "property_package": m.fs.properties_cold_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "pressure_change_type": PressureChangeType.fixed_per_stage,
            "flow_direction": FlowDirection.backward,
        },
    )

    # fully specify system
    hot_ch_flow_mass = 1
    hot_ch_mass_frac_TDS = 0.035
    hot_ch_pressure = 7e5
    membrane_pressure_drop = -0.5e5
    membrane_area = 12
    hot_ch_mass_frac_H2O = 1 - hot_ch_mass_frac_TDS
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_TDS
    )
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_H2O
    )

    m.fs.unit.hot_ch_inlet.pressure[0].fix(hot_ch_pressure)
    m.fs.unit.hot_ch_inlet.temperature[0].fix(273.15 + 90)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.permeability_coef.fix(1e-10)
    m.fs.unit.membrane_thickness.fix(1e-4)
    m.fs.unit.membrane_thermal_conductivity.fix(0.2)
    m.fs.unit.hot_ch.cp_modulus.fix(1.1)

    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(hot_ch_flow_mass)
    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0)
    m.fs.unit.cold_ch_inlet.pressure[0].fix(7e5)
    m.fs.unit.cold_ch_inlet.temperature[0].fix(273.15 + 25)

    m.fs.unit.hot_ch.deltaP_channel.fix(0)
    m.fs.unit.cold_ch.deltaP_channel.fix(0)

    m.fs.unit.length.fix(8)
    m.fs.unit.hot_ch.channel_height.fix(0.0019)
    m.fs.unit.hot_ch.spacer_porosity.fix(0.77)
    m.fs.unit.cold_ch.channel_height.fix(0.0019)
    m.fs.unit.cold_ch.spacer_porosity.fix(0.77)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation1D_temperature_polarization_calculated_concentration_polarization_fixed(
    UnitTestHarness
):
    def configure(self):
        m = build_temperature_polarization_calculated_concentration_polarization_fixed()

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.92185905
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 305.049789
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 351.863122
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.00359508

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"],
                "out": m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"],
            },
        }

        return m


def build_temperature_polarization_calculated_concentration_polarization_calculated_K_fixed():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_hot_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_cold_ch = props_w.WaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation1D(
        hot_ch={
            "property_package": m.fs.properties_hot_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.calculated,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.fixed,
            "pressure_change_type": PressureChangeType.fixed_per_stage,
            "flow_direction": FlowDirection.forward,
        },
        cold_ch={
            "property_package": m.fs.properties_cold_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "pressure_change_type": PressureChangeType.fixed_per_stage,
            "flow_direction": FlowDirection.backward,
        },
    )

    # fully specify system
    hot_ch_flow_mass = 1
    hot_ch_mass_frac_TDS = 0.035
    hot_ch_pressure = 7e5
    membrane_pressure_drop = -0.5e5
    membrane_area = 12
    hot_ch_mass_frac_H2O = 1 - hot_ch_mass_frac_TDS
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_TDS
    )
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_H2O
    )

    m.fs.unit.hot_ch_inlet.pressure[0].fix(hot_ch_pressure)
    m.fs.unit.hot_ch_inlet.temperature[0].fix(273.15 + 90)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.permeability_coef.fix(1e-10)
    m.fs.unit.membrane_thickness.fix(1e-4)
    m.fs.unit.membrane_thermal_conductivity.fix(0.2)
    m.fs.unit.hot_ch.K.fix(3.15e-5)

    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(hot_ch_flow_mass)
    m.fs.unit.cold_ch_inlet.pressure[0].fix(7e5)
    m.fs.unit.cold_ch_inlet.temperature[0].fix(273.15 + 25)

    m.fs.unit.hot_ch.deltaP_channel.fix(0)
    m.fs.unit.cold_ch.deltaP_channel.fix(0)

    m.fs.unit.length.fix(8)
    m.fs.unit.hot_ch.channel_height.fix(0.0019)
    m.fs.unit.hot_ch.spacer_porosity.fix(0.77)
    m.fs.unit.cold_ch.channel_height.fix(0.0019)
    m.fs.unit.cold_ch.spacer_porosity.fix(0.77)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation1D_temperature_polarization_calculated_concentration_polarization_calculated_K_fixed(
    UnitTestHarness
):
    def configure(self):
        m = (
            build_temperature_polarization_calculated_concentration_polarization_calculated_K_fixed()
        )

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.922
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 305.062335
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 351.85863
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.00358303
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.47002591
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.826866

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"],
                "out": m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"],
            },
        }

        return m


def build_temperature_polarization_calculated_concentration_polarization_calculated_K_calculated():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_hot_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_cold_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation1D(
        hot_ch={
            "property_package": m.fs.properties_hot_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.calculated,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.calculated,
            "pressure_change_type": PressureChangeType.fixed_per_stage,
            "flow_direction": FlowDirection.forward,
        },
        cold_ch={
            "property_package": m.fs.properties_cold_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "pressure_change_type": PressureChangeType.fixed_per_stage,
            "flow_direction": FlowDirection.backward,
        },
    )

    # fully specify system
    hot_ch_flow_mass = 1
    hot_ch_mass_frac_TDS = 0.035
    hot_ch_pressure = 7e5
    membrane_pressure_drop = -5e5
    membrane_area = 12
    hot_ch_mass_frac_H2O = 1 - hot_ch_mass_frac_TDS
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_TDS
    )
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_H2O
    )

    m.fs.unit.hot_ch_inlet.pressure[0].fix(hot_ch_pressure)
    m.fs.unit.hot_ch_inlet.temperature[0].fix(273.15 + 90)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.permeability_coef.fix(1e-10)
    m.fs.unit.membrane_thickness.fix(1e-4)
    m.fs.unit.membrane_thermal_conductivity.fix(0.2)

    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(hot_ch_flow_mass)
    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0)
    m.fs.unit.cold_ch_inlet.pressure[0].fix(7e5)
    m.fs.unit.cold_ch_inlet.temperature[0].fix(273.15 + 25)

    m.fs.unit.hot_ch.deltaP_channel.fix(membrane_pressure_drop)
    m.fs.unit.cold_ch.deltaP_channel.fix(membrane_pressure_drop)

    m.fs.unit.length.fix(8)
    m.fs.unit.hot_ch.channel_height.fix(0.0019)
    m.fs.unit.hot_ch.spacer_porosity.fix(0.77)
    m.fs.unit.cold_ch.channel_height.fix(0.0019)
    m.fs.unit.cold_ch.spacer_porosity.fix(0.77)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation1D_temperature_polarization_calculated_concentration_polarization_calculated_K_calculated(
    UnitTestHarness
):
    def configure(self):
        m = (
            build_temperature_polarization_calculated_concentration_polarization_calculated_K_calculated()
        )

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.921668829
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 305.125855
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 351.9792
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.00361093
        self.unit_solutions[m.fs.unit.recovery_mass[0]] = 0.0449028
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.47326646
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.82814152

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"],
                "out": m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"],
            },
        }

        return m


def build_temperature_polarization_calculated_concentration_polarization_calculated_K_calculated_pressure_fixed_per_unit_length():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_hot_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_cold_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation1D(
        hot_ch={
            "property_package": m.fs.properties_hot_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.calculated,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.calculated,
            "pressure_change_type": PressureChangeType.fixed_per_unit_length,
            "flow_direction": FlowDirection.forward,
        },
        cold_ch={
            "property_package": m.fs.properties_cold_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "pressure_change_type": PressureChangeType.fixed_per_unit_length,
            "flow_direction": FlowDirection.backward,
        },
    )

    # fully specify system
    hot_ch_flow_mass = 1
    hot_ch_mass_frac_TDS = 0.035
    hot_ch_pressure = 7e5
    membrane_pressure_drop = -5e5
    membrane_area = 12
    hot_ch_mass_frac_H2O = 1 - hot_ch_mass_frac_TDS
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_TDS
    )
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_H2O
    )

    m.fs.unit.hot_ch_inlet.pressure[0].fix(hot_ch_pressure)
    m.fs.unit.hot_ch_inlet.temperature[0].fix(273.15 + 90)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.permeability_coef.fix(1e-10)
    m.fs.unit.membrane_thickness.fix(1e-4)
    m.fs.unit.membrane_thermal_conductivity.fix(0.2)

    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(hot_ch_flow_mass)
    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0)
    m.fs.unit.cold_ch_inlet.pressure[0].fix(7e5)
    m.fs.unit.cold_ch_inlet.temperature[0].fix(273.15 + 25)

    length = 8
    m.fs.unit.length.fix(length)

    m.fs.unit.hot_ch.dP_dx.fix(membrane_pressure_drop / length)
    m.fs.unit.cold_ch.dP_dx.fix(membrane_pressure_drop / length)

    m.fs.unit.hot_ch.channel_height.fix(0.0019)
    m.fs.unit.hot_ch.spacer_porosity.fix(0.77)
    m.fs.unit.cold_ch.channel_height.fix(0.0019)
    m.fs.unit.cold_ch.spacer_porosity.fix(0.77)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation1D_temperature_polarization_calculated_concentration_polarization_calculated_K_calculated_pressure_fixed_per_unit_length(
    UnitTestHarness
):
    def configure(self):
        m = (
            build_temperature_polarization_calculated_concentration_polarization_calculated_K_calculated_pressure_fixed_per_unit_length()
        )

        self.unit_solutions[m.fs.unit.hot_ch.length] = 8
        self.unit_solutions[m.fs.unit.cold_ch.length] = 8
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9216688
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 305.125855
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 351.9792
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.00361093
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.4732665
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.8281415
        self.unit_solutions[m.fs.unit.hot_ch.deltaP_channel[0]] = -500000.0
        self.unit_solutions[m.fs.unit.cold_ch.deltaP_channel[0]] = -500000.0

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"],
                "out": m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"],
            },
        }

        return m


def build_temperature_polarization_calculated_concentration_polarization_calculated_K_calculated_pressure_calculated():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_hot_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_cold_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation1D(
        hot_ch={
            "property_package": m.fs.properties_hot_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.calculated,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.calculated,
            "pressure_change_type": PressureChangeType.calculated,
            "flow_direction": FlowDirection.forward,
        },
        cold_ch={
            "property_package": m.fs.properties_cold_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "temperature_polarization_type": TemperaturePolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "pressure_change_type": PressureChangeType.calculated,
            "flow_direction": FlowDirection.backward,
        },
    )

    # fully specify system
    hot_ch_flow_mass = 1
    hot_ch_mass_frac_TDS = 0.035
    hot_ch_pressure = 7e5
    membrane_pressure_drop = -5e5
    membrane_area = 12
    hot_ch_mass_frac_H2O = 1 - hot_ch_mass_frac_TDS
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_TDS
    )
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_H2O
    )

    m.fs.unit.hot_ch_inlet.pressure[0].fix(hot_ch_pressure)
    m.fs.unit.hot_ch_inlet.temperature[0].fix(273.15 + 90)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.permeability_coef.fix(1e-10)
    m.fs.unit.membrane_thickness.fix(1e-4)
    m.fs.unit.membrane_thermal_conductivity.fix(0.2)

    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(hot_ch_flow_mass)
    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0)
    m.fs.unit.cold_ch_inlet.pressure[0].fix(7e5)
    m.fs.unit.cold_ch_inlet.temperature[0].fix(273.15 + 25)

    length = 8
    m.fs.unit.length.fix(length)

    m.fs.unit.hot_ch.channel_height.fix(0.0019)
    m.fs.unit.hot_ch.spacer_porosity.fix(0.77)
    m.fs.unit.cold_ch.channel_height.fix(0.0019)
    m.fs.unit.cold_ch.spacer_porosity.fix(0.77)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation1D_temperature_polarization_calculated_concentration_polarization_calculated_K_calculated_pressure_calculated(
    UnitTestHarness
):
    def configure(self):
        m = (
            build_temperature_polarization_calculated_concentration_polarization_calculated_K_calculated_pressure_calculated()
        )

        self.unit_solutions[m.fs.unit.hot_ch.length] = 8
        self.unit_solutions[m.fs.unit.cold_ch.length] = 8
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9216924
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 305.093352
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 351.93729
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.00360897
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.4731037
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.82749678
        self.unit_solutions[m.fs.unit.hot_ch.deltaP_channel[0]] = -295331.803525
        self.unit_solutions[m.fs.unit.cold_ch.deltaP_channel[0]] = -331332.910349

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"],
                "out": m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"],
            },
        }

        return m


def build_temperature_polarization_none_vmd():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_hot_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_cold_ch = props_w.WaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation1D(
        MD_configuration_Type=MDconfigurationType.VMD,
        hot_ch={
            "property_package": m.fs.properties_hot_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": False,
            "temperature_polarization_type": TemperaturePolarizationType.none,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "flow_direction": FlowDirection.forward,
        },
        cold_ch={
            "property_package": m.fs.properties_cold_ch,
            "temperature_polarization_type": TemperaturePolarizationType.none,
            "has_pressure_change": False,
            "flow_direction": FlowDirection.forward,
        },
    )

    m.fs.unit.length.fix(0.1)

    # Fully specify the system
    hot_ch_flow_mass = 1.0
    hot_ch_mass_frac_TDS = 0.035
    hot_ch_pressure = 7e5
    membrane_area = 1.0
    hot_ch_mass_frac_H2O = 1 - hot_ch_mass_frac_TDS

    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_TDS
    )
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_H2O
    )
    m.fs.unit.hot_ch_inlet.pressure[0].fix(hot_ch_pressure)
    m.fs.unit.hot_ch_inlet.temperature[0].fix(273.15 + 90)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.permeability_coef.fix(1e-10)
    m.fs.unit.membrane_thickness.fix(1e-4)

    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    m.fs.unit.cold_ch_inlet.pressure[0].fix(10000)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation1D_temperature_polarization_none_vmd(UnitTestHarness):
    def configure(self):
        m = build_temperature_polarization_none_vmd()

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9333282339579946
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 342.25342
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 342.253422
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.03241786

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Vap", "H2O"],
                "out": m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_outlet.flow_mass_phase_comp[0, "Vap", "H2O"],
            },
        }

        return m


def build_temperature_polarization_fixed_hot_vmd():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_hot_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_cold_ch = props_w.WaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation1D(
        MD_configuration_Type=MDconfigurationType.VMD,
        hot_ch={
            "property_package": m.fs.properties_hot_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": False,
            "temperature_polarization_type": TemperaturePolarizationType.fixed,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "flow_direction": FlowDirection.forward,
        },
        cold_ch={
            "property_package": m.fs.properties_cold_ch,
            "temperature_polarization_type": TemperaturePolarizationType.none,
            "has_pressure_change": False,
            "flow_direction": FlowDirection.forward,
        },
    )

    m.fs.unit.length.fix(0.1)

    # Fully specify the system
    hot_ch_flow_mass = 1.0
    hot_ch_mass_frac_TDS = 0.035
    hot_ch_pressure = 7e5
    membrane_area = 1.0
    hot_ch_mass_frac_H2O = 1 - hot_ch_mass_frac_TDS

    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_TDS
    )
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_H2O
    )
    m.fs.unit.hot_ch_inlet.pressure[0].fix(hot_ch_pressure)
    m.fs.unit.hot_ch_inlet.temperature[0].fix(273.15 + 90)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.permeability_coef.fix(1e-10)
    m.fs.unit.membrane_thickness.fix(1e-4)

    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    m.fs.unit.cold_ch_inlet.pressure[0].fix(10000)
    m.fs.unit.hot_ch.h_conv.fix(2400)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation1D_temperature_polarization_fixed_hot_vmd(UnitTestHarness):
    def configure(self):
        m = build_temperature_polarization_fixed_hot_vmd()

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9482002512081292
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 352.839383
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = (
            337.84769812257065
        )
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.016983716

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Vap", "H2O"],
                "out": m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_outlet.flow_mass_phase_comp[0, "Vap", "H2O"],
            },
        }

        return m


def build_temperature_polarization_calculated_hot_vmd():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_hot_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_cold_ch = props_w.WaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation1D(
        MD_configuration_Type=MDconfigurationType.VMD,
        hot_ch={
            "property_package": m.fs.properties_hot_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": False,
            "temperature_polarization_type": TemperaturePolarizationType.calculated,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "flow_direction": FlowDirection.forward,
        },
        cold_ch={
            "property_package": m.fs.properties_cold_ch,
            "temperature_polarization_type": TemperaturePolarizationType.none,
            "has_pressure_change": False,
            "flow_direction": FlowDirection.forward,
        },
    )

    m.fs.unit.length.fix(0.1)

    # Fully specify the system
    hot_ch_flow_mass = 1.0
    hot_ch_mass_frac_TDS = 0.035
    hot_ch_pressure = 7e5
    membrane_area = 1.0
    hot_ch_mass_frac_H2O = 1 - hot_ch_mass_frac_TDS

    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_TDS
    )
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_H2O
    )
    m.fs.unit.hot_ch_inlet.pressure[0].fix(hot_ch_pressure)
    m.fs.unit.hot_ch_inlet.temperature[0].fix(273.15 + 90)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.permeability_coef.fix(1e-10)
    m.fs.unit.membrane_thickness.fix(1e-4)

    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    m.fs.unit.cold_ch_inlet.pressure[0].fix(10000)

    m.fs.unit.hot_ch.channel_height.fix(0.0019)
    m.fs.unit.hot_ch.spacer_porosity.fix(0.77)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation1D_temperature_polarization_calculated_hot_vmd(
    UnitTestHarness
):
    def configure(self):
        m = build_temperature_polarization_calculated_hot_vmd()

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9450790917786648
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 350.72706
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 339.2750084657155
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.02018963

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Vap", "H2O"],
                "out": m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_outlet.flow_mass_phase_comp[0, "Vap", "H2O"],
            },
        }

        return m


def build_temperature_polarization_concentration_polarization_calculated_hot_vmd():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_hot_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_cold_ch = props_w.WaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation1D(
        MD_configuration_Type=MDconfigurationType.VMD,
        hot_ch={
            "property_package": m.fs.properties_hot_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": False,
            "temperature_polarization_type": TemperaturePolarizationType.calculated,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.calculated,
            "flow_direction": FlowDirection.forward,
        },
        cold_ch={
            "property_package": m.fs.properties_cold_ch,
            "temperature_polarization_type": TemperaturePolarizationType.none,
            "has_pressure_change": False,
            "flow_direction": FlowDirection.forward,
        },
    )

    m.fs.unit.length.fix(0.1)

    # Fully specify the system
    hot_ch_flow_mass = 1.0
    hot_ch_mass_frac_TDS = 0.035
    hot_ch_pressure = 7e5
    membrane_area = 1.0
    hot_ch_mass_frac_H2O = 1 - hot_ch_mass_frac_TDS

    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_TDS
    )
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_H2O
    )
    m.fs.unit.hot_ch_inlet.pressure[0].fix(hot_ch_pressure)
    m.fs.unit.hot_ch_inlet.temperature[0].fix(273.15 + 90)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.permeability_coef.fix(1e-10)
    m.fs.unit.membrane_thickness.fix(1e-4)

    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    m.fs.unit.cold_ch_inlet.pressure[0].fix(10000)

    m.fs.unit.hot_ch.channel_height.fix(0.0019)
    m.fs.unit.hot_ch.spacer_porosity.fix(0.77)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation1D_temperature_polarization_concentration_polarization_calculated_hot_vmd(
    UnitTestHarness
):
    def configure(self):
        m = (
            build_temperature_polarization_concentration_polarization_calculated_hot_vmd()
        )

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9450790917786648
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 350.82795
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 339.763283
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.020033

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Vap", "H2O"],
                "out": m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_outlet.flow_mass_phase_comp[0, "Vap", "H2O"],
            },
        }

        return m


def build_temperature_polarization_concentration_polarization_pressure_calculated_hot_vmd():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_hot_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_cold_ch = props_w.WaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation1D(
        MD_configuration_Type=MDconfigurationType.VMD,
        hot_ch={
            "property_package": m.fs.properties_hot_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "pressure_change_type": PressureChangeType.calculated,
            "temperature_polarization_type": TemperaturePolarizationType.calculated,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.calculated,
            "flow_direction": FlowDirection.forward,
        },
        cold_ch={
            "property_package": m.fs.properties_cold_ch,
            "temperature_polarization_type": TemperaturePolarizationType.none,
            "has_pressure_change": False,
            "flow_direction": FlowDirection.forward,
        },
    )

    m.fs.unit.length.fix(0.1)

    # Fully specify the system
    hot_ch_flow_mass = 1.0
    hot_ch_mass_frac_TDS = 0.035
    hot_ch_pressure = 7e5
    membrane_area = 1.0
    hot_ch_mass_frac_H2O = 1 - hot_ch_mass_frac_TDS

    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_TDS
    )
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_H2O
    )
    m.fs.unit.hot_ch_inlet.pressure[0].fix(hot_ch_pressure)
    m.fs.unit.hot_ch_inlet.temperature[0].fix(273.15 + 90)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.permeability_coef.fix(1e-10)
    m.fs.unit.membrane_thickness.fix(1e-4)

    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    m.fs.unit.cold_ch_inlet.pressure[0].fix(10000)

    m.fs.unit.hot_ch.channel_height.fix(0.0019)
    m.fs.unit.hot_ch.spacer_porosity.fix(0.77)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation1D_temperature_polarization_concentration_polarization_pressure_calculated_hot_vmd(
    UnitTestHarness
):
    def configure(self):
        m = (
            build_temperature_polarization_concentration_polarization_pressure_calculated_hot_vmd()
        )

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9450790917786648
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 350.827972
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 339.763298
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.020033

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Vap", "H2O"],
                "out": m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_outlet.flow_mass_phase_comp[0, "Vap", "H2O"],
            },
        }

        return m


def build_temperature_polarization_none_pgmd():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_hot_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_cold_ch = props_w.WaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation1D(
        MD_configuration_Type=MDconfigurationType.GMD,
        hot_ch={
            "property_package": m.fs.properties_hot_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": False,
            "temperature_polarization_type": TemperaturePolarizationType.none,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "flow_direction": FlowDirection.forward,
        },
        cold_ch={
            "property_package": m.fs.properties_hot_ch,
            "temperature_polarization_type": TemperaturePolarizationType.none,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "has_pressure_change": False,
            "flow_direction": FlowDirection.backward,
        },
        gap_ch={
            "property_package": m.fs.properties_cold_ch,
            "temperature_polarization_type": TemperaturePolarizationType.fixed,
            "has_pressure_change": False,
            "flow_direction": FlowDirection.forward,
        },
    )

    m.fs.unit.length.fix(2)

    # fully specify system
    hot_ch_flow_mass = 1
    hot_ch_mass_frac_TDS = 0.035
    hot_ch_pressure = 7e5
    membrane_area = 3
    hot_ch_mass_frac_H2O = 1 - hot_ch_mass_frac_TDS
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_TDS
    )
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_H2O
    )

    m.fs.unit.hot_ch_inlet.pressure[0].fix(hot_ch_pressure)
    m.fs.unit.hot_ch_inlet.temperature[0].fix(273.15 + 90)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.length.fix(2)

    m.fs.unit.permeability_coef.fix(1e-10)
    m.fs.unit.membrane_thickness.fix(1e-4)
    m.fs.unit.gap_thickness.fix(0.0004)
    m.fs.unit.membrane_thermal_conductivity.fix(0.2)
    m.fs.unit.gap_thermal_conductivity.fix(0.06)
    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_H2O
    )
    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_TDS
    )

    m.fs.unit.gap_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(0)
    m.fs.unit.gap_ch_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)

    m.fs.unit.cold_ch_inlet.pressure[0].fix(100000)

    m.fs.unit.gap_ch_inlet.pressure[0].fix(100000)

    m.fs.unit.cold_ch_inlet.temperature[0].fix(273.15 + 25)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation1D_temperature_polarization_none_pgmd(UnitTestHarness):
    def configure(self):
        m = build_temperature_polarization_none_pgmd()

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9580840582054698
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 356.6334998833635
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = (
            304.69178989433595
        )
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.002317366
        self.unit_solutions[m.fs.unit.gap_ch_outlet.temperature[0]] = 326.62398751842414

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.gap_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.gap_ch_inlet.flow_mass_phase_comp[0, "Vap", "H2O"],
                "out": m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.cold_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.gap_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.gap_ch_outlet.flow_mass_phase_comp[0, "Vap", "H2O"],
            },
        }

        return m


########


def build_temperature_polarization_concentration_polarization_pressure_calculated_pgmd():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_hot_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_cold_ch = props_w.WaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation1D(
        MD_configuration_Type=MDconfigurationType.GMD,
        hot_ch={
            "property_package": m.fs.properties_hot_ch,
            "property_package_vapor": m.fs.properties_vapor,
            "has_pressure_change": True,
            "pressure_change_type": PressureChangeType.calculated,
            "temperature_polarization_type": TemperaturePolarizationType.calculated,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.calculated,
            "flow_direction": FlowDirection.forward,
        },
        cold_ch={
            "property_package": m.fs.properties_hot_ch,
            "temperature_polarization_type": TemperaturePolarizationType.calculated,
            "concentration_polarization_type": ConcentrationPolarizationType.none,
            "mass_transfer_coefficient": MassTransferCoefficient.none,
            "has_pressure_change": True,
            "pressure_change_type": PressureChangeType.calculated,
            "flow_direction": FlowDirection.backward,
        },
        gap_ch={
            "property_package": m.fs.properties_cold_ch,
            "temperature_polarization_type": TemperaturePolarizationType.fixed,
            "has_pressure_change": False,
            "flow_direction": FlowDirection.forward,
        },
    )

    m.fs.unit.length.fix(2)

    # fully specify system
    hot_ch_flow_mass = 1
    hot_ch_mass_frac_TDS = 0.035
    hot_ch_pressure = 7e5
    membrane_area = 3
    hot_ch_mass_frac_H2O = 1 - hot_ch_mass_frac_TDS
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_TDS
    )
    m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_H2O
    )

    m.fs.unit.hot_ch_inlet.pressure[0].fix(hot_ch_pressure)
    m.fs.unit.hot_ch_inlet.temperature[0].fix(273.15 + 90)
    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.length.fix(2)

    m.fs.unit.permeability_coef.fix(1e-10)
    m.fs.unit.membrane_thickness.fix(1e-4)
    m.fs.unit.gap_thickness.fix(0.0004)
    m.fs.unit.membrane_thermal_conductivity.fix(0.2)
    m.fs.unit.gap_thermal_conductivity.fix(0.06)
    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_H2O
    )
    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_TDS
    )

    m.fs.unit.gap_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(0)
    m.fs.unit.gap_ch_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)

    m.fs.unit.cold_ch_inlet.pressure[0].fix(100000)

    m.fs.unit.gap_ch_inlet.pressure[0].fix(100000)

    m.fs.unit.cold_ch_inlet.temperature[0].fix(273.15 + 25)

    m.fs.unit.hot_ch.channel_height.fix(0.0019)
    m.fs.unit.hot_ch.spacer_porosity.fix(0.77)

    m.fs.unit.cold_ch.channel_height.fix(0.0019)
    m.fs.unit.cold_ch.spacer_porosity.fix(0.77)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation1D_temperature_polarization_concentration_polarization_pressure_calculated_pgmd(
    UnitTestHarness
):
    def configure(self):
        m = (
            build_temperature_polarization_concentration_polarization_pressure_calculated_pgmd()
        )

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9580840582054698
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 357.132013
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = (
            304.69178989433595
        )
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.002241485
        self.unit_solutions[m.fs.unit.gap_ch_outlet.temperature[0]] = 327.03213

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.gap_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.gap_ch_inlet.flow_mass_phase_comp[0, "Vap", "H2O"],
                "out": m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.cold_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
                + m.fs.unit.cold_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.gap_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
                + m.fs.unit.gap_ch_outlet.flow_mass_phase_comp[0, "Vap", "H2O"],
            },
        }

        return m
