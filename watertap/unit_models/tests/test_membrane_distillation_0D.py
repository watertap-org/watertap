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
from watertap.unit_models.MD.membrane_distillation_0D import MembraneDistillation0D
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
    m.fs.unit = MembraneDistillation0D(
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
        MD_configuration_Type=MDconfigurationType.DCMD,
    )

    # fully specify system
    hot_ch_flow_mass = 1
    hot_ch_mass_frac_TDS = 0.035
    hot_ch_pressure = 7e5
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
    m.fs.unit.cold_ch_inlet.pressure[0].fix(101325)
    m.fs.unit.cold_ch_inlet.temperature[0].fix(273.15 + 25)

    m.fs.unit.hot_ch.deltaP.fix(0)
    m.fs.unit.cold_ch.deltaP.fix(0)

    m.fs.unit.hot_ch.h_conv.fix(2400)
    m.fs.unit.cold_ch.h_conv.fix(2400)

    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestMembraneDisillation0D(UnitTestHarness):
    def configure(self):
        m = build()

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.92471851
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 312.652577449503
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 343.401843712127
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.003356790
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.5268
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.69618

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
    m.fs.unit = MembraneDistillation0D(
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
    m.fs.unit.permeability_coef.fix(1e-10)
    m.fs.unit.membrane_thickness.fix(1e-4)
    m.fs.unit.membrane_thermal_conductivity.fix(0.2)
    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(hot_ch_flow_mass)
    m.fs.unit.cold_ch_inlet.pressure[0].fix(7e5)
    m.fs.unit.cold_ch_inlet.temperature[0].fix(273.15 + 25)
    m.fs.unit.hot_ch.deltaP.fix(0)
    m.fs.unit.cold_ch.deltaP.fix(0)

    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestMembraneDisillation0D_temperature_polarization_none(UnitTestHarness):
    def configure(self):
        m = build_temperature_polarization_none()

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9016872
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 301.76628
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 351.44918
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.01266
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.687139
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.819987

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
    m.fs.unit = MembraneDistillation0D(
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

    m.fs.unit.hot_ch.deltaP.fix(0)
    m.fs.unit.cold_ch.deltaP.fix(0)

    m.fs.unit.hot_ch.h_conv.fix(2400)
    m.fs.unit.cold_ch.h_conv.fix(2400)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation0D_temperature_polarization_fixed(UnitTestHarness):
    def configure(self):
        m = build_temperature_polarization_fixed()

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.92471851
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 312.652577449503
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 343.401843712127
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.003356790
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.5268
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.69618

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
    m.fs.unit = MembraneDistillation0D(
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

    m.fs.unit.hot_ch.deltaP.fix(0)
    m.fs.unit.cold_ch.deltaP.fix(0)

    m.fs.unit.length.fix(8)
    m.fs.unit.hot_ch.channel_height.fix(0.0019)
    m.fs.unit.hot_ch.spacer_porosity.fix(0.77)
    m.fs.unit.cold_ch.channel_height.fix(0.0019)
    m.fs.unit.cold_ch.spacer_porosity.fix(0.77)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation0D_temperature_polarization_calculated(UnitTestHarness):
    def configure(self):
        m = build_temperature_polarization_calculated()

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9002124090930
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 300.20215
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 352.806
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.0053989
        self.unit_solutions[m.fs.unit.recovery_mass[0]] = 0.06713739
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.68575
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.8408

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
    m.fs.properties_cold_ch = props_w.WaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation0D(
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
    m.fs.unit.cold_ch_inlet.pressure[0].fix(7e5)
    m.fs.unit.cold_ch_inlet.temperature[0].fix(273.15 + 25)

    m.fs.unit.hot_ch.deltaP.fix(0)
    m.fs.unit.cold_ch.deltaP.fix(0)

    m.fs.unit.length.fix(8)
    m.fs.unit.hot_ch.channel_height.fix(0.0019)
    m.fs.unit.hot_ch.spacer_porosity.fix(0.77)
    m.fs.unit.cold_ch.channel_height.fix(0.0019)
    m.fs.unit.cold_ch.spacer_porosity.fix(0.77)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation0D_temperature_polarization_calculated_concentration_polarization_fixed(
    UnitTestHarness
):
    def configure(self):
        m = build_temperature_polarization_calculated_concentration_polarization_fixed()

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9005880857452575
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 300.2498
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 352.785
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.0053676595
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.68237393
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.8405494334

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
    m.fs.unit = MembraneDistillation0D(
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

    m.fs.unit.hot_ch.deltaP.fix(0)
    m.fs.unit.cold_ch.deltaP.fix(0)

    m.fs.unit.length.fix(8)
    m.fs.unit.hot_ch.channel_height.fix(0.0019)
    m.fs.unit.hot_ch.spacer_porosity.fix(0.77)
    m.fs.unit.cold_ch.channel_height.fix(0.0019)
    m.fs.unit.cold_ch.spacer_porosity.fix(0.77)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation0D_temperature_polarization_calculated_concentration_polarization_calculated_K_fixed(
    UnitTestHarness
):
    def configure(self):
        m = (
            build_temperature_polarization_calculated_concentration_polarization_calculated_K_fixed()
        )

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9015774522415636
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 300.377886
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 352.504663
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.00528521
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.6733963
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.83968273

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
    m.fs.properties_cold_ch = props_w.WaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation0D(
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
    m.fs.unit.cold_ch_inlet.pressure[0].fix(7e5)
    m.fs.unit.cold_ch_inlet.temperature[0].fix(273.15 + 25)

    m.fs.unit.hot_ch.deltaP.fix(membrane_pressure_drop)
    m.fs.unit.cold_ch.deltaP.fix(membrane_pressure_drop)

    m.fs.unit.length.fix(8)
    m.fs.unit.hot_ch.channel_height.fix(0.0019)
    m.fs.unit.hot_ch.spacer_porosity.fix(0.77)
    m.fs.unit.cold_ch.channel_height.fix(0.0019)
    m.fs.unit.cold_ch.spacer_porosity.fix(0.77)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation0D_temperature_polarization_calculated_concentration_polarization_calculated_K_calculated(
    UnitTestHarness
):
    def configure(self):
        m = (
            build_temperature_polarization_calculated_concentration_polarization_calculated_K_calculated()
        )

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9008354100672312
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 300.3940
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 352.8690612
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.00534704
        self.unit_solutions[m.fs.unit.recovery_mass[0]] = 0.066491803
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.6800786
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.841831

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
    m.fs.properties_cold_ch = props_w.WaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation0D(
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
    m.fs.unit.cold_ch_inlet.pressure[0].fix(7e5)
    m.fs.unit.cold_ch_inlet.temperature[0].fix(273.15 + 25)

    lenght = 8
    m.fs.unit.length.fix(lenght)

    m.fs.unit.hot_ch.dP_dx.fix(membrane_pressure_drop / lenght)
    m.fs.unit.cold_ch.dP_dx.fix(membrane_pressure_drop / lenght)

    m.fs.unit.hot_ch.channel_height.fix(0.0019)
    m.fs.unit.hot_ch.spacer_porosity.fix(0.77)
    m.fs.unit.cold_ch.channel_height.fix(0.0019)
    m.fs.unit.cold_ch.spacer_porosity.fix(0.77)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation0D_temperature_polarization_calculated_concentration_polarization_calculated_K_calculated_pressure_fixed_per_unit_length(
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
        ] = 0.9008354100672311
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 300.3940665
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 352.86906
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.00534704
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.6800786
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.8418317
        self.unit_solutions[m.fs.unit.hot_ch.deltaP[0]] = -500000.0
        self.unit_solutions[m.fs.unit.cold_ch.deltaP[0]] = -500000.0

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
    m.fs.properties_cold_ch = props_w.WaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation0D(
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
    lenght = 8
    m.fs.unit.length.fix(lenght)
    m.fs.unit.hot_ch.channel_height.fix(0.0019)
    m.fs.unit.hot_ch.spacer_porosity.fix(0.77)
    m.fs.unit.cold_ch.channel_height.fix(0.0019)
    m.fs.unit.cold_ch.spacer_porosity.fix(0.77)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation0D_temperature_polarization_calculated_concentration_polarization_calculated_K_calculated_pressure_calculated(
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
        ] = 0.9007196
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 300.52407
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 352.84163
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.00535669
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.6811818
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.841409
        self.unit_solutions[m.fs.unit.hot_ch.deltaP[0]] = -301454.81094
        self.unit_solutions[m.fs.unit.cold_ch.deltaP[0]] = -343467.7791

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
