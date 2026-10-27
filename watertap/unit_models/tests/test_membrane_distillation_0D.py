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

    m.fs.properties_hot_ch.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.properties_hot_ch.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )
    m.fs.properties_cold_ch.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.properties_cold_ch.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )
    m.fs.properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.properties_vapor.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Vap", "H2O")
    )

    iscale.set_scaling_factor(
        m.fs.unit.cold_ch.properties_in[0].flow_mass_phase_comp["Vap", "H2O"], 1
    )
    iscale.set_scaling_factor(
        m.fs.unit.cold_ch.properties_out[0].flow_mass_phase_comp["Vap", "H2O"], 1
    )
    iscale.set_scaling_factor(
        m.fs.unit.cold_ch.properties_interface[0, 0].flow_mass_phase_comp["Vap", "H2O"],
        1,
    )
    iscale.set_scaling_factor(
        m.fs.unit.cold_ch.properties_interface[0, 1].flow_mass_phase_comp["Vap", "H2O"],
        1,
    )
    iscale.set_scaling_factor(m.fs.unit.area, 1e-1)
    iscale.set_scaling_factor(m.fs.unit.hot_ch.heat, 1e-3)
    iscale.set_scaling_factor(m.fs.unit.cold_ch.heat, 1e-3)
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
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 314.0520934
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 344.102435
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.00334578745
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.517108340
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.7069605503

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
        ] = 0.9033266209853699
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 303.6704079347907
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = (
            352.05089103799094
        )
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.01233
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.6633469
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.8292445

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
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 314.05609
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 344.10965
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.00334521
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.51706
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.70707

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
        ] = 0.9029749
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 302.03889
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 353.54345
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.0051688
        self.unit_solutions[m.fs.unit.recovery_mass[0]] = 0.0642747
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.65
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.852207

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
        ] = 0.903349276
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 302.08527
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 353.522566
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.00513756
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.6466382
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.8518856

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
        ] = 0.9042205
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 302.19553
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 353.471809
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.00506495
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.6387468
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.8511047

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
        ] = 0.903569
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 302.22345
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 353.60605
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.00511924
        self.unit_solutions[m.fs.unit.recovery_mass[0]] = 0.063659
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.6446025
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.85317

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
        ] = 0.903569
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 302.22345
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 353.60605
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.00511924
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.6446025
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.85317
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
        ] = 0.903455
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 302.169
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 353.5788
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.00512875
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.64568
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.85275
        self.unit_solutions[m.fs.unit.hot_ch.deltaP[0]] = -299896.723878
        self.unit_solutions[m.fs.unit.cold_ch.deltaP[0]] = -342521.37806

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
    m.fs.unit = MembraneDistillation0D(
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


class TestMembraneDisillation0D_temperature_polarization_none_vmd(UnitTestHarness):
    def configure(self):
        m = build_temperature_polarization_none_vmd()

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9285175510218185
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 338.6605219221607
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 338.6605219221607
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.037008

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
    m.fs.unit = MembraneDistillation0D(
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


class TestMembraneDisillation0D_temperature_polarization_fixed_hot_vmd(UnitTestHarness):
    def configure(self):
        m = build_temperature_polarization_fixed_hot_vmd()

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9479100277176233
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 352.655289
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 338.0086
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.017254657

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
    m.fs.unit = MembraneDistillation0D(
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
    m.fs.unit.length.fix(1)

    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    m.fs.unit.cold_ch_inlet.pressure[0].fix(10000)

    m.fs.unit.hot_ch.channel_height.fix(0.0019)
    m.fs.unit.hot_ch.spacer_porosity.fix(0.77)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation0D_temperature_polarization_calculated_hot_vmd(
    UnitTestHarness
):
    def configure(self):
        m = build_temperature_polarization_calculated_hot_vmd()

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9348449235681139
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 343.476969
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 340.669024
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.03057997

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
    m.fs.unit = MembraneDistillation0D(
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
    m.fs.unit.length.fix(1)

    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    m.fs.unit.cold_ch_inlet.pressure[0].fix(10000)

    m.fs.unit.hot_ch.channel_height.fix(0.0019)
    m.fs.unit.hot_ch.spacer_porosity.fix(0.77)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation0D_temperature_polarization_concentration_polarization_calculated_hot_vmd(
    UnitTestHarness
):
    def configure(self):
        m = (
            build_temperature_polarization_concentration_polarization_calculated_hot_vmd()
        )

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9348449235681139
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 343.558025
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 340.747587
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.03046524

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


#
def build_temperature_polarization_concentration_polarization_pressure_calculated_hot_vmd():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties_hot_ch = props_sw.SeawaterParameterBlock()
    m.fs.properties_cold_ch = props_w.WaterParameterBlock()
    m.fs.properties_vapor = props_w.WaterParameterBlock()
    m.fs.unit = MembraneDistillation0D(
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
    m.fs.unit.length.fix(1)

    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    m.fs.unit.cold_ch_inlet.pressure[0].fix(10000)

    m.fs.unit.hot_ch.channel_height.fix(0.0019)
    m.fs.unit.hot_ch.spacer_porosity.fix(0.77)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation0D_temperature_polarization_concentration_polarization_pressure_calculated_hot_vmd(
    UnitTestHarness
):
    def configure(self):
        m = (
            build_temperature_polarization_concentration_polarization_pressure_calculated_hot_vmd()
        )

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9348449235681139
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 343.56891
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 340.756717
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.03047074
        self.unit_solutions[m.fs.unit.hot_ch.deltaP[0]] = -73386.0456

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
    m.fs.unit = MembraneDistillation0D(
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

    # Fully specify the system
    hot_ch_flow_mass = 1.0
    hot_ch_mass_frac_TDS = 0.035
    hot_ch_pressure = 7e5
    membrane_area = 3.0
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
    m.fs.unit.gap_thickness.fix(0.0004)
    m.fs.unit.membrane_thermal_conductivity.fix(0.2)
    m.fs.unit.gap_thermal_conductivity.fix(0.06)

    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_H2O
    )
    m.fs.unit.cold_ch_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(
        hot_ch_flow_mass * hot_ch_mass_frac_TDS
    )
    m.fs.unit.cold_ch_inlet.temperature[0].fix(273.15 + 25)
    m.fs.unit.cold_ch_inlet.pressure[0].fix(100000)

    m.fs.unit.gap_ch_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(1e-8)
    m.fs.unit.gap_ch_inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(0)
    m.fs.unit.gap_ch_inlet.temperature[0].set_value(273.15 + 30)
    m.fs.unit.gap_ch_inlet.pressure[0].fix(100000)

    iscale.calculate_scaling_factors(m)

    return m


class TestMembraneDisillation0D_temperature_polarization_none_pgmd(UnitTestHarness):
    def configure(self):
        m = build_temperature_polarization_none_pgmd()

        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.9580840582054698
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 356.6334998833635
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = (
            304.69178989433595
        )
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.002318637
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
