from pyomo.environ import ConcreteModel

from idaes.core.solvers import get_solver


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
        ] = 0.9294923084883953
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 313.4770191978068
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 342.8447987723201
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.0029589742926337237
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.46958331303607553
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.6876122888049245

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
        ] = 0.9244443003141366
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 307.4425600993805
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 347.7273575270859
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.008111139937172647
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.4809369122475431
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.7627285773397838

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
        ] = 0.929497268800918
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 313.4706323418694
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 342.841898880965
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.0029585609332568297
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.46956031608126614
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.6875676750917696

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
        ] = 0.9233923911250498
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 304.337512844963
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 350.4583511214421
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.0034673007395
        self.unit_solutions[m.fs.unit.recovery_mass[0]] = 0.0431166931346
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.467514426654
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.804743863406

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
        ] = 0.923602738825
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 304.3575352876
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 350.4506210584
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.0034497717645

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
        ] = 0.9237376005292054
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 304.3676091577648
        self.unit_solutions[
            m.fs.unit.cold_ch_outlet.temperature[0]
        ] = 350.44794988928794
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.003438533289232888
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.4639450335411154
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.8045838444505841

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
        ] = 0.9234150642000756
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 304.4323257757313
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 350.5701165470125
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.0034654113166603566
        self.unit_solutions[m.fs.unit.recovery_mass[0]] = 0.04309319772012876
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.4671469978106365
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.8064633314925004

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
        ] = 0.923415064200076
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 304.4323257757313
        self.unit_solutions[m.fs.unit.cold_ch_outlet.temperature[0]] = 350.5701165470125
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.0034654113166603523
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.4671469978106365
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.8064633314925004
        self.unit_solutions[m.fs.unit.hot_ch.deltaP_channel[0]] = -500000.0
        self.unit_solutions[m.fs.unit.cold_ch.deltaP_channel[0]] = -500000.0

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
        ] = 0.9234395262038794
        self.unit_solutions[
            m.fs.unit.hot_ch_outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ] = 0.035
        self.unit_solutions[m.fs.unit.hot_ch_outlet.temperature[0]] = 304.40094664041476
        self.unit_solutions[
            m.fs.unit.cold_ch_outlet.temperature[0]
        ] = 350.52705264537326
        self.unit_solutions[m.fs.unit.flux_mass_avg[0]] = 0.0034633728163433928
        self.unit_solutions[m.fs.unit.thermal_efficiency[0]] = 0.4669795564991426
        self.unit_solutions[m.fs.unit.effectiveness[0]] = 0.8058008099288198
        self.unit_solutions[m.fs.unit.hot_ch.deltaP_channel[0]] = -297491.6707887916
        self.unit_solutions[m.fs.unit.cold_ch.deltaP_channel[0]] = -332494.6036866968
        return m
