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

from pyomo.environ import (
    ConcreteModel,
    value,
)
from idaes.core import (
    FlowsheetBlock,
)
from watertap.unit_models.osmotically_assisted_reverse_osmosis_0D import (
    OsmoticallyAssistedReverseOsmosis0D,
)
from watertap.unit_models.reverse_osmosis_base import TransportModel, ModuleType
import watertap.property_models.NaCl_prop_pack as props
from watertap.core.solvers import get_solver
from idaes.core.util.scaling import (
    calculate_scaling_factors,
)

from watertap.core import (
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)

from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = OsmoticallyAssistedReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.fixed,
        mass_transfer_coefficient=MassTransferCoefficient.none,
        has_full_reporting=True,
    )

    # fully specify system
    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.035
    feed_pressure = 50e5
    feed_temperature = 273.15 + 25
    membrane_area = 50
    A = 4.2e-12
    B = 1.3e-8
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

    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        0.16189443280346605
    )
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        0.0024736353420558967
    )
    m.fs.unit.permeate_inlet.pressure[0].fix(101325)
    m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)

    m.fs.unit.permeate_side.cp_modulus.fix(permeate_cp_mod)
    m.fs.unit.feed_side.deltaP.fix(0)
    m.fs.unit.permeate_side.deltaP.fix(0)

    # scaling
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    calculate_scaling_factors(m)

    return m


class TestOsmoticallyAssistedReverseOsmosis(UnitTestHarness):
    def configure(self):
        m = build()
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            0.00671211
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            5.2729e-7
        )
        self.unit_solutions[
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.62939
        self.unit_solutions[
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
        ] = 0.03497
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
            / m.fs.unit.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = value(m.fs.unit.feed_side.cp_modulus[0, 0.0, "NaCl"])
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
            / m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = value(m.fs.unit.feed_side.cp_modulus[0, 1.0, "NaCl"])
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
            / m.fs.unit.permeate_side.properties_out[0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = value(m.fs.unit.permeate_side.cp_modulus[0, 0.0, "NaCl"])
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
            / m.fs.unit.permeate_side.properties_in[0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = value(m.fs.unit.permeate_side.cp_modulus[0, 1.0, "NaCl"])
        self.unit_solutions[m.fs.unit.feed_side.deltaP[0]] = 0
        self.unit_solutions[m.fs.unit.permeate_side.deltaP[0]] = 0

        return m


def build_SKK():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = OsmoticallyAssistedReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.fixed,
        mass_transfer_coefficient=MassTransferCoefficient.none,
        transport_model=TransportModel.SKK,
        has_full_reporting=True,
    )

    # fully specify system
    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.035
    feed_pressure = 50e5
    feed_temperature = 273.15 + 25
    membrane_area = 50
    A = 4.2e-12
    B = 1.3e-8
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
    m.fs.unit.reflect_coeff.fix(0.9)
    m.fs.unit.feed_side.cp_modulus.fix(feed_cp_mod)

    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        0.16189443280346605
    )
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        0.0024736353420558967
    )
    m.fs.unit.permeate_inlet.pressure[0].fix(101325)
    m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)

    m.fs.unit.permeate_side.cp_modulus.fix(permeate_cp_mod)
    m.fs.unit.feed_side.deltaP.fix(0)
    m.fs.unit.permeate_side.deltaP.fix(0)

    # scaling
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    calculate_scaling_factors(m)

    return m


class TestOsmoticallyAssistedReverseOsmosis_SKK(UnitTestHarness):
    def configure(self):
        m = build_SKK()

        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            7.971e-03
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            3.801e-5
        )
        self.unit_solutions[
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.56644
        self.unit_solutions[
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
        ] = 0.03309
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
            / m.fs.unit.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = value(m.fs.unit.feed_side.cp_modulus[0, 0.0, "NaCl"])
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
            / m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = value(m.fs.unit.feed_side.cp_modulus[0, 1.0, "NaCl"])
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
            / m.fs.unit.permeate_side.properties_out[0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = value(m.fs.unit.permeate_side.cp_modulus[0, 0.0, "NaCl"])
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
            / m.fs.unit.permeate_side.properties_in[0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = value(m.fs.unit.permeate_side.cp_modulus[0, 1.0, "NaCl"])
        self.unit_solutions[m.fs.unit.feed_side.deltaP[0]] = 0
        self.unit_solutions[m.fs.unit.permeate_side.deltaP[0]] = 0

        return m


def build_CP_calculation_with_kf_fixed():
    """
    Testing 0D-OARO with ConcentrationPolarizationType.calculated option enabled.
    This option makes use of an alternative constraint for the feed-side, membrane-interface concentration.
    Additionally, feed and permeate-side mass trasnfer coefficients are fixed.
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

    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        0.251513390470481
    )
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        0.0024704183950224
    )
    m.fs.unit.permeate_inlet.pressure[0].fix(151325.0)
    m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)

    m.fs.unit.feed_side.deltaP.fix(membrane_pressure_drop)
    m.fs.unit.permeate_side.deltaP.fix(membrane_pressure_drop)

    m.fs.unit.feed_side.K[0, 0.0, "NaCl"].fix(kf)
    m.fs.unit.feed_side.K[0, 1.0, "NaCl"].fix(kf)
    m.fs.unit.permeate_side.K[0, 0.0, "NaCl"].fix(kp)
    m.fs.unit.permeate_side.K[0, 1.0, "NaCl"].fix(kp)

    # scaling
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    calculate_scaling_factors(m)

    return m


class TestOsmoticallyAssistedReverseOsmosis_CP_calculation_with_kf_fixed(
    UnitTestHarness
):
    def configure(self):
        m = build_CP_calculation_with_kf_fixed()
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            4.9197e-3
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            5.9163e-7
        )
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.2515
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 0.00247
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_out[0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.4975
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_out[0].flow_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 0.0025
        self.unit_solutions[
            m.fs.unit.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 35.7511
        self.unit_solutions[
            m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 47.775
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 43.656
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 53.425
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 9.7496
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_out[0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 4.9939
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_out[0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = value(
            self.unit_solutions[
                m.fs.unit.permeate_side.properties[0, 0].conc_mass_phase_comp[
                    "Liq", "NaCl"
                ]
            ]
        )
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = value(
            self.unit_solutions[
                m.fs.unit.permeate_side.properties[0, 1].conc_mass_phase_comp[
                    "Liq", "NaCl"
                ]
            ]
        )
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 1.3745
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 4.6858
        self.unit_solutions[m.fs.unit.feed_side.cp_modulus[0, 0, "NaCl"]] = 1.221
        self.unit_solutions[m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"]] = 1.118
        self.unit_solutions[m.fs.unit.permeate_side.cp_modulus[0, 0, "NaCl"]] = 0.2752
        self.unit_solutions[m.fs.unit.permeate_side.cp_modulus[0, 1, "NaCl"]] = 0.4806

        return m


def build_CP_calculation_with_kf_calculation():
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
    membrane_pressure_drop = -0.5e5
    membrane_area = 50
    A = 4.2e-12
    B = 1.3e-8
    pressure_atmospheric = 101325

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

    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        0.7428108105680903
    )
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        0.004970619823292316
    )
    m.fs.unit.permeate_inlet.pressure[0].fix(151325.0)
    m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)

    m.fs.unit.structural_parameter.fix(300e-6)

    m.fs.unit.feed_side.deltaP.fix(membrane_pressure_drop)
    m.fs.unit.permeate_side.deltaP.fix(membrane_pressure_drop)
    m.fs.unit.permeate_side.channel_height.fix(0.001)
    m.fs.unit.permeate_side.spacer_porosity.fix(0.75)
    m.fs.unit.feed_side.channel_height.fix(0.002)
    m.fs.unit.feed_side.spacer_porosity.fix(0.75)
    length = 20
    m.fs.unit.length.fix(length)

    # scaling
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    calculate_scaling_factors(m)

    return m


class TestOsmoticallyAssistedReverseOsmosis_CP_calculation_with_kf_calculation(
    UnitTestHarness
):
    def configure(self):
        m = build_CP_calculation_with_kf_calculation()

        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            5.044e-3
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            5.876e-7
        )
        self.unit_solutions[
            m.fs.unit.feed_side.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 35.75
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 42.23
        self.unit_solutions[
            m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 48.185
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 52.783
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 6.647
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 3.364
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_out[0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 4.9939
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 1.2509
        self.unit_solutions[
            m.fs.unit.permeate_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = value(
            self.unit_solutions[
                m.fs.unit.permeate_side.properties_out[0].conc_mass_phase_comp[
                    "Liq", "NaCl"
                ]
            ]
        )
        self.unit_solutions[
            m.fs.unit.permeate_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        ] = value(
            self.unit_solutions[
                m.fs.unit.permeate_side.properties_in[0].conc_mass_phase_comp[
                    "Liq", "NaCl"
                ]
            ]
        )

        return m


def build_Pdrop_calculation():
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
    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.035
    feed_pressure = 50e5
    feed_temperature = 273.15 + 25
    membrane_area = 50
    A = 4.2e-12
    B = 1.3e-8
    pressure_atmospheric = 101325

    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)

    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        0.7719873688757197
    )
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        0.004970483059191647
    )
    m.fs.unit.permeate_inlet.pressure[0].fix(405467.2130804832)
    m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)

    m.fs.unit.area.fix(membrane_area)

    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)

    m.fs.unit.structural_parameter.fix(300e-6)

    m.fs.unit.permeate_side.channel_height.fix(0.001)
    m.fs.unit.permeate_side.spacer_porosity.fix(0.75)
    m.fs.unit.feed_side.channel_height.fix(0.002)
    m.fs.unit.feed_side.spacer_porosity.fix(0.75)
    m.fs.unit.feed_side.velocity[0, 0].fix(0.1)

    # scaling
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    calculate_scaling_factors(m)

    return m


class TestOsmoticallyAssistedReverseOsmosis_Pdrop_calculation(UnitTestHarness):
    def configure(self):
        m = build_Pdrop_calculation()

        self.unit_solutions[m.fs.unit.feed_side.deltaP[0]] = -0.3915e5
        self.unit_solutions[m.fs.unit.permeate_side.deltaP[0]] = -3.0414e5
        self.unit_solutions[m.fs.unit.feed_side.deltaP[0] / m.fs.unit.length] = -5109.74
        self.unit_solutions[m.fs.unit.permeate_side.deltaP[0] / m.fs.unit.length] = (
            -39700.326
        )
        self.unit_solutions[m.fs.unit.feed_side.N_Re[0, 0.0]] = 145.197
        self.unit_solutions[m.fs.unit.feed_side.velocity[0, 0.0]] = 0.1
        self.unit_solutions[m.fs.unit.feed_side.N_Re[0, 1.0]] = 110.557
        self.unit_solutions[m.fs.unit.feed_side.velocity[0, 1.0]] = 0.0771
        self.unit_solutions[m.fs.unit.permeate_side.N_Re[0, 0.0]] = 154.650
        self.unit_solutions[m.fs.unit.permeate_side.velocity[0, 0.0]] = 0.2045
        self.unit_solutions[m.fs.unit.permeate_side.N_Re[0, 1.0]] = 119.793
        self.unit_solutions[m.fs.unit.permeate_side.velocity[0, 1.0]] = 0.1588
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            4.460e-3
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            5.903e-7
        )
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 44.225
        self.unit_solutions[
            m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 46.316
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 51.550
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 3.562
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 1.392
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_out[0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 4.994

        return m


def build_Pdrop_fixed_per_unit_length():
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
    feed_pressure = 50e5
    feed_temperature = 273.15 + 25
    membrane_area = 50
    A = 4.2e-12
    B = 1.3e-8
    pressure_atmospheric = 101325
    feed_pressure_drop = -0.3915e5
    perm_pressure_drop = -3.0414e5
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

    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        0.7710753558388397
    )
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        0.0049705168627404974
    )
    m.fs.unit.permeate_inlet.pressure[0].fix(405465.0)
    m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)

    m.fs.unit.structural_parameter.fix(300e-6)

    m.fs.unit.permeate_side.channel_height.fix(0.001)
    m.fs.unit.permeate_side.spacer_porosity.fix(0.75)
    m.fs.unit.feed_side.channel_height.fix(0.002)
    m.fs.unit.feed_side.spacer_porosity.fix(0.75)
    length = 8
    m.fs.unit.length.fix(length)
    m.fs.unit.feed_side.dP_dx.fix(feed_pressure_drop / length)
    m.fs.unit.permeate_side.dP_dx.fix(perm_pressure_drop / length)

    # scaling
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )
    calculate_scaling_factors(m)

    return m


class TestOsmoticallyAssistedReverseOsmosis_Pdrop_fixed_per_unit_length(
    UnitTestHarness
):
    def configure(self):
        m = build_Pdrop_fixed_per_unit_length()

        self.unit_solutions[m.fs.unit.feed_side.deltaP[0]] = -39150
        self.unit_solutions[m.fs.unit.permeate_side.deltaP[0]] = -304140
        self.unit_solutions[m.fs.unit.feed_side.N_Re[0, 0.0]] = 151.623
        self.unit_solutions[m.fs.unit.feed_side.N_Re[0, 1.0]] = 115.302
        self.unit_solutions[m.fs.unit.permeate_side.N_Re[0, 0.0]] = 161.494
        self.unit_solutions[m.fs.unit.permeate_side.N_Re[0, 1.0]] = 124.946
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            4.478e-3
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            5.900e-7
        )
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 44.128
        self.unit_solutions[
            m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 46.372
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 51.539
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 3.565
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 1.385
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_out[0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 4.994

        return m


def build_friction_factor_spiral_wound():
    """Testing 0D-OARO with FrictionFactor.spiral_wound option."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = OsmoticallyAssistedReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        pressure_change_type=PressureChangeType.calculated,
        friction_factor=FrictionFactor.spiral_wound,
    )

    # fully specify system
    feed_flow_mass = 1
    feed_mass_frac_NaCl = 0.035
    feed_pressure = 50e5
    feed_temperature = 273.15 + 25
    membrane_area = 50
    A = 4.2e-12
    B = 1.3e-8
    pressure_atmospheric = 101325

    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)

    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        0.7719873688757197
    )
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        0.004970483059191647
    )
    m.fs.unit.permeate_inlet.pressure[0].fix(405467.2130804832)
    m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)

    m.fs.unit.area.fix(membrane_area)

    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)

    m.fs.unit.structural_parameter.fix(300e-6)

    m.fs.unit.permeate_side.channel_height.fix(0.001)
    m.fs.unit.permeate_side.spacer_porosity.fix(0.75)
    m.fs.unit.feed_side.channel_height.fix(0.002)
    m.fs.unit.feed_side.spacer_porosity.fix(0.75)
    m.fs.unit.feed_side.velocity[0, 0].fix(0.1)

    # scaling
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    calculate_scaling_factors(m)

    return m


class TestOsmoticallyAssistedReverseOsmosis_friction_factor_spiral_wound(
    UnitTestHarness
):
    def configure(self):
        m = build_friction_factor_spiral_wound()

        self.unit_solutions[m.fs.unit.feed_side.deltaP[0]] = -0.3019e5
        self.unit_solutions[m.fs.unit.permeate_side.deltaP[0]] = -2.4123e5
        self.unit_solutions[m.fs.unit.feed_side.deltaP[0] / m.fs.unit.length] = -3940.35
        self.unit_solutions[m.fs.unit.permeate_side.deltaP[0] / m.fs.unit.length] = (
            -31488.754
        )
        self.unit_solutions[m.fs.unit.feed_side.N_Re[0, 0.0]] = 145.197
        self.unit_solutions[m.fs.unit.feed_side.velocity[0, 0.0]] = 0.1
        self.unit_solutions[m.fs.unit.feed_side.N_Re[0, 1.0]] = 110.973
        self.unit_solutions[m.fs.unit.feed_side.velocity[0, 1.0]] = 0.0774
        self.unit_solutions[m.fs.unit.permeate_side.N_Re[0, 0.0]] = 154.231
        self.unit_solutions[m.fs.unit.permeate_side.velocity[0, 0.0]] = 0.2040
        self.unit_solutions[m.fs.unit.permeate_side.N_Re[0, 1.0]] = 119.793
        self.unit_solutions[m.fs.unit.permeate_side.velocity[0, 1.0]] = 0.1588
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            4.407e-3
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            5.881e-7
        )
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 43.978
        self.unit_solutions[
            m.fs.unit.feed_side.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 46.152
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 51.466
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 3.520
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 0.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 1.443
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_out[0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 5.007

        return m
