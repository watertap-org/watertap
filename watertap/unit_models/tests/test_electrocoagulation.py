#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
import pytest

from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    value,
    units as pyunits,
)

from idaes.core import FlowsheetBlock, MaterialFlowBasis
from idaes.core import UnitModelCostingBlock

from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    set_scaling_factor,
)
from idaes.core.util.testing import initialization_tester

from watertap.costing import WaterTAPCosting
from watertap.core.util.initialization import check_dof
from watertap.core.solvers import get_solver
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models.electrocoagulation import Electrocoagulation
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

solver = get_solver()


def build_ec0():
    """
    Fixed overpotential model.
    """

    flow_vol_phase = 0.25 * pyunits.Mgallons / pyunits.day
    conc_tds = 20 * pyunits.kg / pyunits.m**3

    gap = 30 * pyunits.mm
    electrolysis_time = 5 * pyunits.min
    metal_loading = 200 * pyunits.mg / pyunits.liter
    electrode_thickness = 1 * pyunits.mm
    overpotential = 1.25
    floc_retention_time = 2
    cell_voltage = 8
    current_efficiency = 1.1

    ec_feed = {
        "solute_list": ["TDS"],
        "mw_data": {"TDS": 58.44e-3},  # NaCl
        "material_flow_basis": MaterialFlowBasis.mass,
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = MCASParameterBlock(**ec_feed)
    m.fs.unit = ec = Electrocoagulation(
        property_package=m.fs.properties,
        electrode_material="aluminum",
        overpotential_calculation="fixed",
    )

    set_scaling_factor(ec.properties_in[0].flow_mass_phase_comp["Liq", "H2O"], 0.1)
    set_scaling_factor(ec.properties_in[0].flow_mass_phase_comp["Liq", "TDS"], 10)

    set_scaling_factor(ec.properties_out[0].flow_mass_phase_comp["Liq", "H2O"], 0.1)
    set_scaling_factor(ec.properties_out[0].flow_mass_phase_comp["Liq", "TDS"], 10)

    set_scaling_factor(ec.properties_waste[0].flow_mass_phase_comp["Liq", "H2O"], 10)
    set_scaling_factor(ec.properties_waste[0].flow_mass_phase_comp["Liq", "TDS"], 100)

    set_scaling_factor(ec.applied_current, 1e-4)

    calculate_scaling_factors(m)

    m.fs.unit.properties_in.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): flow_vol_phase,
            ("conc_mass_phase_comp", ("Liq", "TDS")): conc_tds,
            ("temperature", None): 300,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )

    ec.electrode_thickness.fix(electrode_thickness)
    ec.floc_retention_time.fix(floc_retention_time)
    ec.electrode_gap.fix(gap)
    ec.electrolysis_time.fix(electrolysis_time)

    ec.cell_voltage.fix(cell_voltage)
    ec.current_efficiency.fix(current_efficiency)
    ec.coagulant_dose.fix(metal_loading)
    ec.overpotential.fix(overpotential)

    ec.removal_frac_mass_comp["TDS"].set_value(0.1)

    return m


def build_ec1():
    """
    From Dubrawski et al (2014) paper
    Fig 6A
    For E_cell = 10, I_cell = ~4 Amp
    http://dx.doi.org/10.1016/j.electacta.2014.02.089
    """

    flow_vol_phase = 0.1 * pyunits.liter / pyunits.minute
    conc_tds = 0.5 * pyunits.kg / pyunits.m**3

    gap = 2 * pyunits.mm
    electrode_area_total = 38 * pyunits.mm * 240 * pyunits.mm
    electrolysis_time = 1 * pyunits.min
    metal_loading = 87 * pyunits.mg / pyunits.liter
    electrode_thickness = 1 * pyunits.mm

    ec_feed = {
        "solute_list": ["TDS"],
        "mw_data": {"TDS": 58.44e-3},  # NaCl
        "material_flow_basis": MaterialFlowBasis.mass,
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = MCASParameterBlock(**ec_feed)
    m.fs.unit = ec = Electrocoagulation(
        property_package=m.fs.properties,
        electrode_material="iron",
        overpotential_calculation="detailed",
    )
    set_scaling_factor(ec.properties_in[0].flow_mass_phase_comp["Liq", "H2O"], 1e2)
    set_scaling_factor(ec.properties_in[0].flow_mass_phase_comp["Liq", "TDS"], 1e5)

    set_scaling_factor(ec.properties_out[0].flow_mass_phase_comp["Liq", "H2O"], 1e2)
    set_scaling_factor(ec.properties_out[0].flow_mass_phase_comp["Liq", "TDS"], 1e5)

    set_scaling_factor(ec.properties_waste[0].flow_mass_phase_comp["Liq", "H2O"], 1e5)
    set_scaling_factor(ec.properties_waste[0].flow_mass_phase_comp["Liq", "TDS"], 1e7)

    # only because dealing with bench scale
    set_scaling_factor(m.fs.unit.electrode_volume, 1e4)
    set_scaling_factor(m.fs.unit.reactor_volume, 1e4)
    set_scaling_factor(m.fs.unit.floc_basin_vol, 1e4)

    calculate_scaling_factors(m)

    m.fs.unit.properties_in.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): flow_vol_phase,
            ("conc_mass_phase_comp", ("Liq", "TDS")): conc_tds,
            ("temperature", None): 300,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )

    ec.electrode_thickness.fix(electrode_thickness)
    ec.cell_voltage.fix(10)
    ec.electrode_gap.fix(gap)
    ec.coagulant_dose.fix(metal_loading)
    ec.electrolysis_time.fix(electrolysis_time)
    ec.anode_area.fix(electrode_area_total)
    ec.floc_retention_time.fix(2)
    ec.tafel_slope_cathode.fix()
    ec.tafel_slope_anode.fix()

    return m


def build_ec2():
    """
    From Gu et al. (2009) paper
    Fig 6A
    For i = 8 mA/cm2, P_dF = ~15000 uW/cm2
    https://doi.org/10.1021/ie801086c
    """

    flow_vol_phase = 0.1 * pyunits.liter / pyunits.minute
    conc_tds = 1.05 * pyunits.kg / pyunits.m**3
    electrolysis_time = 1 * pyunits.min
    gap = 0.4 * pyunits.cm
    anode_area = 34 * pyunits.cm * 3.2 * pyunits.cm
    current_density = 8 * pyunits.milliamp / pyunits.cm**2
    electrode_thickness = 0.32 * pyunits.cm
    current_efficiency = 1.28
    k1 = 430
    k2 = 1000

    ec_feed = {
        "solute_list": ["TDS"],
        "mw_data": {"TDS": 58.44e-3},  # NaCl
        "material_flow_basis": MaterialFlowBasis.mass,
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = MCASParameterBlock(**ec_feed)
    m.fs.unit = ec = Electrocoagulation(
        property_package=m.fs.properties,
        electrode_material="aluminum",
        overpotential_calculation="regression",
    )
    set_scaling_factor(ec.properties_in[0].flow_mass_phase_comp["Liq", "H2O"], 1e2)
    set_scaling_factor(ec.properties_in[0].flow_mass_phase_comp["Liq", "TDS"], 1e5)

    set_scaling_factor(ec.properties_out[0].flow_mass_phase_comp["Liq", "H2O"], 1e2)
    set_scaling_factor(ec.properties_out[0].flow_mass_phase_comp["Liq", "TDS"], 1e5)

    set_scaling_factor(ec.properties_waste[0].flow_mass_phase_comp["Liq", "H2O"], 1e5)
    set_scaling_factor(ec.properties_waste[0].flow_mass_phase_comp["Liq", "TDS"], 1e7)

    set_scaling_factor(m.fs.unit.electrode_volume, 1e4)
    set_scaling_factor(m.fs.unit.reactor_volume, 1e4)
    set_scaling_factor(m.fs.unit.floc_basin_vol, 1e4)
    calculate_scaling_factors(m)

    m.fs.unit.properties_in.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): flow_vol_phase,
            ("conc_mass_phase_comp", ("Liq", "TDS")): conc_tds,
            ("temperature", None): 300,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )

    ec.electrode_thickness.fix(electrode_thickness)
    ec.electrode_gap.fix(gap)
    ec.floc_retention_time.fix(2)
    ec.electrolysis_time.fix(electrolysis_time)
    ec.current_density.fix(current_density)
    ec.current_efficiency.fix(current_efficiency)
    ec.anode_area.fix(anode_area)
    ec.overpotential_k1.fix(k1)
    ec.overpotential_k2.fix(k2)

    return m


def build_ec3():
    """
    From Gu et al. (2009) paper
    Fig 6B
    For i = 10 mA/cm2, P_dF = ~7400 uW/cm2
    https://doi.org/10.1021/ie801086c
    """

    flow_vol_phase = 0.1 * pyunits.liter / pyunits.minute
    conc_tds = 1.05 * pyunits.kg / pyunits.m**3
    electrolysis_time = 1 * pyunits.min
    gap = 0.4 * pyunits.cm
    anode_area = 34 * pyunits.cm * 3.2 * pyunits.cm
    current_density = 10 * pyunits.milliamp / pyunits.cm**2
    electrode_thickness = 0.32 * pyunits.cm
    current_efficiency = 1
    k1 = 75
    k2 = 600

    ec_feed = {
        "solute_list": ["TDS"],
        "mw_data": {"TDS": 58.44e-3},  # NaCl
        "material_flow_basis": MaterialFlowBasis.mass,
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = MCASParameterBlock(**ec_feed)
    m.fs.unit = ec = Electrocoagulation(
        property_package=m.fs.properties,
        electrode_material="iron",
        overpotential_calculation="regression",
    )
    set_scaling_factor(ec.properties_in[0].flow_mass_phase_comp["Liq", "H2O"], 1e2)
    set_scaling_factor(ec.properties_in[0].flow_mass_phase_comp["Liq", "TDS"], 1e5)

    set_scaling_factor(ec.properties_out[0].flow_mass_phase_comp["Liq", "H2O"], 1e2)
    set_scaling_factor(ec.properties_out[0].flow_mass_phase_comp["Liq", "TDS"], 1e5)

    set_scaling_factor(ec.properties_waste[0].flow_mass_phase_comp["Liq", "H2O"], 1e5)
    set_scaling_factor(ec.properties_waste[0].flow_mass_phase_comp["Liq", "TDS"], 1e7)

    set_scaling_factor(m.fs.unit.electrode_volume, 1e4)
    set_scaling_factor(m.fs.unit.reactor_volume, 1e4)
    set_scaling_factor(m.fs.unit.floc_basin_vol, 1e4)
    calculate_scaling_factors(m)

    m.fs.unit.properties_in.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): flow_vol_phase,
            ("conc_mass_phase_comp", ("Liq", "TDS")): conc_tds,
            ("temperature", None): 300,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )

    ec.electrode_thickness.fix(electrode_thickness)
    ec.electrode_gap.fix(gap)
    ec.floc_retention_time.fix(2)
    ec.electrolysis_time.fix(electrolysis_time)
    ec.current_density.fix(current_density)
    ec.current_efficiency.fix(current_efficiency)
    ec.anode_area.fix(anode_area)
    ec.overpotential_k1.fix(k1)
    ec.overpotential_k2.fix(k2)

    return m


def build_ec_costing():
    """
    Multi-component with fixed overpotential.
    Developed from NMSU Test Case #1
    """
    ec_feed = {
        "solute_list": ["TDS", "Foo_2+", "Bar_-"],
        "mw_data": {"TDS": 31.4038218e-3, "Foo_2+": 100e-3, "Bar_-": 60e-3},
        "material_flow_basis": "mass",
    }

    flow_vol_phase = 0.5 * pyunits.Mgallons / pyunits.day
    conc_tds = 130 * pyunits.kg / pyunits.m**3
    conc_foo = 1.625 * pyunits.kg / pyunits.m**3
    conc_bar = 3.821 * pyunits.kg / pyunits.m**3

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = MCASParameterBlock(**ec_feed)
    m.fs.unit = ec = Electrocoagulation(
        property_package=m.fs.properties,
        overpotential_calculation="fixed",
        reactor_material="stainless_steel",
    )

    set_scaling_factor(ec.properties_in[0].flow_mass_phase_comp["Liq", "H2O"], 1e-2)
    set_scaling_factor(ec.properties_in[0].flow_mass_phase_comp["Liq", "TDS"], 1e2)
    set_scaling_factor(ec.properties_in[0].flow_mass_phase_comp["Liq", "Foo_2+"], 1e2)
    set_scaling_factor(ec.properties_in[0].flow_mass_phase_comp["Liq", "Bar_-"], 1e1)

    set_scaling_factor(ec.properties_out[0].flow_mass_phase_comp["Liq", "H2O"], 1e-2)
    set_scaling_factor(ec.properties_out[0].flow_mass_phase_comp["Liq", "TDS"], 10)
    set_scaling_factor(ec.properties_out[0].flow_mass_phase_comp["Liq", "Foo_2+"], 1e3)
    set_scaling_factor(ec.properties_out[0].flow_mass_phase_comp["Liq", "Bar_-"], 1e3)

    set_scaling_factor(ec.properties_waste[0].flow_mass_phase_comp["Liq", "H2O"], 1)
    set_scaling_factor(ec.properties_waste[0].flow_mass_phase_comp["Liq", "TDS"], 1)
    set_scaling_factor(ec.properties_waste[0].flow_mass_phase_comp["Liq", "Foo_2+"], 1)
    set_scaling_factor(ec.properties_waste[0].flow_mass_phase_comp["Liq", "Bar_-"], 1)

    set_scaling_factor(ec.applied_current, 1e-4)

    calculate_scaling_factors(m)

    m.fs.unit.properties_in.calculate_state(
        var_args={
            ("flow_vol_phase", ("Liq")): flow_vol_phase,
            ("conc_mass_phase_comp", ("Liq", "TDS")): conc_tds,
            ("conc_mass_phase_comp", ("Liq", "Foo_2+")): conc_foo,
            ("conc_mass_phase_comp", ("Liq", "Bar_-")): conc_bar,
            ("temperature", None): 298,
            ("pressure", None): 101325,
        },
        hold_state=True,
    )

    ec.electrode_thickness.fix(0.02)
    ec.current_density.fix(300)
    ec.electrolysis_time.fix(30)
    ec.electrode_gap.fix(0.02)
    ec.current_efficiency.fix(1.5)
    ec.charge_loading_rate.fix(668)
    ec.cell_voltage.fix(10)
    ec.floc_retention_time.fix(12)
    ec.removal_frac_mass_comp["TDS"].set_value(0.1)
    ec.removal_frac_mass_comp["Foo_2+"].set_value(0.98)
    ec.removal_frac_mass_comp["Bar_-"].set_value(0.9)

    return m


class TestEC_noTDS:
    @pytest.mark.unit
    def test_no_tds_in_feed(self):
        ec_feed_no_tds = {
            "solute_list": ["foo", "bar", "baz"],
            "mw_data": {"foo": 10e-3, "bar": 222e-3, "baz": 39e-3},
            "charge": {"foo": 3},
        }
        error_msg = "TDS must be in feed stream for solution conductivity estimation."
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(**ec_feed_no_tds)
        with pytest.raises(ConfigurationError, match=error_msg):
            m.fs.unit = Electrocoagulation(property_package=m.fs.properties)


class TestEC0(UnitTestHarness):
    def configure(self):
        m = build_ec0()

        self.unit_solutions[m.fs.unit.coagulant_dose] = 0.2
        self.unit_solutions[m.fs.unit.electrode_thickness] = 0.001
        self.unit_solutions[m.fs.unit.electrode_mass] = 128.669
        self.unit_solutions[m.fs.unit.electrode_volume] = 0.04747
        self.unit_solutions[m.fs.unit.cathode_area] = 23.739
        self.unit_solutions[m.fs.unit.reactor_volume] = 3.2859
        self.unit_solutions[m.fs.unit.current_density] = 900
        self.unit_solutions[m.fs.unit.applied_current] = 21365.70
        self.unit_solutions[m.fs.unit.ohmic_resistance] = 0.007499
        self.unit_solutions[m.fs.unit.charge_loading_rate] = 1950.6

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.properties_in[0.0].flow_mass_phase_comp["Liq", "H2O"]
                + m.fs.unit.properties_in[0.0].flow_mass_phase_comp["Liq", "TDS"],
                "out": m.fs.unit.properties_out[0.0].flow_mass_phase_comp["Liq", "H2O"]
                + m.fs.unit.properties_out[0.0].flow_mass_phase_comp["Liq", "TDS"]
                + m.fs.unit.properties_waste[0.0].flow_mass_phase_comp["Liq", "H2O"]
                + m.fs.unit.properties_waste[0.0].flow_mass_phase_comp["Liq", "TDS"],
            },
        }
        return m


class TestEC1(UnitTestHarness):
    def configure(self):
        m = build_ec1()

        self.unit_solutions[m.fs.unit.tafel_slope_cathode] = 0.06339
        self.unit_solutions[m.fs.unit.tafel_slope_anode] = 0.04038
        self.unit_solutions[m.fs.unit.electrode_mass] = 0.1433663
        self.unit_solutions[m.fs.unit.electrode_volume] = 1.824e-05
        self.unit_solutions[m.fs.unit.current_density] = 421.9
        self.unit_solutions[m.fs.unit.applied_current] = 3.848  # ~4 amp
        self.unit_solutions[m.fs.unit.ohmic_resistance] = 0.02
        self.unit_solutions[m.fs.unit.charge_loading_rate] = 2309.0
        self.unit_solutions[m.fs.unit.current_efficiency] = 0.13019
        self.unit_solutions[m.fs.unit.cell_voltage] = 10
        self.unit_solutions[m.fs.unit.overpotential] = 1.5604

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.properties_in[0.0].flow_mass_phase_comp["Liq", "H2O"]
                + m.fs.unit.properties_in[0.0].flow_mass_phase_comp["Liq", "TDS"],
                "out": m.fs.unit.properties_out[0.0].flow_mass_phase_comp["Liq", "H2O"]
                + m.fs.unit.properties_out[0.0].flow_mass_phase_comp["Liq", "TDS"]
                + m.fs.unit.properties_waste[0.0].flow_mass_phase_comp["Liq", "H2O"]
                + m.fs.unit.properties_waste[0.0].flow_mass_phase_comp["Liq", "TDS"],
            },
        }
        return m


class TestEC2(UnitTestHarness):
    def configure(self):
        m = build_ec2()

        self.unit_solutions[m.fs.unit.coagulant_dose] = 0.062307
        self.unit_solutions[m.fs.unit.electrode_mass] = 0.1887
        self.unit_solutions[m.fs.unit.electrode_volume] = 6.963e-05
        self.unit_solutions[m.fs.unit.current_density] = 80
        self.unit_solutions[m.fs.unit.applied_current] = 0.8704
        self.unit_solutions[m.fs.unit.ohmic_resistance] = 0.0190476
        self.unit_solutions[m.fs.unit.charge_loading_rate] = 522.24
        self.unit_solutions[m.fs.unit.current_efficiency] = 1.28
        self.unit_solutions[m.fs.unit.cell_voltage] = 3.4179
        self.unit_solutions[m.fs.unit.overpotential] = 1.8941
        self.unit_solutions[m.fs.unit.power_density_faradaic] = (
            15153.27  # ~15000 uW/cm2
        )

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.properties_in[0.0].flow_mass_phase_comp["Liq", "H2O"]
                + m.fs.unit.properties_in[0.0].flow_mass_phase_comp["Liq", "TDS"],
                "out": m.fs.unit.properties_out[0.0].flow_mass_phase_comp["Liq", "H2O"]
                + m.fs.unit.properties_out[0.0].flow_mass_phase_comp["Liq", "TDS"]
                + m.fs.unit.properties_waste[0.0].flow_mass_phase_comp["Liq", "H2O"]
                + m.fs.unit.properties_waste[0.0].flow_mass_phase_comp["Liq", "TDS"],
            },
        }

        return m


class TestEC3(UnitTestHarness):
    def configure(self):
        m = build_ec3()

        self.unit_solutions[m.fs.unit.coagulant_dose] = 0.18891
        self.unit_solutions[m.fs.unit.electrode_mass] = 0.54730
        self.unit_solutions[m.fs.unit.electrode_volume] = 6.963e-05
        self.unit_solutions[m.fs.unit.current_density] = 100
        self.unit_solutions[m.fs.unit.applied_current] = 1.088
        self.unit_solutions[m.fs.unit.ohmic_resistance] = 0.0190476
        self.unit_solutions[m.fs.unit.charge_loading_rate] = 652.8
        self.unit_solutions[m.fs.unit.current_efficiency] = 1
        self.unit_solutions[m.fs.unit.cell_voltage] = 2.67745
        self.unit_solutions[m.fs.unit.overpotential] = 0.77269
        self.unit_solutions[m.fs.unit.power_density_faradaic] = 7726.93  # ~7400 uW/cm2

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.properties_in[0.0].flow_mass_phase_comp["Liq", "H2O"]
                + m.fs.unit.properties_in[0.0].flow_mass_phase_comp["Liq", "TDS"],
                "out": m.fs.unit.properties_out[0.0].flow_mass_phase_comp["Liq", "H2O"]
                + m.fs.unit.properties_out[0.0].flow_mass_phase_comp["Liq", "TDS"]
                + m.fs.unit.properties_waste[0.0].flow_mass_phase_comp["Liq", "H2O"]
                + m.fs.unit.properties_waste[0.0].flow_mass_phase_comp["Liq", "TDS"],
            },
        }

        return m


class TestECCosting(UnitTestHarness):
    def configure(self):
        m = build_ec_costing()

        self.unit_solutions[m.fs.unit.floc_basin_vol] = 15.772
        self.unit_solutions[m.fs.unit.coagulant_dose] = 0.093395
        self.unit_solutions[m.fs.unit.electrode_mass] = 5287.542
        self.unit_solutions[m.fs.unit.electrode_volume] = 1.9511
        self.unit_solutions[m.fs.unit.reactor_volume] = 39.43137
        self.unit_solutions[m.fs.unit.applied_current] = 14633.42
        self.unit_solutions[m.fs.unit.ohmic_resistance] = 0.0007692
        self.unit_solutions[m.fs.unit.overpotential] = 9.7692

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.properties_in[0.0].flow_mass_phase_comp["Liq", "H2O"]
                + m.fs.unit.properties_in[0.0].flow_mass_phase_comp["Liq", "TDS"]
                + m.fs.unit.properties_in[0.0].flow_mass_phase_comp["Liq", "Foo_2+"]
                + m.fs.unit.properties_in[0.0].flow_mass_phase_comp["Liq", "Bar_-"],
                "out": m.fs.unit.properties_out[0.0].flow_mass_phase_comp["Liq", "H2O"]
                + m.fs.unit.properties_out[0.0].flow_mass_phase_comp["Liq", "TDS"]
                + m.fs.unit.properties_out[0.0].flow_mass_phase_comp["Liq", "Foo_2+"]
                + m.fs.unit.properties_out[0.0].flow_mass_phase_comp["Liq", "Bar_-"]
                + m.fs.unit.properties_waste[0.0].flow_mass_phase_comp["Liq", "H2O"]
                + m.fs.unit.properties_waste[0.0].flow_mass_phase_comp["Liq", "TDS"]
                + m.fs.unit.properties_waste[0.0].flow_mass_phase_comp["Liq", "Foo_2+"]
                + m.fs.unit.properties_waste[0.0].flow_mass_phase_comp["Liq", "Bar_-"],
            },
        }
        return m

    @pytest.mark.component
    def test_costing(self):
        m = build_ec_costing()
        m.fs.costing = WaterTAPCosting()
        m.fs.costing.base_currency = pyunits.USD_2023
        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.electrocoagulation.sludge_handling_cost.fix(0.025)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.properties_out[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            m.fs.unit.properties_out[0].flow_vol_phase["Liq"], name="SEC"
        )
        check_dof(m, fail_flag=True)
        initialization_tester(m)
        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_capital_cost": 1168409.07,
            "aggregate_fixed_operating_cost": 311632.55,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 146.33,
            "aggregate_flow_aluminum": 64565.49,
            "aggregate_flow_costs": {"electricity": 118796.73, "aluminum": 162263.39},
            "total_capital_cost": 1168409.07,
            "total_operating_cost": 599638.93,
            "LCOW": 1.1831,
            "SEC": 1.9064,
        }

        for v, r in sys_cost_results.items():
            mv = m.fs.costing.find_component(v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

        ec_cost_results = {
            "capital_cost": 1168409.07,
            "fixed_operating_cost": 311632.55,
            "capital_cost_reactor": 1034159.79,
            "capital_cost_electrodes": 23582.44,
            "capital_cost_power_supply": 101782.18,
            "capital_cost_floc_reactor": 8884.65,
            "annual_sludge_management": 311632.55,
            "direct_capital_cost": 283511.02,
            "annual_electrode_replacement_mass_flow": 64565.49,
        }

        for v, r in ec_cost_results.items():
            mv = m.fs.unit.costing.find_component(v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)
