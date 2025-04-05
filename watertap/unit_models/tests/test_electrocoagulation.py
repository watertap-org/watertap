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

from pyomo.environ import ConcreteModel, units as pyunits

from idaes.core import FlowsheetBlock, MaterialFlowBasis
from idaes.core import UnitModelCostingBlock

from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    set_scaling_factor,
)

from watertap.costing import WaterTAPCosting
from watertap.core.solvers import get_solver
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models.electrocoagulation import Electrocoagulation
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

solver = get_solver()


def build_ec1():
    """
    From Dubrawski (2014) paper
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
        overpotential_calculation="nernst",
    )
    set_scaling_factor(ec.properties_in[0].flow_mass_phase_comp["Liq", "H2O"], 1e2)
    set_scaling_factor(ec.properties_in[0].flow_mass_phase_comp["Liq", "TDS"], 1e5)

    set_scaling_factor(ec.properties_out[0].flow_mass_phase_comp["Liq", "H2O"], 1e2)
    set_scaling_factor(ec.properties_out[0].flow_mass_phase_comp["Liq", "TDS"], 1e5)

    set_scaling_factor(ec.properties_waste[0].flow_mass_phase_comp["Liq", "H2O"], 1e5)
    set_scaling_factor(ec.properties_waste[0].flow_mass_phase_comp["Liq", "TDS"], 1e7)

    # only because dealing with bench scale
    set_scaling_factor(m.fs.unit.electrode_volume, 1e4)
    set_scaling_factor(m.fs.unit.cell_volume, 1e4)
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

    ec.electrode_thick.fix(electrode_thickness)
    ec.cell_voltage.fix(10)
    ec.electrode_gap.fix(gap)
    ec.coagulant_dose.fix(metal_loading)
    ec.electrolysis_time.fix(electrolysis_time)
    ec.anode_area.fix(electrode_area_total)
    ec.floc_retention_time.fix(2)
    ec.tafel_slope_cathode.fix()
    ec.tafel_slope_anode.fix()

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


class TestEC1(UnitTestHarness):
    def configure(self):
        m = build_ec1()

        self.unit_solutions[m.fs.unit.tafel_slope_cathode] = 0.146
        self.unit_solutions[m.fs.unit.tafel_slope_anode] = 0.093
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
                "out": m.fs.unit.properties_in[0.0].flow_mass_phase_comp["Liq", "H2O"]
                + m.fs.unit.properties_in[0.0].flow_mass_phase_comp["Liq", "TDS"]
                + m.fs.unit.properties_waste[0.0].flow_mass_phase_comp["Liq", "TDS"],
            },
        }
        return m
