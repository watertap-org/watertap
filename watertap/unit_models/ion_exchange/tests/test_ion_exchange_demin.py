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
    value,
    assert_optimal_termination,
    units as pyunits,
)

from idaes.core import (
    FlowsheetBlock,
    UnitModelCostingBlock,
    MaterialFlowBasis,
)
import idaes.core.util.scaling as iscale
from idaes.core.util.testing import initialization_tester

from watertap.costing import WaterTAPCosting
from watertap.core.util.initialization import check_dof
from watertap.core.solvers import get_solver
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models import IonExchangeDemin
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

__author__ = "Kurban Sitterley"

solver = get_solver()
zero = 1e-8
relative_tolerance = 1e-3


def build_cation_exchange():

    # example is taken from WBS cation exchange model
    # https://www.epa.gov/sdwa/drinking-water-treatment-technology-unit-cost-models
    # using all defaults except flow rate
    # Q = 1 MGD
    # Cin = 200 mg/L as CaCO3
    # EBCT = 2.5 min
    # loading rate = 9.02 gpm / ft2
    # resin capacity = 27 kilograins / ft3
    #   1 kilograin / ft3 = 0.0458 eq/L
    #   strong acid polystyrenic macroporous

    ion_props = dict(
        solute_list=["Ca_2+", "Inert", "Mg_2+"],
        material_flow_basis=MaterialFlowBasis.mass,
        mw_data={
            "H2O": 0.018,
            "Ca_2+": 0.04,
            "Mg_2+": 0.024,
            "Inert": 0.10,
        },
        charge={"Ca_2+": 2, "Inert": 0, "Mg_2+": 2},
    )

    conc_in = {
        "Ca_2+": 0.04,
        "Mg_2+": 0.0243,
        "Inert": 0.3,
    }

    var_args = dict()
    for k, v in conc_in.items():
        var_args[("conc_mass_phase_comp", ("Liq", k))] = v * pyunits.g / pyunits.liter
    var_args[("flow_vol_phase", ("Liq"))] = 1 * pyunits.Mgallons / pyunits.day
    var_args[("temperature", (None))] = 298
    var_args[("pressure", (None))] = 101325

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**ion_props)
    ix_config = {
        "property_package": m.fs.properties,
    }
    m.fs.unit = IonExchangeDemin(**ix_config)
    pf = m.fs.unit.process_flow

    pf.properties_in.calculate_state(var_args=var_args, hold_state=True)
    pf.properties_in[0].total_hardness

    loading_rate = 9.02 * pyunits.gallon / pyunits.minute / pyunits.ft**2

    m.fs.unit.loading_rate.fix(loading_rate)
    m.fs.unit.ebct.fix(150)
    m.fs.unit.resin_capacity_ax.fix(1.2366)
    m.fs.unit.resin_capacity_cx.fix(1.2366)
    m.fs.unit.resin_diam.fix()
    m.fs.unit.resin_density.fix()
    m.fs.unit.bed_porosity.fix()
    m.fs.unit.number_columns.fix(2)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "Ca_2+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "Mg_2+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "Inert")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 0.1, index=("Liq", "H2O")
    )
    iscale.set_scaling_factor(
        pf.properties_out[0.0].flow_mass_phase_comp["Liq", "Ca_2+"],
        1e6,
    )
    iscale.set_scaling_factor(
        pf.properties_out[0.0].flow_mass_phase_comp["Liq", "Mg_2+"],
        1e6,
    )

    iscale.calculate_scaling_factors(m)

    return m


def build_anion_exchange():

    ion_props = dict(
        solute_list=["Cl_-", "SO4_2-"],
        material_flow_basis=MaterialFlowBasis.mass,
        mw_data={
            "H2O": 0.018,
            "Cl_-": 0.035,
            "SO4_2-": 0.096,
        },
        charge={"Cl_-": -1, "SO4_2-": -2},
    )

    conc_in = {
        "Cl_-": 0.02,
        "SO4_2-": 0.013,
    }

    var_args = dict()
    for k, v in conc_in.items():
        var_args[("conc_mass_phase_comp", ("Liq", k))] = v * pyunits.g / pyunits.liter
    var_args[("flow_vol_phase", ("Liq"))] = 10 * pyunits.Mgallons / pyunits.day
    var_args[("temperature", (None))] = 298
    var_args[("pressure", (None))] = 101325

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**ion_props)
    ix_config = {
        "property_package": m.fs.properties,
    }
    m.fs.unit = IonExchangeDemin(**ix_config)
    pf = m.fs.unit.process_flow

    pf.properties_in.calculate_state(var_args=var_args, hold_state=True)

    loading_rate = 6.5 * pyunits.gallon / pyunits.minute / pyunits.ft**2

    m.fs.unit.loading_rate.fix(loading_rate)
    m.fs.unit.ebct.fix(300)
    m.fs.unit.resin_capacity_ax.fix(1.2)
    m.fs.unit.resin_capacity_cx.fix(1)
    m.fs.unit.resin_diam.fix()
    m.fs.unit.resin_density.fix(800)
    m.fs.unit.bed_porosity.fix()
    m.fs.unit.number_columns.fix(10)
    m.fs.unit.regen_dose.set_value(300)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "Cl_-")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "SO4_2-")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 0.1, index=("Liq", "H2O")
    )
    iscale.set_scaling_factor(
        pf.properties_out[0.0].flow_mass_phase_comp["Liq", "Cl_-"],
        1e6,
    )
    iscale.set_scaling_factor(
        pf.properties_out[0.0].flow_mass_phase_comp["Liq", "SO4_2-"],
        1e6,
    )

    iscale.calculate_scaling_factors(m)
    m.fs.unit.regeneration_stream[0].conc_mass_phase_comp

    return m


def build_mixed_bed():

    ion_props = dict(
        solute_list=["Cl_-", "SO4_2-", "Ca_2+", "Mg_2+", "DOC"],
        material_flow_basis=MaterialFlowBasis.mass,
        mw_data={
            "H2O": 0.018,
            "Cl_-": 0.035,
            "SO4_2-": 0.096,
            "Ca_2+": 0.04,
            "Mg_2+": 0.024,
            "DOC": 0.235,
        },
        charge={"Cl_-": -1, "SO4_2-": -2, "Ca_2+": 2, "Mg_2+": 2, "DOC": 0},
    )

    conc_in = {
        "Cl_-": 0.02,
        "SO4_2-": 0.013,
        "Ca_2+": 0.04,
        "Mg_2+": 0.0243,
        "DOC": 0.3,
    }

    var_args = dict()
    for k, v in conc_in.items():
        var_args[("conc_mass_phase_comp", ("Liq", k))] = v * pyunits.g / pyunits.liter
    var_args[("flow_vol_phase", ("Liq"))] = 2.5 * pyunits.Mgallons / pyunits.day
    var_args[("temperature", (None))] = 298
    var_args[("pressure", (None))] = 101325

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**ion_props)
    ix_config = {
        "property_package": m.fs.properties,
    }
    m.fs.unit = IonExchangeDemin(**ix_config)
    pf = m.fs.unit.process_flow

    pf.properties_in.calculate_state(var_args=var_args, hold_state=True)

    loading_rate = 6.5 * pyunits.gallon / pyunits.minute / pyunits.ft**2

    m.fs.unit.loading_rate.fix(loading_rate)
    m.fs.unit.ebct.fix(240)
    m.fs.unit.resin_capacity_ax.fix(0.6)
    m.fs.unit.resin_capacity_cx.fix(0.8)
    m.fs.unit.resin_diam.fix()
    m.fs.unit.resin_density.fix(650)
    m.fs.unit.bed_porosity.fix()
    m.fs.unit.number_columns.fix(5)
    m.fs.unit.regen_dose.set_value(450)

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "Cl_-")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "SO4_2-")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e-2, index=("Liq", "H2O")
    )

    iscale.set_scaling_factor(
        pf.properties_out[0.0].flow_mass_phase_comp["Liq", "Cl_-"],
        1e6,
    )
    iscale.set_scaling_factor(
        pf.properties_out[0.0].flow_mass_phase_comp["Liq", "SO4_2-"],
        1e6,
    )
    iscale.set_scaling_factor(
        pf.properties_out[0.0].flow_mass_phase_comp["Liq", "Ca_2+"],
        1e6,
    )
    iscale.set_scaling_factor(
        pf.properties_out[0.0].flow_mass_phase_comp["Liq", "Mg_2+"],
        1e6,
    )

    iscale.calculate_scaling_factors(m)

    return m


class TestIXDeminCation(UnitTestHarness):
    def configure(self):
        m = build_cation_exchange()

        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        self.unit_solutions[m.fs.unit.bed_volume] = 3.28594
        self.unit_solutions[m.fs.unit.bed_volume_total] = 6.5718
        self.unit_solutions[m.fs.unit.bed_depth] = 0.9188
        self.unit_solutions[m.fs.unit.column_height] = 3.0172
        self.unit_solutions[m.fs.unit.bed_diameter] = 2.1339
        self.unit_solutions[m.fs.unit.bed_depth_to_diam_ratio] = 0.43058
        self.unit_solutions[m.fs.unit.number_columns] = 2
        self.unit_solutions[m.fs.unit.number_columns_redundant] = 1
        self.unit_solutions[m.fs.unit.breakthrough_time] = 46549
        self.unit_solutions[m.fs.unit.bv] = 310.33
        self.unit_solutions[m.fs.unit.loading_rate] = 0.006125456
        self.unit_solutions[m.fs.unit.N_Re] = 4.287819
        self.unit_solutions[m.fs.unit.N_Pe_particle] = 0.10056423
        self.unit_solutions[m.fs.unit.N_Pe_bed] = 132
        self.unit_solutions[m.fs.unit.resin_capacity_op] = 1.2366
        self.unit_solutions[
            m.fs.unit.process_flow.properties_in[0.0].total_hardness
        ] = 201.42
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "H2O"]
        ] = 0
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Ca_2+"]
        ] = -0.001735
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Inert"]
        ] = 0
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Mg_2+"]
        ] = -0.001054

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "Ca_2+"
                ],
                "out": m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "Ca_2+"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp[
                    "Liq", "Ca_2+"
                ],
            },
            "Check 2": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "Mg_2+"
                ],
                "out": m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "Mg_2+"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp[
                    "Liq", "Mg_2+"
                ],
            },
            "Check 3": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ],
                "out": m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ],
            },
        }

        return m

    @pytest.mark.component
    def test_costing(self):
        m = build_cation_exchange()

        m.fs.costing = WaterTAPCosting()
        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.unit.regen_dose.set_value(240)
        m.fs.costing.nacl.cost.set_value(0.230294)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(
            m.fs.unit.process_flow.properties_out[0].flow_vol_phase["Liq"]
        )
        m.fs.costing.add_specific_energy_consumption(
            m.fs.unit.process_flow.properties_out[0].flow_vol_phase["Liq"], name="SEC"
        )

        sys_cost_results = {
            "aggregate_capital_cost": 588302.62,
            "aggregate_fixed_operating_cost": 2693.98,
            "aggregate_flow_electricity": 1.5706,
            "aggregate_flow_NaCl": 1435400.39,
            "aggregate_flow_costs": {"electricity": 963.76, "NaCl": 334389.81},
            "total_capital_cost": 588302.62,
            "total_operating_cost": 322161.28,
            "LCOW": 0.306194,
            "SEC": 0.009958,
        }

        check_dof(m, fail_flag=True)
        initialization_tester(m)
        results = solver.solve(m)
        assert_optimal_termination(results)

        for v, r in sys_cost_results.items():
            mv = m.fs.costing.find_component(v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

        ix_cost_results = {
            "capital_cost": 588302.62,
            "fixed_operating_cost": 2693.98,
            "capital_cost_vessel": 62472.46,
            "capital_cost_resin": 17959.92,
            "capital_cost_regen_tank": 5826.13,
            "capital_cost_backwash_tank": 47028.0,
            "operating_cost_hazardous": 0.0,
            "flow_mass_regen_soln": 1435400.39,
            "total_pumping_power": 1.5706,
            "regeneration_tank_vol": 560.03,
            "backwash_tank_vol": 21990.39,
            "direct_capital_cost": 294151.31,
        }

        for v, r in ix_cost_results.items():
            mv = m.fs.unit.costing.find_component(v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)


class TestIXDeminAnion(UnitTestHarness):
    def configure(self):
        m = build_anion_exchange()

        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        self.unit_solutions[m.fs.unit.bed_volume] = 13.14379
        self.unit_solutions[m.fs.unit.bed_volume_total] = 131.43
        self.unit_solutions[m.fs.unit.bed_depth] = 1.324
        self.unit_solutions[m.fs.unit.column_height] = 3.9072
        self.unit_solutions[m.fs.unit.bed_diameter] = 3.5549
        self.unit_solutions[m.fs.unit.bed_depth_to_diam_ratio] = 0.3725
        self.unit_solutions[m.fs.unit.number_columns] = 10
        self.unit_solutions[m.fs.unit.number_columns_redundant] = 2.5
        self.unit_solutions[m.fs.unit.breakthrough_time] = 431737
        self.unit_solutions[m.fs.unit.bv] = 1439
        self.unit_solutions[m.fs.unit.loading_rate] = 0.004414
        self.unit_solutions[m.fs.unit.resin_capacity_op] = 1.2
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "H2O"]
        ] = 0
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Cl_-"]
        ] = -0.00867491
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "SO4_2-"]
        ] = -0.005639

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "Cl_-"
                ],
                "out": m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "Cl_-"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp[
                    "Liq", "Cl_-"
                ],
            },
            "Check 2": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "SO4_2-"
                ],
                "out": m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "SO4_2-"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp[
                    "Liq", "SO4_2-"
                ],
            },
            "Check 3": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ],
                "out": m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ],
            },
        }

        return m

    @pytest.mark.component
    def test_costing(self):
        m = build_anion_exchange()

        m.fs.costing = WaterTAPCosting()
        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(
            m.fs.unit.process_flow.properties_out[0].flow_vol_phase["Liq"]
        )
        m.fs.costing.add_specific_energy_consumption(
            m.fs.unit.process_flow.properties_out[0].flow_vol_phase["Liq"], name="SEC"
        )

        sys_cost_results = {
            "aggregate_flow_electricity": 17.69,
            "aggregate_flow_NaCl": 3538788.93,
            "aggregate_flow_costs": {"electricity": 10860.27, "NaCl": 322176.99},
            "total_capital_cost": 5664077.65,
            "total_operating_cost": 529815.76,
            "aggregate_direct_capital_cost": 2832038.82,
            "maintenance_labor_chemical_operating_cost": 169922.32,
            "total_fixed_operating_cost": 230082.21,
            "total_variable_operating_cost": 299733.54,
            "total_annualized_cost": 1096223.52,
            "LCOW": 0.088098,
            "SEC": 0.011221,
        }

        check_dof(m, fail_flag=True)
        initialization_tester(m)
        results = solver.solve(m)
        assert_optimal_termination(results)

        for v, r in sys_cost_results.items():
            mv = m.fs.costing.find_component(v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

        ix_cost_results = {
            "capital_cost": 5664077.65,
            "fixed_operating_cost": 60159.88,
            "capital_cost_vessel": 112455.69,
            "capital_cost_resin": 96255.81,
            "capital_cost_regen_tank": 18843.29,
            "capital_cost_backwash_tank": 204301.7,
            "flow_mass_regen_soln": 3538788.93,
            "total_pumping_power": 17.69,
            "regeneration_tank_vol": 2800.17,
            "backwash_tank_vol": 411462.44,
            "direct_capital_cost": 2832038.82,
        }

        for v, r in ix_cost_results.items():
            mv = m.fs.unit.costing.find_component(v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)


class TestIXDeminMixed(UnitTestHarness):
    def configure(self):
        m = build_mixed_bed()

        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        self.unit_solutions[m.fs.unit.bed_volume] = 5.25751
        self.unit_solutions[m.fs.unit.bed_volume_total] = 26.2875
        self.unit_solutions[m.fs.unit.bed_depth] = 1.0594
        self.unit_solutions[m.fs.unit.bed_porosity] = 0.4
        self.unit_solutions[m.fs.unit.column_height] = 3.3258
        self.unit_solutions[m.fs.unit.bed_diameter] = 2.5137
        self.unit_solutions[m.fs.unit.bed_depth_to_diam_ratio] = 0.4214
        self.unit_solutions[m.fs.unit.number_columns] = 5
        self.unit_solutions[m.fs.unit.number_columns_redundant] = 1.25
        self.unit_solutions[m.fs.unit.breakthrough_time] = 38121
        self.unit_solutions[m.fs.unit.bv] = 158.84
        self.unit_solutions[m.fs.unit.resin_capacity_ax] = 0.6
        self.unit_solutions[m.fs.unit.resin_capacity_cx] = 0.8
        self.unit_solutions[m.fs.unit.resin_capacity_op] = 0.7654
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "H2O"]
        ] = 0
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Cl_-"]
        ] = -0.0021687
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "SO4_2-"]
        ] = -0.001409
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Ca_2+"]
        ] = -0.004337
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Mg_2+"]
        ] = -0.002635

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "Cl_-"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                    "Liq", "Cl_-"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp[
                    "Liq", "Cl_-"
                ],
            },
            "Check 2": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "SO4_2-"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                    "Liq", "SO4_2-"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp[
                    "Liq", "SO4_2-"
                ],
            },
            "Check 3": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "Ca_2+"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                    "Liq", "Ca_2+"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp[
                    "Liq", "Ca_2+"
                ],
            },
            "Check 4": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "Mg_2+"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                    "Liq", "Mg_2+"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp[
                    "Liq", "Mg_2+"
                ],
            },
            "Check 5": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ],
            },
        }

        return m

    @pytest.mark.component
    def test_costing(self):
        m = build_mixed_bed()

        m.fs.costing = WaterTAPCosting()
        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(
            m.fs.unit.process_flow.properties_out[0].flow_vol_phase["Liq"]
        )
        m.fs.costing.add_specific_energy_consumption(
            m.fs.unit.process_flow.properties_out[0].flow_vol_phase["Liq"], name="SEC"
        )

        sys_cost_results = {
            "aggregate_flow_electricity": 3.1654,
            "aggregate_flow_NaCl": 9838326.68,
            "aggregate_flow_costs": {"electricity": 1942.39, "NaCl": 895696.97},
            "total_capital_cost": 1541452.33,
            "total_operating_cost": 863627.1,
            "total_fixed_operating_cost": 55751.67,
            "total_variable_operating_cost": 807875.42,
            "total_annualized_cost": 1017772.33,
            "LCOW": 0.327194,
            "SEC": 0.008028,
        }

        check_dof(m, fail_flag=True)
        initialization_tester(m)
        results = solver.solve(m)
        assert_optimal_termination(results)

        for v, r in sys_cost_results.items():
            mv = m.fs.costing.find_component(v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

        ix_cost_results = {
            "capital_cost": 1541452.33,
            "fixed_operating_cost": 9508.1,
            "capital_cost_vessel": 75946.17,
            "capital_cost_resin": 30425.92,
            "capital_cost_regen_tank": 12982.45,
            "capital_cost_backwash_tank": 92918.03,
            "operating_cost_hazardous": 0.0,
            "flow_mass_regen_soln": 9838326.68,
            "total_pumping_power": 3.1654,
            "regeneration_tank_vol": 1680.1,
            "backwash_tank_vol": 85504.49,
            "direct_capital_cost": 770726.16,
        }

        for v, r in ix_cost_results.items():
            mv = m.fs.unit.costing.find_component(v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)
