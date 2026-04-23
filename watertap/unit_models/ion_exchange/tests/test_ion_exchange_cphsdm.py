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
)

from idaes.core import (
    FlowsheetBlock,
    UnitModelCostingBlock,
)
import idaes.core.util.scaling as iscale
from idaes.core.util.testing import initialization_tester

from watertap.costing import WaterTAPCosting
from watertap.core.util.initialization import check_dof
from watertap.core.solvers import get_solver
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models import IonExchangeCPHSDM
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

__author__ = "Kurban Sitterley"

solver = get_solver()
zero = 1e-8
relative_tolerance = 1e-3


def build_ix_hand():
    """
    Duplicate test from GAC model build_hand
    to ensure results are identical
    """

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(
        solute_list=["DCE"],
        mw_data={"H2O": 0.018, "DCE": 0.09896},
        charge={"DCE": 1},
        ignore_neutral_charge=True,
    )
    ix_config = {"property_package": m.fs.properties, "target_component": "DCE"}
    m.fs.unit = IonExchangeCPHSDM(**ix_config)

    pf = m.fs.unit.process_flow
    pf.properties_in[0].pressure.fix(101325)
    pf.properties_in[0].temperature.fix(273.15 + 25)
    pf.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix(55555.55426666667)
    pf.properties_in[0].flow_mol_phase_comp["Liq", "DCE"].fix(0.0002344381568310428)

    # adjust variable bounds to use GAC inputs
    m.fs.unit.bed_depth.setub(6)
    m.fs.unit.bed_diameter.setub(10)
    m.fs.unit.column_height.setub(100)
    m.fs.unit.resin_density.setlb(0)
    m.fs.unit.loading_rate.setub(None)

    m.fs.unit.freundlich_k.fix(37.9e-6 * (1e6**0.8316))
    m.fs.unit.freundlich_ninv.fix(0.8316)
    m.fs.unit.resin_density_app.fix(722)
    m.fs.unit.resin_diam.fix(0.00106)
    m.fs.unit.ebct.fix(300)
    m.fs.unit.bed_porosity.fix(0.449)
    m.fs.unit.bed_depth.fix(6)
    m.fs.unit.c_norm.fix(0.5)
    m.fs.unit.film_mass_transfer_coeff.fix(3.29e-5)
    m.fs.unit.surf_diff_coeff.fix(1.77e-13)
    m.fs.unit.a0.fix(3.68421)
    m.fs.unit.a1.fix(13.1579)
    m.fs.unit.b0.fix(0.784576)
    m.fs.unit.b1.fix(0.239663)
    m.fs.unit.b2.fix(0.484422)
    m.fs.unit.b3.fix(0.003206)
    m.fs.unit.b4.fix(0.134987)
    m.fs.unit.number_columns.fix(1)

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e4, index=("Liq", "DCE")
    )
    iscale.calculate_scaling_factors(m)

    return m


class TestIXCPHSDMHand(UnitTestHarness):
    def configure(self):
        m = build_ix_hand()

        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance
        self.unit_solutions[m.fs.unit.resin_density] = 397.82
        self.unit_solutions[m.fs.unit.bed_porosity] = 0.449
        self.unit_solutions[m.fs.unit.column_height] = 14.1722
        self.unit_solutions[m.fs.unit.bed_diameter] = 7.9788
        self.unit_solutions[m.fs.unit.bed_depth_to_diam_ratio] = 0.75198
        self.unit_solutions[m.fs.unit.number_columns] = 1
        self.unit_solutions[m.fs.unit.number_columns_redundant] = 1
        self.unit_solutions[m.fs.unit.breakthrough_time] = 2554274
        self.unit_solutions[m.fs.unit.bv] = 8514
        self.unit_solutions[m.fs.unit.loading_rate] = 0.02
        self.unit_solutions[m.fs.unit.c_eq["DCE"]] = 0.0005178
        self.unit_solutions[m.fs.unit.solute_dist_param] = 19775
        self.unit_solutions[m.fs.unit.N_Bi] = 6.113
        self.unit_solutions[m.fs.unit.min_N_St] = 35.68
        self.unit_solutions[m.fs.unit.min_ebct] = 1043.1
        self.unit_solutions[m.fs.unit.throughput] = 0.9882
        self.unit_solutions[m.fs.unit.min_t_contact] = 468.4
        self.unit_solutions[m.fs.unit.min_breakthrough_time] = 9153491
        self.unit_solutions[m.fs.unit.tb_traps["DCE", 0]] = 0.0
        self.unit_solutions[m.fs.unit.tb_traps["DCE", 1]] = 969736.0476
        self.unit_solutions[m.fs.unit.tb_traps["DCE", 2]] = 1621753.556
        self.unit_solutions[m.fs.unit.tb_traps["DCE", 3]] = 1980033.2921
        self.unit_solutions[m.fs.unit.tb_traps["DCE", 4]] = 2276213.8754
        self.unit_solutions[m.fs.unit.tb_traps["DCE", 5]] = 2554274.2291
        self.unit_solutions[m.fs.unit.traps["DCE", 1]] = 0.001898261424
        self.unit_solutions[m.fs.unit.traps["DCE", 2]] = 0.018187650702
        self.unit_solutions[m.fs.unit.traps["DCE", 3]] = 0.027176682157
        self.unit_solutions[m.fs.unit.traps["DCE", 4]] = 0.036670733463
        self.unit_solutions[m.fs.unit.traps["DCE", 5]] = 0.047762679049
        self.unit_solutions[m.fs.unit.c_norm_avg["DCE"]] = 0.131696
        self.unit_solutions[m.fs.unit.min_tb_traps[1]] = 7568953
        self.unit_solutions[m.fs.unit.min_tb_traps[2]] = 8220971
        self.unit_solutions[m.fs.unit.min_tb_traps[3]] = 8579250
        self.unit_solutions[m.fs.unit.min_tb_traps[4]] = 8875431
        self.unit_solutions[m.fs.unit.min_tb_traps[5]] = 9153491
        self.unit_solutions[m.fs.unit.throughput_traps[1]] = 0.81710
        self.unit_solutions[m.fs.unit.throughput_traps[2]] = 0.88749
        self.unit_solutions[m.fs.unit.throughput_traps[3]] = 0.92617
        self.unit_solutions[m.fs.unit.throughput_traps[4]] = 0.95814
        self.unit_solutions[m.fs.unit.throughput_traps[5]] = 0.98816
        self.unit_solutions[m.fs.unit.mass_adsorbed] = 51.454
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "H2O"]
        ] = 85.945
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "DCE"]
        ] = 0.00020356
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 1.5470
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp["Liq", "DCE"]
        ] = 2.0144e-05
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_vol_phase["Liq"]
        ] = 0.001547

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "DCE"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "DCE"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "DCE"],
            },
            "Check 2": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "DCE"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "DCE"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "DCE"],
            },
        }

        return m

    @pytest.mark.component
    def test_costing(self):
        m = build_ix_hand()

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
            "aggregate_capital_cost": 9216178.73,
            "aggregate_fixed_operating_cost": 163970.28,
            "aggregate_flow_electricity": 896.21,
            "aggregate_flow_NaCl": 2217092.5,
            "aggregate_flow_costs": {"electricity": 549935.15, "NaCl": 201847.64},
            "total_capital_cost": 9216178.73,
            "total_operating_cost": 1117060.16,
            "LCOW": 0.071779,
            "SEC": 0.248948,
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
            "capital_cost": 9216178.73,
            "fixed_operating_cost": 163970.28,
            "capital_cost_vessel": 427335.95,
            "capital_cost_resin": 1639702.86,
            "capital_cost_regen_tank": 184447.54,
            "capital_cost_backwash_tank": 289564.18,
            "flow_mass_regen_soln": 2217092.5,
            "total_pumping_power": 896.21,
            "regeneration_tank_vol": 63912.59,
            "backwash_tank_vol": 824877.23,
            "direct_capital_cost": 4608089.36,
        }

        for v, r in ix_cost_results.items():
            mv = m.fs.unit.costing.find_component(v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)


def build_ix_hand_surrogate():
    """
    Duplicate test from GAC model build_hand
    to ensure results are identical
    """

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(
        solute_list=["DCE"],
        mw_data={"H2O": 0.018, "DCE": 0.09896},
        charge={"DCE": 1},
        ignore_neutral_charge=True,
    )
    ix_config = {
        "property_package": m.fs.properties,
        "target_component": "DCE",
        "cphsdm_calculation_method": "surrogate",
    }
    m.fs.unit = IonExchangeCPHSDM(**ix_config)

    pf = m.fs.unit.process_flow
    pf.properties_in[0].pressure.fix(101325)
    pf.properties_in[0].temperature.fix(273.15 + 25)
    pf.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix(55555.55426666667)
    pf.properties_in[0].flow_mol_phase_comp["Liq", "DCE"].fix(0.0002344381568310428)

    # adjust variable bounds to use GAC inputs
    m.fs.unit.bed_depth.setub(6)
    m.fs.unit.bed_diameter.setub(10)
    m.fs.unit.column_height.setub(100)
    m.fs.unit.resin_density.setlb(0)
    m.fs.unit.loading_rate.setub(None)

    m.fs.unit.freundlich_k.fix(37.9e-6 * (1e6**0.8316))
    m.fs.unit.freundlich_ninv.fix(0.8316)
    m.fs.unit.resin_density_app.fix(722)
    m.fs.unit.resin_diam.fix(0.00106)
    m.fs.unit.ebct.fix(300)
    m.fs.unit.bed_porosity.fix(0.449)
    m.fs.unit.bed_depth.fix(6)
    m.fs.unit.c_norm.fix(0.5)
    m.fs.unit.film_mass_transfer_coeff.fix(3.29e-5)
    m.fs.unit.surf_diff_coeff.fix(1.77e-13)
    m.fs.unit.number_columns.fix(1)

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e4, index=("Liq", "DCE")
    )
    iscale.calculate_scaling_factors(m)

    return m


class TestIXCPHSDMSurrogateHand(UnitTestHarness):
    def configure(self):
        m = build_ix_hand_surrogate()

        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance
        self.unit_solutions[m.fs.unit.resin_density] = 397.82
        self.unit_solutions[m.fs.unit.bed_porosity] = 0.449
        self.unit_solutions[m.fs.unit.column_height] = 14.1722
        self.unit_solutions[m.fs.unit.bed_diameter] = 7.9788
        self.unit_solutions[m.fs.unit.bed_depth_to_diam_ratio] = 0.75198
        self.unit_solutions[m.fs.unit.number_columns] = 1
        self.unit_solutions[m.fs.unit.number_columns_redundant] = 1
        self.unit_solutions[m.fs.unit.breakthrough_time] = 2530674
        self.unit_solutions[m.fs.unit.bv] = 8435
        self.unit_solutions[m.fs.unit.loading_rate] = 0.02
        self.unit_solutions[m.fs.unit.c_eq["DCE"]] = 0.0005178
        self.unit_solutions[m.fs.unit.solute_dist_param] = 19775
        self.unit_solutions[m.fs.unit.N_Bi] = 6.113
        self.unit_solutions[m.fs.unit.min_N_St] = 47.76
        self.unit_solutions[m.fs.unit.min_ebct] = 1396
        self.unit_solutions[m.fs.unit.throughput] = 0.989252
        self.unit_solutions[m.fs.unit.min_t_contact] = 626.96
        self.unit_solutions[m.fs.unit.min_breakthrough_time] = 12266003
        self.unit_solutions[m.fs.unit.tb_traps["DCE", 1]] = 722271
        self.unit_solutions[m.fs.unit.tb_traps["DCE", 2]] = 1452770
        self.unit_solutions[m.fs.unit.tb_traps["DCE", 3]] = 1856540
        self.unit_solutions[m.fs.unit.tb_traps["DCE", 4]] = 2202334
        self.unit_solutions[m.fs.unit.tb_traps["DCE", 5]] = 2530674
        self.unit_solutions[m.fs.unit.traps["DCE", 1]] = 0.0014270
        self.unit_solutions[m.fs.unit.traps["DCE", 2]] = 0.0205669
        self.unit_solutions[m.fs.unit.traps["DCE", 3]] = 0.0309129
        self.unit_solutions[m.fs.unit.traps["DCE", 4]] = 0.0432127
        self.unit_solutions[m.fs.unit.traps["DCE", 5]] = 0.0569253
        self.unit_solutions[m.fs.unit.c_norm_avg["DCE"]] = 0.1530447
        self.unit_solutions[m.fs.unit.min_tb_traps[1]] = 10457601
        self.unit_solutions[m.fs.unit.min_tb_traps[2]] = 11188099
        self.unit_solutions[m.fs.unit.min_tb_traps[3]] = 11591869
        self.unit_solutions[m.fs.unit.min_tb_traps[4]] = 11937663
        self.unit_solutions[m.fs.unit.min_tb_traps[5]] = 12266003
        self.unit_solutions[m.fs.unit.throughput_traps[1]] = 0.84340
        self.unit_solutions[m.fs.unit.throughput_traps[2]] = 0.90232
        self.unit_solutions[m.fs.unit.throughput_traps[3]] = 0.93488
        self.unit_solutions[m.fs.unit.throughput_traps[4]] = 0.96277
        self.unit_solutions[m.fs.unit.throughput_traps[5]] = 0.98925
        self.unit_solutions[m.fs.unit.mass_adsorbed] = 49.726

        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "H2O"]
        ] = 86.74
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "DCE"]
        ] = 0.0001985
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 1.5613
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp["Liq", "DCE"]
        ] = 1.964e-05
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_vol_phase["Liq"]
        ] = 0.0015614

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "DCE"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "DCE"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "DCE"],
            },
            "Check 2": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "DCE"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "DCE"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "DCE"],
            },
        }

        return m

    @pytest.mark.component
    def test_costing(self):
        m = build_ix_hand_surrogate()

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
            "aggregate_flow_electricity": 896.2,
            "aggregate_flow_NaCl": 2237704.39,
            "aggregate_flow_costs": {"electricity": 549927.47, "NaCl": 203724.18},
            "total_capital_cost": 9216178.73,
            "total_operating_cost": 1118742.14,
            "total_fixed_operating_cost": 440455.64,
            "total_variable_operating_cost": 678286.49,
            "total_annualized_cost": 2040360.01,
            "LCOW": 0.071839,
            "SEC": 0.248945,
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
            "capital_cost": 9216178.73,
            "fixed_operating_cost": 163970.28,
            "capital_cost_vessel": 427335.95,
            "capital_cost_resin": 1639702.86,
            "capital_cost_regen_tank": 184447.54,
            "capital_cost_backwash_tank": 289564.18,
            "flow_mass_regen_soln": 2237704.39,
            "total_pumping_power": 896.2,
            "regeneration_tank_vol": 63912.59,
            "backwash_tank_vol": 824877.23,
            "direct_capital_cost": 4608089.36,
        }

        for v, r in ix_cost_results.items():
            mv = m.fs.unit.costing.find_component(v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)


def build_ix_crittenden():
    # trial problem from Crittenden, 2012 for removal of TCE
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = MCASParameterBlock(
        solute_list=["TCE"],
        mw_data={"H2O": 0.018, "TCE": 0.1314},
        diffus_calculation="HaydukLaudie",
        molar_volume_data={("Liq", "TCE"): 9.81e-5},
        charge={"TCE": -1},
        ignore_neutral_charge=True,
    )
    m.fs.properties.visc_d_phase["Liq"] = 1.3097e-3
    m.fs.properties.dens_mass_const = 999.7
    ix_config = {"property_package": m.fs.properties, "target_component": "TCE"}
    m.fs.unit = IonExchangeCPHSDM(**ix_config)

    unit_feed = m.fs.unit.process_flow.properties_in[0]
    unit_feed.pressure.fix(101325)
    unit_feed.temperature.fix(273.15 + 25)
    unit_feed.flow_mol_phase_comp["Liq", "H2O"].fix(823.8)
    unit_feed.flow_mol_phase_comp["Liq", "TCE"].fix(5.6444e-05)

    m.fs.unit.resin_density.setlb(0)

    m.fs.unit.freundlich_k.fix(1062e-6 * (1e6**0.48))
    m.fs.unit.freundlich_ninv.fix(0.48)
    m.fs.unit.resin_density_app.fix(803.4)
    m.fs.unit.resin_diam.fix(0.001026)
    m.fs.unit.ebct.fix(600)
    m.fs.unit.bed_porosity.fix(0.44)
    m.fs.unit.loading_rate.fix((5 / 3600))
    m.fs.unit.c_norm["TCE"].fix(0.80)
    m.fs.unit.surf_diff_coeff.fix(1.24e-14)
    m.fs.unit.film_mass_transfer_coeff.fix(3.73e-05)
    m.fs.unit.a0.fix(0.8)
    m.fs.unit.a1.fix(0)
    m.fs.unit.b0.fix(0.023)
    m.fs.unit.b1.fix(0.793673)
    m.fs.unit.b2.fix(0.039324)
    m.fs.unit.b3.fix(0.009326)
    m.fs.unit.b4.fix(0.08275)
    m.fs.unit.number_columns.fix(1)

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-2, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e5, index=("Liq", "TCE")
    )
    iscale.calculate_scaling_factors(m)

    return m


class TestIXCPHSDMCrittenden(UnitTestHarness):
    def configure(self):
        m = build_ix_crittenden()

        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance
        self.unit_solutions[m.fs.unit.bed_volume] = 8.899
        self.unit_solutions[m.fs.unit.bed_depth] = 0.833
        self.unit_solutions[m.fs.unit.column_height] = 2.829
        self.unit_solutions[m.fs.unit.bed_diameter] = 3.688
        self.unit_solutions[m.fs.unit.bed_depth_to_diam_ratio] = 0.226
        self.unit_solutions[m.fs.unit.breakthrough_time] = 13687913
        self.unit_solutions[m.fs.unit.bv] = 22813
        self.unit_solutions[m.fs.unit.c_eq["TCE"]] = 0.020971
        self.unit_solutions[m.fs.unit.mass_adsorbed] = 78.24
        self.unit_solutions[m.fs.unit.solute_dist_param] = 42886
        self.unit_solutions[m.fs.unit.N_Bi] = 45.79
        self.unit_solutions[m.fs.unit.resin_density_app] = 803.4
        self.unit_solutions[m.fs.unit.min_N_St] = 36.63
        self.unit_solutions[m.fs.unit.min_ebct] = 899.7
        self.unit_solutions[m.fs.unit.throughput] = 1.139
        self.unit_solutions[m.fs.unit.min_t_contact] = 395.9
        self.unit_solutions[m.fs.unit.min_breakthrough_time] = 19344728
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "TCE"]
        ] = -4.350108e-05
        self.unit_solutions[
            m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp["Liq", "TCE"]
        ] = 1.2942e-05
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "H2O"]
        ] = 0.567751
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "TCE"]
        ] = 4.3501e-05
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.010219
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp["Liq", "TCE"]
        ] = 5.71e-06
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_vol_phase["Liq"]
        ] = 1.022e-05
        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "TCE"],
            },
            "Check 2": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "TCE"],
            },
        }

        return m

    @pytest.mark.component
    def test_costing(self):
        m = build_ix_crittenden()

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
            "aggregate_flow_electricity": 0.159214,
            "aggregate_flow_NaCl": 12301.3,
            "aggregate_flow_costs": {"electricity": 97.69, "NaCl": 1119.93},
            "total_capital_cost": 800840.15,
            "total_operating_cost": 31638.59,
            "LCOW": 0.265197,
            "SEC": 0.002981,
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
            "capital_cost": 800840.15,
            "fixed_operating_cost": 6517.52,
            "capital_cost_vessel": 100274.86,
            "capital_cost_resin": 65175.2,
            "capital_cost_regen_tank": 14179.12,
            "capital_cost_backwash_tank": 55340.8,
            "flow_mass_regen_soln": 12301.3,
            "total_pumping_power": 0.159214,
            "regeneration_tank_vol": 1896.01,
            "backwash_tank_vol": 30422.66,
            "direct_capital_cost": 400420.07,
        }

        for v, r in ix_cost_results.items():
            mv = m.fs.unit.costing.find_component(v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)


def build_ix_crittenden_surrogate():
    # trial problem from Crittenden, 2012 for removal of TCE
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = MCASParameterBlock(
        solute_list=["TCE"],
        mw_data={"H2O": 0.018, "TCE": 0.1314},
        diffus_calculation="HaydukLaudie",
        molar_volume_data={("Liq", "TCE"): 9.81e-5},
        charge={"TCE": -1},
        ignore_neutral_charge=True,
    )
    m.fs.properties.visc_d_phase["Liq"] = 1.3097e-3
    m.fs.properties.dens_mass_const = 999.7
    ix_config = {
        "property_package": m.fs.properties,
        "target_component": "TCE",
        "cphsdm_calculation_method": "surrogate",
    }
    m.fs.unit = IonExchangeCPHSDM(**ix_config)

    unit_feed = m.fs.unit.process_flow.properties_in[0]
    unit_feed.pressure.fix(101325)
    unit_feed.temperature.fix(273.15 + 25)
    unit_feed.flow_mol_phase_comp["Liq", "H2O"].fix(823.8)
    unit_feed.flow_mol_phase_comp["Liq", "TCE"].fix(5.6444e-05)

    m.fs.unit.resin_density.setlb(0)

    m.fs.unit.freundlich_k.fix(1062e-6 * (1e6**0.48))
    m.fs.unit.freundlich_ninv.fix(0.48)
    m.fs.unit.resin_density_app.fix(803.4)
    m.fs.unit.resin_diam.fix(0.001026)
    m.fs.unit.ebct.fix(600)
    m.fs.unit.bed_porosity.fix(0.44)
    m.fs.unit.loading_rate.fix((5 / 3600))
    m.fs.unit.c_norm["TCE"].fix(0.80)
    m.fs.unit.surf_diff_coeff.fix(1.24e-14)
    m.fs.unit.film_mass_transfer_coeff.fix(3.73e-05)
    m.fs.unit.number_columns.fix(1)

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-2, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e5, index=("Liq", "TCE")
    )
    iscale.calculate_scaling_factors(m)

    return m


class TestIXCPHSDMSurrogateCrittenden(UnitTestHarness):
    def configure(self):
        m = build_ix_crittenden_surrogate()

        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance
        self.unit_solutions[m.fs.unit.bed_volume] = 8.899
        self.unit_solutions[m.fs.unit.bed_depth] = 0.833
        self.unit_solutions[m.fs.unit.column_height] = 2.829
        self.unit_solutions[m.fs.unit.bed_diameter] = 3.688
        self.unit_solutions[m.fs.unit.bed_depth_to_diam_ratio] = 0.226
        self.unit_solutions[m.fs.unit.breakthrough_time] = 12573768
        self.unit_solutions[m.fs.unit.bv] = 20956
        self.unit_solutions[m.fs.unit.c_eq["TCE"]] = 0.020971
        self.unit_solutions[m.fs.unit.mass_adsorbed] = 74.4781
        self.unit_solutions[m.fs.unit.min_N_St] = 32.0981
        self.unit_solutions[m.fs.unit.min_breakthrough_time] = 16127386
        self.unit_solutions[m.fs.unit.min_ebct] = 788.3
        self.unit_solutions[m.fs.unit.min_t_contact] = 346.85
        self.unit_solutions[m.fs.unit.min_tb_traps[1]] = 11384169
        self.unit_solutions[m.fs.unit.min_tb_traps[2]] = 12149408
        self.unit_solutions[m.fs.unit.min_tb_traps[3]] = 12667838
        self.unit_solutions[m.fs.unit.min_tb_traps[4]] = 13555300
        self.unit_solutions[m.fs.unit.min_tb_traps[5]] = 16127386
        self.unit_solutions[m.fs.unit.traps["TCE", 1]] = 0.003113
        self.unit_solutions[m.fs.unit.traps["TCE", 2]] = 0.006618
        self.unit_solutions[m.fs.unit.traps["TCE", 3]] = 0.012627
        self.unit_solutions[m.fs.unit.traps["TCE", 4]] = 0.035554
        self.unit_solutions[m.fs.unit.traps["TCE", 5]] = 0.143447
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "TCE"]
        ] = -4.350e-05
        self.unit_solutions[
            m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp["Liq", "TCE"]
        ] = 1.2942e-05
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "H2O"]
        ] = 0.618032
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.0111245
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp["Liq", "TCE"]
        ] = 5.92e-06
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_vol_phase["Liq"]
        ] = 1.113e-05
        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "TCE"],
            },
            "Check 2": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "TCE"],
            },
        }

        return m


def build_ix_cphsdm_inert():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(
        solute_list=["TCE", "BGSOL", "BGCAT", "BGAN"],
        mw_data={
            "H2O": 0.018,
            "TCE": 0.1314,
            "BGSOL": 0.1,
            "BGCAT": 0.1,
            "BGAN": 0.1,
        },
        charge={"BGCAT": 1, "BGAN": -2, "TCE": -1},
        diffus_calculation="HaydukLaudie",
        molar_volume_data={("Liq", "TCE"): 9.81e-5},
        diffusivity_data={
            ("Liq", "BGSOL"): 1e-5,
            ("Liq", "BGCAT"): 1e-5,
            ("Liq", "BGAN"): 1e-5,
        },
        ignore_neutral_charge=True,
    )
    m.fs.properties.visc_d_phase["Liq"] = 1.3097e-3
    m.fs.properties.dens_mass_const = 1000
    ix_config = {
        "property_package": m.fs.properties,
        "target_component": "TCE",
        "film_transfer_coefficient_type": "calculated",
        "surface_diffusion_coefficient_type": "calculated",
    }
    m.fs.unit = IonExchangeCPHSDM(**ix_config)

    unit_feed = m.fs.unit.process_flow.properties_in[0]
    unit_feed.pressure.fix(101325)
    unit_feed.temperature.fix(273.15 + 25)
    unit_feed.flow_mol_phase_comp["Liq", "H2O"].fix(824.0736620370348)
    unit_feed.flow_mol_phase_comp["Liq", "TCE"].fix(5.644342973110135e-05)
    unit_feed.flow_mol_phase_comp["Liq", "BGSOL"].fix(5e-05)
    unit_feed.flow_mol_phase_comp["Liq", "BGCAT"].fix(2e-05)
    unit_feed.flow_mol_phase_comp["Liq", "BGAN"].fix(1e-05)

    m.fs.unit.resin_density.setlb(0)

    m.fs.unit.freundlich_k.fix(1062e-6 * (1e6**0.48))
    m.fs.unit.freundlich_ninv.fix(0.48)
    m.fs.unit.resin_density_app.fix(803.4)
    m.fs.unit.resin_diam.fix(0.001026)
    m.fs.unit.ebct.fix(10 * 60)
    m.fs.unit.bed_porosity.fix(0.44)
    m.fs.unit.loading_rate.fix(5 / 3600)
    m.fs.unit.c_norm["TCE"].fix(0.80)
    m.fs.unit.resin_porosity.fix(0.641)
    m.fs.unit.tortuosity.fix(1)
    m.fs.unit.spdfr.fix(1)
    m.fs.unit.shape_correction_factor.fix(1.5)
    m.fs.unit.a0.fix(0.8)
    m.fs.unit.a1.fix(0)
    m.fs.unit.b0.fix(0.023)
    m.fs.unit.b1.fix(0.793673)
    m.fs.unit.b2.fix(0.039324)
    m.fs.unit.b3.fix(0.009326)
    m.fs.unit.b4.fix(0.08275)
    m.fs.unit.number_columns.fix(1)

    # scaling
    prop = m.fs.properties
    prop.set_default_scaling("flow_mol_phase_comp", 1e-2, index=("Liq", "H2O"))
    for j in prop.ion_set | prop.solute_set:
        prop.set_default_scaling("flow_mol_phase_comp", 1e5, index=("Liq", j))

    iscale.calculate_scaling_factors(m)

    return m


class TestIXCPHSDMInert(UnitTestHarness):
    def configure(self):
        m = build_ix_cphsdm_inert()

        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        self.unit_solutions[m.fs.unit.N_Re] = 2.473
        self.unit_solutions[m.fs.unit.N_Pe_particle] = 0.0772
        self.unit_solutions[m.fs.unit.N_Pe_bed] = 62.71
        self.unit_solutions[m.fs.unit.N_Sc["TCE"]] = 2001.2
        self.unit_solutions[m.fs.unit.film_mass_transfer_coeff] = 2.5988e-05
        self.unit_solutions[m.fs.unit.surf_diff_coeff] = 1.2449e-14

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ]
                + m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "BGSOL"
                ]
                + m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "BGCAT"
                ]
                + m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "BGAN"
                ],
                "out": m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ]
                + m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "BGSOL"
                ]
                + m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "BGCAT"
                ]
                + m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "BGAN"
                ],
            },
            "Check 2": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "TCE"],
            },
            "Check 3": {
                "in": m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp[
                    "Liq", "BGSOL"
                ],
                "out": 0,
            },
            "Check 4": {
                "in": m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp[
                    "Liq", "BGCAT"
                ],
                "out": 0,
            },
        }

        return m
