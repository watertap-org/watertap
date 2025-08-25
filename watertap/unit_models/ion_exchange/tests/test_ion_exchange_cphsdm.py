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
from watertap.unit_models import IonExchangeCPHSDM
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

__author__ = "Kurban Sitterley"

solver = get_solver()
zero = 1e-8
relative_tolerance = 1e-3


def build_cphsdm():
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


class TestIXCPHSDM(UnitTestHarness):
    def configure(self):
        m = build_cphsdm()

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
        self.unit_solutions[m.fs.unit.tb_traps[0]] = 0
        self.unit_solutions[m.fs.unit.tb_traps[1]] = 969736
        self.unit_solutions[m.fs.unit.tb_traps[2]] = 1621753
        self.unit_solutions[m.fs.unit.tb_traps[3]] = 1980033
        self.unit_solutions[m.fs.unit.tb_traps[4]] = 2276213
        self.unit_solutions[m.fs.unit.tb_traps[5]] = 2554274
        self.unit_solutions[m.fs.unit.traps[1]] = 0.00189826
        self.unit_solutions[m.fs.unit.traps[2]] = 0.01818765
        self.unit_solutions[m.fs.unit.traps[3]] = 0.02717668
        self.unit_solutions[m.fs.unit.traps[4]] = 0.03667073
        self.unit_solutions[m.fs.unit.traps[5]] = 0.04776267
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
        m = build_cphsdm()

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
