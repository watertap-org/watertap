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
from watertap.unit_models.ion_exchange import IonExchangeBaseData
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
        # },
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


class TestIXDeminCation(UnitTestHarness):
    def configure(self):
        m = build_cation_exchange()

        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        self.unit_solutions[m.fs.unit.bed_volume] = 3.28594
        self.unit_solutions[m.fs.unit.bed_volume_total] = 6.5718
        self.unit_solutions[m.fs.unit.bed_depth] = 0.9188
        self.unit_solutions[m.fs.unit.column_height] = 2.3451
        self.unit_solutions[m.fs.unit.bed_diameter] = 2.1339
        self.unit_solutions[m.fs.unit.col_height_to_diam_ratio] = 1.0990
        self.unit_solutions[m.fs.unit.number_columns] = 2
        self.unit_solutions[m.fs.unit.number_columns_redundant] = 1
        self.unit_solutions[m.fs.unit.breakthrough_time] = 46549
        self.unit_solutions[m.fs.unit.bv] = 310.33
        self.unit_solutions[m.fs.unit.ebct] = 150
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
        ix = m.fs.unit

        m.fs.costing = WaterTAPCosting()
        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.unit.costing.regen_dose.set_value(240)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.process_flow.properties_out[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            m.fs.unit.process_flow.properties_out[0].flow_vol_phase["Liq"], name="SEC"
        )

        sys_cost_results = {
            "aggregate_capital_cost": 483892.9,
            "aggregate_fixed_operating_cost": 1993.8,
            "aggregate_variable_operating_cost": 0.0,
            "aggregate_flow_electricity": 1.2152,
            "aggregate_flow_NaCl": 1116422.5,
            "aggregate_flow_costs": {"electricity": 745.72, "NaCl": 101640.8},
            "total_capital_cost": 483892.9,
            "total_operating_cost": 108658.5,
            "LCOW": 0.170524,
        }

        check_dof(m, fail_flag=True)
        initialization_tester(m)
        results = solver.solve(m)
        assert_optimal_termination(results)

        for v, r in sys_cost_results.items():
            mv = getattr(m.fs.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

        ix_cost_results = {
            "capital_cost": 483892.941,
            "capital_cost_vessel": 48385.243,
            "capital_cost_resin": 13292.384,
            "capital_cost_regen_tank": 29343.804,
            "capital_cost_backwash_tank": 27569.78,
            "flow_mass_regen_soln": 1116422.518,
            "total_pumping_power": 1.215287,
            "regeneration_tank_vol": 5139.677,
            "backwash_tank_vol": 7581.307,
            "direct_capital_cost": 241946.47,
        }

        for v, r in ix_cost_results.items():
            mv = getattr(m.fs.unit.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)
