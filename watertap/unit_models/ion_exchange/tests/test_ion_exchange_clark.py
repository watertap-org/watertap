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
)
import idaes.core.util.scaling as iscale
from idaes.core.util.testing import initialization_tester

from watertap.costing import WaterTAPCosting
from watertap.core.util.initialization import check_dof
from watertap.core.solvers import get_solver
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models import IonExchangeClark
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

__author__ = "Kurban Sitterley"

solver = get_solver()
zero = 1e-12
relative_tolerance = 1e-3


def build_clark1():
    """
    Test case adapted from experiment 2 in
    https://doi.org/10.3390/w17030329
    """

    target_component = "CrO4_2-"
    solute_list = [target_component, "NO3_-", "HCO3_-", "SO4_2-"]

    ion_props = {
        "solute_list": solute_list,
        "mw_data": {
            "H2O": 0.018,
            target_component: 0.116,
            "NO3_-": 0.062,
            "HCO3_-": 0.061,
            "SO4_2-": 0.096,
        },
        "molar_volume_data": {
            ("Liq", target_component): 4.3e-5,
            ("Liq", "NO3_-"): 3.45e-5,
            ("Liq", "HCO3_-"): 2.89e-5,
            ("Liq", "SO4_2-"): 2.5e-5,
        },
        "diffus_calculation": "HaydukLaudie",
        "charge": {target_component: -2, "NO3_-": -1, "HCO3_-": -1, "SO4_2-": -2},
    }

    bed_volume = 5712 * pyunits.milliliter
    bed_diameter = 21 * pyunits.centimeter
    bed_depth = pyunits.convert(
        bed_volume / (3.14159 * (bed_diameter / 2) ** 2), to_units=pyunits.m
    )
    loading_rate = pyunits.convert(
        22 * pyunits.milliliter / pyunits.minute / pyunits.cm**2,
        to_units=pyunits.m / pyunits.s,
    )

    tb50 = 327 * pyunits.hours
    bv50 = pyunits.convert(
        (loading_rate * tb50) / bed_depth, to_units=pyunits.dimensionless
    )

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**ion_props)
    ix_config = {
        "property_package": m.fs.properties,
        "target_component": target_component,
        "add_steady_state_approximation": False,
    }
    m.fs.unit = ix = IonExchangeClark(**ix_config)

    ix.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix(7.0513)
    ix.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "CrO4_2-"].fix(
        1.1264e-08
    )
    ix.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "NO3_-"].fix(2.8665e-06)
    ix.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "HCO3_-"].fix(
        3.4275e-04
    )
    ix.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "SO4_2-"].fix(
        8.3307e-06
    )
    ix.process_flow.properties_in[0].pressure.fix(101325)
    ix.process_flow.properties_in[0].temperature.fix(298)

    ix.bed_depth.setlb(0)
    ix.bed_diameter.setlb(0)
    ix.ebct.setlb(0)
    ix.freundlich_n.fix(2)
    ix.bv_50.fix(bv50)
    ix.mass_transfer_coeff.fix(0.15934630)
    ix.resin_density.fix(720)
    ix.resin_diam.fix(7.25e-4)
    ix.number_columns.fix(1)
    ix.c_norm.fix(0.15)
    ix.bed_depth.fix(bed_depth)
    ix.ebct.fix(45)

    m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e6, index=("Liq", target_component)
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e6, index=("Liq", "NO3_-")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e4, index=("Liq", "HCO3_-")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e6, index=("Liq", "SO4_2-")
    )

    iscale.calculate_scaling_factors(m)

    return m


class TestIXClark1(UnitTestHarness):
    def configure(self):
        m = build_clark1()

        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        self.unit_solutions[m.fs.unit.bed_volume] = 0.0057125
        self.unit_solutions[m.fs.unit.bed_depth] = 0.164914
        self.unit_solutions[m.fs.unit.bed_diameter] = 0.210
        self.unit_solutions[m.fs.unit.breakthrough_time] = 892889
        self.unit_solutions[m.fs.unit.bv] = 19841
        self.unit_solutions[m.fs.unit.bv_50] = 26173

        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "H2O"]
        ] = 0
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "CrO4_2-"]
        ] = -9.5744e-09
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "NO3_-"]
        ] = 0
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "HCO3_-"]
        ] = 0
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "SO4_2-"]
        ] = 0

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "CrO4_2-"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "CrO4_2-"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp[
                    "Liq", "CrO4_2-"
                ],
            },
            "Check 2": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "CrO4_2-"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "CrO4_2-"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp[
                    "Liq", "CrO4_2-"
                ],
            },
        }

        return m


def build_clark2():

    target_component = "Cl_-"

    ion_props = {
        "solute_list": [target_component],
        "diffusivity_data": {("Liq", target_component): 1e-9},
        "mw_data": {"H2O": 0.018, target_component: 35.45e-3},
        "charge": {target_component: -1},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**ion_props)
    ix_config = {
        "property_package": m.fs.properties,
        "target_component": target_component,
        "regenerant": "NaCl",
    }
    m.fs.unit = ix = IonExchangeClark(**ix_config)

    m.fs.unit.process_flow.properties_in.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 0.5,
            ("conc_mass_phase_comp", ("Liq", target_component)): 1e-6,
            ("pressure", None): 101325,
            ("temperature", None): 298,
        },
        hold_state=True,
    )

    m.fs.unit.freundlich_n.fix(1.2)
    m.fs.unit.bv_50.fix(20000)
    ix.mass_transfer_coeff.fix(0.15934630)
    m.fs.unit.resin_density.fix(720)
    m.fs.unit.resin_diam.fix(6.75e-4)
    m.fs.unit.c_norm.fix(0.25)

    m.fs.unit.number_columns.fix(16)
    m.fs.unit.bed_depth.fix(1.476)
    m.fs.unit.bed_porosity.fix(0.5)
    m.fs.unit.ebct.fix(240)

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e6, index=("Liq", "Cl_-")
    )
    iscale.calculate_scaling_factors(m)

    return m


class TestIXClark2(UnitTestHarness):
    def configure(self):
        m = build_clark2()

        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance
        self.unit_solutions[m.fs.unit.bed_volume] = 7.5
        self.unit_solutions[m.fs.unit.bed_volume_total] = 120
        self.unit_solutions[m.fs.unit.bed_depth] = 1.476
        self.unit_solutions[m.fs.unit.bed_porosity] = 0.5
        self.unit_solutions[m.fs.unit.column_height] = 4.240
        self.unit_solutions[m.fs.unit.bed_diameter] = 2.543
        self.unit_solutions[m.fs.unit.breakthrough_time] = 4320000
        self.unit_solutions[m.fs.unit.cycle_time] = 4327205
        self.unit_solutions[m.fs.unit.loading_rate] = 0.00615
        self.unit_solutions[m.fs.unit.bv_50] = 20000
        self.unit_solutions[m.fs.unit.c_traps[0]] = 0
        self.unit_solutions[m.fs.unit.c_traps[1]] = 0.01
        self.unit_solutions[m.fs.unit.c_traps[2]] = 0.07
        self.unit_solutions[m.fs.unit.c_traps[3]] = 0.13
        self.unit_solutions[m.fs.unit.c_traps[4]] = 0.19
        self.unit_solutions[m.fs.unit.c_traps[5]] = 0.25
        self.unit_solutions[m.fs.unit.tb_traps[0]] = 0
        self.unit_solutions[m.fs.unit.tb_traps[1]] = 3344557
        self.unit_solutions[m.fs.unit.tb_traps[2]] = 3825939
        self.unit_solutions[m.fs.unit.tb_traps[3]] = 4034117
        self.unit_solutions[m.fs.unit.tb_traps[4]] = 4188551
        self.unit_solutions[m.fs.unit.tb_traps[5]] = 4320000
        self.unit_solutions[m.fs.unit.traps[1]] = 0.0038710
        self.unit_solutions[m.fs.unit.traps[2]] = 0.0044572
        self.unit_solutions[m.fs.unit.traps[3]] = 0.0048189
        self.unit_solutions[m.fs.unit.traps[4]] = 0.0057197
        self.unit_solutions[m.fs.unit.traps[5]] = 0.0066941
        self.unit_solutions[m.fs.unit.c_norm_avg["Cl_-"]] = 0.02556
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "H2O"]
        ] = 0
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Cl_-"]
        ] = -1.37438e-05

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "Cl_-"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "Cl_-"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "Cl_-"],
            },
            "Check 2": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "Cl_-"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "Cl_-"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "Cl_-"],
            },
        }
        return m

    @pytest.mark.component
    def test_solution(self):
        m = build_clark2()

        m.fs.costing = WaterTAPCosting()
        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(
            m.fs.unit.process_flow.properties_in[0].flow_vol_phase["Liq"]
        )
        m.fs.costing.add_specific_energy_consumption(
            m.fs.unit.process_flow.properties_in[0].flow_vol_phase["Liq"], name="SEC"
        )

        check_dof(m, fail_flag=True)
        initialization_tester(m)
        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_flow_electricity": 30.79,
            "aggregate_flow_costs": {"electricity": 18893.46, "NaCl": 29877.81},
            "total_capital_cost": 6042905.25,
            "total_operating_cost": 280105.99,
            "total_variable_operating_cost": 43894.15,
            "total_annualized_cost": 884396.52,
            "LCOW": 0.0623,
            "SEC": 0.0171,
        }

        for v, r in sys_cost_results.items():
            mv = getattr(m.fs.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

        ix_cost_results = {
            "capital_cost": 6042905.25,
            "fixed_operating_cost": 54924.68,
            "capital_cost_vessel": 85841.84,
            "capital_cost_resin": 54924.68,
            "capital_cost_regen_tank": 12515.54,
            "capital_cost_backwash_tank": 193606.47,
            "flow_mass_regen_soln": 328177.61,
            "total_pumping_power": 30.79,
            "regeneration_tank_vol": 1597.81,
            "backwash_tank_vol": 369626.09,
            "direct_capital_cost": 3021452.62,
        }

        for v, r in ix_cost_results.items():
            mv = getattr(m.fs.unit.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)


def build_inert():

    inert = "Cl_-"
    target_component = "Ca_2+"

    ion_props = {
        "solute_list": [target_component, inert],
        "diffusivity_data": {("Liq", target_component): 1e-9, ("Liq", inert): 9.2e-10},
        "mw_data": {"H2O": 0.018, target_component: 35.45e-3, inert: 40e-3},
        "charge": {target_component: -1, inert: 2},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**ion_props)
    ix_config = {
        "property_package": m.fs.properties,
        "target_component": target_component,
        "regenerant": "single_use",
    }
    m.fs.unit = ix = IonExchangeClark(**ix_config)

    m.fs.unit.process_flow.properties_in.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 0.5,
            ("conc_mass_phase_comp", ("Liq", target_component)): 1e-6,
            ("conc_mass_phase_comp", ("Liq", inert)): 1e-7,
            ("pressure", None): 101325,
            ("temperature", None): 298,
        },
        hold_state=True,
    )
    m.fs.unit.freundlich_n.fix(1.2)
    m.fs.unit.bv_50.fix(20000)
    ix.mass_transfer_coeff.fix(0.15934630)
    m.fs.unit.resin_density.fix(720)
    m.fs.unit.resin_diam.fix(6.75e-4)
    m.fs.unit.c_norm.fix(0.345)

    m.fs.unit.number_columns.fix(16)
    m.fs.unit.bed_depth.fix(1.476)
    m.fs.unit.bed_porosity.fix(0.5)
    m.fs.unit.ebct.fix(240)

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e5, index=("Liq", "Cl_-")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e5, index=("Liq", "Ca_2+")
    )
    iscale.calculate_scaling_factors(m)

    return m


class TestIXClarkWithInert(UnitTestHarness):
    def configure(self):
        m = build_inert()

        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        self.unit_solutions[m.fs.unit.bed_porosity] = 0.5
        self.unit_solutions[m.fs.unit.breakthrough_time] = 4506965
        self.unit_solutions[m.fs.unit.bv] = 18779
        self.unit_solutions[m.fs.unit.c_traps[0]] = 0.0
        self.unit_solutions[m.fs.unit.c_traps[1]] = 0.01
        self.unit_solutions[m.fs.unit.c_traps[2]] = 0.09375
        self.unit_solutions[m.fs.unit.c_traps[3]] = 0.1775
        self.unit_solutions[m.fs.unit.c_traps[4]] = 0.26125
        self.unit_solutions[m.fs.unit.c_traps[5]] = 0.345
        self.unit_solutions[m.fs.unit.column_height] = 4.24
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Ca_2+"]
        ] = -1.4033e-05
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Cl_-"]
        ] = 0.0
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "H2O"]
        ] = 0.0

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "Cl_-"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "Ca_2+"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "Cl_-"
                ]
                + m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "Ca_2+"
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "H2O"]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "Cl_-"]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp[
                    "Liq", "Ca_2+"
                ],
            },
        }
        return m

    @pytest.mark.component
    def test_costing(self):

        m = build_inert()

        m.fs.costing = WaterTAPCosting()
        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(
            m.fs.unit.process_flow.properties_out[0].flow_vol_phase["Liq"]
        )
        m.fs.costing.add_specific_energy_consumption(
            m.fs.unit.process_flow.properties_out[0].flow_vol_phase["Liq"], name="SEC"
        )

        check_dof(m, fail_flag=True)
        initialization_tester(m)
        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        sys_cost_results = {
            "aggregate_flow_electricity": 30.78,
            "aggregate_flow_costs": {"electricity": 18891.32},
            "total_capital_cost": 6017874.15,
            "total_operating_cost": 6350826.88,
            "total_fixed_operating_cost": 6333824.69,
            "LCOW": 0.489589,
            "SEC": 0.017103,
        }

        for v, r in sys_cost_results.items():
            mv = getattr(m.fs.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)

        ix_cost_results = {
            "capital_cost": 6017874.15,
            "fixed_operating_cost": 6153288.46,
            "capital_cost_vessel": 85841.84,
            "capital_cost_resin": 54924.68,
            "capital_cost_backwash_tank": 193606.47,
            "total_pumping_power": 30.78,
            "flow_vol_resin": 840.23,
            "single_use_resin_replacement_cost": 6153288.46,
            "backwash_tank_vol": 369626.09,
            "direct_capital_cost": 3008937.07,
        }

        for v, r in ix_cost_results.items():
            mv = getattr(m.fs.unit.costing, v)
            if mv.is_indexed():
                for i, s in r.items():
                    assert pytest.approx(s, rel=1e-3) == value(mv[i])
            else:
                assert pytest.approx(r, rel=1e-3) == value(mv)
