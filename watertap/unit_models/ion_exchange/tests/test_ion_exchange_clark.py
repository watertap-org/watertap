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
    10.1061/(asce)ee.1943-7870.0001997
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
    ix.freundlich_n[target_component].fix(2)
    ix.bv_50[target_component].fix(bv50)
    ix.mass_transfer_coeff[target_component].fix(0.15934630)
    ix.resin_density.fix(720)
    ix.resin_diam.fix(7.25e-4)
    ix.number_columns.fix(1)
    ix.c_norm[target_component].fix(0.15)
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
        self.unit_solutions[m.fs.unit.bv_50["CrO4_2-"]] = 26173

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
        self.unit_solutions[m.fs.unit.bv_50["Cl_-"]] = 20000
        self.unit_solutions[m.fs.unit.c_traps["Cl_-", 0]] = 0
        self.unit_solutions[m.fs.unit.c_traps["Cl_-", 1]] = 0.01
        self.unit_solutions[m.fs.unit.c_traps["Cl_-", 2]] = 0.07
        self.unit_solutions[m.fs.unit.c_traps["Cl_-", 3]] = 0.13
        self.unit_solutions[m.fs.unit.c_traps["Cl_-", 4]] = 0.19
        self.unit_solutions[m.fs.unit.c_traps["Cl_-", 5]] = 0.25
        self.unit_solutions[m.fs.unit.tb_traps["Cl_-", 0]] = 0
        self.unit_solutions[m.fs.unit.tb_traps["Cl_-", 1]] = 3344557
        self.unit_solutions[m.fs.unit.tb_traps["Cl_-", 2]] = 3825939
        self.unit_solutions[m.fs.unit.tb_traps["Cl_-", 3]] = 4034117
        self.unit_solutions[m.fs.unit.tb_traps["Cl_-", 4]] = 4188551
        self.unit_solutions[m.fs.unit.tb_traps["Cl_-", 5]] = 4320000
        self.unit_solutions[m.fs.unit.traps["Cl_-", 1]] = 0.0038710
        self.unit_solutions[m.fs.unit.traps["Cl_-", 2]] = 0.0044572
        self.unit_solutions[m.fs.unit.traps["Cl_-", 3]] = 0.0048189
        self.unit_solutions[m.fs.unit.traps["Cl_-", 4]] = 0.0057197
        self.unit_solutions[m.fs.unit.traps["Cl_-", 5]] = 0.0066941
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
        self.unit_solutions[m.fs.unit.c_traps["Ca_2+", 0]] = 0.0
        self.unit_solutions[m.fs.unit.c_traps["Ca_2+", 1]] = 0.01
        self.unit_solutions[m.fs.unit.c_traps["Ca_2+", 2]] = 0.09375
        self.unit_solutions[m.fs.unit.c_traps["Ca_2+", 3]] = 0.1775
        self.unit_solutions[m.fs.unit.c_traps["Ca_2+", 4]] = 0.26125
        self.unit_solutions[m.fs.unit.c_traps["Ca_2+", 5]] = 0.345
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


def build_clark_mc():
    """
    Clark parameters fit from Franke Case Study
    """
    target_component = "PFOA"
    solute_list = ["PFOA", "PFHxA", "PFBS", "PFHxS", "PFPeS"]

    ion_props = {
        "solute_list": solute_list,
        "mw_data": {
            "PFOA": 0.41407,
            "PFHxA": 0.31405,
            "PFBS": 0.3001,
            "PFHxS": 0.40012,
            "PFPeS": 0.350107,
            "H2O": 0.018,
        },
        "diffus_calculation": "HaydukLaudie",
        "charge": {"PFOA": -1, "PFHxA": -1, "PFBS": -1, "PFHxS": -1, "PFPeS": -1},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**ion_props)
    c0_dict = {
        "PFOA": 2.2e-08,
        "PFHxA": 2.5e-08,
        "PFBS": 3.3e-08,
        "PFHxS": 3e-07,
        "PFPeS": 5.9e-08,
    }

    ix_config = {
        "property_package": m.fs.properties,
        "target_component": target_component,
        "reactive_components": solute_list,
        "add_steady_state_approximation": True,
    }

    m.fs.unit = ix = IonExchangeClark(**ix_config)

    loading_rate = 0.005 * pyunits.m / pyunits.s

    freund_n_dict = {
        "PFOA": 1.05,
        "PFHxA": 1.05,
        "PFBS": 1.05,
        "PFHxS": 1.05,
        "PFPeS": 1.3734,
    }
    k_dict = {
        "PFOA": 0.171560619,
        "PFHxA": 0.223008672,
        "PFBS": 0.13932198,
        "PFHxS": 0.183103632,
        "PFPeS": 0.050370161,
    }
    bv50_dict = {
        "PFOA": 26067.25,
        "PFHxA": 17909.56,
        "PFBS": 52860.16,
        "PFHxS": 77802.27,
        "PFPeS": 53248.92,
    }

    for s in solute_list:
        m.fs.unit.freundlich_n[s].fix(freund_n_dict[s])
        m.fs.unit.bv_50[s].fix(bv50_dict[s])
        m.fs.unit.mass_transfer_coeff[s].fix(k_dict[s])

    ix.resin_density.fix(720)
    ix.resin_diam.fix(7.25e-4)
    ix.bed_porosity.fix(0.4)
    ix.loading_rate.fix(loading_rate)
    ix.number_columns.fix(2)
    ix.number_columns_redundant.fix(1)
    ix.c_norm[target_component].fix(0.65)

    var_dict = {
        ("flow_vol_phase", "Liq"): 0.04381,  # m3/s
        ("pressure", None): 101325,
        ("temperature", None): 298,
    }
    ix.process_flow.properties_in[0].conc_mass_phase_comp
    ix.process_flow.properties_in[0].flow_vol_phase
    ix.process_flow.properties_in[0].flow_mol_phase_comp

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-1, index=("Liq", "H2O")
    )
    for s, c in c0_dict.items():
        var_dict[("conc_mass_phase_comp", ("Liq", s))] = c * pyunits.g / pyunits.liter
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e9, index=("Liq", s)
        )
        m.fs.unit.process_flow.properties_in[0].conc_mass_phase_comp["Liq", s].setlb(
            c * 0.9999
        )
    iscale.set_scaling_factor(ix.traps["PFBS", 1], 1e6)
    iscale.calculate_scaling_factors(m)
    ix.process_flow.properties_in.calculate_state(
        var_args=var_dict,
        hold_state=True,
    )

    return m


@pytest.mark.skip
# BUG: Will always end up with tb_traps[*,1] very small for one of the components
# removing problematic component from stream will shift error to another component.
# Will solve to an acceptable level but terminate with:
#   Cannot recompute multipliers for feasibility problem.  Error in eq_mult_calculator
# Possibly need to rethink tb_traps constraint for multi-component Clark model
# or could be a mass balance issue.
class TestIXClarkMultiComponent(UnitTestHarness):
    def configure(self):
        m = build_clark_mc()

        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance
        self.default_small = 1e-6
        self.unit_solutions[m.fs.unit.N_Pe_bed] = 183.2591
        self.unit_solutions[m.fs.unit.N_Pe_particle] = 0.092776471132
        self.unit_solutions[m.fs.unit.N_Re] = 3.625
        self.unit_solutions[m.fs.unit.bed_depth] = 1.4320751
        self.unit_solutions[m.fs.unit.bed_depth_to_diam_ratio] = 0.606351016427
        self.unit_solutions[m.fs.unit.bed_diameter] = 2.36179221
        self.unit_solutions[m.fs.unit.bed_porosity] = 0.4
        self.unit_solutions[m.fs.unit.bed_volume] = 6.27392104
        self.unit_solutions[m.fs.unit.bed_volume_total] = 12.5478
        self.unit_solutions[m.fs.unit.breakthrough_time] = 8931464.6158
        self.unit_solutions[m.fs.unit.bv] = 31183.6459
        self.unit_solutions[m.fs.unit.bv_50["PFBS"]] = 52860.16
        self.unit_solutions[m.fs.unit.bv_50["PFHxA"]] = 17909.56
        self.unit_solutions[m.fs.unit.bv_50["PFHxS"]] = 77802.27
        self.unit_solutions[m.fs.unit.bv_50["PFOA"]] = 26067.25
        self.unit_solutions[m.fs.unit.bv_50["PFPeS"]] = 53248.92
        self.unit_solutions[m.fs.unit.c_norm["PFBS"]] = 0.214854203634
        self.unit_solutions[m.fs.unit.c_norm["PFHxA"]] = 0.936114402191
        self.unit_solutions[m.fs.unit.c_norm["PFHxS"]] = 0.04349564889
        self.unit_solutions[m.fs.unit.c_norm["PFOA"]] = 0.65
        self.unit_solutions[m.fs.unit.c_norm["PFPeS"]] = 0.028949706047
        self.unit_solutions[m.fs.unit.c_norm_avg["PFBS"]] = 0.086331035991
        self.unit_solutions[m.fs.unit.c_norm_avg["PFHxA"]] = 0.400484932106
        self.unit_solutions[m.fs.unit.c_norm_avg["PFHxS"]] = 0.012622265084
        self.unit_solutions[m.fs.unit.c_norm_avg["PFOA"]] = 0.235666590501
        self.unit_solutions[m.fs.unit.c_norm_avg["PFPeS"]] = 0.007115781791
        self.unit_solutions[m.fs.unit.c_traps["PFBS", 0]] = 0.0
        self.unit_solutions[m.fs.unit.c_traps["PFBS", 1]] = 0.010000250243
        self.unit_solutions[m.fs.unit.c_traps["PFBS", 2]] = 0.061213591755
        self.unit_solutions[m.fs.unit.c_traps["PFBS", 3]] = 0.112427124055
        self.unit_solutions[m.fs.unit.c_traps["PFBS", 4]] = 0.163640668003
        self.unit_solutions[m.fs.unit.c_traps["PFBS", 5]] = 0.21485421527
        self.unit_solutions[m.fs.unit.c_traps["PFHxA", 0]] = 0.0
        self.unit_solutions[m.fs.unit.c_traps["PFHxA", 1]] = 0.010000250243
        self.unit_solutions[m.fs.unit.c_traps["PFHxA", 2]] = 0.241528610898
        self.unit_solutions[m.fs.unit.c_traps["PFHxA", 3]] = 0.47305720638
        self.unit_solutions[m.fs.unit.c_traps["PFHxA", 4]] = 0.704585805191
        self.unit_solutions[m.fs.unit.c_traps["PFHxA", 5]] = 0.936114404861
        self.unit_solutions[m.fs.unit.c_traps["PFHxS", 0]] = 0.0
        self.unit_solutions[m.fs.unit.c_traps["PFHxS", 1]] = 0.010000250243
        self.unit_solutions[m.fs.unit.c_traps["PFHxS", 2]] = 0.018374048358
        self.unit_solutions[m.fs.unit.c_traps["PFHxS", 3]] = 0.026747917945
        self.unit_solutions[m.fs.unit.c_traps["PFHxS", 4]] = 0.035121807869
        self.unit_solutions[m.fs.unit.c_traps["PFHxS", 5]] = 0.04349570638
        self.unit_solutions[m.fs.unit.c_traps["PFOA", 0]] = 0.0
        self.unit_solutions[m.fs.unit.c_traps["PFOA", 1]] = 0.010000250243
        self.unit_solutions[m.fs.unit.c_traps["PFOA", 2]] = 0.170000014706
        self.unit_solutions[m.fs.unit.c_traps["PFOA", 3]] = 0.330000007575
        self.unit_solutions[m.fs.unit.c_traps["PFOA", 4]] = 0.490000005102
        self.unit_solutions[m.fs.unit.c_traps["PFOA", 5]] = 0.650000003846
        self.unit_solutions[m.fs.unit.c_traps["PFPeS", 0]] = 0.0
        self.unit_solutions[m.fs.unit.c_traps["PFPeS", 1]] = 0.010000250243
        self.unit_solutions[m.fs.unit.c_traps["PFPeS", 2]] = 0.014737596261
        self.unit_solutions[m.fs.unit.c_traps["PFPeS", 3]] = 0.019474981459
        self.unit_solutions[m.fs.unit.c_traps["PFPeS", 4]] = 0.024212382831
        self.unit_solutions[m.fs.unit.c_traps["PFPeS", 5]] = 0.028949792433
        self.unit_solutions[m.fs.unit.column_height] = 4.1439469
        self.unit_solutions[m.fs.unit.ebct] = 286.415
        self.unit_solutions[m.fs.unit.freundlich_n["PFBS"]] = 1.05
        self.unit_solutions[m.fs.unit.freundlich_n["PFHxA"]] = 1.05
        self.unit_solutions[m.fs.unit.freundlich_n["PFHxS"]] = 1.05
        self.unit_solutions[m.fs.unit.freundlich_n["PFOA"]] = 1.05
        self.unit_solutions[m.fs.unit.freundlich_n["PFPeS"]] = 1.3734
        self.unit_solutions[m.fs.unit.loading_rate] = 0.005
        self.unit_solutions[m.fs.unit.mass_transfer_coeff["PFBS"]] = 0.13932198
        self.unit_solutions[m.fs.unit.mass_transfer_coeff["PFHxA"]] = 0.223008672
        self.unit_solutions[m.fs.unit.mass_transfer_coeff["PFHxS"]] = 0.183103632
        self.unit_solutions[m.fs.unit.mass_transfer_coeff["PFOA"]] = 0.171560619
        self.unit_solutions[m.fs.unit.mass_transfer_coeff["PFPeS"]] = 0.050370161
        self.unit_solutions[m.fs.unit.number_columns] = 2.0
        self.unit_solutions[m.fs.unit.number_columns_redundant] = 1.0
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "H2O"]
        ] = 0.0
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "PFBS"]
        ] = -4.401e-09
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "PFHxA"]
        ] = -2.09e-09
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "PFHxS"]
        ] = -3.2433e-08
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "PFOA"]
        ] = -1.779e-09
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "PFPeS"]
        ] = -7.33e-09
        self.unit_solutions[m.fs.unit.resin_density] = 720.0
        self.unit_solutions[m.fs.unit.resin_diam] = 0.000725
        self.unit_solutions[m.fs.unit.service_flow_rate] = 12.5691
        self.unit_solutions[m.fs.unit.tb_traps["PFBS", 0]] = 0.0
        self.unit_solutions[m.fs.unit.tb_traps["PFBS", 1]] = 11763.6423
        self.unit_solutions[m.fs.unit.tb_traps["PFBS", 2]] = 4159543.5962
        self.unit_solutions[m.fs.unit.tb_traps["PFBS", 3]] = 6139656.7296
        self.unit_solutions[m.fs.unit.tb_traps["PFBS", 4]] = 7642043.6123
        self.unit_solutions[m.fs.unit.tb_traps["PFBS", 5]] = 8931464.8934
        self.unit_solutions[m.fs.unit.tb_traps["PFHxA", 0]] = 0.0
        self.unit_solutions[m.fs.unit.tb_traps["PFHxA", 1]] = 1927421.8274
        self.unit_solutions[m.fs.unit.tb_traps["PFHxA", 2]] = 3947323.8328
        self.unit_solutions[m.fs.unit.tb_traps["PFHxA", 3]] = 5003845.1843
        self.unit_solutions[m.fs.unit.tb_traps["PFHxA", 4]] = 6240248.8316
        self.unit_solutions[m.fs.unit.tb_traps["PFHxA", 5]] = 8931464.6853
        self.unit_solutions[m.fs.unit.tb_traps["PFHxS", 0]] = 0.0
        self.unit_solutions[m.fs.unit.tb_traps["PFHxS", 1]] = 5341407.7834
        self.unit_solutions[m.fs.unit.tb_traps["PFHxS", 2]] = 6679266.1388
        self.unit_solutions[m.fs.unit.tb_traps["PFHxS", 3]] = 7600051.9337
        self.unit_solutions[m.fs.unit.tb_traps["PFHxS", 4]] = 8324076.2888
        self.unit_solutions[m.fs.unit.tb_traps["PFHxS", 5]] = 8931468.4867
        self.unit_solutions[m.fs.unit.tb_traps["PFOA", 0]] = 0.0
        self.unit_solutions[m.fs.unit.tb_traps["PFOA", 1]] = 1407686.5249
        self.unit_solutions[m.fs.unit.tb_traps["PFOA", 2]] = 4531006.2772
        self.unit_solutions[m.fs.unit.tb_traps["PFOA", 3]] = 6007004.7055
        self.unit_solutions[m.fs.unit.tb_traps["PFOA", 4]] = 7377203.3815
        self.unit_solutions[m.fs.unit.tb_traps["PFOA", 5]] = 8931464.658
        self.unit_solutions[m.fs.unit.tb_traps["PFPeS", 0]] = 0.0
        self.unit_solutions[m.fs.unit.tb_traps["PFPeS", 1]] = 7489488.4335
        self.unit_solutions[m.fs.unit.tb_traps["PFPeS", 2]] = 7997389.2178
        self.unit_solutions[m.fs.unit.tb_traps["PFPeS", 3]] = 8374321.3446
        self.unit_solutions[m.fs.unit.tb_traps["PFPeS", 4]] = 8676951.0166
        self.unit_solutions[m.fs.unit.tb_traps["PFPeS", 5]] = 8931468.916
        self.unit_solutions[m.fs.unit.traps["PFBS", 1]] = 6.58567e-06
        self.unit_solutions[m.fs.unit.traps["PFBS", 2]] = 0.016535884639
        self.unit_solutions[m.fs.unit.traps["PFBS", 3]] = 0.019248144956
        self.unit_solutions[m.fs.unit.traps["PFBS", 4]] = 0.023219070694
        self.unit_solutions[m.fs.unit.traps["PFBS", 5]] = 0.02732135003
        self.unit_solutions[m.fs.unit.traps["PFHxA", 1]] = 0.001079033578
        self.unit_solutions[m.fs.unit.traps["PFHxA", 2]] = 0.028442347863
        self.unit_solutions[m.fs.unit.traps["PFHxA", 3]] = 0.042264913989
        self.unit_solutions[m.fs.unit.traps["PFHxA", 4]] = 0.081511944904
        self.unit_solutions[m.fs.unit.traps["PFHxA", 5]] = 0.247186691769
        self.unit_solutions[m.fs.unit.traps["PFHxS", 1]] = 0.002990292949
        self.unit_solutions[m.fs.unit.traps["PFHxS", 2]] = 0.00212511484
        self.unit_solutions[m.fs.unit.traps["PFHxS", 3]] = 0.002325914583
        self.unit_solutions[m.fs.unit.traps["PFHxS", 4]] = 0.002507716866
        self.unit_solutions[m.fs.unit.traps["PFHxS", 5]] = 0.002673225844
        self.unit_solutions[m.fs.unit.traps["PFOA", 1]] = 0.000788068813
        self.unit_solutions[m.fs.unit.traps["PFOA", 2]] = 0.031472910909
        self.unit_solutions[m.fs.unit.traps["PFOA", 3]] = 0.041314570191
        self.unit_solutions[m.fs.unit.traps["PFOA", 4]] = 0.062899142228
        self.unit_solutions[m.fs.unit.traps["PFOA", 5]] = 0.099191898357
        self.unit_solutions[m.fs.unit.traps["PFPeS", 1]] = 0.004192857817
        self.unit_solutions[m.fs.unit.traps["PFPeS", 2]] = 0.000703376553
        self.unit_solutions[m.fs.unit.traps["PFPeS", 3]] = 0.000721931622
        self.unit_solutions[m.fs.unit.traps["PFPeS", 4]] = 0.000740141003
        self.unit_solutions[m.fs.unit.traps["PFPeS", 5]] = 0.000757474795

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "PFOA"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "PFBS"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "PFHxA"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "PFHxS"
                ]
                + m.fs.unit.process_flow.properties_in[0.0].flow_mol_phase_comp[
                    "Liq", "PFPeS"
                ],
                "out": m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "PFOA"
                ]
                + m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "PFBS"
                ]
                + m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "PFHxA"
                ]
                + m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "PFHxS"
                ]
                + m.fs.unit.process_flow.properties_out[0.0].flow_mol_phase_comp[
                    "Liq", "PFPeS"  ###
                ]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "PFOA"]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "PFBS"]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "PFHxA"]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp["Liq", "PFHxS"]
                + m.fs.unit.regeneration_stream[0.0].flow_mol_phase_comp[
                    "Liq", "PFPeS"
                ],
            },
        }
        return m
