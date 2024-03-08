#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
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
    value,
    ConcreteModel,
    Var,
    units as pyunits,
)
from idaes.core import (
    FlowsheetBlock,
    UnitModelCostingBlock,
)
from idaes.core.solvers import get_solver
from watertap.costing import WaterTAPCosting
from watertap.costing.unit_models.gac import cost_gac
from watertap.costing.unit_models.tests.unit_costing_test_harness import (
    UnitCostingTestHarness,
    DummyUnitModel,
)

__author__ = "Hunter Barber"

solver = get_solver()
zero = 1e-8
relative_tolerance = 1e-3


def build_pressure():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.unit = DummyUnitModel()
    m.fs.unit.default_costing_method = cost_gac

    m.fs.unit.bed_volume = Var(units=pyunits.meter**3)
    m.fs.unit.bed_volume.fix(8.900)
    m.fs.unit.bed_mass_gac = Var(units=pyunits.kg)
    m.fs.unit.bed_mass_gac.fix(4004)
    m.fs.unit.gac_usage_rate = Var(units=pyunits.kg / pyunits.second)
    m.fs.unit.gac_usage_rate.fix(0.0002925)

    m.fs.costing = WaterTAPCosting()
    m.fs.costing.base_currency = pyunits.USD_2020

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.costing.cost_process()

    return m


class TestGACPressureCosting(UnitCostingTestHarness):
    def configure(self):
        m = build_pressure()

        # arguments for UnitTestHarness
        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        cost_blk = m.fs.unit.costing
        self.cost_solutions[cost_blk.contactor_cost] = 56900
        self.cost_solutions[cost_blk.adsorbent_unit_cost] = 4.359
        self.cost_solutions[cost_blk.adsorbent_cost] = 17450
        self.cost_solutions[cost_blk.other_process_cost] = 81690
        self.cost_solutions[cost_blk.capital_cost] = 2 * 156000
        self.cost_solutions[cost_blk.gac_makeup_cost] = 12680
        self.cost_solutions[cost_blk.gac_regen_cost] = 27660
        self.cost_solutions[cost_blk.energy_consumption] = 0.01631
        self.cost_solutions[cost_blk.fixed_operating_cost] = 40370

        return m


def build_gravity():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.unit = DummyUnitModel()
    m.fs.unit.default_costing_method = cost_gac

    m.fs.unit.bed_volume = Var(units=pyunits.meter**3)
    m.fs.unit.bed_volume.fix(8.900)
    m.fs.unit.bed_mass_gac = Var(units=pyunits.kg)
    m.fs.unit.bed_mass_gac.fix(4004)
    m.fs.unit.gac_usage_rate = Var(units=pyunits.kg / pyunits.second)
    m.fs.unit.gac_usage_rate.fix(0.0002925)

    m.fs.costing = WaterTAPCosting()
    m.fs.costing.base_currency = pyunits.USD_2020

    m.fs.unit.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method_arguments={"contactor_type": "gravity"},
    )
    m.fs.costing.cost_process()

    return m


class TestGACGravityCosting(UnitCostingTestHarness):
    def configure(self):
        m = build_gravity()

        # arguments for UnitTestHarness
        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        cost_blk = m.fs.unit.costing
        self.cost_solutions[cost_blk.contactor_cost] = 163200
        self.cost_solutions[cost_blk.adsorbent_unit_cost] = 4.359
        self.cost_solutions[cost_blk.adsorbent_cost] = 17450
        self.cost_solutions[cost_blk.other_process_cost] = 159500
        self.cost_solutions[cost_blk.capital_cost] = 2 * 340200
        self.cost_solutions[cost_blk.gac_makeup_cost] = 12680
        self.cost_solutions[cost_blk.gac_regen_cost] = 27660
        self.cost_solutions[cost_blk.energy_consumption] = 2.476
        self.cost_solutions[cost_blk.fixed_operating_cost] = 40370

        return m


def build_modular():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.unit = DummyUnitModel()
    m.fs.unit.default_costing_method = cost_gac

    m.fs.unit.bed_volume = Var(units=pyunits.meter**3)
    m.fs.unit.bed_volume.fix(8.900)
    m.fs.unit.bed_mass_gac = Var(units=pyunits.kg)
    m.fs.unit.bed_mass_gac.fix(4004)
    m.fs.unit.gac_usage_rate = Var(units=pyunits.kg / pyunits.second)
    m.fs.unit.gac_usage_rate.fix(0.0002925)

    m.fs.costing = WaterTAPCosting()
    m.fs.costing.base_currency = pyunits.USD_2020

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.costing.cost_process()

    m.fs.costing.gac_pressure.num_contactors_op.fix(4)
    m.fs.costing.gac_pressure.num_contactors_redundant.fix(2)

    return m


class TestGACModularCosting(UnitCostingTestHarness):
    def configure(self):
        m = build_modular()

        # arguments for UnitTestHarness
        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        cost_blk = m.fs.unit.costing
        self.cost_solutions[cost_blk.contactor_cost] = 89040
        self.cost_solutions[cost_blk.other_process_cost] = 69690
        self.cost_solutions[cost_blk.capital_cost] = 2 * 176200

        return m


def build_max():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.unit = DummyUnitModel()
    m.fs.unit.default_costing_method = cost_gac

    m.fs.unit.bed_volume = Var(units=pyunits.meter**3)
    m.fs.unit.bed_volume.fix(10 * 8.900)
    m.fs.unit.bed_mass_gac = Var(units=pyunits.kg)
    m.fs.unit.bed_mass_gac.fix(10 * 4004)
    m.fs.unit.gac_usage_rate = Var(units=pyunits.kg / pyunits.second)
    m.fs.unit.gac_usage_rate.fix(10 * 0.0002925)

    m.fs.costing = WaterTAPCosting()
    m.fs.costing.base_currency = pyunits.USD_2020

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.costing.cost_process()

    return m


class TestGACMaxCosting(UnitCostingTestHarness):
    def configure(self):
        m = build_max()

        # arguments for UnitTestHarness
        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        cost_blk = m.fs.unit.costing
        assert value(m.fs.unit.bed_mass_gac) > value(
            m.fs.costing.gac_pressure.bed_mass_max_ref
        )
        self.cost_solutions[cost_blk.bed_mass_gac_ref] = value(
            m.fs.costing.gac_pressure.bed_mass_max_ref
        )
        self.cost_solutions[cost_blk.adsorbent_unit_cost] = 3.651
        self.cost_solutions[cost_blk.adsorbent_cost] = 146200

        return m
