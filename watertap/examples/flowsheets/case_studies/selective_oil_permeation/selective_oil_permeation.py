###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

import os

from pyomo.environ import (
    ConcreteModel,
    units as pyunits,
    assert_optimal_termination,
    Expression,
    value,
    TransformationFactory,
)

from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.models.unit_models import Product
import idaes.core.util.scaling as iscale
from idaes.core import UnitModelCostingBlock

from watertap.core.util.initialization import assert_degrees_of_freedom

from watertap.core.wt_database import Database
import watertap.core.zero_order_properties as prop_ZO
from watertap.unit_models.zero_order import (
    FeedZO,
    SelectiveOilPermeationZO,
    PumpElectricityZO,
)
from watertap.core.zero_order_costing import ZeroOrderCosting


def main():
    m = build()

    set_operating_conditions(m)
    assert_degrees_of_freedom(m, 0)
    assert_units_consistent(m)

    initialize_system(m)

    results = solve(m)
    assert_optimal_termination(results)
    display_results(m)

    add_costing(m)
    m.fs.costing.initialize()

    assert_degrees_of_freedom(m, 0)
    assert_units_consistent(m)
    results = solve(m)
    assert_optimal_termination(results)
    display_costing(m)

    return m, results


def build():
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.prop = prop_ZO.WaterParameterBlock(
        solute_list=[
            "oil",
        ]
    )

    # unit models
    # feed
    m.fs.feed = FeedZO(property_package=m.fs.prop)

    # pump
    m.fs.pump = PumpElectricityZO(property_package=m.fs.prop, database=m.db)

    # selective oil permeation
    m.fs.sop = SelectiveOilPermeationZO(property_package=m.fs.prop, database=m.db)

    # product streams
    m.fs.byproduct_oil = Product(property_package=m.fs.prop)
    m.fs.treated_water = Product(property_package=m.fs.prop)

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.pump.inlet)
    m.fs.s02 = Arc(source=m.fs.pump.outlet, destination=m.fs.sop.inlet)
    m.fs.s03 = Arc(source=m.fs.sop.byproduct, destination=m.fs.byproduct_oil.inlet)
    m.fs.s04 = Arc(source=m.fs.sop.treated, destination=m.fs.treated_water.inlet)

    # expand arcs
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):
    # ---specifications---
    # feed
    flow_vol = 3.8 * pyunits.L / pyunits.min
    conc_mass_oil = 500 * pyunits.mg / pyunits.liter

    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "oil"].fix(conc_mass_oil)
    solve(m.fs.feed)

    m.fs.pump.load_parameters_from_database(use_default_removal=True)
    m.fs.pump.lift_height.fix(2.86)  # head of 2.86 m (equivalent to 0.28 bar)

    m.fs.sop.load_parameters_from_database(use_default_removal=True)


def initialize_system(m):
    seq = SequentialDecomposition()
    seq.options.tear_set = []
    seq.options.iterLim = 1
    seq.run(m, lambda u: u.initialize())


def solve(blk, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results


def display_results(m):
    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++ DISPLAY RESULTS ++++++++++++++++++++")
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

    m.fs.sop.report()

    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


def add_costing(m):
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "selective_oil_permeation_case_study.yaml",
    )

    m.fs.costing = ZeroOrderCosting(case_study_definition=source_file)

    m.fs.pump.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.sop.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.costing.cost_process()

    m.fs.costing.add_electricity_intensity(m.fs.feed.properties[0].flow_vol)

    m.fs.costing.byproduct_oil_volume = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.byproduct_oil.properties[0].flow_vol,
            to_units=pyunits.bbl / m.fs.costing.base_period,
        ),
        doc="Volume of oil byproduct per year",
    )

    m.fs.costing.value_byproduct_oil = Expression(
        expr=pyunits.convert(
            m.fs.costing.byproduct_oil_volume * m.fs.costing.oil_cost,
            to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
        ),
        doc="Value of oil byproduct per year",
    )

    m.fs.costing.water_recovery_volume = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.treated_water.properties[0].flow_vol,
            to_units=pyunits.m**3 / m.fs.costing.base_period,
        ),
        doc="Volume of water recovered per year",
    )

    m.fs.costing.value_water_recovery = Expression(
        expr=pyunits.convert(
            m.fs.costing.water_recovery_volume * m.fs.costing.water_cost,
            to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
        ),
        doc="Value of water recovered per year",
    )

    m.fs.costing.LCOT = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
        )
        / (
            pyunits.convert(
                m.fs.feed.properties[0].flow_vol,
                to_units=pyunits.m**3 / m.fs.costing.base_period,
            )
            * m.fs.costing.utilization_factor
        ),
        doc="Levelized Cost of Treatment with respect to influent flowrate",
    )

    m.fs.costing.LCOT_with_revenue = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
            - m.fs.costing.value_byproduct_oil
            - m.fs.costing.value_water_recovery
        )
        / (
            pyunits.convert(
                m.fs.feed.properties[0].flow_vol,
                to_units=pyunits.m**3 / m.fs.costing.base_period,
            )
            * m.fs.costing.utilization_factor
        ),
        doc="Levelized Cost of Treatment with respect to influent flowrate including revenue from water recovered and oil byproduct",
    )


def display_costing(m):
    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++ DISPLAY COSTING ++++++++++++++++++++")
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

    print("\n---------- Levelized costs ----------")
    LCOT = value(
        pyunits.convert(m.fs.costing.LCOT, to_units=pyunits.USD_2020 / pyunits.m**3)
    )
    print(f"Levelized Cost of Treatment: {LCOT:.3f} $/m^3 feed water")

    LCOT_with_revenue = value(
        pyunits.convert(
            m.fs.costing.LCOT_with_revenue, to_units=pyunits.USD_2020 / pyunits.m**3
        )
    )
    print(
        f"Levelized Cost of Treatment with respect to influent flowrate including revenue: {LCOT_with_revenue:.3f} $/m^3 feed water"
    )

    print("\n------------- Capital costs -------------")

    total_capital_cost = value(
        pyunits.convert(m.fs.costing.total_capital_cost, to_units=pyunits.USD_2020)
    )
    print(f"Total Capital Costs: {total_capital_cost:.3f} $")

    print(
        f"Pump capital cost: {value(pyunits.convert(m.fs.pump.costing.capital_cost, to_units=pyunits.USD_2020)):.2f} $"
    )

    print(
        f"SOP capital cost: {value(pyunits.convert(m.fs.sop.costing.capital_cost, to_units=pyunits.USD_2020)):.2f} $"
    )

    print("\n------------- Operating costs -------------")

    electricity_intensity = value(
        pyunits.convert(
            m.fs.costing.electricity_intensity, to_units=pyunits.kWh / pyunits.m**3
        )
    )
    print(
        f"Electricity Intensity with respect to influent flowrate: {electricity_intensity:.3f} kWh/m^3"
    )

    total_operating_costs = value(
        pyunits.convert(
            m.fs.costing.total_operating_cost, to_units=pyunits.USD_2020 / pyunits.year
        )
    )
    print(f"Total operating costs: {total_operating_costs:.2f} $/year")

    fixed_operating_costs = value(
        pyunits.convert(
            m.fs.costing.total_fixed_operating_cost,
            to_units=pyunits.USD_2020 / pyunits.year,
        )
    )
    print(f"Fixed operating costs: {fixed_operating_costs:.2f} $/year")

    variable_operating_costs = value(
        pyunits.convert(
            m.fs.costing.total_variable_operating_cost,
            to_units=pyunits.USD_2020 / pyunits.year,
        )
    )
    print(f"Variable operating costs: {variable_operating_costs:.2f} $/year")

    print("\n------------- Revenue -------------")

    oil_revenue = value(
        pyunits.convert(
            m.fs.costing.value_byproduct_oil,
            to_units=pyunits.USD_2020 / pyunits.year,
        )
    )
    print(f"Oil revenue: {oil_revenue:.2f} $/year")

    water_revenue = value(
        pyunits.convert(
            m.fs.costing.value_water_recovery,
            to_units=pyunits.USD_2020 / pyunits.year,
        )
    )
    print(f"Water revenue: {water_revenue:.2f} $/year")

    total_revenue = value(
        pyunits.convert(
            m.fs.costing.value_byproduct_oil + m.fs.costing.value_water_recovery,
            to_units=pyunits.USD_2020 / pyunits.year,
        )
    )
    print(f"Total revenue: {total_revenue:.2f} $/year")

    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


if __name__ == "__main__":
    m, results = main()
