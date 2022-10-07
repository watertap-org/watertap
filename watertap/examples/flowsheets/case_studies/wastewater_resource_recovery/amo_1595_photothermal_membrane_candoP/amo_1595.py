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
    PhotothermalMembraneZO,
    CANDOPZO,
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
            "nitrogen",
            "phosphates",
            "bioconcentrated_phosphorous",
            "nitrous_oxide",
        ]
    )

    # unit models
    # feed
    m.fs.feed = FeedZO(property_package=m.fs.prop)

    # pump
    m.fs.pump = PumpElectricityZO(property_package=m.fs.prop, database=m.db)

    # photothermal membrane
    m.fs.photothermal = PhotothermalMembraneZO(
        property_package=m.fs.prop, database=m.db
    )

    # CANDO+P reactor
    # Note that CANDOPZO model electricity costing includes pumping costs for
    # the unit, so a pump has not been included in this flowsheet between the
    # photothermal membrane and CANDO+P units
    m.fs.candop = CANDOPZO(property_package=m.fs.prop, database=m.db)

    # product streams
    m.fs.photothermal_water = Product(property_package=m.fs.prop)
    m.fs.candop_treated = Product(property_package=m.fs.prop)
    m.fs.candop_byproduct = Product(property_package=m.fs.prop)

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.pump.inlet)
    m.fs.s02 = Arc(source=m.fs.pump.outlet, destination=m.fs.photothermal.inlet)
    m.fs.s03 = Arc(
        source=m.fs.photothermal.byproduct, destination=m.fs.photothermal_water.inlet
    )
    m.fs.s04 = Arc(source=m.fs.photothermal.treated, destination=m.fs.candop.inlet)
    m.fs.s05 = Arc(
        source=m.fs.candop.byproduct, destination=m.fs.candop_byproduct.inlet
    )
    m.fs.s06 = Arc(source=m.fs.candop.treated, destination=m.fs.candop_treated.inlet)

    # expand arcs
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):
    # ---specifications---
    # feed
    flow_vol = 800 * pyunits.m**3 / pyunits.hr
    conc_mass_nitrogen = 35 * pyunits.mg / pyunits.liter
    conc_mass_phosphates = 6 * pyunits.mg / pyunits.liter
    conc_mass_bcp = 1e-6 * pyunits.mg / pyunits.liter
    conc_mass_no2 = 1e-6 * pyunits.mg / pyunits.liter

    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "nitrogen"].fix(conc_mass_nitrogen)
    m.fs.feed.conc_mass_comp[0, "phosphates"].fix(conc_mass_phosphates)
    m.fs.feed.conc_mass_comp[0, "bioconcentrated_phosphorous"].fix(conc_mass_bcp)
    m.fs.feed.conc_mass_comp[0, "nitrous_oxide"].fix(conc_mass_no2)
    solve(m.fs.feed)

    # pump
    m.fs.pump.load_parameters_from_database(use_default_removal=True)
    m.fs.pump.lift_height.fix(4.2)  # head of 4.2 m

    # photothermal membrane
    m.fs.photothermal.load_parameters_from_database(use_default_removal=True)

    # CANDO+P reactor
    m.fs.candop.load_parameters_from_database(use_default_removal=True)


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

    unit_list = ["pump", "photothermal", "candop"]
    for u in unit_list:
        m.fs.component(u).report()

    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


def add_costing(m):
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "amo_1595_case_study.yaml",
    )

    m.fs.costing = ZeroOrderCosting(case_study_definition=source_file)

    costing_kwargs = {"flowsheet_costing_block": m.fs.costing}
    m.fs.pump.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.photothermal.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.candop.costing = UnitModelCostingBlock(**costing_kwargs)

    m.fs.costing.cost_process()

    m.fs.costing.add_electricity_intensity(m.fs.feed.properties[0].flow_vol)

    m.fs.costing.annual_water_inlet = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.feed.properties[0].flow_vol,
            to_units=pyunits.m**3 / m.fs.costing.base_period,
        )
    )

    m.fs.costing.bcp_recovery_mass = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.candop_treated.flow_mass_comp[0, "bioconcentrated_phosphorous"],
            to_units=pyunits.kg / m.fs.costing.base_period,
        ),
        doc="Mass of BCP generated per year",
    )

    m.fs.costing.value_bcp_recovery = Expression(
        expr=pyunits.convert(
            m.fs.costing.bcp_recovery_mass * m.fs.costing.bcp_cost,
            to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
        ),
        doc="Value of BCP generated per year",
    )

    m.fs.costing.water_byproduct_volume = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.photothermal_water.properties[0].flow_vol,
            to_units=pyunits.m**3 / m.fs.costing.base_period,
        ),
        doc="Volume of byproduct water generated per year",
    )

    m.fs.costing.value_water_byproduct = Expression(
        expr=pyunits.convert(
            m.fs.costing.water_byproduct_volume * m.fs.costing.water_cost,
            to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
        ),
        doc="Value of byproduct water generated per year",
    )

    m.fs.costing.N2O_byproduct_mass = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.candop_byproduct.flow_mass_comp[0, "nitrous_oxide"],
            to_units=pyunits.kg / m.fs.costing.base_period,
        ),
        doc="Mass of byproduct nitrous oxide generated per year",
    )

    m.fs.costing.value_N2O_byproduct = Expression(
        expr=pyunits.convert(
            m.fs.costing.N2O_byproduct_mass * m.fs.costing.nitrous_oxide_cost,
            to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
        ),
        doc="Value of byproduct nitrous oxide generated per year",
    )

    m.fs.costing.LCOW = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
        )
        / (
            pyunits.convert(
                m.fs.photothermal_water.properties[0].flow_vol,
                to_units=pyunits.m**3 / m.fs.costing.base_period,
            )
            * m.fs.costing.utilization_factor
        ),
        doc="Levelized Cost of Water",
    )

    m.fs.costing.LCOW_with_revenue = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
            - m.fs.costing.value_bcp_recovery
            - m.fs.costing.value_N2O_byproduct
        )
        / (
            pyunits.convert(
                m.fs.photothermal_water.properties[0].flow_vol,
                to_units=pyunits.m**3 / m.fs.costing.base_period,
            )
            * m.fs.costing.utilization_factor
        ),
        doc="Levelized Cost of Water including revenue from BCP and N20",
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
            - m.fs.costing.value_bcp_recovery
            - m.fs.costing.value_N2O_byproduct
            - m.fs.costing.value_water_byproduct
        )
        / (
            pyunits.convert(
                m.fs.feed.properties[0].flow_vol,
                to_units=pyunits.m**3 / m.fs.costing.base_period,
            )
            * m.fs.costing.utilization_factor
        ),
        doc="Levelized Cost of Treatment with respect to influent flowrate including revenue from BCP, N2O, and byproduct water",
    )

    m.fs.costing.LC_BCP = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
        )
        / (
            pyunits.convert(
                m.fs.candop_treated.properties[0].flow_mass_comp[
                    "bioconcentrated_phosphorous"
                ],
                to_units=pyunits.kg / m.fs.costing.base_period,
            )
            * m.fs.costing.utilization_factor
        ),
        doc="Levelized Cost of BCP (not accounting for any revenue)",
    )

    m.fs.costing.LC_N2O = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
        )
        / (
            pyunits.convert(
                m.fs.candop_byproduct.properties[0].flow_mass_comp["nitrous_oxide"],
                to_units=pyunits.kg / m.fs.costing.base_period,
            )
            * m.fs.costing.utilization_factor
        ),
        doc="Levelized Cost of N2O (not accounting for any revenue)",
    )


def display_costing(m):
    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++ DISPLAY COSTING ++++++++++++++++++++")
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

    print("\n---------- Levelized costs ----------")
    LCOW = value(
        pyunits.convert(m.fs.costing.LCOW, to_units=pyunits.USD_2020 / pyunits.m**3)
    )
    print(f"Levelized Cost of Water: {LCOW:.3f} $/m^3 product water")

    LCOW_with_revenue = value(
        pyunits.convert(
            m.fs.costing.LCOW_with_revenue, to_units=pyunits.USD_2020 / pyunits.m**3
        )
    )
    print(
        f"Levelized Cost of Water including revenue: {LCOW_with_revenue:.3f} $/m^3 product water"
    )

    LCOT = value(
        pyunits.convert(m.fs.costing.LCOT, to_units=pyunits.USD_2020 / pyunits.m**3)
    )
    print(
        f"Levelized Cost of Treatment with respect to influent flowrate: {LCOT:.3f} $/m^3 feed water"
    )

    LCOT_with_revenue = value(
        pyunits.convert(
            m.fs.costing.LCOT_with_revenue, to_units=pyunits.USD_2020 / pyunits.m**3
        )
    )
    print(
        f"Levelized Cost of Treatment with respect to influent flowrate including revenue: {LCOT_with_revenue:.3f} $/m^3 feed water"
    )

    LC_BCP = value(
        pyunits.convert(m.fs.costing.LC_BCP, to_units=pyunits.USD_2020 / pyunits.kg)
    )
    print(f"Levelized Cost of BCP (not including any revenue): {LC_BCP:.2f} $/kg")

    LC_N2O = value(
        pyunits.convert(m.fs.costing.LC_N2O, to_units=pyunits.USD_2020 / pyunits.kg)
    )
    print(f"Levelized Cost of N2O (not including any revenue): {LC_N2O:.2f} $/kg")

    print("\n------------- Capital costs -------------")
    DCC_normalized = value(
        pyunits.convert(
            (
                m.fs.pump.costing.direct_capital_cost
                + m.fs.photothermal.costing.direct_capital_cost
                + m.fs.candop.costing.direct_capital_cost
            )
            / m.fs.feed.properties[0].flow_vol,
            to_units=pyunits.kUSD_2020 / (pyunits.m**3 / pyunits.hr),
        )
    )
    print(f"Normalized direct capital costs: {DCC_normalized:.2f} k$/(m^3/hr)")

    ICC_normalized = value(
        pyunits.convert(
            m.fs.costing.total_capital_cost / m.fs.feed.properties[0].flow_vol,
            to_units=pyunits.kUSD_2020 / (pyunits.m**3 / pyunits.hr),
        )
    )
    print(f"Normalized total capital costs: {ICC_normalized:.2f} k$/(m^3/hr)")

    total_capital_cost = value(
        pyunits.convert(m.fs.costing.total_capital_cost, to_units=pyunits.MUSD_2020)
    )
    print(f"Total Capital Costs: {total_capital_cost:.3f} M$")

    print(
        f"Pump capital cost: {value(pyunits.convert(m.fs.pump.costing.capital_cost, to_units=pyunits.kUSD_2020)):.2f} k$"
    )

    print(
        f"Photothermal membrane capital cost: {value(pyunits.convert(m.fs.photothermal.costing.capital_cost, to_units=pyunits.MUSD_2020)):.2f} M$"
    )

    print(
        f"CANDO+P capital cost: {value(pyunits.convert(m.fs.candop.costing.capital_cost, to_units=pyunits.MUSD_2020)):.2f} M$"
    )

    print("\n------------- Operating costs -------------")
    FMC_normalized = value(
        pyunits.convert(
            m.fs.costing.maintenance_cost / m.fs.costing.total_capital_cost,
            to_units=1 / pyunits.year,
        )
    )
    print(f"Normalized maintenance costs: {FMC_normalized:.2f} 1/year")

    OFOC_normalized = value(
        pyunits.convert(
            m.fs.costing.aggregate_fixed_operating_cost
            / m.fs.costing.total_capital_cost,
            to_units=1 / pyunits.year,
        )
    )
    print(f"Normalized other fixed operating cost: {OFOC_normalized:.2f} 1/year")

    EC_normalized = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_costs["electricity"]
            * m.fs.costing.utilization_factor
            / m.fs.costing.annual_water_inlet,
            to_units=pyunits.USD_2020 / pyunits.m**3,
        )
    )
    print(f"Normalized electricity cost: {EC_normalized:.3f} $/m^3 of feed")

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
            m.fs.costing.total_operating_cost, to_units=pyunits.kUSD_2020 / pyunits.year
        )
    )
    print(f"Total operating costs: {total_operating_costs:.2f} k$/year")

    fixed_operating_costs = value(
        pyunits.convert(
            m.fs.costing.total_fixed_operating_cost,
            to_units=pyunits.kUSD_2020 / pyunits.year,
        )
    )
    print(f"Fixed operating costs: {fixed_operating_costs:.2f} k$/year")

    variable_operating_costs = value(
        pyunits.convert(
            m.fs.costing.total_variable_operating_cost,
            to_units=pyunits.kUSD_2020 / pyunits.year,
        )
    )
    print(f"Variable operating costs: {variable_operating_costs:.2f} k$/year")

    electricity_operating_costs = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_costs["electricity"]
            * m.fs.costing.utilization_factor,
            to_units=pyunits.kUSD_2020 / pyunits.year,
        )
    )
    print(f"Electricity operating costs: {electricity_operating_costs:.2f} k$/year")

    print("\n------------- Revenue -------------")

    BCP_revenue = value(
        pyunits.convert(
            m.fs.costing.value_bcp_recovery,
            to_units=pyunits.kUSD_2020 / pyunits.year,
        )
    )
    print(f"BCP revenue: {BCP_revenue:.2f} k$/year")

    N2O_revenue = value(
        pyunits.convert(
            m.fs.costing.value_N2O_byproduct,
            to_units=pyunits.kUSD_2020 / pyunits.year,
        )
    )
    print(f"N2O revenue: {N2O_revenue:.2f} k$/year")

    water_revenue = value(
        pyunits.convert(
            m.fs.costing.value_water_byproduct,
            to_units=pyunits.kUSD_2020 / pyunits.year,
        )
    )
    print(f"Water revenue: {water_revenue:.2f} k$/year")

    total_revenue = value(
        pyunits.convert(
            m.fs.costing.value_bcp_recovery
            + m.fs.costing.value_N2O_byproduct
            + m.fs.costing.value_water_byproduct,
            to_units=pyunits.kUSD_2020 / pyunits.year,
        )
    )
    print(f"Total revenue: {total_revenue:.2f} k$/year")

    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


if __name__ == "__main__":
    m, results = main()
