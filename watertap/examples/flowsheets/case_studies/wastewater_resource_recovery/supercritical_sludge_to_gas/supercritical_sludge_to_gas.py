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
    Set,
    Expression,
    value,
    TransformationFactory,
    units as pyunits,
    assert_optimal_termination,
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
    ATHTLZO,
    SaltPrecipitationZO,
    HTGZO,
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
    # display_reports(m)

    add_costing(m)
    initialize_costing(m)
    assert_degrees_of_freedom(m, 0)
    assert_units_consistent(m)

    results = solve(m)
    assert_optimal_termination(results)

    display_metrics_results(m)
    display_additional_results(m)

    return m, results


def build():
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.prop = prop_ZO.WaterParameterBlock(
        solute_list=[
            "organic_solid",
            "organic_liquid",
            "inorganic_solid",
            "carbon_dioxide",
        ]
    )

    # unit models
    m.fs.feed = FeedZO(property_package=m.fs.prop)
    m.fs.ATHTL = ATHTLZO(property_package=m.fs.prop, database=m.db)
    m.fs.salt_precipitation = SaltPrecipitationZO(
        property_package=m.fs.prop, database=m.db
    )
    m.fs.HTG = HTGZO(property_package=m.fs.prop, database=m.db)

    m.fs.product_H2O = Product(property_package=m.fs.prop)
    # CO2 from ATHTL to sulfer conversion unit, final product should be solid sulfer, H2 and CO2
    m.fs.product_CO2 = Product(property_package=m.fs.prop)
    # sulfates, nitrates & phosphates (organics & inorganics) from salt precipitation unit
    m.fs.product_salts = Product(property_package=m.fs.prop)
    # CO2 from HTG for renewable natural gas product
    m.fs.product_natural_gas = Product(property_package=m.fs.prop)

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.ATHTL.inlet)
    m.fs.s02 = Arc(source=m.fs.ATHTL.treated, destination=m.fs.salt_precipitation.inlet)
    m.fs.s03 = Arc(source=m.fs.salt_precipitation.treated, destination=m.fs.HTG.inlet)
    m.fs.s04 = Arc(source=m.fs.HTG.treated, destination=m.fs.product_H2O.inlet)
    m.fs.s05 = Arc(
        source=m.fs.HTG.byproduct, destination=m.fs.product_natural_gas.inlet
    )
    m.fs.s06 = Arc(source=m.fs.ATHTL.byproduct, destination=m.fs.product_CO2.inlet)
    m.fs.s07 = Arc(
        source=m.fs.salt_precipitation.byproduct, destination=m.fs.product_salts.inlet
    )
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):
    # ---specifications---
    # feed
    flow_mass_H2O = 400 * pyunits.metric_ton / pyunits.day
    flow_mass_inorganic_solid = 28.8 * pyunits.metric_ton / pyunits.day
    flow_mass_organic_solid = 71.2 * pyunits.metric_ton / pyunits.day

    m.fs.feed.flow_mass_comp[0, "H2O"].fix(flow_mass_H2O)
    m.fs.feed.flow_mass_comp[0, "inorganic_solid"].fix(flow_mass_inorganic_solid)
    m.fs.feed.flow_mass_comp[0, "organic_solid"].fix(flow_mass_organic_solid)
    m.fs.feed.conc_mass_comp[0, "organic_liquid"].fix(1e-5)
    m.fs.feed.conc_mass_comp[0, "carbon_dioxide"].fix(1e-5)
    solve(m.fs.feed)

    # ATHTL
    m.fs.ATHTL.load_parameters_from_database(use_default_removal=True)

    # Salt precipitation
    m.fs.salt_precipitation.load_parameters_from_database(use_default_removal=True)

    # HTG
    m.fs.HTG.load_parameters_from_database(use_default_removal=True)


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
    unit_list = ["feed", "ATHTL", "salt_precipitation", "HTG"]
    for u in unit_list:
        m.fs.component(u).report()


def add_costing(m):
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "supercritical_sludge_to_gas_global_costing.yaml",
    )
    m.fs.costing = ZeroOrderCosting(case_study_definition=source_file)
    # typing aid
    costing_kwargs = {"flowsheet_costing_block": m.fs.costing}
    m.fs.ATHTL.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.salt_precipitation.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.HTG.costing = UnitModelCostingBlock(**costing_kwargs)

    m.fs.costing.cost_process()
    m.fs.costing.add_electricity_intensity(m.fs.product_H2O.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product_H2O.properties[0].flow_vol)

    m.fs.costing.LCOG = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
        )
        / (
            m.fs.costing.utilization_factor
            * pyunits.convert(
                m.fs.product_natural_gas.properties[0].flow_mass_comp["carbon_dioxide"],
                to_units=pyunits.kg / m.fs.costing.base_period,
            )
        ),
        doc="Levelized cost of CO2 for production of renewable natural gas",
    )

    m.fs.costing.LCOC = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
        )
        / (
            m.fs.costing.utilization_factor
            * pyunits.convert(
                m.fs.product_CO2.properties[0].flow_mass_comp["carbon_dioxide"],
                to_units=pyunits.kg / m.fs.costing.base_period,
            )
        ),
        doc="Levelized cost of CO2 from AT-HTL",
    )

    m.fs.costing.LCOS = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
        )
        / (
            m.fs.costing.utilization_factor
            * pyunits.convert(
                (
                    m.fs.product_salts.properties[0].flow_mass_comp["organic_solid"]
                    + m.fs.product_salts.properties[0].flow_mass_comp["inorganic_solid"]
                ),
                to_units=pyunits.kg / m.fs.costing.base_period,
            )
        ),
        doc="Levelized cost of Sulfates, Nitrates and Phosphates",
    )

    # other levelized costs
    m.fs.costing.annual_water_inlet = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.feed.properties[0].flow_vol,
            to_units=pyunits.m**3 / m.fs.costing.base_period,
        )
    )

    m.fs.costing.annual_water_production = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.product_H2O.properties[0].flow_vol,
            to_units=pyunits.m**3 / m.fs.costing.base_period,
        )
    )

    m.fs.costing.annual_CO2_renewable_production = Expression(
        expr=(
            m.fs.costing.utilization_factor
            * pyunits.convert(
                m.fs.HTG.byproduct.flow_mass_comp[0, "carbon_dioxide"],
                to_units=pyunits.kg / m.fs.costing.base_period,
            )
        )
    )

    m.fs.costing.annual_CO2_sulfur_production = Expression(
        expr=(
            m.fs.costing.utilization_factor
            * pyunits.convert(
                m.fs.ATHTL.byproduct.flow_mass_comp[0, "carbon_dioxide"],
                to_units=pyunits.kg / m.fs.costing.base_period,
            )
        )
    )

    m.fs.costing.annual_sulphate_nitrate_phosphate_production = Expression(
        expr=(
            m.fs.costing.utilization_factor
            * pyunits.convert(
                (
                    m.fs.salt_precipitation.byproduct.flow_mass_comp[0, "organic_solid"]
                    + m.fs.salt_precipitation.byproduct.flow_mass_comp[
                        0, "inorganic_solid"
                    ]
                ),
                to_units=pyunits.kg / m.fs.costing.base_period,
            )
        )
    )

    m.fs.costing.total_annualized_cost = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
        )
    )

    m.fs.costing.LCOT = Expression(
        expr=(m.fs.costing.total_annualized_cost / m.fs.costing.annual_water_inlet),
        doc="Levelized Cost of Treatment",
    )


def initialize_costing(m):
    m.fs.costing.initialize()


def display_costing(m):
    m.fs.costing.total_capital_cost.display()
    m.fs.costing.total_operating_cost.display()
    m.fs.costing.LCOW.display()

    print("\nUnit Capital Costs\n")
    for u in m.fs.costing._registered_unit_costing:
        print(
            u.name,
            " :   ",
            value(pyunits.convert(u.capital_cost, to_units=pyunits.USD_2020)),
        )

    print("\nUtility Costs\n")
    for f in m.fs.costing.flow_types:
        print(
            f,
            " :   ",
            value(
                pyunits.convert(
                    m.fs.costing.aggregate_flow_costs[f],
                    to_units=pyunits.USD_2020 / pyunits.year,
                )
            ),
        )

    print("")
    total_capital_cost = value(
        pyunits.convert(m.fs.costing.total_capital_cost, to_units=pyunits.MUSD_2020)
    )
    print(f"Total Capital Costs: {total_capital_cost:.2f} M$")
    total_operating_cost = value(
        pyunits.convert(
            m.fs.costing.total_operating_cost, to_units=pyunits.MUSD_2020 / pyunits.year
        )
    )
    print(f"Total Operating Costs: {total_operating_cost:.2f} M$/year")
    electricity_intensity = value(
        pyunits.convert(
            m.fs.costing.electricity_intensity, to_units=pyunits.kWh / pyunits.m**3
        )
    )
    print(f"Electricity Intensity: {electricity_intensity:.4f} kWh/m^3")
    LCOW = value(
        pyunits.convert(m.fs.costing.LCOW, to_units=pyunits.USD_2020 / pyunits.m**3)
    )
    print(f"Levelized Cost of Water: {LCOW:.4f} $/m^3")
    LCOG = value(
        pyunits.convert(m.fs.costing.LCOG, to_units=pyunits.USD_2020 / pyunits.kg)
    )
    print(f"Levelized Cost of Natural Gas Based on Carbon Dioxide: {LCOG:.3f} $/kg")
    LCOC = value(
        pyunits.convert(m.fs.costing.LCOC, to_units=pyunits.USD_2020 / pyunits.kg)
    )
    print(f"Levelized Cost of Carbon Dioxide from AT-HTL: {LCOC:.3f} $/kg")
    LCOS = value(
        pyunits.convert(m.fs.costing.LCOS, to_units=pyunits.USD_2020 / pyunits.kg)
    )
    print(f"Levelized Cost of Sulfates, Nitrates and Phosphates: {LCOS:.3f} $/kg")


def display_metrics_results(m):
    print("----------Levelized costs----------")
    LCOT = value(
        pyunits.convert(
            m.fs.costing.LCOT, to_units=m.fs.costing.base_currency / pyunits.m**3
        )
    )
    print(f"Levelized Cost of Treatment: {LCOT:.4f} $/m3 of feed")
    LCOW = value(
        pyunits.convert(
            m.fs.costing.LCOW, to_units=m.fs.costing.base_currency / pyunits.m**3
        )
    )
    print(f"Levelized Cost of Water: {LCOW:.4f} $/m3 of product")
    LCOG = value(
        pyunits.convert(m.fs.costing.LCOG, to_units=pyunits.USD_2020 / pyunits.kg)
    )
    print(f"Levelized Cost of Natural Gas Based on Carbon Dioxide: {LCOG:.4f} $/kg")
    LCOC = value(
        pyunits.convert(m.fs.costing.LCOC, to_units=pyunits.USD_2020 / pyunits.kg)
    )
    print(f"Levelized Cost of Carbon Dioxide from AT-HTL: {LCOC:.4f} $/kg")
    LCOS = value(
        pyunits.convert(m.fs.costing.LCOS, to_units=pyunits.USD_2020 / pyunits.kg)
    )
    print(f"Levelized Cost of Sulfates, Nitrates and Phosphates: {LCOS:.4f} $/kg")

    print("----------Capital costs----------")
    DCC_normalized = value(
        pyunits.convert(
            (
                m.fs.ATHTL.costing.direct_capital_cost
                + m.fs.salt_precipitation.costing.direct_capital_cost
                + m.fs.HTG.costing.direct_capital_cost
            )
            / m.fs.feed.properties[0].flow_vol,
            to_units=m.fs.costing.base_currency / (pyunits.m**3 / pyunits.day),
        )
    )
    print(f"Normalized direct capital costs: {DCC_normalized:.4f} $/(m3/day)")
    ICC_normalized = value(
        pyunits.convert(
            m.fs.costing.total_capital_cost / m.fs.feed.properties[0].flow_vol,
            to_units=m.fs.costing.base_currency / (pyunits.m**3 / pyunits.day),
        )
    )
    print(f"Normalized total capital costs: {ICC_normalized:.4f} $/(m3/day)")

    print("----------Operating costs----------")
    FMC_normalized = value(
        pyunits.convert(
            m.fs.costing.maintenance_cost / m.fs.costing.total_capital_cost,
            to_units=1 / pyunits.a,
        )
    )
    print(f"Normalized maintenance costs: {FMC_normalized:.4f} 1/year")
    EC_normalized = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_costs["electricity"]
            / m.fs.costing.annual_water_inlet,
            to_units=m.fs.costing.base_currency / pyunits.m**3,
        )
    )
    print(f"Normalized electricity cost: {EC_normalized:.4f} $/m3 of feed")
    Catalyst_ATHTL_normalized = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_costs["catalyst_ATHTL"]
            / m.fs.costing.annual_water_inlet,
            to_units=m.fs.costing.base_currency / pyunits.m**3,
        )
    )
    print(
        f"Normalized catalyst for AT-HTL cost: {Catalyst_ATHTL_normalized:.4f} $/m3 of feed"
    )
    Catalyst_HTG_normalized = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_costs["catalyst_HTG"]
            / m.fs.costing.annual_water_inlet,
            to_units=m.fs.costing.base_currency / pyunits.m**3,
        )
    )
    print(
        f"Normalized catalyst for HTG cost: {Catalyst_HTG_normalized:.4f} $/m3 of feed"
    )

    print("----------Performance metrics----------")
    volumetric_recovery = value(
        m.fs.product_H2O.properties[0].flow_vol / m.fs.feed.properties[0].flow_vol
    )
    print(f"Water recovery: {volumetric_recovery:.4f} m3 of product/m3 of feed")
    OrganicR_normalized = value(
        pyunits.convert(
            1
            - (
                m.fs.product_H2O.properties[0].flow_mass_comp["organic_solid"]
                + m.fs.product_H2O.properties[0].flow_mass_comp["organic_liquid"]
            )
            / (
                m.fs.feed.properties[0].flow_mass_comp["organic_solid"]
                + m.fs.feed.properties[0].flow_mass_comp["organic_liquid"]
            ),
            to_units=pyunits.dimensionless,
        )
    )
    print(f"Organic removal: {OrganicR_normalized:.4f} dimensionless")
    InorganicR_normalized = value(
        pyunits.convert(
            1
            - m.fs.product_H2O.properties[0].flow_mass_comp["inorganic_solid"]
            / m.fs.feed.properties[0].flow_mass_comp["inorganic_solid"],
            to_units=pyunits.dimensionless,
        )
    )
    print(f"Inorganic removal: {InorganicR_normalized:.4f} dimensionless")

    print("----------Energy intensity----------")
    SEC = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_electricity / m.fs.feed.properties[0].flow_vol,
            to_units=pyunits.kWh / pyunits.m**3,
        )
    )
    print(f"Specific electricity consumption: {SEC:.4f} kWh/m3 of feed")


def display_additional_results(m):
    print("----------Outlets----------")
    product_H2O_flow = value(
        pyunits.convert(
            m.fs.product_H2O.properties[0].flow_vol,
            to_units=pyunits.m**3 / pyunits.hr,
        )
    )
    print(f"H2O outlet flow: {product_H2O_flow:.4f} m3/h")
    product_H2O_Organic_Solid = value(
        pyunits.convert(
            m.fs.product_H2O.properties[0].conc_mass_comp["organic_solid"],
            to_units=pyunits.g / pyunits.L,
        )
    )
    print(f"H2O outlet organic(s) conc: {product_H2O_Organic_Solid:.4f} g/L")
    product_H2O_Organic_Liquid = value(
        pyunits.convert(
            m.fs.product_H2O.properties[0].conc_mass_comp["organic_liquid"],
            to_units=pyunits.g / pyunits.L,
        )
    )
    print(f"H2O outlet Organic(l) conc: {product_H2O_Organic_Liquid:.4f} g/L")
    product_H2O_Inorganic_Solid = value(
        pyunits.convert(
            m.fs.product_H2O.properties[0].conc_mass_comp["inorganic_solid"],
            to_units=pyunits.g / pyunits.L,
        )
    )
    print(f"H2O outlet inorganic(s) conc: {product_H2O_Inorganic_Solid:.4f} g/L")

    print("----------Capital costs----------")
    total_capital_costs = value(m.fs.costing.total_capital_cost) / 1e6
    print(f"Total capital costs: {total_capital_costs:.4f} $M")
    ATHTL_capital_costs = value(m.fs.ATHTL.costing.capital_cost) / 1e6
    print(f"AT-HTL capital costs: {ATHTL_capital_costs:.4f} $M")
    salt_precipitation_capital_costs = (
        value(m.fs.salt_precipitation.costing.capital_cost) / 1e6
    )
    print(
        f"Salt precipitation capital costs: {salt_precipitation_capital_costs:.4f} $M"
    )
    HTG_capital_costs = value(m.fs.HTG.costing.capital_cost) / 1e6
    print(f"HTG capital costs: {HTG_capital_costs:.4f} $M")

    print("----------Operating costs----------")
    total_operating_costs = value(m.fs.costing.total_operating_cost) / 1e6
    print(f"Total operating costs: {total_operating_costs:.4f} $M/year")
    fixed_operating_costs = value(m.fs.costing.total_fixed_operating_cost) / 1e6
    print(f"Fixed operating costs: {fixed_operating_costs:.4f} $M/year")
    electricity_operating_costs = (
        value(m.fs.costing.aggregate_flow_costs["electricity"]) / 1e3
    )
    print(f"Electricity operating costs: {electricity_operating_costs:.4f} $k/year")
    catalyst_ATHTL_operating_costs = (
        value(m.fs.costing.aggregate_flow_costs["catalyst_ATHTL"]) / 1e6
    )
    print(
        f"Catalyst in AT-HTL operating costs: {catalyst_ATHTL_operating_costs:.4f} $M/year"
    )
    catalyst_HTG_operating_costs = (
        value(m.fs.costing.aggregate_flow_costs["catalyst_HTG"]) / 1e6
    )
    print(
        f"Catalyst in HTG operating costs: {catalyst_HTG_operating_costs:.4f} $M/year"
    )


if __name__ == "__main__":
    m, results = main()
