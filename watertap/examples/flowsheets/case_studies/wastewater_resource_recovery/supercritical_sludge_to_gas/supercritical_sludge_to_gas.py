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
from idaes.core.util import get_solver
from idaes.generic_models.unit_models import Product
import idaes.core.util.scaling as iscale
from idaes.generic_models.costing import UnitModelCostingBlock

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
    display_results(m)

    add_costing(m)
    initialize_costing(m)
    assert_degrees_of_freedom(m, 0)
    assert_units_consistent(m)

    results = solve(m)
    display_costing(m)

    return m, results


def build():
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.prop = prop_ZO.WaterParameterBlock(
        default={
            "solute_list": [
                "organic_solid",
                "organic_liquid",
                "inorganic_solid",
                "carbon_dioxide",
            ]
        }
    )

    # unit models
    m.fs.feed = FeedZO(default={"property_package": m.fs.prop})
    m.fs.ATHTL = ATHTLZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
        },
    )
    m.fs.salt_precipitation = SaltPrecipitationZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
        },
    )
    m.fs.HTG = HTGZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
        },
    )

    m.fs.product_H2O = Product(default={"property_package": m.fs.prop})
    # CO2 from ATHTL to sulfer conversion unit, final product should be solid sulfer, H2 and CO2
    m.fs.product_CO2 = Product(default={"property_package": m.fs.prop})
    # sulfates, nitrates & phosphates (organics & inorganics) from salt precipitation unit
    m.fs.product_salts = Product(default={"property_package": m.fs.prop})
    # CO2 from HTG for renewable natural gas product
    m.fs.product_natural_gas = Product(default={"property_package": m.fs.prop})

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
    m.fs.costing = ZeroOrderCosting(default={"case_study_definition": source_file})
    # typing aid
    costing_kwargs = {"default": {"flowsheet_costing_block": m.fs.costing}}
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


if __name__ == "__main__":
    m, results = main()
