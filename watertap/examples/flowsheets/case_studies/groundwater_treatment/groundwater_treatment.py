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
from idaes.core.util import get_solver
from idaes.generic_models.unit_models import Product
import idaes.core.util.scaling as iscale
from idaes.generic_models.costing import UnitModelCostingBlock

from watertap.core.util.initialization import assert_degrees_of_freedom

from watertap.core.wt_database import Database
import watertap.core.zero_order_properties as prop_ZO
from watertap.unit_models.zero_order import (
    FeedZO,
    MicrobialBatteryZO,
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

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.prop = prop_ZO.WaterParameterBlock(
        default={
            "solute_list": [
                "arsenic",
                "uranium",
                "nitrate",
                "phosphates",
                "iron",
                "filtration_media",
            ]
        }
    )

    # unit models
    # feed
    m.fs.feed = FeedZO(default={"property_package": m.fs.prop})

    # pump
    m.fs.pump = PumpElectricityZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
        },
    )

    # microbial battery
    m.fs.micbatt = MicrobialBatteryZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
        },
    )

    # product streams
    m.fs.filtered_water = Product(default={"property_package": m.fs.prop})
    m.fs.byproduct = Product(default={"property_package": m.fs.prop})

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.pump.inlet)
    m.fs.s02 = Arc(source=m.fs.pump.outlet, destination=m.fs.micbatt.inlet)
    m.fs.s03 = Arc(source=m.fs.micbatt.treated, destination=m.fs.filtered_water.inlet)
    m.fs.s04 = Arc(source=m.fs.micbatt.byproduct, destination=m.fs.byproduct.inlet)

    # expand arcs
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    iscale.calculate_scaling_factors(m)

    return


def set_operating_conditions(m):
    # ---specifications---
    # feed
    flow_vol = 1 * pyunits.m**3 / pyunits.day
    conc_mass_arsenic = 0.04 * pyunits.mg / pyunits.liter
    conc_mass_uranium = 0.06 * pyunits.mg / pyunits.liter
    conc_mass_nitrate = 10 * pyunits.mg / pyunits.liter
    conc_mass_phosphates = 0.1 * pyunits.mg / pyunits.liter
    conc_mass_iron = 0.5 * pyunits.mg / pyunits.liter
    conc_mass_filtration_media = 0.1 * pyunits.mg / pyunits.liter

    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "arsenic"].fix(conc_mass_arsenic)
    m.fs.feed.conc_mass_comp[0, "uranium"].fix(conc_mass_uranium)
    m.fs.feed.conc_mass_comp[0, "nitrate"].fix(conc_mass_nitrate)
    m.fs.feed.conc_mass_comp[0, "phosphates"].fix(conc_mass_phosphates)
    m.fs.feed.conc_mass_comp[0, "iron"].fix(conc_mass_iron)
    m.fs.feed.conc_mass_comp[0, "filtration_media"].fix(conc_mass_filtration_media)

    # pump
    m.fs.pump.load_parameters_from_database(use_default_removal=True)
    m.fs.pump.lift_height.fix(30)  # head of 30 m to pump groundwater

    # microbial battery
    m.fs.micbatt.load_parameters_from_database(use_default_removal=True)


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

    unit_list = ["feed", "pump", "micbatt"]
    for u in unit_list:
        m.fs.component(u).report()

    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


def add_costing(m):
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "groundwater_treatment_case_study.yaml",
    )

    m.fs.costing = ZeroOrderCosting(default={"case_study_definition": source_file})

    costing_kwargs = {"default": {"flowsheet_costing_block": m.fs.costing}}
    m.fs.pump.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.micbatt.costing = UnitModelCostingBlock(**costing_kwargs)

    m.fs.costing.cost_process()

    m.fs.costing.add_electricity_intensity(m.fs.feed.properties[0].flow_vol)

    m.fs.costing.spent_media_recovery_mass = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.byproduct.flow_mass_comp[0, "filtration_media"],
            to_units=pyunits.kg / m.fs.costing.base_period,
        ),
        doc="Mass of spent filtration media generated per year (accounting for utilization factor)",
    )

    m.fs.costing.cost_spent_media_disposal = Expression(
        expr=pyunits.convert(
            m.fs.costing.spent_media_recovery_mass
            * m.fs.costing.filtration_media_disposal_cost,
            to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
        ),
        doc="Disposal cost of spent filtration media per year (accounting for utilization factor)",
    )

    m.fs.costing.fresh_media_mass = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.feed.flow_mass_comp[0, "filtration_media"],
            to_units=pyunits.kg / m.fs.costing.base_period,
        ),
        doc="Mass of fresh filtration media used per year (accounting for utilization factor)",
    )

    m.fs.costing.cost_fresh_media = Expression(
        expr=pyunits.convert(
            m.fs.costing.fresh_media_mass * m.fs.costing.filtration_media_cost,
            to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
        ),
        doc="Cost of fresh filtration media used per year (accounting for utilization factor)",
    )

    m.fs.costing.water_product_volume = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.filtered_water.properties[0].flow_vol,
            to_units=pyunits.m**3 / m.fs.costing.base_period,
        ),
        doc="Volume of product water generated per year (accounting for utilization factor)",
    )

    m.fs.costing.value_water_product = Expression(
        expr=pyunits.convert(
            m.fs.costing.water_product_volume * m.fs.costing.water_cost,
            to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
        ),
        doc="Dollar value of product water generated per year (accounting for utilization factor)",
    )

    m.fs.costing.LCOW = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
            + m.fs.costing.cost_spent_media_disposal
            + m.fs.costing.cost_fresh_media
        )
        / (
            pyunits.convert(
                m.fs.filtered_water.properties[0].flow_vol,
                to_units=pyunits.m**3 / m.fs.costing.base_period,
            )
            * m.fs.costing.utilization_factor
        ),
        doc="Levelized Cost of Water when accounting for seed media acquisition and spent media disposal costs",
    )

    m.fs.costing.LCOT = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
            + m.fs.costing.cost_spent_media_disposal
            + m.fs.costing.cost_fresh_media
            - m.fs.costing.value_water_product
        )
        / (
            pyunits.convert(
                m.fs.feed.properties[0].flow_vol,
                to_units=pyunits.m**3 / m.fs.costing.base_period,
            )
            * m.fs.costing.utilization_factor
        ),
        doc="Levelized Cost of Treatment with respect to influent flowrate (accounts for seed media acquisition and spent media disposal costs and sale of treated water)",
    )


def display_costing(m):
    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++ DISPLAY COSTING ++++++++++++++++++++")
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

    print("\n------------- Unit Capital Costs -------------")
    for u in m.fs.costing._registered_unit_costing:
        print(
            f"{u.name} :   {value(pyunits.convert(u.capital_cost, to_units=pyunits.USD_2020)):.3f} $"
        )

    print("\n------------- Utility Costs -------------")
    for f in m.fs.costing.flow_types:
        print(
            f"{f} :   {value(pyunits.convert(m.fs.costing.aggregate_flow_costs[f], to_units=pyunits.USD_2020 / pyunits.year)):.3f} $/year"
        )

    print("\n---------------------------------------")

    total_capital_cost = value(
        pyunits.convert(m.fs.costing.total_capital_cost, to_units=pyunits.MUSD_2020)
    )
    print(f"Total Capital Costs: {total_capital_cost:.3f} M$")

    total_operating_cost = value(
        pyunits.convert(
            m.fs.costing.total_operating_cost, to_units=pyunits.MUSD_2020 / pyunits.year
        )
    )
    print(f"Total Operating Costs: {total_operating_cost:.3f} M$/year")

    electricity_intensity = value(
        pyunits.convert(
            m.fs.costing.electricity_intensity, to_units=pyunits.kWh / pyunits.m**3
        )
    )
    print(
        f"Electricity Intensity with respect to influent flowrate: {electricity_intensity:.3f} kWh/m^3"
    )

    LCOW = value(
        pyunits.convert(m.fs.costing.LCOW, to_units=pyunits.USD_2020 / pyunits.m**3)
    )
    print(f"Levelized Cost of Water: {LCOW:.3f} $/m^3")

    LCOT = value(
        pyunits.convert(m.fs.costing.LCOT, to_units=pyunits.USD_2020 / pyunits.m**3)
    )
    print(
        f"Levelized Cost of Treatment with respect to influent flowrate: {LCOT:.3f} $/m^3"
    )

    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


if __name__ == "__main__":
    m, results = main()
