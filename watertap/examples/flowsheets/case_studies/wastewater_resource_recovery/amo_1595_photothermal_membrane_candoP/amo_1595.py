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

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.prop = prop_ZO.WaterParameterBlock(
        default={
            "solute_list": [
                "nitrogen",
                "phosphates",
                "bioconcentrated_phosphorous",
                "nitrous_oxide",
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

    # photothermal membrane
    m.fs.photothermal = PhotothermalMembraneZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
        },
    )

    # CANDO+P reactor
    # Note that CANDOPZO model electricity costing includes pumping costs for
    # the unit, so a pump has not been included in this flowsheet between the
    # photothermal membrane and CANDO+P units
    m.fs.candop = CANDOPZO(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
        },
    )

    # product streams
    m.fs.photothermal_water = Product(default={"property_package": m.fs.prop})
    m.fs.candop_treated = Product(default={"property_package": m.fs.prop})
    m.fs.candop_byproduct = Product(default={"property_package": m.fs.prop})

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

    unit_list = ["feed", "pump", "photothermal", "candop"]
    for u in unit_list:
        m.fs.component(u).report()

    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


def add_costing(m):
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "amo_1595_case_study.yaml",
    )

    m.fs.costing = ZeroOrderCosting(default={"case_study_definition": source_file})

    costing_kwargs = {"default": {"flowsheet_costing_block": m.fs.costing}}
    m.fs.pump.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.photothermal.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.candop.costing = UnitModelCostingBlock(**costing_kwargs)

    m.fs.costing.cost_process()

    m.fs.costing.add_electricity_intensity(m.fs.feed.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.photothermal_water.properties[0].flow_vol)

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
        doc="Dollar value of BCP generated per year",
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
        doc="Dollar value of byproduct water generated per year",
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
        doc="Dollar value of byproduct nitrous oxide generated per year",
    )

    m.fs.costing.LCOW_bcp = Expression(
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
        doc="Levelized Cost of Water when accounting for salable BCP",
    )

    m.fs.costing.LCOT = Expression(
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
        doc="Levelized Cost of Treatment with respect to influent flowrate (accounts for salable BCP)",
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

    LCOW_with_bcp = value(
        pyunits.convert(
            m.fs.costing.LCOW_bcp, to_units=pyunits.USD_2020 / pyunits.m**3
        )
    )
    print(
        f"Levelized Cost of Water (accounting for BCP sale): {LCOW_with_bcp:.3f} $/m^3"
    )

    LCOT = value(
        pyunits.convert(m.fs.costing.LCOT, to_units=pyunits.USD_2020 / pyunits.m**3)
    )
    print(
        f"Levelized Cost of Treatment with respect to influent flowrate: {LCOW:.3f} $/m^3"
    )

    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


if __name__ == "__main__":
    m, results = main()
