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

from idaes.core.util.model_statistics import degrees_of_freedom
from watertap.core.util.initialization import assert_degrees_of_freedom

from watertap.core.wt_database import Database
import watertap.core.zero_order_properties as prop_ZO
from watertap.unit_models.zero_order import (
    FeedZO,
    MagprexZO,
    CentrifugeZO,
    StruviteClassifierZO,
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
    m.fs.prop = prop_ZO.WaterParameterBlock(solute_list=["phosphates", "struvite"])

    # unit models
    # feed
    m.fs.feed = FeedZO(property_package=m.fs.prop)

    # Magprex reactor treatment for struvite precipitation
    m.fs.magprex = MagprexZO(property_package=m.fs.prop, database=m.db)

    # centrifuge treatment
    m.fs.centrifuge = CentrifugeZO(property_package=m.fs.prop, database=m.db)

    # struvite classifier treatment
    m.fs.classifier = StruviteClassifierZO(property_package=m.fs.prop, database=m.db)

    # product and waste streams
    m.fs.struvite_product = Product(property_package=m.fs.prop)
    m.fs.biosolid_product = Product(property_package=m.fs.prop)
    m.fs.centrate = Product(property_package=m.fs.prop)

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.magprex.inlet)
    m.fs.s02 = Arc(source=m.fs.magprex.treated, destination=m.fs.centrifuge.inlet)
    m.fs.s03 = Arc(source=m.fs.centrifuge.treated, destination=m.fs.centrate.inlet)
    m.fs.s04 = Arc(
        source=m.fs.centrifuge.byproduct, destination=m.fs.biosolid_product.inlet
    )
    m.fs.s05 = Arc(source=m.fs.magprex.byproduct, destination=m.fs.classifier.inlet)
    m.fs.s06 = Arc(
        source=m.fs.classifier.outlet, destination=m.fs.struvite_product.inlet
    )

    # expand arcs
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):
    # ---specifications---
    # feed
    flow_vol = 126.2 * pyunits.m**3 / pyunits.hr
    conc_mass_phosphates = 350 * pyunits.g / pyunits.m**3
    conc_mass_struvite = 1e-10 * pyunits.g / pyunits.m**3
    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "phosphates"].fix(conc_mass_phosphates)
    m.fs.feed.conc_mass_comp[0, "struvite"].fix(conc_mass_struvite)

    solve(m.fs.feed)

    # Magprex reactor
    m.fs.magprex.load_parameters_from_database(use_default_removal=True)

    # clarifier
    m.fs.centrifuge.load_parameters_from_database(use_default_removal=True)

    # struvite classifier
    m.fs.classifier.load_parameters_from_database(use_default_removal=True)


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

    unit_list = ["feed", "magprex", "centrifuge", "classifier"]
    for u in unit_list:
        m.fs.component(u).report()

    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


def add_costing(m):
    # yaml file path
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "magprex_case_1575.yaml",
    )

    # add costing block
    m.fs.costing = ZeroOrderCosting(case_study_definition=source_file)
    costing_kwargs = {"default": {"flowsheet_costing_block": m.fs.costing}}
    m.fs.magprex.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.centrifuge.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.classifier.costing = UnitModelCostingBlock(**costing_kwargs)

    m.fs.costing.cost_process()

    m.fs.costing.add_electricity_intensity(m.fs.feed.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.feed.properties[0].flow_vol)

    m.fs.costing.LCOS = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
        )
        / (
            m.fs.costing.utilization_factor
            * pyunits.convert(
                m.fs.struvite_product.properties[0].flow_mass_comp["struvite"],
                to_units=pyunits.kg / m.fs.costing.base_period,
            )
        ),
        doc="Levelized cost of struvite",
    )


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
    unit_list = ["feed", "magprex", "centrifuge", "classifier"]
    for u in unit_list:
        m.fs.component(u).report()


def display_costing(m):
    print("\nUnit Capital Costs\n")
    for u in m.fs.costing._registered_unit_costing:
        print(
            u.name,
            " : {price:0.3f} $".format(
                price=value(pyunits.convert(u.capital_cost, to_units=pyunits.USD_2018))
            ),
        )
    print("\nUtility Costs\n")
    for f in m.fs.costing.flow_types:
        print(
            f,
            " :    {price:0.3f} M$/year".format(
                price=value(
                    pyunits.convert(
                        m.fs.costing.aggregate_flow_costs[f],
                        to_units=pyunits.MUSD_2018 / pyunits.year,
                    )
                )
            ),
        )

    print("")
    total_capital_cost = value(
        pyunits.convert(m.fs.costing.total_capital_cost, to_units=pyunits.MUSD_2018)
    )
    print(f"Total Capital Costs: {total_capital_cost:.3f} M$")
    total_operating_cost = value(
        pyunits.convert(
            m.fs.costing.total_operating_cost, to_units=pyunits.MUSD_2018 / pyunits.year
        )
    )
    print(f"Total Operating Costs: {total_operating_cost:.3f} M$/year")
    electricity_intensity = value(
        pyunits.convert(
            m.fs.costing.electricity_intensity, to_units=pyunits.kWh / pyunits.m**3
        )
    )
    print(f"Electricity Intensity: {electricity_intensity:.3f} kWh/m^3")
    LCOW = value(
        pyunits.convert(m.fs.costing.LCOW, to_units=pyunits.USD_2018 / pyunits.m**3)
    )
    print(f"Levelized Cost of Water: {LCOW:.3f} $/m^3")

    LCOS = value(
        pyunits.convert(m.fs.costing.LCOS, to_units=pyunits.USD_2018 / pyunits.kg)
    )
    print(f"Levelized Cost of Struvite: {LCOS:.3f} $/kg")


if __name__ == "__main__":
    m, results = main()
