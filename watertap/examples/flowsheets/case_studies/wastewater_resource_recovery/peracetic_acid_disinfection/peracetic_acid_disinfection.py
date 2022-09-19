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
    PeraceticAcidDisinfectionZO,
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
            "peracetic_acid",
            "total_coliforms_fecal_ecoli",
        ]
    )

    # unit models
    # feed
    m.fs.feed = FeedZO(property_package=m.fs.prop)

    # PAA treatment
    m.fs.PAA = PeraceticAcidDisinfectionZO(
        property_package=m.fs.prop,
        database=m.db,
    )

    # product stream
    m.fs.treated_water = Product(property_package=m.fs.prop)

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.PAA.inlet)
    m.fs.s02 = Arc(source=m.fs.PAA.treated, destination=m.fs.treated_water.inlet)

    # expand arcs
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):
    # ---specifications---
    # feed
    flow_vol = 240000 * pyunits.m**3 / pyunits.day
    conc_mass_PAA = 1.8 * pyunits.mg / pyunits.liter
    conc_mass_ecoli = 2.0e-4 * pyunits.mg / pyunits.liter
    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "peracetic_acid"].fix(conc_mass_PAA)
    m.fs.feed.conc_mass_comp[0, "total_coliforms_fecal_ecoli"].fix(conc_mass_ecoli)

    solve(m.fs.feed)

    m.fs.PAA.load_parameters_from_database(use_default_removal=True)


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

    unit_list = ["feed", "PAA"]
    for u in unit_list:
        m.fs.component(u).report()

    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


def add_costing(m):
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "peracetic_acid_case_study.yaml",
    )

    m.fs.costing = ZeroOrderCosting(case_study_definition=source_file)

    m.fs.PAA.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.costing.cost_process()

    m.fs.costing.add_electricity_intensity(m.fs.feed.properties[0].flow_vol)

    m.fs.costing.PAA_mass = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.feed.flow_mass_comp[0, "peracetic_acid"],
            to_units=pyunits.kg / m.fs.costing.base_period,
        ),
        doc="Mass of peracetic acid consumed per year (accounting for utilization factor)",
    )

    m.fs.costing.disinfection_solution_volume = Expression(
        expr=pyunits.convert(
            m.fs.costing.PAA_mass
            / m.fs.PAA.disinfection_solution_wt_frac_PAA
            / m.fs.PAA.disinfection_solution_density,
            to_units=pyunits.gallon / m.fs.costing.base_period,
        ),
        doc="Volume of disinfection solution used per year (accounting for utilization factor)",
    )

    m.fs.costing.cost_disinfection_solution = Expression(
        expr=pyunits.convert(
            m.fs.costing.disinfection_solution_volume
            * m.fs.costing.disinfection_solution_cost,
            to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
        ),
        doc="Cost of disinfection solution used per year (accounting for utilization factor)",
    )

    m.fs.costing.LCOT = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
            + m.fs.costing.cost_disinfection_solution
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
    print(f"\nTotal Capital Costs: {total_capital_cost:.3f} M$")

    total_operating_cost = value(
        pyunits.convert(
            m.fs.costing.total_operating_cost, to_units=pyunits.MUSD_2020 / pyunits.year
        )
    )
    print(f"\nTotal Operating Costs: {total_operating_cost:.3f} M$/year")

    electricity_intensity = value(
        pyunits.convert(
            m.fs.costing.electricity_intensity, to_units=pyunits.kWh / pyunits.m**3
        )
    )
    print(
        f"\nElectricity Intensity with respect to influent flowrate: {electricity_intensity:.3f} kWh/m^3"
    )

    PAA_mass = value(
        pyunits.convert(
            m.fs.costing.PAA_mass, to_units=pyunits.kg / m.fs.costing.base_period
        )
    )
    print(
        f"\nMass of peracetic acid consumed per year (accounting for utilization factor): {PAA_mass:.3f} kg/year"
    )

    disinfection_solution_volume = value(
        pyunits.convert(
            m.fs.costing.disinfection_solution_volume,
            to_units=pyunits.gallon / m.fs.costing.base_period,
        )
    )
    print(
        f"\nVolume of disinfection solution used per year (accounting for utilization factor): {disinfection_solution_volume:.3f} gal/year"
    )

    cost_disinfection_solution = value(
        pyunits.convert(
            m.fs.costing.cost_disinfection_solution,
            to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
        )
    )
    print(
        f"\nCost of disinfection solution used per year (accounting for utilization factor): {cost_disinfection_solution:.2f} $/year"
    )

    LCOT = value(
        pyunits.convert(m.fs.costing.LCOT, to_units=pyunits.USD_2020 / pyunits.m**3)
    )
    print(
        f"\nLevelized Cost of Treatment with respect to influent flowrate: {LCOT:.3f} $/m^3"
    )

    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


if __name__ == "__main__":
    m, results = main()
