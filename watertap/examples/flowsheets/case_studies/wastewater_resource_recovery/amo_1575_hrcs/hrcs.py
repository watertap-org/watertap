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

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import Product, Mixer, MomentumMixingType, MixingType
import idaes.core.util.scaling as iscale

from idaes.core.util.model_statistics import degrees_of_freedom
from watertap.core.util.initialization import assert_degrees_of_freedom

from watertap.core.wt_database import Database
import watertap.core.zero_order_properties as prop_ZO
from watertap.unit_models.zero_order import (
    FeedZO,
    ClarifierZO,
    HRCSZO,
)
from watertap.core.zero_order_costing import ZeroOrderCosting
from watertap.costing import WaterTAPCosting


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
        solute_list=["tss", "cod", "oxygen", "carbon_dioxide"]
    )

    # unit models
    # feed
    m.fs.feed = FeedZO(property_package=m.fs.prop)

    # mixer for recycle stream
    m.fs.mixer = Mixer(
        property_package=m.fs.prop,
        inlet_list=["inlet1", "inlet2"],
        momentum_mixing_type=MomentumMixingType.none,
        energy_mixing_type=MixingType.none,
    )

    # HR-CS treatment
    m.fs.HRCS = HRCSZO(property_package=m.fs.prop, database=m.db)

    # clarifier treatment
    m.fs.clarifier = ClarifierZO(property_package=m.fs.prop, database=m.db)

    # product and waste streams
    m.fs.product = Product(property_package=m.fs.prop)
    m.fs.disposal = Product(property_package=m.fs.prop)

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.mixer.inlet1)
    m.fs.s02 = Arc(source=m.fs.mixer.outlet, destination=m.fs.HRCS.inlet)
    m.fs.s03 = Arc(source=m.fs.HRCS.treated, destination=m.fs.clarifier.inlet)
    m.fs.s04 = Arc(source=m.fs.HRCS.byproduct, destination=m.fs.disposal.inlet)
    m.fs.s05 = Arc(source=m.fs.clarifier.treated, destination=m.fs.product.inlet)
    m.fs.s06 = Arc(source=m.fs.clarifier.byproduct, destination=m.fs.mixer.inlet2)

    # expand arcs
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):
    # ---specifications---
    # feed
    flow_vol = 47317.6 * pyunits.m**3 / pyunits.hr
    conc_mass_tss = 243 * pyunits.mg / pyunits.liter
    conc_mass_co2 = 0.001 * pyunits.mg / pyunits.liter
    mass_flow_cod = 336 * pyunits.ton / pyunits.day
    mass_flow_oxygen = (
        895.15 * pyunits.ton / pyunits.day
    )  # equivalent to the number of moles of carbon being added
    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "tss"].fix(conc_mass_tss)
    m.fs.feed.conc_mass_comp[0, "carbon_dioxide"].fix(conc_mass_co2)
    m.fs.feed.properties[0].flow_mass_comp["cod"].fix(mass_flow_cod)
    m.fs.feed.properties[0].flow_mass_comp["oxygen"].fix(mass_flow_oxygen)

    solve(m.fs.feed)

    # HR-CS
    m.fs.HRCS.load_parameters_from_database(use_default_removal=True)

    # clarifier
    m.fs.clarifier.load_parameters_from_database(use_default_removal=True)
    m.fs.clarifier.removal_frac_mass_comp[0, "tss"].fix(0.914)
    m.fs.clarifier.removal_frac_mass_comp[0, "cod"].fix(0.72)
    m.fs.clarifier.recovery_frac_mass_H2O[0].fix(1)

    assert_degrees_of_freedom(m, 0)


def initialize_system(m, solver=None):

    # set up solver
    if solver is None:
        solver = get_solver()
    optarg = solver.options

    # populate initial properties throughout the system
    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.s01)
    propagate_state(m.fs.s02)
    m.fs.HRCS.initialize(optarg=optarg)
    propagate_state(m.fs.s03)
    m.fs.clarifier.initialize(optarg=optarg)
    propagate_state(m.fs.s04)
    m.fs.disposal.initialize()
    propagate_state(m.fs.s05)
    m.fs.product.initialize()
    propagate_state(m.fs.s06)
    m.fs.mixer.initialize(optarg=optarg)


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

    unit_list = ["feed", "HRCS", "clarifier"]
    for u in unit_list:
        m.fs.component(u).report()

    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


def add_costing(m):
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "hrcs_case_1575.yaml",
    )

    m.fs.costing = ZeroOrderCosting(case_study_definition=source_file)
    m.fs.watertap_costing = WaterTAPCosting()

    m.fs.mixer.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.watertap_costing
    )
    costing_kwargs = {"flowsheet_costing_block": m.fs.costing}

    m.fs.HRCS.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.clarifier.costing = UnitModelCostingBlock(**costing_kwargs)

    m.fs.costing.cost_process()
    m.fs.watertap_costing.cost_process()

    m.fs.costing.add_electricity_intensity(m.fs.feed.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.feed.properties[0].flow_vol)

    # Combine results from costing packages and calculate overall metrics
    @m.Expression()
    def total_capital_cost(b):
        return pyunits.convert(
            m.fs.costing.total_capital_cost, to_units=m.fs.costing.base_currency
        ) + pyunits.convert(
            m.fs.watertap_costing.total_investment_cost,
            to_units=m.fs.costing.base_currency,
        )

    @m.Expression()
    def total_operating_cost(b):
        return (
            pyunits.convert(
                m.fs.costing.total_fixed_operating_cost,
                to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
            )
            + pyunits.convert(
                m.fs.costing.total_variable_operating_cost,
                to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
            )
            + pyunits.convert(
                m.fs.watertap_costing.total_operating_cost,
                to_units=m.fs.costing.base_currency / m.fs.costing.base_period,
            )
        )


def display_results(m):
    unit_list = ["feed", "mixer", "HRCS", "clarifier"]
    for u in unit_list:
        m.fs.component(u).report()


def display_costing(m):

    print("\n----------Unit Capital Costs----------\n")
    for u in m.fs.costing._registered_unit_costing:
        print(
            u.name,
            " : {price:0.3f} $".format(
                price=value(pyunits.convert(u.capital_cost, to_units=pyunits.USD_2018))
            ),
        )
    print("\n----------Utility Costs----------\n")
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

    print("\n----------Capital costs----------\n")
    total_capital_costs = value(m.fs.costing.total_capital_cost) / 1e3
    print(f"Total capital costs: {total_capital_costs:.3f} $k")

    print("\n----------Operating costs----------\n")
    total_operating_costs = value(m.fs.costing.total_operating_cost) / 1e6
    print(f"Total operating costs: {total_operating_costs:.5f} $M/year")

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


if __name__ == "__main__":
    m, results = main()
