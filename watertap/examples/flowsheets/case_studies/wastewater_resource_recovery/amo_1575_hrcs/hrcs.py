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
    PrimarySeparatorZO,
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
    display_additional_results(m)

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
    m.fs.clarifier = ClarifierZO(
        property_package=m.fs.prop,
        database=m.db,
        process_subtype="HRCS_clarifier",
    )

    # Separator for recycle loop
    m.fs.sep = PrimarySeparatorZO(property_package=m.fs.prop, database=m.db)

    # product and waste streams
    m.fs.product = Product(property_package=m.fs.prop)
    m.fs.WAS_product = Product(property_package=m.fs.prop)
    m.fs.disposal = Product(property_package=m.fs.prop)

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.mixer.inlet1)
    m.fs.s02 = Arc(source=m.fs.mixer.outlet, destination=m.fs.HRCS.inlet)
    m.fs.s03 = Arc(source=m.fs.HRCS.treated, destination=m.fs.clarifier.inlet)
    m.fs.s04 = Arc(source=m.fs.HRCS.byproduct, destination=m.fs.disposal.inlet)
    m.fs.s05 = Arc(source=m.fs.clarifier.treated, destination=m.fs.product.inlet)
    m.fs.s06 = Arc(source=m.fs.clarifier.byproduct, destination=m.fs.sep.inlet)
    m.fs.s07 = Arc(source=m.fs.sep.treated, destination=m.fs.WAS_product.inlet)
    m.fs.s08 = Arc(source=m.fs.sep.byproduct, destination=m.fs.mixer.inlet2)

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
    mass_flow_oxygen = 895.15 * pyunits.ton / pyunits.day
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

    # separator
    m.fs.sep.load_parameters_from_database(use_default_removal=True)
    m.fs.sep.removal_frac_mass_comp[0, "tss"].fix(0.95)
    m.fs.sep.removal_frac_mass_comp[0, "cod"].fix(0.95)
    m.fs.sep.recovery_frac_mass_H2O[0].fix(0.05)

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
    m.fs.sep.initialize(optarg=optarg)
    propagate_state(m.fs.s07)
    m.fs.WAS_product.initialize()
    propagate_state(m.fs.s08)
    m.fs.mixer.initialize(optarg=optarg)


def solve(blk, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results


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

    m.fs.HRCS.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.clarifier.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.sep.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

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

    m.fs.costing.annual_water_production = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.product.properties[0].flow_vol,
            to_units=pyunits.m**3 / m.fs.costing.base_period,
        )
    )

    m.fs.costing.annual_water_inlet = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.feed.properties[0].flow_vol,
            to_units=pyunits.m**3 / m.fs.costing.base_period,
        )
    )

    m.fs.costing.total_annualized_cost = Expression(
        expr=(
            m.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.total_operating_cost
        )
    )

    m.fs.costing.LCOW_with_revenue = Expression(
        expr=(
            (
                m.fs.costing.total_annualized_cost
                + m.fs.costing.aggregate_flow_costs["ferric_chloride"]
            )
            / m.fs.costing.annual_water_production
        ),
        doc="Levelized Cost of Water With Revenue",
    )

    m.fs.costing.LCOT = Expression(
        expr=(m.fs.costing.total_annualized_cost / m.fs.costing.annual_water_inlet),
        doc="Levelized Cost of Treatment Without Revenue",
    )

    m.fs.costing.LCOT_with_revenue = Expression(
        expr=(
            (
                m.fs.costing.total_annualized_cost
                + m.fs.costing.aggregate_flow_costs["ferric_chloride"]
            )
            / m.fs.costing.annual_water_inlet
        ),
        doc="Levelized Cost of Treatment With Revenue",
    )


def display_results(m):
    unit_list = ["feed", "HRCS", "clarifier", "sep", "mixer"]
    for u in unit_list:
        m.fs.component(u).report()

    print("\n----------Performance Metrics----------\n")
    volumetric_recovery = value(
        m.fs.product.properties[0].flow_vol / m.fs.feed.properties[0].flow_vol
    )
    print(f"Volumetric recovery: {volumetric_recovery:.4f} m3 of product/m3 of feed")

    sys_water_recovery = (
        m.fs.product.flow_mass_comp[0, "H2O"]() / m.fs.feed.flow_mass_comp[0, "H2O"]()
    )
    print(f"Water Recovery: {sys_water_recovery*100 : .3f} %")

    carbon_capture = value(
        (
            m.fs.WAS_product.flow_mass_comp[0, "cod"]()
            + m.fs.product.flow_mass_comp[0, "cod"]()
        )
        / m.fs.feed.flow_mass_comp[0, "cod"]
    )
    print(f"Carbon Capture (includes non-biomass COD): {carbon_capture*100 : .3f} %")


def display_costing(m):

    print("\n----------Unit Capital Costs----------\n")
    for u in m.fs.costing._registered_unit_costing:
        print(
            u.name,
            " : {price:0.3f} $".format(
                price=value(pyunits.convert(u.capital_cost, to_units=pyunits.USD_2018))
            ),
        )
    for u in m.fs.watertap_costing._registered_unit_costing:
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
            " :    {price:0.3f} $M/year".format(
                price=value(
                    pyunits.convert(
                        m.fs.costing.aggregate_flow_costs[f],
                        to_units=pyunits.MUSD_2018 / pyunits.year,
                    )
                )
            ),
        )

    ferric_chloride_cost = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_costs["ferric_chloride"],
            to_units=pyunits.USD_2018 / pyunits.year,
        )
        / pyunits.convert(
            m.fs.product.properties[0].flow_vol, to_units=pyunits.m**3 / pyunits.year
        )
    )
    print(f"ferric chloride: { ferric_chloride_cost: .4f} $/m3 of product")

    print("\n----------Total Costs----------\n")

    total_annualized_cost = value(
        pyunits.convert(
            m.fs.costing.total_annualized_cost,
            to_units=pyunits.MUSD_2018 / pyunits.year,
        )
    )
    print(f"Total Annual Cost: { total_annualized_cost: .4f} $M/year")

    total_capital_cost = value(
        pyunits.convert(m.total_capital_cost, to_units=pyunits.MUSD_2018)
    )
    print(f"Total Capital Costs: {total_capital_cost:.3f} $M")

    total_operating_cost = value(
        pyunits.convert(
            m.total_operating_cost, to_units=pyunits.MUSD_2018 / pyunits.year
        )
    )
    print(f"Total Operating Costs: {total_operating_cost:.3f} $M/year")

    opex_fraction = 100 * value(total_operating_cost / total_annualized_cost)
    print(f"Opex Fraction of Annual Cost: {opex_fraction:.4f} %")

    SEC = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_electricity
            / m.fs.product.properties[0].flow_vol,
            to_units=pyunits.kWh / pyunits.m**3,
        )
    )
    print(f"Specific electricity consumption: {SEC:.4f} kWh/m3 of product")

    print("\n----------Levelized Costs----------\n")

    LCOW_with_revenue = value(
        pyunits.convert(
            m.fs.costing.LCOW_with_revenue,
            to_units=m.fs.costing.base_currency / pyunits.m**3,
        )
    )
    print(
        f"Levelized Cost of Water with Revenue: {LCOW_with_revenue:.4f} $/m3 of product"
    )
    LCOW = value(
        pyunits.convert(m.fs.costing.LCOW, to_units=pyunits.USD_2018 / pyunits.m**3)
    )
    print(f"Levelized Cost of Water: {LCOW:.3f} $/m^3 of product")

    LCOT_with_revenue = value(
        pyunits.convert(
            m.fs.costing.LCOT_with_revenue, to_units=pyunits.USD_2018 / pyunits.m**3
        )
    )
    print(f"Levelized Cost of Treatment With Revenue: {LCOT_with_revenue:.3f} $/m^3")
    LCOT = value(
        pyunits.convert(m.fs.costing.LCOT, to_units=pyunits.USD_2018 / pyunits.m**3)
    )
    print(f"Levelized Cost of Treatment: {LCOT:.3f} $/m^3")


def display_additional_results(m):
    print("\n----------Outlets----------\n")
    product_H2O_flow = value(
        pyunits.convert(
            m.fs.product.properties[0].flow_vol,
            to_units=pyunits.m**3 / pyunits.hr,
        )
    )
    print(f"H2O outlet flow: {product_H2O_flow:.4f} m3/h")

    product_H2O_tss = value(
        pyunits.convert(
            m.fs.product.properties[0].conc_mass_comp["tss"],
            to_units=pyunits.g / pyunits.L,
        )
    )
    print(f"H2O outlet tss conc: {product_H2O_tss:.4f} g/L")

    product_H2O_cod = value(
        pyunits.convert(
            m.fs.product.properties[0].conc_mass_comp["cod"],
            to_units=pyunits.g / pyunits.L,
        )
    )
    print(f"H2O outlet cod conc: {product_H2O_cod:.4f} g/L")

    product_H2O_oxgyen = value(
        pyunits.convert(
            m.fs.product.properties[0].conc_mass_comp["oxygen"],
            to_units=pyunits.g / pyunits.L,
        )
    )
    print(f"H2O outlet oxygen conc: {product_H2O_oxgyen:.4f} g/L")

    product_H2O_carbon_dioxide = value(
        pyunits.convert(
            m.fs.product.properties[0].conc_mass_comp["carbon_dioxide"],
            to_units=pyunits.g / pyunits.L,
        )
    )
    print(f"H2O outlet CO2 conc: {product_H2O_carbon_dioxide:.4f} g/L")


if __name__ == "__main__":
    m, results = main()
