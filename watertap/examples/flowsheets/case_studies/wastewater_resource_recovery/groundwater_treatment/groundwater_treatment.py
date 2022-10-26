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

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.prop = prop_ZO.WaterParameterBlock(
        solute_list=[
            "arsenic",
            "uranium",
            "nitrate",
            "phosphates",
            "iron",
            "filtration_media",
        ]
    )

    # unit models
    # feed
    m.fs.feed = FeedZO(property_package=m.fs.prop)

    # pump
    m.fs.pump = PumpElectricityZO(property_package=m.fs.prop, database=m.db)

    # microbial battery
    m.fs.micbatt = MicrobialBatteryZO(property_package=m.fs.prop, database=m.db)

    # product streams
    m.fs.filtered_water = Product(property_package=m.fs.prop)
    m.fs.byproduct = Product(property_package=m.fs.prop)

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.pump.inlet)
    m.fs.s02 = Arc(source=m.fs.pump.outlet, destination=m.fs.micbatt.inlet)
    m.fs.s03 = Arc(source=m.fs.micbatt.treated, destination=m.fs.filtered_water.inlet)
    m.fs.s04 = Arc(source=m.fs.micbatt.byproduct, destination=m.fs.byproduct.inlet)

    # expand arcs
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):
    # ---specifications---
    # feed
    flow_vol = 1 * pyunits.m**3 / pyunits.day
    conc_mass_arsenic = 0.04 * pyunits.mg / pyunits.L
    conc_mass_uranium = 0.06 * pyunits.mg / pyunits.L
    conc_mass_nitrate = 10 * pyunits.mg / pyunits.L
    conc_mass_phosphates = 0.1 * pyunits.mg / pyunits.L
    conc_mass_iron = 0.5 * pyunits.mg / pyunits.L
    conc_mass_filtration_media = 0.1 * pyunits.mg / pyunits.L

    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "arsenic"].fix(conc_mass_arsenic)
    m.fs.feed.conc_mass_comp[0, "uranium"].fix(conc_mass_uranium)
    m.fs.feed.conc_mass_comp[0, "nitrate"].fix(conc_mass_nitrate)
    m.fs.feed.conc_mass_comp[0, "phosphates"].fix(conc_mass_phosphates)
    m.fs.feed.conc_mass_comp[0, "iron"].fix(conc_mass_iron)
    m.fs.feed.conc_mass_comp[0, "filtration_media"].fix(conc_mass_filtration_media)

    solve(m.fs.feed)

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

    unit_list = ["pump", "micbatt"]
    for u in unit_list:
        m.fs.component(u).report()

    print("\n---------- Feed volumetric flowrate ----------")
    feed_vol_flow = value(
        pyunits.convert(m.fs.feed.flow_vol[0], to_units=pyunits.m**3 / pyunits.day)
    )
    print(f"{feed_vol_flow} m^3/day")

    print("\n---------- Feed water concentrations ----------")
    as_feed_conc = value(
        pyunits.convert(
            m.fs.feed.conc_mass_comp[0, "arsenic"], to_units=pyunits.mg / pyunits.L
        )
    )
    print(f"Arsenic: {as_feed_conc:.6f} mg/L")

    u_feed_conc = value(
        pyunits.convert(
            m.fs.feed.conc_mass_comp[0, "uranium"], to_units=pyunits.mg / pyunits.L
        )
    )
    print(f"Uranium: {u_feed_conc:.6f} mg/L")

    nitrates_feed_conc = value(
        pyunits.convert(
            m.fs.feed.conc_mass_comp[0, "nitrate"], to_units=pyunits.mg / pyunits.L
        )
    )
    print(f"Nitrates: {nitrates_feed_conc:.6f} mg/L")

    phosphates_feed_conc = value(
        pyunits.convert(
            m.fs.feed.conc_mass_comp[0, "phosphates"], to_units=pyunits.mg / pyunits.L
        )
    )
    print(f"Phosphates: {phosphates_feed_conc:.6f} mg/L")

    iron_feed_conc = value(
        pyunits.convert(
            m.fs.feed.conc_mass_comp[0, "iron"], to_units=pyunits.mg / pyunits.L
        )
    )
    print(f"Iron: {iron_feed_conc:.6f} mg/L")

    media_feed_conc = value(
        pyunits.convert(
            m.fs.feed.conc_mass_comp[0, "filtration_media"],
            to_units=pyunits.mg / pyunits.L,
        )
    )
    print(f"Filtration media: {media_feed_conc:.6f} mg/L")

    print("\n---------- Treated water flowrate ----------")
    treated_vol_flow = value(
        pyunits.convert(
            m.fs.filtered_water.properties[0].flow_vol,
            to_units=pyunits.m**3 / pyunits.day,
        )
    )
    print(f"{treated_vol_flow} m^3/day")

    print("\n---------- Treated water concentrations ----------")
    as_treated_conc = value(
        pyunits.convert(
            m.fs.micbatt.properties_treated[0].conc_mass_comp["arsenic"],
            to_units=pyunits.mg / pyunits.L,
        )
    )
    print(f"Arsenic: {as_treated_conc:.6f} mg/L")

    u_treated_conc = value(
        pyunits.convert(
            m.fs.micbatt.properties_treated[0].conc_mass_comp["uranium"],
            to_units=pyunits.mg / pyunits.L,
        )
    )
    print(f"Uranium: {u_treated_conc:.6f} mg/L")

    nitrates_treated_conc = value(
        pyunits.convert(
            m.fs.micbatt.properties_treated[0].conc_mass_comp["nitrate"],
            to_units=pyunits.mg / pyunits.L,
        )
    )
    print(f"Nitrates: {nitrates_treated_conc:.6f} mg/L")

    phosphates_treated_conc = value(
        pyunits.convert(
            m.fs.micbatt.properties_treated[0].conc_mass_comp["phosphates"],
            to_units=pyunits.mg / pyunits.L,
        )
    )
    print(f"Phosphates: {phosphates_treated_conc:.6f} mg/L")

    iron_treated_conc = value(
        pyunits.convert(
            m.fs.micbatt.properties_treated[0].conc_mass_comp["iron"],
            to_units=pyunits.mg / pyunits.L,
        )
    )
    print(f"Iron: {iron_treated_conc:.6f} mg/L")

    print("\n---------- Byproduct mass flow rates ----------")
    as_byproduct_mass_flow = value(
        pyunits.convert(
            m.fs.micbatt.byproduct.flow_mass_comp[0, "arsenic"],
            to_units=pyunits.mg / pyunits.day,
        )
    )
    print(f"Arsenic: {as_byproduct_mass_flow:.2f} mg/day")

    u_byproduct_mass_flow = value(
        pyunits.convert(
            m.fs.micbatt.byproduct.flow_mass_comp[0, "uranium"],
            to_units=pyunits.mg / pyunits.day,
        )
    )
    print(f"Uranium: {u_byproduct_mass_flow:.2f} mg/day")

    n_byproduct_mass_flow = value(
        pyunits.convert(
            m.fs.micbatt.byproduct.flow_mass_comp[0, "nitrate"],
            to_units=pyunits.mg / pyunits.day,
        )
    )
    print(f"Nitrates: {n_byproduct_mass_flow:.2f} mg/day")

    p_byproduct_mass_flow = value(
        pyunits.convert(
            m.fs.micbatt.byproduct.flow_mass_comp[0, "phosphates"],
            to_units=pyunits.mg / pyunits.day,
        )
    )
    print(f"Phosphates: {p_byproduct_mass_flow:.2f} mg/day")

    iron_byproduct_mass_flow = value(
        pyunits.convert(
            m.fs.micbatt.byproduct.flow_mass_comp[0, "iron"],
            to_units=pyunits.mg / pyunits.day,
        )
    )
    print(f"Iron: {iron_byproduct_mass_flow:.2f} mg/day")

    media_byproduct_mass_flow = value(
        pyunits.convert(
            m.fs.micbatt.byproduct.flow_mass_comp[0, "filtration_media"],
            to_units=pyunits.mg / pyunits.day,
        )
    )
    print(f"Filtration media: {media_byproduct_mass_flow:.2f} mg/day")

    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


def add_costing(m):
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "groundwater_treatment_case_study.yaml",
    )

    m.fs.costing = ZeroOrderCosting(case_study_definition=source_file)

    costing_kwargs = {"flowsheet_costing_block": m.fs.costing}
    m.fs.pump.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.micbatt.costing = UnitModelCostingBlock(**costing_kwargs)

    m.fs.costing.cost_process()

    m.fs.costing.add_electricity_intensity(m.fs.feed.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.feed.properties[0].flow_vol)

    m.fs.costing.annual_water_inlet = Expression(
        expr=m.fs.costing.utilization_factor
        * pyunits.convert(
            m.fs.feed.properties[0].flow_vol,
            to_units=pyunits.m**3 / m.fs.costing.base_period,
        )
    )

    m.fs.costing.LCOT = Expression(
        expr=(
            m.fs.costing.total_capital_cost * m.fs.costing.capital_recovery_factor
            + m.fs.costing.total_operating_cost
        )
        / (
            pyunits.convert(
                m.fs.feed.flow_vol[0],
                to_units=pyunits.m**3 / m.fs.costing.base_period,
            )
            * m.fs.costing.utilization_factor
        ),
        doc="Levelized Cost of Treatment (with respect to feed flowrate)",
    )


def display_costing(m):
    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++ DISPLAY COSTING ++++++++++++++++++++")
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

    print("\n---------- Levelized costs ----------")
    LCOT = value(
        pyunits.convert(m.fs.costing.LCOT, to_units=pyunits.USD_2020 / pyunits.m**3)
    )
    print(f"Levelized cost of treatment: {LCOT:.2f} $/m^3")

    LCOW = value(
        pyunits.convert(m.fs.costing.LCOW, to_units=pyunits.USD_2020 / pyunits.m**3)
    )
    print(f"Levelized cost of water: {LCOW:.2f} $/m^3")

    print("\n------------- Capital costs -------------")
    DCC_normalized = value(
        pyunits.convert(
            (m.fs.pump.costing.capital_cost + m.fs.micbatt.costing.capital_cost)
            / m.fs.costing.TIC
            / m.fs.feed.properties[0].flow_vol,
            to_units=m.fs.costing.base_currency / (pyunits.m**3 / pyunits.day),
        )
    )
    print(f"Normalized direct capital costs: {DCC_normalized:.2f} $/(m^3/day)")

    ICC_normalized = value(
        pyunits.convert(
            m.fs.costing.total_capital_cost / m.fs.feed.properties[0].flow_vol,
            to_units=m.fs.costing.base_currency / (pyunits.m**3 / pyunits.day),
        )
    )
    print(f"Normalized total capital costs: {ICC_normalized:.2f} $/(m^3/day)")

    total_capital_costs = value(
        pyunits.convert(m.fs.costing.total_capital_cost, to_units=pyunits.USD_2020)
    )
    print(f"Total capital costs: {total_capital_costs:.2f} $")

    print(
        f"Pump capital cost: {value(pyunits.convert(m.fs.pump.costing.capital_cost, to_units=pyunits.USD_2020)):.2f} $"
    )

    print(
        f"Microbial battery capital cost: {value(pyunits.convert(m.fs.micbatt.costing.capital_cost, to_units=pyunits.USD_2020)):.2f} $"
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
            to_units=m.fs.costing.base_currency / pyunits.m**3,
        )
    )
    print(f"Normalized electricity cost: {EC_normalized:.3f} $/m^3 of feed")

    fresh_media_cost_normalized = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_costs["filtration_media"]
            * m.fs.costing.utilization_factor
            / m.fs.costing.annual_water_inlet,
            to_units=m.fs.costing.base_currency / pyunits.m**3,
        )
    )
    print(
        f"Normalized fresh media cost: {fresh_media_cost_normalized:.3f} $/m^3 of feed"
    )

    spent_media_disposal_cost_normalized = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_costs["filtration_media_disposal"]
            * m.fs.costing.utilization_factor
            / m.fs.costing.annual_water_inlet,
            to_units=m.fs.costing.base_currency / pyunits.m**3,
        )
    )
    print(
        f"Normalized spent media disposal cost: {spent_media_disposal_cost_normalized:.3f} $/m^3 of feed"
    )

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

    electricity_operating_costs = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_costs["electricity"]
            * m.fs.costing.utilization_factor,
            to_units=pyunits.USD_2020 / pyunits.year,
        )
    )
    print(f"Electricity operating costs: {electricity_operating_costs:.2f} $/year")

    fresh_media_costs = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_costs["filtration_media"]
            * m.fs.costing.utilization_factor,
            to_units=pyunits.USD_2020 / pyunits.year,
        )
    )
    print(f"Fresh media operating costs: {fresh_media_costs:.2f} $/year")

    media_disposal_costs = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_costs["filtration_media_disposal"]
            * m.fs.costing.utilization_factor,
            to_units=pyunits.USD_2020 / pyunits.year,
        )
    )
    print(f"Spent media disposal operating costs: {media_disposal_costs:.2f} $/year")

    print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


if __name__ == "__main__":
    m, results = main()
