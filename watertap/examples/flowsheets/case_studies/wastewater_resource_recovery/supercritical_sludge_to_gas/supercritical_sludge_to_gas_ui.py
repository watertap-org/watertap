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
from watertap.ui.fsapi import FlowsheetInterface
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.supercritical_sludge_to_gas.supercritical_sludge_to_gas import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    add_costing,
)
from pyomo.environ import units as pyunits, assert_optimal_termination
from pyomo.util.check_units import assert_units_consistent


def export_to_ui():
    return FlowsheetInterface(
        name="Supercritical sludge to gas",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
    )


def export_variables(flowsheet=None, exports=None):
    fs = flowsheet
    # --- Input data ---
    # Feed conditions
    flow_mass = (
        fs.feed.flow_mass_comp[0, "H2O"]
        + fs.feed.flow_mass_comp[0, "organic_solid"]
        + fs.feed.flow_mass_comp[0, "inorganic_solid"]
        + fs.feed.flow_mass_comp[0, "organic_liquid"]
        + fs.feed.flow_mass_comp[0, "carbon_dioxide"]
    )
    exports.add(
        obj=flow_mass,
        name="Mass flow rate",
        ui_units=pyunits.metric_ton / pyunits.day,
        display_units="ton/day",
        rounding=1,
        description="Inlet mass flow rate",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.flow_mass_comp[0, "organic_solid"],
        name="Organics(s) mass flow",
        ui_units=pyunits.metric_ton / pyunits.day,
        display_units="ton/day",
        rounding=1,
        description="Inlet organics(solid) mass flow",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.flow_mass_comp[0, "inorganic_solid"],
        name="Inorganics(s) mass flow",
        ui_units=pyunits.metric_ton / pyunits.day,
        display_units="ton/day",
        rounding=1,
        description="Inlet inorganics(solid) mass flow",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )

    # Unit model data, AT-HTL
    exports.add(
        obj=fs.ATHTL.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="AT-HTL",
        is_output=False,
    )
    exports.add(
        obj=fs.ATHTL.reaction_conversion[0, "hydrothermal_liquefaction"],
        name="Organics(s) conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Organics(s) conversion [g-organics(s) reacted/g-organics(s) inlet]",
        is_input=True,
        input_category="AT-HTL",
        is_output=False,
    )
    exports.add(
        obj=fs.ATHTL.generation_ratio["hydrothermal_liquefaction", "organic_liquid"],
        name="Organics(l) conversion ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=4,
        description="Organics(l) mass conversion ratio with respect to organics(s) [g-organics(l) produced"
        "/g-organics(s) reacted]",
        is_input=True,
        input_category="AT-HTL",
        is_output=False,
    )
    exports.add(
        obj=fs.ATHTL.generation_ratio["hydrothermal_liquefaction", "carbon_dioxide"],
        name="CO2 conversion ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=4,
        description="CO2 mass conversion ratio with respect to organics(s) [g-CO2 produced"
        "/g-organics(s) reacted]",
        is_input=True,
        input_category="AT-HTL",
        is_output=False,
    )
    exports.add(
        obj=fs.ATHTL.energy_electric_flow_mass,
        name="Electricity specific power",
        ui_units=pyunits.kWh / pyunits.metric_ton,
        display_units="kWh/ton of inlet flow rate",
        rounding=2,
        description="Specific energy consumption with respect to influent mass",
        is_input=True,
        input_category="AT-HTL",
        is_output=False,
    )
    exports.add(
        obj=fs.ATHTL.catalyst_dosage,
        name="Catalyst dosing",
        ui_units=pyunits.pound / pyunits.metric_ton,
        display_units="lb/ton of inlet flow rate",
        rounding=2,
        description="Catalyst dosing",
        is_input=True,
        input_category="AT-HTL",
        is_output=False,
    )

    # Unit model data, supercritical salt precipitation
    exports.add(
        obj=fs.salt_precipitation.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Supercritical salt precipitation",
        is_output=False,
    )
    exports.add(
        obj=fs.salt_precipitation.removal_frac_mass_comp[0, "organic_solid"],
        name="Organics(s) removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=0,
        description="Organics(s) removal",
        is_input=True,
        input_category="Supercritical salt precipitation",
        is_output=False,
    )
    exports.add(
        obj=fs.salt_precipitation.removal_frac_mass_comp[0, "inorganic_solid"],
        name="Inorganics(s) removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=4,
        description="Inorganics(s) removal",
        is_input=True,
        input_category="Supercritical salt precipitation",
        is_output=False,
    )
    exports.add(
        obj=fs.salt_precipitation.removal_frac_mass_comp[0, "organic_liquid"],
        name="Organics(l) removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=0,
        description="Organics(l) removal",
        is_input=True,
        input_category="Supercritical salt precipitation",
        is_output=False,
    )
    exports.add(
        obj=fs.salt_precipitation.energy_electric_flow_mass,
        name="Electricity specific power",
        ui_units=pyunits.kWh / pyunits.metric_ton,
        display_units="kWh/ton of inlet flow rate",
        rounding=2,
        description="Specific power relating the power to inlet mass flow rate",
        is_input=True,
        input_category="Supercritical salt precipitation",
        is_output=False,
    )

    # Unit model data, HTG & WWT units
    exports.add(
        obj=fs.HTG.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="HTG",
        is_output=False,
    )
    exports.add(
        obj=fs.HTG.reaction_conversion[0, "hydrothermal_gasification"],
        name="Organics(l) conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Organics(l) conversion [g-organics(l) reacted/g-organics(l) inlet]",
        is_input=True,
        input_category="HTG",
        is_output=False,
    )
    exports.add(
        obj=fs.HTG.generation_ratio["hydrothermal_gasification", "carbon_dioxide"],
        name="CO2 conversion ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=4,
        description="CO2 mass conversion ratio with respect to organics(l) [g-CO2 produced"
        "/g-organics(l) reacted]",
        is_input=True,
        input_category="HTG",
        is_output=False,
    )
    exports.add(
        obj=fs.HTG.energy_electric_flow_mass,
        name="Electricity specific power",
        ui_units=pyunits.kWh / pyunits.metric_ton,
        display_units="kWh/ton of inlet flow rate",
        rounding=2,
        description="Specific power relating the power to inlet mass flow rate",
        is_input=True,
        input_category="HTG",
        is_output=False,
    )
    exports.add(
        obj=fs.HTG.catalyst_dosage,
        name="Catalyst dosing",
        ui_units=pyunits.pound / pyunits.metric_ton,
        display_units="lb/ton of inlet flow rate",
        rounding=2,
        description="Catalyst dosing",
        is_input=True,
        input_category="HTG",
        is_output=False,
    )

    # System costing
    exports.add(
        obj=fs.costing.utilization_factor,
        name="Utilization factor",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Utilization factor - [annual use hours/total hours in year]",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.TIC,
        name="Practical investment factor",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="Practical investment factor - [total investment cost/direct "
        "capital costs]",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.plant_lifetime,
        name="Plant lifetime",
        ui_units=pyunits.year,
        display_units="years",
        rounding=1,
        description="Plant lifetime",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.wacc,
        name="Discount rate",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Discount rate used in calculating the capital annualization",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.maintenance_costs_percent_FCI,
        name="Fixed operating cost factor",
        ui_units=1 / pyunits.year,
        display_units="fraction/year",
        rounding=2,
        description="Fixed operating cost factor - [annual fixed operating cost/total "
        "investment cost]",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.electricity_cost,
        name="Electricity cost",
        ui_units=fs.costing.base_currency / pyunits.kWh,
        display_units="$/kWh",
        rounding=3,
        description="Electricity cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )

    # Outlets
    exports.add(
        obj=fs.product_H2O.properties[0].flow_vol,
        name="Product water flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Outlet product water flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_H2O.properties[0].conc_mass_comp["organic_solid"],
        name="Product water organics(s) concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=2,
        description="Outlet product water organics(s) concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_H2O.properties[0].conc_mass_comp["organic_liquid"],
        name="Product water organics(l) concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=2,
        description="Outlet product water organics(l) concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_H2O.properties[0].conc_mass_comp["inorganic_solid"],
        name="Product water inorganics(s) concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=2,
        description="Outlet product water inorganics(s) concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )

    # System metrics
    exports.add(
        obj=fs.costing.LCOT,
        name="Levelized cost of treatment",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of feed",
        rounding=2,
        description="Levelized cost of treatment including operating and capital costs",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.LCOW,
        name="Levelized cost of water",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of product",
        rounding=2,
        description="Levelized cost of water including operating and capital costs",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )

    # Normalized metrics
    total_capital_norm = fs.costing.total_capital_cost / fs.feed.properties[0].flow_vol
    exports.add(
        obj=total_capital_norm,
        name="Total capital",
        ui_units=fs.costing.base_currency / (pyunits.m**3 / pyunits.day),
        display_units="$/(m3/day)",
        rounding=1,
        description="Normalized total capital costs accounting for indirect "
        "capital and installation - [total capital costs/feed flow rate]",
        is_input=False,
        is_output=True,
        output_category="Normalized cost metrics",
    )
    direct_capital_norm = (
        fs.ATHTL.costing.direct_capital_cost
        + fs.salt_precipitation.costing.direct_capital_cost
        + fs.HTG.costing.direct_capital_cost
    ) / fs.feed.properties[0].flow_vol
    exports.add(
        obj=direct_capital_norm,
        name="Direct capital",
        ui_units=fs.costing.base_currency / (pyunits.m**3 / pyunits.day),
        display_units="$/(m3/day)",
        rounding=1,
        description="Normalized direct capital costs - [total direct capital "
        "costs/feed flow rate] ",
        is_input=False,
        is_output=True,
        output_category="Normalized cost metrics",
    )
    elec_operating_norm = (
        fs.costing.aggregate_flow_costs["electricity"] / fs.costing.annual_water_inlet
    )
    exports.add(
        obj=elec_operating_norm,
        name="Electricity",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of feed",
        rounding=2,
        description="Normalized electricity cost - [annual electricity costs/annual "
        "feed flow rate]",
        is_input=False,
        is_output=True,
        output_category="Normalized cost metrics",
    )

    # performance metrics
    recovery_vol = (
        fs.product_H2O.properties[0].flow_vol / fs.feed.properties[0].flow_vol
    )
    exports.add(
        obj=recovery_vol,
        name="Volumetric recovery",
        ui_units=pyunits.dimensionless,
        display_units="m3 of product/m3 of feed",
        rounding=3,
        description="Normalized volumetric recovery " "flow rate]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_organics = 1 - (
        fs.product_H2O.properties[0].flow_mass_comp["organic_solid"]
        + fs.product_H2O.properties[0].flow_mass_comp["organic_liquid"]
    ) / (
        fs.feed.properties[0].flow_mass_comp["organic_solid"]
        + fs.feed.properties[0].flow_mass_comp["organic_liquid"]
    )
    exports.add(
        obj=removal_organics,
        name="Organics removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Organics removal fraction [1 - outlet organics flow/inlet organics flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_inorganics = (
        1
        - fs.product_H2O.properties[0].flow_mass_comp["inorganic_solid"]
        / fs.feed.properties[0].flow_mass_comp["inorganic_solid"]
    )
    exports.add(
        obj=removal_inorganics,
        name="Inorganics removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Inorganics removal fraction [1 - outlet inorganics flow/inlet inorganics flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )

    # Capital costs
    exports.add(
        obj=fs.costing.total_capital_cost,
        name="Total",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Total capital costs - including investment factor to account "
        "for indirect capital",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.ATHTL.costing.capital_cost,
        name="AT-HTL",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Autothermal hydrothermal liquefaction (AT-HTL)",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.salt_precipitation.costing.capital_cost,
        name="Supercritical salt precipitation",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Supercritical salt precipitation",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.HTG.costing.capital_cost,
        name="HTG",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Hydrothermal gasification (HTG)",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )

    # Operating costs
    exports.add(
        obj=fs.costing.total_operating_cost,
        name="Total",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Total annual operating costs - including electricity, heating, "
        "and fixed",
        is_input=False,
        is_output=True,
        output_category="Operating costs",
    )
    exports.add(
        obj=fs.costing.aggregate_flow_costs["electricity"],
        name="Electricity",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Annual electricity costs ",
        is_input=False,
        is_output=True,
        output_category="Operating costs",
    )
    exports.add(
        obj=fs.costing.aggregate_flow_costs["catalyst_ATHTL"],
        name="Catalyst for AT-HTL",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Annual catalyst cost for AT-HTL",
        is_input=False,
        is_output=True,
        output_category="Operating costs",
    )
    exports.add(
        obj=fs.costing.aggregate_flow_costs["catalyst_HTG"],
        name="Catalyst for HTG",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Annual catalyst cost for HTG",
        is_input=False,
        is_output=True,
        output_category="Operating costs",
    )
    exports.add(
        obj=fs.costing.total_fixed_operating_cost,
        name="Fixed",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Annual fixed operating costs - these costs include material "
        "replacement, maintenance, and labor",
        is_input=False,
        is_output=True,
        output_category="Operating costs",
    )


def build_flowsheet():
    # build and solve initial flowsheet
    m = build()

    set_operating_conditions(m)
    assert_degrees_of_freedom(m, 0)
    assert_units_consistent(m)

    initialize_system(m)

    results = solve(m)
    assert_optimal_termination(results)

    add_costing(m)
    assert_degrees_of_freedom(m, 0)
    m.fs.costing.initialize()

    results = solve(m)
    assert_optimal_termination(results)
    return m.fs


def solve_flowsheet(flowsheet=None):
    fs = flowsheet
    results = solve(fs)
    return results
