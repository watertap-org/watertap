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
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.electrochemical_nutrient_removal.electrochemical_nutrient_removal import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    add_costing,
)
from idaes.core.solvers import get_solver
from pyomo.environ import units as pyunits, assert_optimal_termination
from pyomo.util.check_units import assert_units_consistent


def export_to_ui():
    return FlowsheetInterface(
        name="Electrochemical nutrient removal",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
    )


def export_variables(flowsheet=None, exports=None):
    fs = flowsheet
    # --- Input data ---
    # Feed conditions
    exports.add(
        obj=fs.feed.flow_vol[0],
        name="Volumetric flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Inlet volumetric flow rate",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "nitrogen"],
        name="Nitrogen concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=2,
        description="Inlet nitrogen concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "struvite"],
        name="Struvite concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=2,
        description="Inlet struvite concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "phosphorus"],
        name="Phosphorus concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=2,
        description="Inlet phosphorus concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )

    # Unit model data, pump
    exports.add(
        obj=fs.pump.lift_height,
        name="Lift height",
        ui_units=pyunits.m,
        display_units="m",
        rounding=2,
        description="Lift height for pump",
        is_input=True,
        input_category="Pump",
        is_output=False,
    )
    exports.add(
        obj=fs.pump.eta_pump,
        name="Pump efficiency",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Efficiency of pump",
        is_input=True,
        input_category="Pump",
        is_output=False,
    )
    exports.add(
        obj=fs.pump.eta_motor,
        name="Motor efficiency",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Efficiency of motor",
        is_input=True,
        input_category="Pump",
        is_output=False,
    )

    # Unit cost data, pump
    exports.add(
        obj=fs.costing.pump_electricity.pump_cost["default"],
        name="Pump cost",
        ui_units=fs.costing.base_currency / (pyunits.m**3 / pyunits.hr),
        display_units="$/(m^3/hr)",
        rounding=0,
        description="Pump capital cost parameter",
        is_input=True,
        input_category="Pump costing",
        is_output=False,
    )
    # Unit model data, electroNP
    exports.add(
        obj=fs.electroNP.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="ElectroNP",
        is_output=False,
    )
    exports.add(
        obj=fs.electroNP.reaction_conversion[0, "extract_N_P"],
        name="Phosphorus conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Phosphorus conversion [g-P reacted/g-P inlet]",
        is_input=True,
        input_category="ElectroNP",
        is_output=False,
    )
    exports.add(
        obj=fs.electroNP.energy_electric_flow_mass,
        name="Specific energy",
        ui_units=pyunits.kWh / pyunits.kg,
        display_units="kWh/kg of struvite",
        rounding=2,
        description="ElectroNP specific energy relating the energy to the mass of struvite product",
        is_input=True,
        input_category="ElectroNP",
        is_output=False,
    )
    exports.add(
        obj=fs.electroNP.magnesium_chloride_dosage,
        name="MgCl2 Dosage",
        ui_units=pyunits.dimensionless,
        display_units="g MgCl2/g struvite",
        rounding=2,
        description="MgCl2 dosage per g of struvite product",
        is_input=True,
        input_category="ElectroNP",
        is_output=False,
    )

    # Unit cost data, electroNP
    exports.add(
        obj=fs.costing.electrochemical_nutrient_removal.HRT[None],
        name="HRT",
        ui_units=pyunits.hr,
        display_units="h",
        rounding=1,
        description="Hydraulic retention time",
        is_input=True,
        input_category="ElectroNP costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.electrochemical_nutrient_removal.sizing_cost[None],
        name="ElectroNP cost",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of electroNP",
        rounding=0,
        description="ElectroNP capital cost parameter",
        is_input=True,
        input_category="ElectroNP costing",
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
    exports.add(
        obj=fs.costing.magnesium_chloride_cost,
        name="MgCl2 cost",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=3,
        description="Magnesium chloride cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.struvite_product_cost,
        name="Struvite cost",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=3,
        description="Struvite cost is negative because it is sold",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )

    # Outlets
    exports.add(
        obj=fs.product_H2O.properties[0].flow_vol,
        name="Product H2O flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Outlet product H2O flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_H2O.properties[0].conc_mass_comp["nitrogen"],
        name="Product water N2 concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=3,
        description="Outlet product water nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_H2O.properties[0].conc_mass_comp["phosphorus"],
        name="Product water P concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=3,
        description="Outlet product water phosphorus concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_struvite.properties[0].flow_mass_comp["struvite"],
        name="Product struvite flow rate",
        ui_units=pyunits.kg / pyunits.hr,
        display_units="kg-struvite/h",
        rounding=3,
        description="Outlet product struvite flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )

    # System metrics
    exports.add(
        obj=fs.costing.LCOT,
        name="LCOT",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of feed",
        rounding=3,
        description="Levelized cost of treatment including operating and capital costs",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.LCOT_with_revenue,
        name="LCOT with revenue",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of feed",
        rounding=3,
        description="Levelized cost of treatment including revenue of struvite and consumption costs for MgCl2",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.LCOW,
        name="LCOW",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of centrate",
        rounding=3,
        description="Levelized cost of water including operating and capital costs",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.LCOW_with_revenue,
        name="LCOW with revenue",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of centrate",
        rounding=3,
        description="Levelized cost of water including revenue of struvite and consumption costs for MgCl2",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.LCOS,
        name="Levelized cost of struvite",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg-struvite",
        rounding=3,
        description="Levelized cost of struvite including operating and capital costs",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )  # Normalized metrics
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
        (fs.pump.costing.capital_cost + fs.electroNP.costing.capital_cost)
        / fs.costing.TIC
        / fs.feed.properties[0].flow_vol
    )
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
        description="Normalized volumetric recovery",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_nitrogen = (
        1
        - fs.product_H2O.properties[0].flow_mass_comp["nitrogen"]
        / fs.feed.properties[0].flow_mass_comp["nitrogen"]
    )
    exports.add(
        obj=removal_nitrogen,
        name="Nitrogen removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Nitrogen removal fraction [1 - outlet N2 flow/inlet N2 flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_phosphorus = (
        1
        - fs.product_H2O.properties[0].flow_mass_comp["phosphorus"]
        / fs.feed.properties[0].flow_mass_comp["phosphorus"]
    )
    exports.add(
        obj=removal_phosphorus,
        name="Phosphorus removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Phosphorus removal fraction [1 - outlet P flow/inlet P flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    struvite_production = (
        fs.product_struvite.properties[0].flow_mass_comp["struvite"]
        / fs.feed.properties[0].flow_vol
    )
    exports.add(
        obj=struvite_production,
        name="Struvite production",
        ui_units=pyunits.kg / pyunits.m**3,
        display_units="kg-struvite/m3 of feed",
        rounding=3,
        description="Struvite production [Struvite product flow rate/feed flow rate]",
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
        obj=fs.pump.costing.capital_cost,
        name="Pump",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Pump",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.electroNP.costing.capital_cost,
        name="ElectroNP",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="ElectroNP",
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

    # Revenue
    total_revenue = -(
        fs.costing.aggregate_flow_costs["struvite_product"]
        + fs.costing.aggregate_flow_costs["magnesium_chloride"]
    )
    exports.add(
        obj=total_revenue,
        name="Total",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Total revenue - only includes the purchase of MgCl2 (not struvite revenue)",
        is_input=False,
        is_output=True,
        output_category="Revenue",
    )
    exports.add(
        obj=-fs.costing.aggregate_flow_costs["struvite_product"],
        name="Struvite",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Revenue from selling struvite",
        is_input=False,
        is_output=True,
        output_category="Revenue",
    )
    exports.add(
        obj=-fs.costing.aggregate_flow_costs["magnesium_chloride"],
        name="Magnesium chloride",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Cost from buying magnesium chloride",
        is_input=False,
        is_output=True,
        output_category="Revenue",
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
