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
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.groundwater_treatment.groundwater_treatment import (
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
        name="Groundwater treatment",
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
        obj=fs.feed.conc_mass_comp[0, "arsenic"],
        name="Arsenic concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=2,
        description="Inlet arsenic concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "uranium"],
        name="Uranium concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=2,
        description="Inlet uranium concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "nitrate"],
        name="Nitrate concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=2,
        description="Inlet nitrate concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "phosphates"],
        name="Phosphates concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=2,
        description="Inlet phosphates concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "iron"],
        name="Iron concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=2,
        description="Inlet iron concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "filtration_media"],
        name="Filtration media concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=2,
        description="Inlet filtration media concentration",
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
        obj=fs.costing.pump_electricity.pump_cost[None],
        name="Pump cost",
        ui_units=fs.costing.base_currency / (pyunits.m**3 / pyunits.hr),
        display_units="$/(m^3/hr)",
        rounding=0,
        description="Pump capital cost parameter",
        is_input=True,
        input_category="Pump costing",
        is_output=False,
    )

    # Unit model data, microbial battery
    exports.add(
        obj=fs.micbatt.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Microbial battery",
        is_output=False,
    )
    exports.add(
        obj=fs.micbatt.reaction_conversion[0, "iron_precipitation"],
        name="Iron conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Iron conversion [g-Fe reacted/g-Fe inlet]",
        is_input=True,
        input_category="Microbial battery",
        is_output=False,
    )
    exports.add(
        obj=fs.micbatt.removal_frac_mass_comp[0, "arsenic"],
        name="Arsenic removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Arsenic removal [g-As removed/g-As inlet]",
        is_input=True,
        input_category="Microbial battery",
        is_output=False,
    )
    exports.add(
        obj=fs.micbatt.removal_frac_mass_comp[0, "uranium"],
        name="Uranium removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Uranium removal [g-U removed/g-U inlet]",
        is_input=True,
        input_category="Microbial battery",
        is_output=False,
    )
    exports.add(
        obj=fs.micbatt.removal_frac_mass_comp[0, "nitrate"],
        name="Nitrate removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Nitrate removal [g-nitrate removed/g-nitrate inlet]",
        is_input=True,
        input_category="Microbial battery",
        is_output=False,
    )
    exports.add(
        obj=fs.micbatt.removal_frac_mass_comp[0, "phosphates"],
        name="Phosphates removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Phosphates removal [g-P removed/g-P inlet]",
        is_input=True,
        input_category="Microbial battery",
        is_output=False,
    )
    exports.add(
        obj=fs.micbatt.energy_electric_flow_vol_inlet,
        name="Electricity intensity",
        ui_units=pyunits.kWh / pyunits.m**3,
        display_units="kWh/m3",
        rounding=2,
        description="Electricity intensity relating the intensity to the product volume",
        is_input=True,
        input_category="Microbial battery",
        is_output=False,
    )
    exports.add(
        obj=fs.micbatt.HRT,
        name="HRT",
        ui_units=pyunits.hr,
        display_units="h",
        rounding=1,
        description="Hydraulic retention time",
        is_input=True,
        input_category="Microbial battery",
        is_output=False,
    )
    # Unit cost data, microbial battery
    exports.add(
        obj=fs.costing.microbial_battery.sizing_cost[None],
        name="Microbial battery cost",
        ui_units=(fs.costing.base_currency * pyunits.day) / pyunits.m**3,
        display_units="($ * day)/m3 of microbial battery",
        rounding=0,
        description="Microbial battery capital cost parameter",
        is_input=True,
        input_category="Microbial battery costing",
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
        obj=fs.costing.filtration_media_cost,
        name="Filtration media cost",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=3,
        description="Filtration media cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.filtration_media_disposal_cost,
        name="Filtration media disposal cost",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=3,
        description="Filtration media disposal cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    # Outlets
    exports.add(
        obj=fs.filtered_water.properties[0].flow_vol,
        name="Filtered water flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Outlet filtered water flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.filtered_water.properties[0].conc_mass_comp["arsenic"],
        name="Filtered water arsenic concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=3,
        description="Outlet filtered water arsenic concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.filtered_water.properties[0].conc_mass_comp["uranium"],
        name="Filtered water uranium concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=3,
        description="Outlet filtered water uranium concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.filtered_water.properties[0].conc_mass_comp["nitrate"],
        name="Filtered water nitrate concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=3,
        description="Outlet filtered water nitrate concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.filtered_water.properties[0].conc_mass_comp["phosphates"],
        name="Filtered water phosphates concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=3,
        description="Outlet filtered water phosphates concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.byproduct.properties[0].flow_mass_comp["arsenic"],
        name="Byproduct arsenic flow rate",
        ui_units=pyunits.mg / pyunits.day,
        display_units="mg-arsenic/day",
        rounding=5,
        description="Outlet byproduct arsenic flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.byproduct.properties[0].flow_mass_comp["uranium"],
        name="Byproduct uranium flow rate",
        ui_units=pyunits.mg / pyunits.day,
        display_units="mg-uranium/day",
        rounding=5,
        description="Outlet byproduct uranium flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.byproduct.properties[0].flow_mass_comp["nitrate"],
        name="Byproduct nitrate flow rate",
        ui_units=pyunits.mg / pyunits.day,
        display_units="mg-nitrate/day",
        rounding=5,
        description="Outlet byproduct nitrate flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.byproduct.properties[0].flow_mass_comp["phosphates"],
        name="Byproduct phosphates flow rate",
        ui_units=pyunits.mg / pyunits.day,
        display_units="mg-phosphates/day",
        rounding=5,
        description="Outlet byproduct phosphates flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.byproduct.properties[0].flow_mass_comp["filtration_media"],
        name="Byproduct filtration media flow rate",
        ui_units=pyunits.mg / pyunits.day,
        display_units="mg-media/day",
        rounding=5,
        description="Outlet byproduct filtration media flow rate",
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
        obj=fs.costing.LCOW,
        name="LCOW",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of treated water",
        rounding=3,
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
        (fs.pump.costing.capital_cost + fs.micbatt.costing.capital_cost)
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
        fs.filtered_water.properties[0].flow_vol / fs.feed.properties[0].flow_vol
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
    removal_arsenic = (
        1
        - fs.filtered_water.properties[0].flow_mass_comp["arsenic"]
        / fs.feed.properties[0].flow_mass_comp["arsenic"]
    )
    exports.add(
        obj=removal_arsenic,
        name="Arsenic removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Arsenic removal fraction [1 - outlet As flow/inlet As flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_uranium = (
        1
        - fs.filtered_water.properties[0].flow_mass_comp["uranium"]
        / fs.feed.properties[0].flow_mass_comp["uranium"]
    )
    exports.add(
        obj=removal_uranium,
        name="Uranium removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Uranium removal fraction [1 - outlet U flow/inlet U flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_nitrate = (
        1
        - fs.filtered_water.properties[0].flow_mass_comp["nitrate"]
        / fs.feed.properties[0].flow_mass_comp["nitrate"]
    )
    exports.add(
        obj=removal_nitrate,
        name="Nitrate removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Nitrate removal fraction [1 - outlet NO3 flow/inlet NO3 flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_phosphates = (
        1
        - fs.filtered_water.properties[0].flow_mass_comp["phosphates"]
        / fs.feed.properties[0].flow_mass_comp["phosphates"]
    )
    exports.add(
        obj=removal_phosphates,
        name="Phosphates removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Phosphates removal fraction [1 - outlet P flow/inlet P flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_iron = (
        1
        - fs.filtered_water.properties[0].flow_mass_comp["iron"]
        / fs.feed.properties[0].flow_mass_comp["iron"]
    )
    exports.add(
        obj=removal_iron,
        name="Iron removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Iron removal fraction [1 - outlet Fe flow/inlet Fe flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_filtration_media = (
        1
        - fs.filtered_water.properties[0].flow_mass_comp["filtration_media"]
        / fs.feed.properties[0].flow_mass_comp["filtration_media"]
    )
    exports.add(
        obj=removal_filtration_media,
        name="Filtration media removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Filtration media removal fraction [1 - outlet media flow/inlet media flow]",
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
        obj=fs.micbatt.costing.capital_cost,
        name="Microbial battery",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Microbial battery",
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
