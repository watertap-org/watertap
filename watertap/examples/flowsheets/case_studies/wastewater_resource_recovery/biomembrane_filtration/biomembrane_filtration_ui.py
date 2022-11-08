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
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.biomembrane_filtration.biomembrane_filtration import (
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
        name="Biomembrane filtration",
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
        rounding=1,
        description="Inlet volumetric flow rate",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "bod"],
        name="BOD concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=1,
        description="Inlet biochemical oxygen demand (BOD) concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "tss"],
        name="TSS concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=1,
        description="Inlet total suspended solids (TSS) concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "ammonium_as_nitrogen"],
        name="NH4-N concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=1,
        description="Inlet ammonium as nitrogen (NH4-N) concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "nitrate"],
        name="NO3-N concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=1,
        description="Inlet nitrate (NO3-N) concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )

    # Unit model data, mabr
    exports.add(
        obj=fs.mabr.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="MABR",
        is_output=False,
    )
    exports.add(
        obj=fs.mabr.removal_frac_mass_comp[0, "tss"],
        name="TSS removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="TSS removal [g-TSS byproduct/g-TSS inlet]",
        is_input=True,
        input_category="MABR",
        is_output=False,
    )
    exports.add(
        obj=fs.mabr.reaction_conversion[0, "ammonium_to_nitrate"],
        name="NH4-N conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="NH4-N conversion [g-NH4-N reacted/g-NH4-N inlet]",
        is_input=True,
        input_category="MABR",
        is_output=False,
    )
    exports.add(
        obj=fs.mabr.generation_ratio["ammonium_to_nitrate", "bod"],
        name="BOD conversion ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="BOD mass conversion ratio with respect to NH4-N [g-BOD produced/g-NH4-N"
        " reacted]",
        is_input=True,
        input_category="MABR",
        is_output=False,
    )
    exports.add(
        obj=fs.mabr.generation_ratio["ammonium_to_nitrate", "nitrate"],
        name="NO3-N conversion ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="NO3-N mass conversion ratio with respect to NH4-N [g-NO3-N produced/g-NH4-N"
        " reacted]",
        is_input=True,
        input_category="MABR",
        is_output=False,
    )
    exports.add(
        obj=fs.mabr.nitrogen_removal_rate,
        name="Nitrogen removal rate",
        ui_units=pyunits.g / pyunits.m**2 / pyunits.day,
        display_units="g/m2/day",
        rounding=1,
        description="Nitrogen removal rate",
        is_input=True,
        input_category="MABR",
        is_output=False,
    )
    exports.add(
        obj=fs.mabr.air_flow_rate[0],
        name="Air flow rate",
        ui_units=pyunits.m**3 / pyunits.hr / pyunits.m**2,
        display_units="m3 of air/h/m2",
        rounding=3,
        description="Air flow rate",
        is_input=True,
        input_category="MABR",
        is_output=False,
    )
    exports.add(
        obj=fs.mabr.energy_electric_flow_vol_inlet,
        name="Blower specific power",
        ui_units=pyunits.kWh / pyunits.m**3,
        display_units="kWh/m3 of air",
        rounding=5,
        description="Blower specific power relating the power to the volume of the "
        "air",
        is_input=True,
        input_category="MABR",
        is_output=False,
    )

    # Unit cost data, mabr
    exports.add(
        obj=fs.costing.mabr.reactor_cost[None],
        name="Reactor cost",
        ui_units=fs.costing.base_currency / pyunits.m**2,
        display_units="$/m2 of reactor",
        rounding=0,
        description="Reactor capital cost parameter",
        is_input=True,
        input_category="MABR costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.mabr.blower_cost[None],
        name="Blower cost",
        ui_units=fs.costing.base_currency / (pyunits.m**3 / pyunits.hr),
        display_units="$/(m3/h)",
        rounding=0,
        description="Blower capital cost parameter",
        is_input=True,
        input_category="MABR costing",
        is_output=False,
    )

    # Unit model data, dmbr
    exports.add(
        obj=fs.dmbr.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="DMBR",
        is_output=False,
    )
    exports.add(
        obj=fs.dmbr.removal_frac_mass_comp[0, "tss"],
        name="TSS removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="TSS removal [g-TSS byproduct/g-TSS inlet]",
        is_input=True,
        input_category="DMBR",
        is_output=False,
    )
    exports.add(
        obj=fs.dmbr.reaction_conversion[0, "nitrate_to_nitrogen"],
        name="NO3-N conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="NO3-N conversion [g-NO3-N reacted/g-NO3-N inlet]",
        is_input=True,
        input_category="DMBR",
        is_output=False,
    )
    exports.add(
        obj=fs.dmbr.generation_ratio["nitrate_to_nitrogen", "nitrogen"],
        name="N2 conversion ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="N2 mass conversion ratio with respect to NO3-N [g-N2 "
        "produced/g-NO3-N reacted]",
        is_input=True,
        input_category="DMBR",
        is_output=False,
    )
    exports.add(
        obj=fs.dmbr.reaction_conversion[0, "BOD_usage"],
        name="BOD conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="BOD conversion [g-BOD reacted/g-BOD inlet]",
        is_input=True,
        input_category="DMBR",
        is_output=False,
    )
    exports.add(
        obj=fs.dmbr.energy_electric_flow_vol_inlet,
        name="Electricity specific power",
        ui_units=pyunits.kWh / pyunits.m**3,
        display_units="kWh/m3 of reactor",
        rounding=2,
        description="Electricity specific power relating the power to the volume of the "
        "reactor",
        is_input=True,
        input_category="DMBR",
        is_output=False,
    )

    # Unit cost data, dmbr
    exports.add(
        obj=fs.costing.dmbr.water_flux[None],
        name="Water flux",
        ui_units=pyunits.L / pyunits.m**2 / pyunits.hr,
        display_units="$/m2/h",
        rounding=0,
        description="Reactor sizing parameter - water flux",
        is_input=True,
        input_category="DMBR costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.dmbr.reactor_cost[None],
        name="Reactor cost",
        ui_units=fs.costing.base_currency / pyunits.m**2,
        display_units="$/m2 of reactor",
        rounding=1,
        description="Reactor capital cost parameter",
        is_input=True,
        input_category="DMBR costing",
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
        obj=fs.product_H2O.properties[0].conc_mass_comp["bod"],
        name="Product water BOD concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=4,
        description="Outlet product water biological oxygen demand concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_H2O.properties[0].conc_mass_comp["tss"],
        name="Product water TSS concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=4,
        description="Outlet product water total suspended solids concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_H2O.properties[0].conc_mass_comp["ammonium_as_nitrogen"],
        name="Product water NH4-N concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=4,
        description="Outlet product water NH4-N concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_H2O.properties[0].conc_mass_comp["nitrate"],
        name="Product water NO3-N concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=4,
        description="Outlet product water NO3-N concentration",
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
        fs.mabr.costing.direct_capital_cost
        + fs.dmbr.costing.direct_capital_cost
        + fs.P1.costing.direct_capital_cost
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
        description="Normalized heating cost - [annual heating costs/annual feed "
        "flow rate]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_bod = (
        1
        - fs.product_H2O.properties[0].flow_mass_comp["bod"]
        / fs.feed.properties[0].flow_mass_comp["bod"]
    )
    exports.add(
        obj=removal_bod,
        name="BOD removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="BOD removal fraction [1 - outlet BOD flow/inlet BOD flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_tss = (
        1
        - fs.product_H2O.properties[0].flow_mass_comp["tss"]
        / fs.feed.properties[0].flow_mass_comp["tss"]
    )
    exports.add(
        obj=removal_tss,
        name="TSS removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="TSS removal fraction [1 - outlet TSS flow/inlet TSS flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_TN = 1 - (
        fs.product_H2O.properties[0].flow_mass_comp["ammonium_as_nitrogen"]
        + fs.product_H2O.properties[0].flow_mass_comp["nitrate"]
    ) / (
        fs.feed.properties[0].flow_mass_comp["ammonium_as_nitrogen"]
        + fs.feed.properties[0].flow_mass_comp["nitrate"]
    )
    exports.add(
        obj=removal_TN,
        name="Total N removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Total N removal fraction [1 - outlet total N flow/inlet total N flow]",
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
        obj=fs.mabr.costing.capital_cost,
        name="MABR",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="MABR",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.dmbr.costing.capital_cost,
        name="DMBR",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="DMBR",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.P1.costing.capital_cost,
        name="Pump",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Pump",
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
