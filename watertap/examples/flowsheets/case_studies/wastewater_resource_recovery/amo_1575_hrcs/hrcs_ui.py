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
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.amo_1575_hrcs.hrcs import (
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
        name="HRCS",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
    )


def export_variables(flowsheet=None, exports=None):
    fs = flowsheet

    def _base_curr(x):
        return pyunits.convert(x, to_units=fs.costing.base_currency)

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
        obj=fs.feed.conc_mass_comp[0, "tss"],
        name="TSS concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=2,
        description="Inlet total suspended solids concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.properties[0].flow_mass_comp["cod"],
        name="COD mass flow",
        ui_units=pyunits.ton / pyunits.day,
        display_units="ton/day",
        rounding=2,
        description="Inlet COD mass flow",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.properties[0].flow_mass_comp["oxygen"],
        name="O2 mass flow",
        ui_units=pyunits.ton / pyunits.day,
        display_units="ton/day",
        rounding=2,
        description="Inlet oxygen mass flow",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "carbon_dioxide"],
        name="CO2 concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=2,
        description="Inlet carbon dioxide concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    # Unit model data, HRCS
    exports.add(
        obj=fs.HRCS.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="HRCS",
        is_output=False,
    )
    exports.add(
        obj=fs.HRCS.reaction_conversion[0, "oxidation"],
        name="COD conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="COD conversion [g-COD reacted/g-COD inlet]",
        is_input=True,
        input_category="HRCS",
        is_output=False,
    )
    exports.add(
        obj=fs.HRCS.removal_frac_mass_comp[0, "tss"],
        name="TSS removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Total suspended solids removal [g-TSS removed/g-TSS inlet]",
        is_input=True,
        input_category="HRCS",
        is_output=False,
    )
    exports.add(
        obj=fs.HRCS.energy_electric_flow_vol_inlet,
        name="Aeration energy",
        ui_units=pyunits.kWh / pyunits.m**3,
        display_units="kWh/m3",
        rounding=2,
        description="Specific aeration energy relating the energy to the volume of treated water",
        is_input=True,
        input_category="HRCS",
        is_output=False,
    )

    # Unit cost data, HRCS
    exports.add(
        obj=fs.costing.hrcs.SRT[None],
        name="SRT",
        ui_units=pyunits.hr,
        display_units="h",
        rounding=1,
        description="Solids retention time",
        is_input=True,
        input_category="HRCS costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.hrcs.sizing_cost[None],
        name="HRCS cost",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of HRCS",
        rounding=0,
        description="HRCS capital cost parameter",
        is_input=True,
        input_category="HRCS costing",
        is_output=False,
    )

    # Unit model data, clarifier
    exports.add(
        obj=fs.clarifier.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.clarifier.removal_frac_mass_comp[0, "tss"],
        name="TSS removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=3,
        description="Total suspended solids removal [g-TSS removed/g-TSS inlet]",
        is_input=True,
        input_category="Clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.clarifier.removal_frac_mass_comp[0, "cod"],
        name="COD removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=3,
        description="Total COD [g-COD removed/g-COD inlet]",
        is_input=True,
        input_category="Clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.clarifier.energy_electric_flow_vol_inlet,
        name="Electricity intensity",
        ui_units=pyunits.kWh / pyunits.m**3,
        display_units="kWh/m3",
        rounding=2,
        description="Clarifier electricity intensity relating the energy to the volume of treated water",
        is_input=True,
        input_category="Clarifier",
        is_output=False,
    )
    exports.add(
        obj=fs.clarifier.ferric_chloride_dose[0],
        name="FeCl3 dosage",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg FeCl3/L sludge treated",
        rounding=2,
        description="Ferric chloride dosage per L of sludge treated",
        is_input=True,
        input_category="Clarifier",
        is_output=False,
    )
    # Unit cost data, clarifier
    exports.add(
        obj=fs.costing.clarifier.HRT["HRCS_clarifier"],
        name="HRT",
        ui_units=pyunits.hr,
        display_units="h",
        rounding=1,
        description="Hydraulic retention time",
        is_input=True,
        input_category="Clarifier costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.clarifier.sizing_cost["HRCS_clarifier"],
        name="Clarifier cost",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of clarifier",
        rounding=0,
        description="Clarifier capital cost parameter",
        is_input=True,
        input_category="Clarifier costing",
        is_output=False,
    )

    # Unit model data, primary separator
    exports.add(
        obj=fs.sep.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Water recovery in purge stream [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Separator",
        is_output=False,
    )
    exports.add(
        obj=fs.sep.removal_frac_mass_comp[0, "tss"],
        name="TSS removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=3,
        description="Total suspended solids removal to recycle [g-TSS removed/g-TSS inlet]",
        is_input=True,
        input_category="Separator",
        is_output=False,
    )
    exports.add(
        obj=fs.sep.removal_frac_mass_comp[0, "cod"],
        name="COD removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=3,
        description="COD removal to recycle [g-COD removed/g-COD inlet]",
        is_input=True,
        input_category="Separator",
        is_output=False,
    )
    exports.add(
        obj=fs.sep.energy_electric_flow_vol_inlet,
        name="Electricity intensity",
        ui_units=pyunits.kWh / pyunits.m**3,
        display_units="kWh/m3",
        rounding=2,
        description="Separator electricity intensity relating the energy to the volume of treated water",
        is_input=True,
        input_category="Separator",
        is_output=False,
    )

    # Unit cost data, separator - uses default costing parameters

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
        obj=fs.costing.ferric_chloride_cost,
        name="FeCl3 cost",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=3,
        description="Ferric chloride cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )

    # Outlets
    exports.add(
        obj=fs.product.properties[0].flow_vol,
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
        obj=fs.product.properties[0].conc_mass_comp["tss"],
        name="Product water TSS concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=5,
        description="Outlet product water total suspended solids concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product.properties[0].flow_mass_comp["cod"],
        name="Product water COD mass flow",
        ui_units=pyunits.ton / pyunits.day,
        display_units="ton/day",
        rounding=5,
        description="Outlet product water COD mass flow",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.WAS_product.properties[0].flow_vol,
        name="WAS flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Outlet waste activated sludge flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.WAS_product.properties[0].conc_mass_comp["H2O"],
        name="WAS H2O concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=5,
        description="Outlet waste activated sludge water concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.WAS_product.properties[0].conc_mass_comp["tss"],
        name="WAS TSS concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=5,
        description="Outlet waste activated sludge total suspended solids concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.WAS_product.properties[0].flow_mass_comp["cod"],
        name="WAS COD mass flow",
        ui_units=pyunits.ton / pyunits.day,
        display_units="ton/day",
        rounding=5,
        description="Outlet waste activated sludge COD mass flow",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.disposal.properties[0].flow_vol,
        name="Gaseous release flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Outlet gaseous release flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.disposal.properties[0].flow_mass_comp["carbon_dioxide"],
        name="Gaseous CO2 flow rate",
        ui_units=pyunits.ton / pyunits.day,
        display_units="ton/day",
        rounding=5,
        description="Outlet gaseous carbon dioxide release mass flow",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.disposal.properties[0].flow_mass_comp["oxygen"],
        name="Gaseous O2 flow rate",
        ui_units=pyunits.ton / pyunits.day,
        display_units="ton/day",
        rounding=5,
        description="Outlet gaseous oxygen release mass flow",
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
        description="Levelized cost of treatment including consumption costs for FeCl3",
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
        description="Levelized cost of water including consumption costs for FeCl3",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )

    # Normalized metrics
    total_capital_norm = (
        fs.costing.total_capital_cost
        + _base_curr(fs.watertap_costing.total_capital_cost)
    ) / fs.feed.properties[0].flow_vol
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
        fs.HRCS.costing.direct_capital_cost
        + fs.clarifier.costing.direct_capital_cost
        + fs.sep.costing.direct_capital_cost
        # In WaterTAPCosting package `capital_cost` doesn't have any adders
        + _base_curr(fs.mixer.costing.capital_cost)
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
    recovery_vol = fs.product.properties[0].flow_vol / fs.feed.properties[0].flow_vol
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
    carbon_capture = (
        100
        * (
            fs.WAS_product.flow_mass_comp[0, "cod"]
            + fs.product.flow_mass_comp[0, "cod"]
        )
        / fs.feed.flow_mass_comp[0, "cod"]
    )
    exports.add(
        obj=carbon_capture,
        name="Carbon capture",
        ui_units=pyunits.dimensionless,
        display_units="%",
        rounding=3,
        description="Normalized carbon capture (WAS captured COD + effluent biomass COD + effluent non-biomass COD",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_COD = (
        1
        - fs.product.properties[0].flow_mass_comp["cod"]
        / fs.feed.properties[0].flow_mass_comp["cod"]
    )
    exports.add(
        obj=removal_COD,
        name="COD removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="COD removal fraction [1 - outlet COD flow/inlet COD flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_TSS = (
        1
        - fs.product.properties[0].flow_mass_comp["tss"]
        / fs.feed.properties[0].flow_mass_comp["tss"]
    )
    exports.add(
        obj=removal_TSS,
        name="TSS removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="TSS removal fraction [1 - outlet TSS flow/inlet TSS flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )

    # Capital costs
    total_capital = fs.costing.total_capital_cost + _base_curr(
        fs.watertap_costing.total_capital_cost
    )
    exports.add(
        obj=total_capital,
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
        obj=fs.HRCS.costing.capital_cost,
        name="HRCS",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="HRCS",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.clarifier.costing.capital_cost,
        name="Clarifier",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Clarifier",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.sep.costing.capital_cost,
        name="Separator",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Separator",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=_base_curr(fs.mixer.costing.capital_cost),
        name="Mixer",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Mixer",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )

    # Operating costs
    total_operating = (
        fs.costing.total_fixed_operating_cost
        + fs.costing.total_variable_operating_cost
        + pyunits.convert(
            fs.watertap_costing.total_operating_cost,
            to_units=fs.costing.base_currency / pyunits.year,
        )
    )
    exports.add(
        obj=total_operating,
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
