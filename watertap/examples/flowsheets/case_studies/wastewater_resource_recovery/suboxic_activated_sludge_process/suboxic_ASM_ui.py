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
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.suboxic_activated_sludge_process.suboxic_activated_sludge_process import (
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
        name="Suboxic ASM",
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
        description="Inlet biological oxygen demand (BOD) concentration",
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
        obj=fs.feed.conc_mass_comp[0, "tkn"],
        name="TKN concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=1,
        description="Inlet total Kjeldahl nitrogen (TKN) concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "phosphorus"],
        name="Total P concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=1,
        description="Inlet total phosphorus (TP) concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )

    # Unit model data, suboxic ASM
    exports.add(
        obj=fs.suboxicASM.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Suboxic ASM",
        is_output=False,
    )
    exports.add(
        obj=fs.suboxicASM.removal_frac_mass_comp[0, "bod"],
        name="BOD removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=4,
        description="BOD removal [g-BOD reacted/g-BOD byproduct]",
        is_input=True,
        input_category="Suboxic ASM",
        is_output=False,
    )
    exports.add(
        obj=fs.suboxicASM.removal_frac_mass_comp[0, "tss"],
        name="TSS removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=4,
        description="TSS removal [g-TSS reacted/g-TSS byproduct]",
        is_input=True,
        input_category="Suboxic ASM",
        is_output=False,
    )
    exports.add(
        obj=fs.suboxicASM.removal_frac_mass_comp[0, "tkn"],
        name="TKN removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=4,
        description="TKN removal [g-TKN reacted/g-TKN byproduct]",
        is_input=True,
        input_category="Suboxic ASM",
        is_output=False,
    )
    exports.add(
        obj=fs.suboxicASM.removal_frac_mass_comp[0, "phosphorus"],
        name="TP removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=4,
        description="TP removal [g-TP reacted/g-TP byproduct]",
        is_input=True,
        input_category="Suboxic ASM",
        is_output=False,
    )
    exports.add(
        obj=fs.suboxicASM.energy_electric_flow_vol_inlet,
        name="Electricity specific power",
        ui_units=pyunits.kWh / pyunits.m**3,
        display_units="kWh/m3",
        rounding=1,
        description="Specific power relating the power to the inlet flow rate",
        is_input=True,
        input_category="Suboxic ASM",
        is_output=False,
    )

    # Unit cost data, suboxic ASM
    exports.add(
        obj=fs.costing.suboxic_activated_sludge_process.aeration_basin_cost[None],
        name="Aeration basin cost",
        ui_units=fs.costing.base_currency / (pyunits.m**3 / pyunits.hr),
        display_units="$/(m3/h)",
        rounding=1,
        description="Aeration basin capital cost parameter",
        is_input=True,
        input_category="Suboxic ASM costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.suboxic_activated_sludge_process.other_equipment_cost[None],
        name="Other equipment cost",
        ui_units=fs.costing.base_currency / (pyunits.m**3 / pyunits.hr),
        display_units="$/(m3/h)",
        rounding=1,
        description="Other equipment capital cost parameter",
        is_input=True,
        input_category="Suboxic ASM costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.suboxic_activated_sludge_process.control_system_cost[None],
        name="Control system cost",
        ui_units=fs.costing.base_currency / (pyunits.m**3 / pyunits.hr),
        display_units="$/(m3/h)",
        rounding=1,
        description="Control system capital cost parameter",
        is_input=True,
        input_category="Suboxic ASM costing",
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
        rounding=2,
        description="Outlet product water biological oxygen demand (BOD) concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_H2O.properties[0].conc_mass_comp["tss"],
        name="Product water TSS concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=2,
        description="Outlet product water total suspended solids (TSS) concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_H2O.properties[0].conc_mass_comp["tkn"],
        name="Product water TKN concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=2,
        description="Outlet product water total Kjeldahl nitrogen (TKN) concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_H2O.properties[0].conc_mass_comp["phosphorus"],
        name="Product water TP concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=2,
        description="Outlet product water total phosphorus (TP) concentration",
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
        fs.suboxicASM.costing.direct_capital_cost / fs.feed.properties[0].flow_vol
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
        rounding=4,
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
        rounding=4,
        description="Normalized volumetric recovery" "flow rate]",
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
        rounding=4,
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
        rounding=4,
        description="TSS removal fraction [1 - outlet TSS flow/inlet TSS flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_tkn = (
        1
        - fs.product_H2O.properties[0].flow_mass_comp["tkn"]
        / fs.feed.properties[0].flow_mass_comp["tkn"]
    )
    exports.add(
        obj=removal_tkn,
        name="TKN removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=4,
        description="TKN removal fraction [1 - outlet TKN flow/inlet TKN flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_tp = (
        1
        - fs.product_H2O.properties[0].flow_mass_comp["phosphorus"]
        / fs.feed.properties[0].flow_mass_comp["phosphorus"]
    )
    exports.add(
        obj=removal_tp,
        name="TP removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=4,
        description="TP removal fraction [1 - outlet TP flow/inlet TP flow]",
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
        obj=fs.suboxicASM.costing.capital_cost,
        name="Suboxic ASM",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Suboxic ASM",
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
