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
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.peracetic_acid_disinfection.peracetic_acid_disinfection import (
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
        name="Peracetic acid disinfection",
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
        obj=fs.feed.conc_mass_comp[0, "peracetic_acid"],
        name="Peracetic acid concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=2,
        description="Inlet peracetic acid concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "total_coliforms_fecal_ecoli"],
        name="Ecoli concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=7,
        description="Inlet total coliforms fecal ecoli concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )

    # Unit model data, peracetic acid disinfection
    exports.add(
        obj=fs.PAA.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Peracetic acid disinfection",
        is_output=False,
    )
    exports.add(
        obj=fs.PAA.reaction_conversion[0, "paa_decomposition"],
        name="PAA decomposition",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Peracetic acid decomposition [g-PAA decomposed/g-PAA inlet]",
        is_input=True,
        input_category="Peracetic acid disinfection",
        is_output=False,
    )
    exports.add(
        obj=fs.PAA.reaction_conversion[0, "ecoli_inactivation"],
        name="Ecoli inactivation",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Ecoli (total coliforms fecal ecoli) inactivation [g-ecoli inactivated/g-ecoli inlet]",
        is_input=True,
        input_category="Peracetic acid disinfection",
        is_output=False,
    )
    exports.add(
        obj=fs.PAA.energy_electric_flow_vol_inlet,
        name="Electricity intensity",
        ui_units=pyunits.kWh / pyunits.m**3,
        display_units="fraction",
        rounding=2,
        description="Electricity intensity with respect to inlet flow rate of unit",
        is_input=True,
        input_category="Peracetic acid disinfection",
        is_output=False,
    )
    exports.add(
        obj=fs.PAA.ecoli_cell_mass,
        name="Ecoli mass",
        ui_units=pyunits.kg,
        display_units="kg",
        rounding=16,
        description="Ecoli cell mass",
        is_input=True,
        input_category="Peracetic acid disinfection",
        is_output=False,
    )
    exports.add(
        obj=fs.PAA.disinfection_solution_wt_frac_PAA,
        name="PAA weight fraction of disinfection solution",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Disinfection solution weight fraction of peracetic acid",
        is_input=True,
        input_category="Peracetic acid disinfection",
        is_output=False,
    )
    exports.add(
        obj=fs.PAA.disinfection_solution_density,
        name="Disinfection solution density",
        ui_units=pyunits.kg / pyunits.L,
        display_units="kg/L",
        rounding=3,
        description="Disinfection solution density",
        is_input=True,
        input_category="Peracetic acid disinfection",
        is_output=False,
    )
    exports.add(
        obj=fs.PAA.HRT,
        name="HRT",
        ui_units=pyunits.hr,
        display_units="h",
        rounding=0,
        description="Hydraulic retention time",
        is_input=True,
        input_category="Peracetic acid disinfection",
        is_output=False,
    )

    # Unit cost data, peracetic acid disinfection
    exports.add(
        obj=fs.costing.peracetic_acid_disinfection.sizing_cost[None],
        name="Peracetic acid disinfection cost",
        ui_units=(fs.costing.base_currency * pyunits.day) / pyunits.m**3,
        display_units="($*day)/m3 of PAA",
        rounding=0,
        description="Peracetic acid disinfection capital cost parameter",
        is_input=True,
        input_category="Peracetic acid disinfection costing",
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
        obj=fs.costing.disinfection_solution_cost,
        name="Disinfection solution cost",
        ui_units=fs.costing.base_currency / pyunits.gal,
        display_units="$/gal",
        rounding=3,
        description="Disinfection cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )

    # Outlets
    exports.add(
        obj=fs.treated_water.properties[0].flow_vol,
        name="Treated water flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Outlet treated water flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.treated_water.properties[0].conc_mass_comp["peracetic_acid"],
        name="Treated water PAA concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=5,
        description="Outlet treated water peracetic acid concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.treated_water.properties[0].conc_mass_comp[
            "total_coliforms_fecal_ecoli"
        ],
        name="Treated water Ecoli concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=7,
        description="Outlet treated water ecoli concentration",
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
    direct_capital_norm = (fs.PAA.costing.direct_capital_cost) / fs.feed.properties[
        0
    ].flow_vol
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
    #
    # performance metrics
    recovery_vol = (
        fs.treated_water.properties[0].flow_vol / fs.feed.properties[0].flow_vol
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
    removal_PAA = (
        1
        - fs.treated_water.properties[0].flow_mass_comp["peracetic_acid"]
        / fs.feed.properties[0].flow_mass_comp["peracetic_acid"]
    )
    exports.add(
        obj=removal_PAA,
        name="PAA removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Peracetic acid removal fraction [1 - outlet PAA flow/inlet PAA flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_Ecoli = (
        1
        - fs.treated_water.properties[0].flow_mass_comp["total_coliforms_fecal_ecoli"]
        / fs.feed.properties[0].flow_mass_comp["total_coliforms_fecal_ecoli"]
    )
    exports.add(
        obj=removal_Ecoli,
        name="Ecoli removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Ecoli removal fraction [1 - outlet Ecoli flow/inlet Ecoli flow]",
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
        obj=fs.PAA.costing.capital_cost,
        name="Peracetic acid disinfection",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Peracetic acid disinfection",
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
