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
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.amo_1595_photothermal_membrane_candoP.amo_1595 import (
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
        name="Photothermal Membrane CANDO_P",
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
        rounding=3,
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
        rounding=3,
        description="Inlet nitrogen concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "phosphates"],
        name="Phosphates concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=3,
        description="Inlet phosphates concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "bioconcentrated_phosphorous"],
        name="Bioconcentrated phosphorous concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=3,
        description="Inlet bioconcentrated phosphorous concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "nitrous_oxide"],
        name="Nitrous oxide concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=2,
        description="Inlet nitrous oxide concentration",
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

    # Unit model data, photothermal membrane
    exports.add(
        obj=fs.photothermal.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Photothermal membrane",
        is_output=False,
    )
    exports.add(
        obj=fs.photothermal.removal_frac_mass_comp[0, "phosphates"],
        name="Phosphates removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Phosphates removal [g-P removed/g-P inlet]",
        is_input=True,
        input_category="Photothermal membrane",
        is_output=False,
    )
    exports.add(
        obj=fs.photothermal.removal_frac_mass_comp[0, "nitrogen"],
        name="Nitrogen removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Nitrogen removal [g-N2 removed/g-N2 inlet]",
        is_input=True,
        input_category="Photothermal membrane",
        is_output=False,
    )
    exports.add(
        obj=fs.photothermal.water_flux,
        name="Water flux",
        ui_units=pyunits.kg / pyunits.m**2 / pyunits.hr,
        display_units="kg/(m2/h)",
        rounding=2,
        description="Water flux through membrane",
        is_input=True,
        input_category="Photothermal membrane",
        is_output=False,
    )
    # Unit cost data, photothermal membrane
    exports.add(
        obj=fs.costing.photothermal_membrane.membrane_cost[None],
        name="Membrane",
        ui_units=fs.costing.base_currency / pyunits.m**2,
        display_units="$/m2 of photothermal membrane",
        rounding=0,
        description="Membrane cost",
        is_input=True,
        input_category="Photothermal membrane costing",
        is_output=False,
    )

    # Unit model data, CANDO+P reactor
    exports.add(
        obj=fs.candop.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="CANDO+P reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.candop.reaction_conversion[0, "n_reaction"],
        name="Nitrogen conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Nitrogen conversion [g-N2 reacted/g-N2 inlet]",
        is_input=True,
        input_category="CANDO+P reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.candop.reaction_conversion[0, "p_reaction"],
        name="Phosphorus conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Phosphorus conversion [g-P reacted/g-P inlet]",
        is_input=True,
        input_category="CANDO+P reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.candop.electricity_intensity_N,
        name="Specific energy",
        ui_units=pyunits.kWh / pyunits.kg,
        display_units="kWh/kg",
        rounding=2,
        description="Specific energy relating the energy to the mass of product",
        is_input=True,
        input_category="CANDO+P reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.candop.oxygen_nitrogen_ratio,
        name="O2:N2 ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Oxygen:Nitrogen ratio",
        is_input=True,
        input_category="CANDO+P reactor",
        is_output=False,
    )
    # Unit cost data, CANDO+P reactor
    exports.add(
        obj=fs.costing.CANDO_P.sizing_parameter[None],
        name="HRT",
        ui_units=pyunits.hr,
        display_units="h",
        rounding=1,
        description="Hydraulic retention time",
        is_input=True,
        input_category="CANDO+P reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.CANDO_P.sizing_cost[None],
        name="CANDO+P reactor cost",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of reactor",
        rounding=0,
        description="CANDO+P reactor capital cost parameter",
        is_input=True,
        input_category="CANDO+P reactor costing",
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
        obj=fs.costing.bcp_cost,
        name="Bioconcentrated phosphorous cost",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=3,
        description="Bioconcentrated phosphorous cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.water_cost,
        name="Water cost",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3",
        rounding=13,
        description="Water cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.nitrous_oxide_cost,
        name="Nitrous oxide cost",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=3,
        description="Nitrous oxide cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )

    # Outlets
    exports.add(
        obj=fs.photothermal_water.properties[0].flow_vol,
        name="Photothermal membrane product water flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Photothermal membrane outlet water flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.candop_byproduct.properties[0].flow_mass_comp["nitrous_oxide"],
        name="CANDO_P N2O byproduct flow rate",
        ui_units=pyunits.kg / pyunits.hr,
        display_units="kg-N2O/h",
        rounding=5,
        description="CANDO+P outlet nitrous oxide byproduct flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.candop_treated.properties[0].flow_vol,
        name="CANDO+P product flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="CANDO+P outlet product flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.candop_treated.properties[0].conc_mass_comp["H2O"],
        name="CANDO+P product H2O concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=2,
        description="Outlet product water concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.candop_treated.properties[0].conc_mass_comp["nitrogen"],
        name="CANDO+P product N2 concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=2,
        description="Outlet product nitrogen concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.candop_treated.properties[0].conc_mass_comp["phosphates"],
        name="CANDO+P product phosphates concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=2,
        description="Outlet product phosphates concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.candop_treated.properties[0].conc_mass_comp[
            "bioconcentrated_phosphorous"
        ],
        name="CANDO+P product BCP concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=2,
        description="Outlet product bioconcentrated_phosphorous concentration",
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
        description="Levelized cost of treatment including revenue of bcp recovery, N2O byproduct and treated water ",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.LCOW,
        name="LCOW",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of photothermal membrane water",
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
        display_units="$/m3 of photothermal membrane water",
        rounding=3,
        description="Levelized cost of water including revenue of revenue of BCP recovery and N2O byproduct",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.LC_BCP,
        name="Levelized cost of BCP",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg-BCP",
        rounding=3,
        description="Levelized cost of bioconcentrated phosphorous including operating and capital costs",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.LC_N2O,
        name="Levelized cost of nitrous oxide",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg-nitrous oxide",
        rounding=3,
        description="Levelized cost of nitrous oxide including operating and capital costs",
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
        fs.pump.costing.direct_capital_cost
        + fs.photothermal.costing.direct_capital_cost
        + fs.candop.costing.direct_capital_cost
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
        fs.candop_treated.properties[0].flow_vol / fs.feed.properties[0].flow_vol
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
    removal_H2O = (
        1
        - fs.candop_treated.properties[0].flow_mass_comp["H2O"]
        / fs.feed.properties[0].flow_mass_comp["H2O"]
    )
    exports.add(
        obj=removal_H2O,
        name="H2O removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Water removal fraction [1 - outlet H2O flow/inlet H2O flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_N2 = (
        1
        - fs.candop_treated.properties[0].flow_mass_comp["nitrogen"]
        / fs.feed.properties[0].flow_mass_comp["nitrogen"]
    )
    exports.add(
        obj=removal_N2,
        name="N2 removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Nitrogen removal fraction [1 - outlet N2 flow/inlet N2 flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_P = (
        1
        - fs.candop_treated.properties[0].flow_mass_comp["phosphates"]
        / fs.feed.properties[0].flow_mass_comp["phosphates"]
    )
    exports.add(
        obj=removal_P,
        name="Phosphates removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Phosphates removal fraction [1 - outlet P flow/inlet P flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_BCP = fs.candop.removal_frac_mass_comp[0, "bioconcentrated_phosphorous"]
    exports.add(
        obj=removal_BCP,
        name="BCP removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Bioconcentrated phosphorous removal fraction",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    N2O_production = (
        fs.candop_byproduct.properties[0].flow_mass_comp["nitrous_oxide"]
        / fs.feed.properties[0].flow_vol
    )
    exports.add(
        obj=N2O_production,
        name="N2O production",
        ui_units=pyunits.kg / pyunits.m**3,
        display_units="kg-N2O/m3 of feed",
        rounding=3,
        description="Nitrous oxide production [N2O product flow rate/feed flow rate]",
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
        obj=fs.photothermal.costing.capital_cost,
        name="Photothermal membrane",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Photothermal membrane",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.candop.costing.capital_cost,
        name="CANDO_P",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="CANDO_P",
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
    total_revenue = (
        fs.costing.value_bcp_recovery
        + fs.costing.value_N2O_byproduct
        + fs.costing.value_water_byproduct
    )
    exports.add(
        obj=total_revenue,
        name="Total",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Total revenue - including the sale of BCP, N2O, and water",
        is_input=False,
        is_output=True,
        output_category="Revenue",
    )

    exports.add(
        obj=fs.costing.value_bcp_recovery,
        name="BCP",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Revenue from selling bioconcentrated phosphorous",
        is_input=False,
        is_output=True,
        output_category="Revenue",
    )
    exports.add(
        obj=fs.costing.value_N2O_byproduct,
        name="Nitrous oxide",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Revenue from selling nitrous oxide",
        is_input=False,
        is_output=True,
        output_category="Revenue",
    )
    exports.add(
        obj=fs.costing.value_water_byproduct,
        name="Water",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Revenue from selling treated water",
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
