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
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.amo_1690.amo_1690 import (
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
        name="AMO 1690",
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
        ui_units=pyunits.gal / pyunits.day,
        display_units="gal/day",
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
        obj=fs.feed.conc_mass_comp[0, "cod"],
        name="COD concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=2,
        description="Inlet COD concentration",
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
        rounding=2,
        description="Inlet TKN concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "acetic_acid"],
        name="Acetic acid concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=2,
        description="Inlet acetic acid concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "ammonium_as_nitrogen"],
        name="NH4 concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=2,
        description="Inlet ammonium concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )

    # Unit model data, cloth media filtration
    exports.add(
        obj=fs.cmf.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Cloth media filtration",
        is_output=False,
    )
    exports.add(
        obj=fs.cmf.energy_electric_flow_vol_inlet,
        name="Electricity intensity",
        ui_units=pyunits.kWh / pyunits.m**3,
        display_units="kWh/m3",
        rounding=2,
        description="Specific electricity intensity relating the intensity to the volume of product",
        is_input=True,
        input_category="Cloth media filtration",
        is_output=False,
    )
    exports.add(
        obj=fs.cmf.removal_frac_mass_comp[0, "tss"],
        name="TSS removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Total suspended solids removal [g-TSS removed/g-TSS inlet]",
        is_input=True,
        input_category="Cloth media filtration",
        is_output=False,
    )
    exports.add(
        obj=fs.cmf.removal_frac_mass_comp[0, "cod"],
        name="COD removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="COD removal [g-COD removed/g-COD inlet]",
        is_input=True,
        input_category="Cloth media filtration",
        is_output=False,
    )
    exports.add(
        obj=fs.cmf.removal_frac_mass_comp[0, "tkn"],
        name="TKN removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="TKN removal [g-TKN removed/g-TKN inlet]",
        is_input=True,
        input_category="Cloth media filtration",
        is_output=False,
    )

    # Unit cost data, cloth media filtration
    exports.add(
        obj=fs.costing.cloth_media_filtration.sizing_cost[None],
        name="Cloth media filtration cost",
        ui_units=(fs.costing.base_currency * pyunits.day) / pyunits.gal,
        display_units="($ * day)/gal",
        rounding=0,
        description="Cloth media filtration capital cost parameter",
        is_input=True,
        input_category="Cloth media filtration costing",
        is_output=False,
    )

    # Unit model data, anaerobic digestion
    exports.add(
        obj=fs.ad.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Anaerobic digestion",
        is_output=False,
    )
    exports.add(
        obj=fs.ad.reaction_conversion[0, "tss_reaction"],
        name="TSS conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Total suspended solids conversion [g-TSS reacted/g-TSS inlet]",
        is_input=True,
        input_category="Anaerobic digestion",
        is_output=False,
    )
    exports.add(
        obj=fs.ad.reaction_conversion[0, "cod_reaction"],
        name="COD conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="COD conversion [g-COD reacted/g-COD inlet]",
        is_input=True,
        input_category="Anaerobic digestion",
        is_output=False,
    )
    exports.add(
        obj=fs.ad.reaction_conversion[0, "tkn_reaction"],
        name="TKN conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="TKN conversion [g-TKN reacted/g-TKN inlet]",
        is_input=True,
        input_category="Anaerobic digestion",
        is_output=False,
    )
    exports.add(
        obj=fs.ad.energy_electric_flow_vol_inlet,
        name="Electricity intensity",
        ui_units=pyunits.kWh / pyunits.m**3,
        display_units="kWh/m3",
        rounding=2,
        description="Specific electricity intensity relating the intensity to the volume of product",
        is_input=True,
        input_category="Anaerobic digestion",
        is_output=False,
    )
    exports.add(
        obj=fs.ad.biogas_tss_ratio,
        name="Biogas:TSS ratio",
        ui_units=pyunits.m**3 / pyunits.kg,
        display_units="m3/kg",
        rounding=2,
        description="Ratio of biogas volume to total suspended solids mass",
        is_input=True,
        input_category="Anaerobic digestion",
        is_output=False,
    )

    # Unit model data, membrane evaporator
    exports.add(
        obj=fs.me.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Membrane evaporator",
        is_output=False,
    )
    exports.add(
        obj=fs.me.removal_frac_mass_comp[0, "ammonium_as_nitrogen"],
        name="NH4 removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Ammonium removal [g-NH4 removed/g-NH4 inlet]",
        is_input=True,
        input_category="Membrane evaporator",
        is_output=False,
    )
    exports.add(
        obj=fs.me.removal_frac_mass_comp[0, "acetic_acid"],
        name="Acetic acid removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Acetic acid removal [g-CH3COOH removed/g-CH3COOH inlet]",
        is_input=True,
        input_category="Membrane evaporator",
        is_output=False,
    )
    exports.add(
        obj=fs.me.energy_electric_flow_vol_inlet,
        name="Electricity intensity",
        ui_units=pyunits.kWh / pyunits.m**3,
        display_units="kWh/m3",
        rounding=2,
        description="Specific electricity intensity relating the intensity to the volume of product",
        is_input=True,
        input_category="Membrane evaporator",
        is_output=False,
    )
    exports.add(
        obj=fs.me.water_flux,
        name="Water flux",
        ui_units=pyunits.m / pyunits.hr,
        display_units="m/h",
        rounding=2,
        description="Water flux",
        is_input=True,
        input_category="Membrane evaporator",
        is_output=False,
    )

    # Unit cost data, membrane evaporator
    exports.add(
        obj=fs.costing.membrane_evaporator.membrane_cost[None],
        name="Membrane cost",
        ui_units=fs.costing.base_currency / pyunits.m**2,
        display_units="$/m2 of membrane",
        rounding=0,
        description="Membrane evaporator capital cost parameter",
        is_input=True,
        input_category="Membrane evaporator costing",
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
        obj=fs.costing.biogas_cost,
        name="Biogas cost",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3",
        rounding=3,
        description="Biogas cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.fertilizer_cost,
        name="Fertilizer cost",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=3,
        description="Fertilizer cost",
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
        obj=fs.filtered_water.properties[0].conc_mass_comp["tss"],
        name="Filtered water TSS concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=5,
        description="Outlet filtered water total suspended solids concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.filtered_water.properties[0].conc_mass_comp["cod"],
        name="Filtered water COD concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=5,
        description="Outlet filtered water COD concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.filtered_water.properties[0].conc_mass_comp["tkn"],
        name="Filtered water TKN concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=5,
        description="Outlet filtered water TKN concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.ad_byproduct.properties[0].flow_vol,
        name="AD byproduct flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Outlet anaerobic digestion byproduct flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.ad_byproduct.properties[0].conc_mass_comp["H2O"],
        name="AD byproduct water concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=5,
        description="Outlet AD byproduct water concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.ad_byproduct.properties[0].conc_mass_comp["tss"],
        name="AD byproduct TSS concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=5,
        description="Outlet AD byproduct total suspended solids concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.ad_byproduct.properties[0].conc_mass_comp["cod"],
        name="AD byproduct COD concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=5,
        description="Outlet AD byproduct COD concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.ad_byproduct.properties[0].conc_mass_comp["tkn"],
        name="AD byproduct TKN concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=5,
        description="Outlet AD byproduct TKN concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.me_byproduct.properties[0].flow_vol,
        name="Membrane evaporator byproduct flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Outlet membrane evaporator byproduct flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.me_byproduct.properties[0].conc_mass_comp["acetic_acid"],
        name="Membrane evaporator byproduct acetic acid concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=5,
        description="Outlet membrane evaporator byproduct acetic acid concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.me_byproduct.properties[0].conc_mass_comp["ammonium_as_nitrogen"],
        name="Membrane evaporator byproduct ammonium concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=5,
        description="Outlet membrane evaporator byproduct ammonium concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.me_treated.properties[0].flow_vol,
        name="Membrane evaporator treated flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=2,
        description="Outlet membrane evaporator treated flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.me_treated.properties[0].conc_mass_comp["acetic_acid"],
        name="Membrane evaporator treated acetic_acid concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=5,
        description="Outlet membrane evaporator treated acetic acid concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.me_treated.properties[0].conc_mass_comp["ammonium_as_nitrogen"],
        name="Membrane evaporator treated ammonium concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=5,
        description="Outlet membrane evaporator treated ammonium concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )

    # System metrics
    exports.add(
        obj=fs.costing.LCOT,
        name="LCOT",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of feed water",
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
        display_units="$/m3 of feed water",
        rounding=3,
        description="Levelized cost of treatment including revenue of biogas and fertilizer",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.LC_biogas,
        name="LCOB",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of biogas",
        rounding=3,
        description="Levelized cost of biogas including operating and capital costs",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.LC_biogas_with_revenue,
        name="LCOB with revenue",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of biogas",
        rounding=3,
        description="Levelized cost of biogas including revenue of fertilizer",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.LC_fertilizer,
        name="LCOF",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg of fertilizer",
        rounding=3,
        description="Levelized cost of fertilizer including operating and capital costs",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.LC_fertilizer_with_revenue,
        name="LCOF with revenue",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg of fertilizer",
        rounding=3,
        description="Levelized cost of fertilizer including revenue of biogas",
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
        fs.cmf.costing.direct_capital_cost
        + fs.ad.costing.direct_capital_cost
        + fs.me.costing.direct_capital_cost
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
        fs.filtered_water.properties[0].flow_vol / fs.feed.properties[0].flow_vol
    )
    exports.add(
        obj=recovery_vol,
        name="Volumetric recovery",
        ui_units=pyunits.dimensionless,
        display_units="m3 of filtered water/m3 of feed",
        rounding=3,
        description="Normalized volumetric recovery",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_TSS = (
        1
        - fs.filtered_water.properties[0].flow_mass_comp["tss"]
        / fs.feed.properties[0].flow_mass_comp["tss"]
    )
    exports.add(
        obj=removal_TSS,
        name="TSS removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Total suspended solids removal fraction [1 - outlet TSS flow/inlet TSS flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_COD = (
        1
        - fs.filtered_water.properties[0].flow_mass_comp["cod"]
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
    removal_TKN = (
        1
        - fs.filtered_water.properties[0].flow_mass_comp["tkn"]
        / fs.feed.properties[0].flow_mass_comp["tkn"]
    )
    exports.add(
        obj=removal_TKN,
        name="TKN removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="TKN removal fraction [1 - outlet TKN flow/inlet TKN flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    biogas_recovery_volume = fs.costing.utilization_factor * fs.ad.biogas_production[0]
    exports.add(
        obj=biogas_recovery_volume,
        name="Biogas recovery",
        ui_units=pyunits.m**3 / pyunits.year,
        display_units="m3-biogas/year",
        rounding=0,
        description="Biogas recovery volume",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    fertilizer_recovery_mass = fs.costing.utilization_factor * (
        fs.me_byproduct.flow_mass_comp[0, "ammonium_as_nitrogen"]
        + fs.me_byproduct.flow_mass_comp[0, "acetic_acid"]
    )
    exports.add(
        obj=fertilizer_recovery_mass,
        name="Fertilizer recovery",
        ui_units=pyunits.kg / pyunits.year,
        display_units="kg-fertilizer/year",
        rounding=0,
        description="Fertilizer recovery mass",
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
        obj=fs.cmf.costing.capital_cost,
        name="Cloth media filtration",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Cloth media filtration",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.ad.costing.capital_cost,
        name="Anaerobic digestion",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Anaerobic digestion",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.me.costing.capital_cost,
        name="Membrane evaporator",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Membrane evaporator",
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
    total_revenue = (biogas_recovery_volume * fs.costing.biogas_cost) + (
        fertilizer_recovery_mass * fs.costing.fertilizer_cost
    )
    exports.add(
        obj=total_revenue,
        name="Total",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Total revenue - including the sale of struvite and purchase of MgCl2 and dry polymer",
        is_input=False,
        is_output=True,
        output_category="Revenue",
    )
    biogas_revenue = biogas_recovery_volume * fs.costing.biogas_cost
    exports.add(
        obj=biogas_revenue,
        name="Biogas",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Revenue from selling biogas",
        is_input=False,
        is_output=True,
        output_category="Revenue",
    )
    fertilizer_revenue = fertilizer_recovery_mass * fs.costing.fertilizer_cost
    exports.add(
        obj=fertilizer_revenue,
        name="Fertilizer",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Revenue from selling fertilizer",
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
