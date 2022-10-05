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
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.metab.metab import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    add_costing,
    adjust_default_parameters,
)
from idaes.core.solvers import get_solver
from pyomo.environ import units as pyunits, assert_optimal_termination
from pyomo.util.check_units import assert_units_consistent


def export_to_ui():
    return FlowsheetInterface(
        name="Metab",
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
        obj=fs.feed.conc_mass_comp[0, "cod"],
        name="COD concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",  # can this be done by default?
        rounding=2,
        description="Inlet chemical oxygen demand (COD) concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )

    # Unit model data, hydrogen reactor
    exports.add(
        obj=fs.metab_hydrogen.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Hydrogen reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_hydrogen.reaction_conversion[0, "cod_to_hydrogen"],
        name="COD conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="COD conversion [g-COD reacted/g-COD inlet]",
        is_input=True,
        input_category="Hydrogen reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_hydrogen.generation_ratio["cod_to_hydrogen", "hydrogen"],
        name="H2 conversion ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="H2 mass conversion ratio with respect to COD [g-H2 produced/g-COD"
        " reacted]",
        is_input=True,
        input_category="Hydrogen reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_hydrogen.hydraulic_retention_time,
        name="HRT",
        ui_units=pyunits.hr,
        display_units="h",
        rounding=1,
        description="Hydraulic retention time",
        is_input=True,
        input_category="Hydrogen reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_hydrogen.energy_electric_mixer_vol,
        name="Mixer specific power",
        ui_units=pyunits.kW / pyunits.m**3,
        display_units="kW/m3 of reactor",
        rounding=3,
        description="Mixer specific power relating the power to the volume of the "
        "reactor",
        is_input=True,
        input_category="Hydrogen reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_hydrogen.energy_electric_vacuum_flow_vol_byproduct,
        name="Vacuum specific power",
        ui_units=pyunits.kW / (pyunits.kg / pyunits.hr),
        display_units="kW/(kg-H2/h)",
        rounding=2,
        description="Vacuum specific power relating the power to the production of H2",
        is_input=True,
        input_category="Hydrogen reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_hydrogen.energy_thermal_flow_vol_inlet,
        name="Specific heating",
        ui_units=pyunits.MJ / pyunits.m**3,
        display_units="MJ/m3 of water",
        rounding=2,
        description="Specific heating relating the thermal energy input to water "
        "volume",
        is_input=True,
        input_category="Hydrogen reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.bead_bulk_density["hydrogen"],
        name="Bead bulk density",
        ui_units=pyunits.kg / pyunits.m**3,
        display_units="kg/m3 of reactor",
        rounding=1,
        description="Bead bulk density [kg-beads/m3 of reactor]",
        is_input=True,
        input_category="Hydrogen reactor",
        is_output=False,
    )

    # Unit cost data, hydrogen reactor
    exports.add(
        obj=fs.costing.metab.reactor_cost["hydrogen"],
        name="Reactor cost",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of reactor",
        rounding=0,
        description="Reactor capital cost parameter",
        is_input=True,
        input_category="Hydrogen reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.mixer_cost["hydrogen"],
        name="Mixer cost",
        ui_units=fs.costing.base_currency / pyunits.kW,
        display_units="$/kW of mixer",
        rounding=0,
        description="Mixer capital cost parameter",
        is_input=True,
        input_category="Hydrogen reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.bead_cost["hydrogen"],
        name="Bead cost",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=1,
        description="Bead cost parameter",
        is_input=True,
        input_category="Hydrogen reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.bead_replacement_factor["hydrogen"],
        name="Bead replacement factor",
        ui_units=1 / pyunits.year,
        display_units="1/year",
        rounding=2,
        description="Bead replacement factor - amount of initial beads replaced per "
        "year",
        is_input=True,
        input_category="Hydrogen reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.membrane_sidestream_fraction["hydrogen"],
        name="Membrane side stream",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Fraction of reactor volumetric flow that recirculates in membrane "
        "side stream",
        is_input=True,
        input_category="Hydrogen reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.membrane_specific_size["hydrogen"],
        name="Membrane specific size",
        ui_units=pyunits.m**2 / (pyunits.m**3 / pyunits.hr),
        display_units="m2/(m3/h of side stream)",
        rounding=2,
        description="Membrane specific size relating membrane area to side stream "
        "flow rate",
        is_input=True,
        input_category="Hydrogen reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.membrane_cost["hydrogen"],
        name="Membrane cost",
        ui_units=fs.costing.base_currency / pyunits.m**2,
        display_units="$/m2",
        rounding=0,
        description="Membrane cost parameter",
        is_input=True,
        input_category="Hydrogen reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.vacuum_cost["hydrogen"],
        name="Vacuum cost",
        ui_units=fs.costing.base_currency / (pyunits.kg / pyunits.hr),
        display_units="$/(kg-H2/h)",
        rounding=0,
        description="Vacuum cost parameter",
        is_input=True,
        input_category="Hydrogen reactor costing",
        is_output=False,
    )

    # Unit model data, methane reactor
    exports.add(
        obj=fs.metab_methane.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Methane reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_methane.reaction_conversion[0, "cod_to_methane"],
        name="COD conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",  # we should change to %
        rounding=2,
        description="COD conversion [g-COD reacted/g-COD inlet]",
        is_input=True,
        input_category="Methane reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_methane.generation_ratio["cod_to_methane", "methane"],
        name="CH4 conversion ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="CH4 mass conversion ratio with respect to COD [g-CH4 "
        "produced/g-COD reacted]",
        is_input=True,
        input_category="Methane reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_methane.hydraulic_retention_time,
        name="HRT",
        ui_units=pyunits.hr,
        display_units="h",
        rounding=1,
        description="Hydraulic retention time",
        is_input=True,
        input_category="Methane reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_methane.energy_electric_mixer_vol,
        name="Mixer specific power",
        ui_units=pyunits.kW / pyunits.m**3,
        display_units="kW/m3 of reactor",
        rounding=3,
        description="Mixer specific power relating the power to the volume of the "
        "reactor",
        is_input=True,
        input_category="Methane reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_methane.energy_electric_vacuum_flow_vol_byproduct,
        name="Vacuum specific power",
        ui_units=pyunits.kW / (pyunits.kg / pyunits.hr),
        display_units="kW/(kg-H2/h)",
        rounding=2,
        description="Vacuum specific power relating the power to the production of H2",
        is_input=True,
        input_category="Methane reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.metab_methane.energy_thermal_flow_vol_inlet,
        name="Specific heating",
        ui_units=pyunits.MJ / pyunits.m**3,
        display_units="MJ/m3 of water",
        rounding=2,
        description="Specific heating relating the thermal energy input to water "
        "volume",
        is_input=True,
        input_category="Methane reactor",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.bead_bulk_density["methane"],
        name="Bead bulk density",
        ui_units=pyunits.kg / pyunits.m**3,
        display_units="kg/m3 of reactor",
        rounding=1,
        description="Bead bulk density [kg-beads/m3 of reactor]",
        is_input=True,
        input_category="Methane reactor",
        is_output=False,
    )

    # Unit cost data, methane reactor
    exports.add(
        obj=fs.costing.metab.reactor_cost["methane"],
        name="Reactor cost",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of reactor",
        rounding=0,
        description="Reactor capital cost parameter",
        is_input=True,
        input_category="Methane reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.mixer_cost["methane"],
        name="Mixer cost",
        ui_units=fs.costing.base_currency / pyunits.kW,
        display_units="$/kW of mixer",
        rounding=0,
        description="Mixer capital cost parameter",
        is_input=True,
        input_category="Methane reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.bead_cost["methane"],
        name="Bead cost",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=1,
        description="Bead cost parameter",
        is_input=True,
        input_category="Methane reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.bead_replacement_factor["methane"],
        name="Bead replacement factor",
        ui_units=1 / pyunits.year,
        display_units="1/year",
        rounding=2,
        description="Bead replacement factor - amount of initial beads replaced per "
        "year",
        is_input=True,
        input_category="Methane reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.membrane_sidestream_fraction["methane"],
        name="Membrane side stream",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Fraction of reactor volumetric flow that recirculates in membrane "
        "side stream",
        is_input=True,
        input_category="Methane reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.membrane_specific_size["methane"],
        name="Membrane specific size",
        ui_units=pyunits.m**2 / (pyunits.m**3 / pyunits.hr),
        display_units="m2/(m3/h of side stream)",
        rounding=2,
        description="Membrane specific size relating membrane area to side stream "
        "flow rate",
        is_input=True,
        input_category="Methane reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.membrane_cost["methane"],
        name="Membrane cost",
        ui_units=fs.costing.base_currency / pyunits.m**2,
        display_units="$/m2",
        rounding=0,
        description="Membrane cost parameter",
        is_input=True,
        input_category="Methane reactor costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.metab.vacuum_cost["methane"],
        name="Vacuum cost",
        ui_units=fs.costing.base_currency / (pyunits.kg / pyunits.hr),
        display_units="$/(kg-CH4/h)",
        rounding=0,
        description="Vacuum cost parameter",
        is_input=True,
        input_category="Methane reactor costing",
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
        obj=fs.costing.heat_cost,
        name="Heating cost",
        ui_units=fs.costing.base_currency / pyunits.MJ,
        display_units="$/MJ",
        rounding=3,
        description="Heating cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.hydrogen_product_cost,
        name="Hydrogen cost",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=3,
        description="Hydrogen cost is negative because it is sold",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.methane_product_cost,
        name="Methane cost",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=3,
        description="Methane cost is negative because it is sold",
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
        obj=fs.product_H2O.properties[0].conc_mass_comp["cod"],
        name="Product water COD concentration",
        ui_units=pyunits.g / pyunits.L,
        display_units="g/L",
        rounding=2,
        description="Outlet product water chemical oxygen demand concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_hydrogen.properties[0].flow_mass_comp["hydrogen"],
        name="Hydrogen product flow rate",
        ui_units=pyunits.kg / pyunits.hr,
        display_units="kg-H2/h",
        rounding=3,
        description="Outlet hydrogen product flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_methane.properties[0].flow_mass_comp["methane"],
        name="Methane product flow rate",
        ui_units=pyunits.kg / pyunits.hr,
        display_units="kg-CH4/h",
        rounding=3,
        description="Outlet methane product flow rate",
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
    exports.add(
        obj=fs.costing.LCOH,
        name="Levelized cost of hydrogen",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg-H2",
        rounding=2,
        description="Levelized cost of hydrogen including operating and capital costs",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.LCOM,
        name="Levelized cost of methane",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg-CH4",
        rounding=2,
        description="Levelized cost of methane including operating and capital costs",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.LCOCR,
        name="Levelized cost of COD removal",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg-COD removed",
        rounding=2,
        description="Levelized cost of chemical oxygen demand removal including "
        "operating and capital costs",
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
        fs.metab_hydrogen.costing.direct_capital_cost
        + fs.metab_methane.costing.direct_capital_cost
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
    heat_operating_norm = (
        fs.costing.aggregate_flow_costs["heat"] / fs.costing.annual_water_inlet
    )
    exports.add(
        obj=heat_operating_norm,
        name="Heating",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of feed",
        rounding=2,
        description="Normalized heating cost - [annual heating costs/annual "
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
        description="Volumetric recovery of product water",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_cod = (
        1
        - fs.product_H2O.properties[0].flow_mass_comp["cod"]
        / fs.feed.properties[0].flow_mass_comp["cod"]
    )
    exports.add(
        obj=removal_cod,
        name="COD removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="COD removal fraction [1 - outlet COD flow/inlet COD flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    hydrogen_production = (
        fs.product_hydrogen.properties[0].flow_mass_comp["hydrogen"]
        / fs.feed.properties[0].flow_vol
    )
    exports.add(
        obj=hydrogen_production,
        name="Hydrogen production",
        ui_units=pyunits.kg / pyunits.m**3,
        display_units="kg-H2/m3 of feed",
        rounding=3,
        description="Hydrogen production [hydrogen product flow rate/feed flow rate]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )

    methane_production = (
        fs.product_methane.properties[0].flow_mass_comp["methane"]
        / fs.feed.properties[0].flow_vol
    )
    exports.add(
        obj=methane_production,
        name="Methane production",
        ui_units=pyunits.kg / pyunits.m**3,
        display_units="kg-CH4/m3 of feed",
        rounding=3,
        description="Methane production [methane product flow rate/feed flow rate]",
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
        obj=fs.metab_hydrogen.costing.capital_cost,
        name="Hydrogen reactor",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Hydrogen reactor",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.metab_methane.costing.capital_cost,
        name="Methane reactor",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Methane reactor",
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
        obj=fs.costing.aggregate_flow_costs["heat"],
        name="Heating",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Annual heating costs ",
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
        fs.costing.aggregate_flow_costs["hydrogen_product"]
        + fs.costing.aggregate_flow_costs["methane_product"]
    )
    exports.add(
        obj=total_revenue,
        name="Total",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Total revenue - including the sale of hydrogen and methane",
        is_input=False,
        is_output=True,
        output_category="Revenue",
    )
    exports.add(
        obj=-fs.costing.aggregate_flow_costs["hydrogen_product"],
        name="Hydrogen",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Revenue from selling hydrogen",
        is_input=False,
        is_output=True,
        output_category="Revenue",
    )
    exports.add(
        obj=-fs.costing.aggregate_flow_costs["methane_product"],
        name="Methane",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Revenue from selling methane",
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

    adjust_default_parameters(m)

    results = solve(m)
    assert_optimal_termination(results)
    return m.fs


def solve_flowsheet(flowsheet=None):
    fs = flowsheet
    results = solve(fs)
    return results
