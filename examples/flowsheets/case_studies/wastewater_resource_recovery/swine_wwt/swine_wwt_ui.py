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
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.swine_wwt.swine_wwt import (
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
        name="Swine wastewater treatment",
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
        ui_units=pyunits.L / pyunits.day,
        display_units="L/day",
        rounding=0,
        description="Inlet volumetric flow rate",
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
        rounding=0,
        description="Inlet chemical oxygen demand (COD) concentration",
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
        rounding=0,
        description="Inlet NH4-N concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.feed.conc_mass_comp[0, "phosphates"],
        name="PO4 concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=0,
        description="Inlet PO4 concentration",
        is_input=True,
        input_category="Feed",
        is_output=True,
        output_category="Feed",
    )
    exports.add(
        obj=fs.food_waste.conc_mass_comp[0, "cod"],
        name="COD concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=0,
        description="Food waste chemical oxygen demand (COD) concentration",
        is_input=True,
        input_category="Food waste",
        is_output=True,
        output_category="Food waste",
    )

    # Unit model data, anaerobic MBR-MEC
    exports.add(
        obj=fs.mbr_mec.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Anaerobic MBR-MEC",
        is_output=False,
    )
    exports.add(
        obj=fs.mbr_mec.reaction_conversion[0, "cod_to_nonbiodegradable_cod"],
        name="COD conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="COD conversion [g-COD reacted/g-COD inlet]",
        is_input=True,
        input_category="Anaerobic MBR-MEC",
        is_output=False,
    )
    exports.add(
        obj=fs.mbr_mec.generation_ratio[
            "cod_to_nonbiodegradable_cod", "nonbiodegradable_cod"
        ],
        name="Nonbiodigradable COD conversion ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Nonbiodigradable COD mass conversion ratio with respect to COD "
        "[g-Nonbiodigradable COD produced/g-COD reacted]",
        is_input=True,
        input_category="Anaerobic MBR-MEC",
        is_output=False,
    )

    # Unit cost data, anaerobic MBR-MEC
    exports.add(
        obj=fs.costing.anaerobic_mbr_mec.unit_capex[None],
        name="Unit capital expenditure",
        ui_units=fs.costing.base_currency / (pyunits.L / pyunits.day),
        display_units="$/(L/day)",
        rounding=2,
        description="Unit capital expenditure cost parameter",
        is_input=True,
        input_category="Anaerobic MBR-MEC costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.anaerobic_mbr_mec.unit_opex[None],
        name="Unit operational expenditure",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3",
        rounding=2,
        description="Unit operational expenditure cost parameter",
        is_input=True,
        input_category="Anaerobic MBR-MEC costing",
        is_output=False,
    )

    # Unit model data, gas-sparged membrane
    exports.add(
        obj=fs.gas_sparged_membrane.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Gas-sparged membrane",
        is_output=False,
    )
    exports.add(
        obj=fs.gas_sparged_membrane.gas_mass_influent_ratio[0],
        name="Gas mass influent ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=8,
        description="Gas mass influent ratio",
        is_input=True,
        input_category="Gas-sparged membrane",
        is_output=False,
    )

    # Unit model data, VFA separation
    exports.add(
        obj=fs.vfa_recovery.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="VFA separation",
        is_output=False,
    )
    exports.add(
        obj=fs.vfa_recovery.removal_frac_mass_comp[0, "nonbiodegradable_cod"],
        name="Nonbiodegradable COD removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Nonbiodegradable COD removal [g-nonbiodegradable COD byproduct/ g-nonbiodegradable COD inlet]",
        is_input=True,
        input_category="VFA separation",
        is_output=False,
    )
    exports.add(
        obj=fs.vfa_recovery.heat_required_per_vfa_mass[0],
        name="Heat demand",
        ui_units=pyunits.kJ / pyunits.kg,
        display_units="kJ/kg",
        rounding=9,
        description="Heat demand per VFA mass",
        is_input=True,
        input_category="VFA separation",
        is_output=False,
    )

    # Unit cost data, VFA separation
    exports.add(
        obj=fs.costing.vfa_recovery.unit_capex[None],
        name="Unit capital expenditure",
        ui_units=fs.costing.base_currency / (pyunits.L / pyunits.day),
        display_units="$/(L/day)",
        rounding=2,
        description="Unit capital expenditure cost parameter",
        is_input=True,
        input_category="VFA separation costing",
        is_output=False,
    )

    # Unit model data, co-fermentation
    exports.add(
        obj=fs.cofermentation.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Co-fermentation",
        is_output=False,
    )
    exports.add(
        obj=fs.cofermentation.reaction_conversion[0, "cod_to_nonbiodegradable_cod"],
        name="COD conversion",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="COD conversion [g-COD reacted/g-COD inlet]",
        is_input=True,
        input_category="Co-fermentation",
        is_output=False,
    )
    exports.add(
        obj=fs.cofermentation.generation_ratio[
            "cod_to_nonbiodegradable_cod", "nonbiodegradable_cod"
        ],
        name="Nonbiodigradable COD conversion ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Nonbiodigradable COD mass conversion ratio with respect to COD "
        "[g-Nonbiodigradable COD produced/g-COD reacted]",
        is_input=True,
        input_category="Co-fermentation",
        is_output=False,
    )

    # Unit cost data, co-fermentation
    exports.add(
        obj=fs.costing.cofermentation.unit_capex[None],
        name="Unit capital expenditure",
        ui_units=fs.costing.base_currency / (pyunits.L / pyunits.day),
        display_units="$/(L/day)",
        rounding=2,
        description="Unit capital expenditure cost parameter",
        is_input=True,
        input_category="Co-fermentation costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.cofermentation.unit_opex[None],
        name="Unit operational expenditure",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3",
        rounding=2,
        description="Unit operational expenditure cost parameter",
        is_input=True,
        input_category="Co-fermentation costing",
        is_output=False,
    )

    # Unit model data, ion exchange
    exports.add(
        obj=fs.ion_exchange.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Ion exchange",
        is_output=False,
    )
    exports.add(
        obj=fs.ion_exchange.removal_frac_mass_comp[0, "ammonium_as_nitrogen"],
        name="NH4-N removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="NH4-N removal [g-NH4-N byproduct/ g-NH4-N inlet]",
        is_input=True,
        input_category="Ion exchange",
        is_output=False,
    )
    exports.add(
        obj=fs.ion_exchange.nitrogen_clay_ratio[0],
        name="Nitrogen clay ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Mass fraction of nitrogen in clay mixture",
        is_input=True,
        input_category="Ion exchange",
        is_output=False,
    )
    exports.add(
        obj=fs.ion_exchange.NaCl_dose,
        name="Dosage of NaCl addition",
        ui_units=pyunits.kg / pyunits.m**3,
        display_units="kg/m3",
        rounding=2,
        description="Dosage of NaCl addition",
        is_input=True,
        input_category="Ion exchange",
        is_output=False,
    )
    exports.add(
        obj=fs.ion_exchange.resin_replacement,
        name="Resin replacement as a function of flow",
        ui_units=pyunits.kg / pyunits.m**3,
        display_units="kg/m3",
        rounding=2,
        description="Resin replacement as a function of flow",
        is_input=True,
        input_category="Ion exchange",
        is_output=False,
    )

    # Unit cost data, ion exchange
    exports.add(
        obj=fs.costing.ion_exchange.unit_capex["clinoptilolite"],
        name="Unit capital expenditure",
        ui_units=fs.costing.base_currency / (pyunits.L / pyunits.day),
        display_units="$/(L/day)",
        rounding=2,
        description="Unit capital expenditure cost parameter",
        is_input=True,
        input_category="Ion exchange costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.ion_exchange.unit_opex["clinoptilolite"],
        name="Unit operational expenditure",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3",
        rounding=2,
        description="Unit operational expenditure cost parameter",
        is_input=True,
        input_category="Ion exchange costing",
        is_output=False,
    )

    # Unit model data, sedimentation
    exports.add(
        obj=fs.sedimentation.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Sedimentation",
        is_output=False,
    )
    exports.add(
        obj=fs.sedimentation.removal_frac_mass_comp[0, "phosphates"],
        name="PO4 removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="PO4 removal [g-PO4 byproduct/ g-PO4 inlet]",
        is_input=True,
        input_category="Sedimentation",
        is_output=False,
    )
    exports.add(
        obj=fs.sedimentation.settling_velocity[0],
        name="Particle settling velocity",
        ui_units=pyunits.m / pyunits.s,
        display_units="m/s",
        rounding=3,
        description="Particle settling velocity",
        is_input=True,
        input_category="Sedimentation",
        is_output=False,
    )
    exports.add(
        obj=fs.sedimentation.phosphorus_solids_ratio[0],
        name="Phosphorus solids ratio",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=2,
        description="Mass fraction of phosphorus in settleable solids",
        is_input=True,
        input_category="Sedimentation",
        is_output=False,
    )
    exports.add(
        obj=fs.sedimentation.energy_electric_flow_vol_inlet,
        name="Electricity intensity",
        ui_units=pyunits.kWh / pyunits.m**3,
        display_units="fraction",
        rounding=2,
        description="Electricity intensity with respect to inlet flowrate of unit",
        is_input=True,
        input_category="Sedimentation",
        is_output=False,
    )

    # Unit cost data, sedimentation
    exports.add(
        obj=fs.costing.sedimentation.unit_capex["phosphorus_capture"],
        name="Unit capital expenditure",
        ui_units=fs.costing.base_currency / (pyunits.L / pyunits.day),
        display_units="$/(L/day)",
        rounding=2,
        description="Unit capital expenditure cost parameter",
        is_input=True,
        input_category="Sedimentation costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.sedimentation.unit_opex["phosphorus_capture"],
        name="Unit operational expenditure",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3",
        rounding=2,
        description="Unit operational expenditure cost parameter",
        is_input=True,
        input_category="Sedimentation costing",
        is_output=False,
    )

    # Unit model data, constructed wetlands
    exports.add(
        obj=fs.constructed_wetlands.recovery_frac_mass_H2O[0],
        name="Water recovery",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="Water recovery [g-H2O treated/g-H2O inlet]",
        is_input=True,
        input_category="Constructed wetlands",
        is_output=False,
    )
    exports.add(
        obj=fs.constructed_wetlands.removal_frac_mass_comp[0, "ammonium_as_nitrogen"],
        name="NH4-N removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=3,
        description="NH4-N removal [g-NH4-N byproduct/ g-NH4-N inlet]",
        is_input=True,
        input_category="Constructed wetlands",
        is_output=False,
    )

    # Unit cost data, constructed wetlands
    exports.add(
        obj=fs.costing.constructed_wetlands.unit_capex[None],
        name="Unit capital expenditure",
        ui_units=fs.costing.base_currency / (pyunits.L / pyunits.day),
        display_units="$/(L/day)",
        rounding=2,
        description="Unit capital expenditure cost parameter",
        is_input=True,
        input_category="Constructed wetlands costing",
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
        ui_units=fs.costing.base_currency / pyunits.kWh,
        display_units="$/kWh",
        rounding=3,
        description="Heating cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=fs.costing.waste_disposal_cost,
        name="Waste disposal cost",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3",
        rounding=1,
        description="Waste disposal cost",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=-fs.costing.hydrogen_product_cost,
        name="Hydrogen cost",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=0,
        description="Hydrogen cost is negative because it is sold",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=-fs.costing.ammonia_product_cost,
        name="Ammonia cost",
        ui_units=fs.costing.base_currency / pyunits.short_ton,
        display_units="$/tn",
        rounding=0,
        description="Ammonia cost is negative because it is sold",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=-fs.costing.phosphorus_product_cost,
        name="Phosphorus cost",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=0,
        description="Phosphorus cost is negative because it is sold",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=-fs.costing.vfa_product_cost,
        name="VFA cost",
        ui_units=fs.costing.base_currency / pyunits.metric_ton,
        display_units="$/t",
        rounding=0,
        description="VFA cost is negative because it is sold",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=-fs.costing.water_product_cost,
        name="Water cost",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3",
        rounding=0,
        description="Water cost is negative because it is sold",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )
    exports.add(
        obj=-fs.costing.food_waste_tipping_fee_cost,
        name="Food waste tipping fee",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3",
        rounding=0,
        description="Food waste tipping fee cost is negative because it is revenue",
        is_input=True,
        input_category="System costing",
        is_output=False,
    )

    # Outlets
    exports.add(
        obj=fs.product_water.properties[0].flow_vol,
        name="Product water flow rate",
        ui_units=pyunits.m**3 / pyunits.hr,
        display_units="m3/h",
        rounding=4,
        description="Outlet product water flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_water.properties[0].conc_mass_comp["cod"],
        name="Product water COD concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=4,
        description="Outlet product water chemical oxygen demand (COD) concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_water.properties[0].conc_mass_comp["phosphates"],
        name="Product water total P concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=4,
        description="Outlet product water total P concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_water.properties[0].conc_mass_comp["ammonium_as_nitrogen"],
        name="Product water total N concentration",
        ui_units=pyunits.mg / pyunits.L,
        display_units="mg/L",
        rounding=4,
        description="Outlet product water total N concentration",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_ammonia.properties[0].flow_mass_comp["ammonium_as_nitrogen"],
        name="Ammonia product flow rate",
        ui_units=pyunits.kg / pyunits.hr,
        display_units="kg-N/h",
        rounding=4,
        description="Outlet ammonia product flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )
    exports.add(
        obj=fs.product_phosphate.properties[0].flow_mass_comp["phosphates"],
        name="Phosphate product flow rate",
        ui_units=pyunits.kg / pyunits.hr,
        display_units="kg-PO4/h",
        rounding=4,
        description="Outlet phosphate product flow rate",
        is_input=False,
        is_output=True,
        output_category="Outlets",
    )

    # System metrics
    exports.add(
        obj=fs.costing.levelized_costs.LCOT,
        name="Levelized cost of treatment",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of feed",
        rounding=2,
        description="Levelized cost of treatment with revenue from products",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.levelized_costs.LCOW,
        name="Levelized cost of water",
        ui_units=fs.costing.base_currency / pyunits.m**3,
        display_units="$/m3 of feed",
        rounding=2,
        description="Levelized cost of product water with revenue from products",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.levelized_costs.LCOH2,
        name="Levelized cost of hydrogen",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=2,
        description="Levelized cost of hydrogen with revenue from products",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.levelized_costs.LCON,
        name="Levelized cost of nitrogen",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=2,
        description="Levelized cost of nitrogen with revenue from products",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.levelized_costs.LCOVFA,
        name="Levelized cost of VFA",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=2,
        description="Levelized cost of VFA with revenue from products",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.levelized_costs.LCOP,
        name="Levelized cost of phosphate",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=2,
        description="Levelized cost of phosphate with revenue from products",
        is_input=False,
        is_output=True,
        output_category="Levelized cost metrics",
    )
    exports.add(
        obj=fs.costing.levelized_costs.LCOCOD,
        name="Levelized cost of COD removal",
        ui_units=fs.costing.base_currency / pyunits.kg,
        display_units="$/kg",
        rounding=2,
        description="Levelized cost of COD removal with revenue from products",
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
        fs.mbr_mec.costing.direct_capital_cost
        + fs.vfa_recovery.costing.direct_capital_cost
        + fs.ion_exchange.costing.direct_capital_cost
        + fs.sedimentation.costing.direct_capital_cost
        + fs.cofermentation.costing.direct_capital_cost
        + fs.constructed_wetlands.costing.direct_capital_cost
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
        fs.costing.aggregate_flow_costs["electricity"]
        / fs.costing.annual_production.annual_water_inlet
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
        fs.costing.aggregate_flow_costs["heat"]
        / fs.costing.annual_production.annual_water_inlet
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
        fs.product_water.properties[0].flow_vol / fs.feed.properties[0].flow_vol
    )
    exports.add(
        obj=recovery_vol,
        name="Volumetric recovery",
        ui_units=pyunits.dimensionless,
        display_units="m3 of product/m3 of feed",
        rounding=3,
        description="Volumetric recovery",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_cod = (
        1
        - fs.product_water.properties[0].flow_mass_comp["cod"]
        / fs.feed.properties[0].flow_mass_comp["cod"]
    )
    exports.add(
        obj=removal_cod,
        name="COD removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=1,
        description="COD removal fraction [1 - outlet COD flow/inlet COD flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_TP = (
        1
        - fs.product_water.properties[0].flow_mass_comp["phosphates"]
        / fs.feed.properties[0].flow_mass_comp["phosphates"]
    )
    exports.add(
        obj=removal_TP,
        name="Total P removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=4,
        description="Total P (TP) removal fraction [1 - outlet TP flow/inlet TP flow]",
        is_input=False,
        is_output=True,
        output_category="Normalized performance metrics",
    )
    removal_TN = (
        1
        - fs.product_water.properties[0].flow_mass_comp["ammonium_as_nitrogen"]
        / fs.feed.properties[0].flow_mass_comp["ammonium_as_nitrogen"]
    )
    exports.add(
        obj=removal_TN,
        name="Total N removal",
        ui_units=pyunits.dimensionless,
        display_units="fraction",
        rounding=4,
        description="Total N (TN) removal fraction [1 - outlet TN flow/inlet TN flow]",
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
        obj=fs.mbr_mec.costing.capital_cost,
        name="Anaerobic MBR-MEC",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Anaerobic MBR-MEC",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.vfa_recovery.costing.capital_cost,
        name="VFA separation",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="VFA separation",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.ion_exchange.costing.capital_cost,
        name="Ion exchange",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Ion exchange",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.sedimentation.costing.capital_cost,
        name="Sedimentation",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Sedimentation",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.cofermentation.costing.capital_cost,
        name="Co-fermentation",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Co-fermentation",
        is_input=False,
        is_output=True,
        output_category="Capital costs",
    )
    exports.add(
        obj=fs.constructed_wetlands.costing.capital_cost,
        name="Constructed Wetlands",
        ui_units=fs.costing.base_currency,
        display_units="$",
        rounding=0,
        description="Constructed Wetlands",
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
        obj=fs.costing.aggregate_flow_costs["waste_disposal"],
        name="Waste disposal",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Annual waste disposal costs ",
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
        fs.costing.annual_costs_revenues.annual_hydrogen_revenue
        + fs.costing.annual_costs_revenues.annual_water_revenue
        + fs.costing.annual_costs_revenues.annual_ammonia_revenue
        + fs.costing.annual_costs_revenues.annual_phosphorus_revenue
        + fs.costing.annual_costs_revenues.annual_vfa_revenue
        + fs.costing.annual_costs_revenues.annual_food_waste_revenue
    )
    exports.add(
        obj=total_revenue,
        name="Total",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Total revenue",
        is_input=False,
        is_output=True,
        output_category="Revenue",
    )
    exports.add(
        obj=fs.costing.annual_costs_revenues.annual_hydrogen_revenue,
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
        obj=fs.costing.annual_costs_revenues.annual_water_revenue,
        name="Water",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Revenue from selling water",
        is_input=False,
        is_output=True,
        output_category="Revenue",
    )
    exports.add(
        obj=fs.costing.annual_costs_revenues.annual_vfa_revenue,
        name="VFA",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Revenue from selling VFA",
        is_input=False,
        is_output=True,
        output_category="Revenue",
    )
    exports.add(
        obj=fs.costing.annual_costs_revenues.annual_ammonia_revenue,
        name="Ammonia",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Revenue from selling ammonia",
        is_input=False,
        is_output=True,
        output_category="Revenue",
    )
    exports.add(
        obj=fs.costing.annual_costs_revenues.annual_phosphorus_revenue,
        name="Phosphate",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Revenue from selling phosphate",
        is_input=False,
        is_output=True,
        output_category="Revenue",
    )
    exports.add(
        obj=fs.costing.annual_costs_revenues.annual_food_waste_revenue,
        name="Food waste",
        ui_units=fs.costing.base_currency / pyunits.year,
        display_units="$/year",
        rounding=0,
        description="Revenue from food waste tipping fee",
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
